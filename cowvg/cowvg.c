#include "cowvg.h"
#include <stdlib.h>
#include <float.h>

#define NANOSVG_ALL_COLOR_KEYWORDS    // Include full list of color keywords.
#define NANOSVG_IMPLEMENTATION        // Expands implementation

#include "nanosvg.h"

typedef struct {
    unsigned char r;
    unsigned char g;
    unsigned char b;
    unsigned char a;
} CowColor;

typedef struct {
    float *points;
    int point_counts;
} CowPath;

typedef struct {
    CowPath *paths;
    int path_counts;
    int stroke;
    CowColor color;
    struct CowShape *next;
} CowShape;

typedef struct {
    float x, y;
} CowPoint;

static inline int is_point_inside_edge(CowPoint point, CowPoint vertex1, CowPoint vertex2) {
    return ((vertex1.y > point.y) != (vertex2.y > point.y)) &&
           (point.x < (vertex2.x - vertex1.x) * (point.y - vertex1.y) / (vertex2.y - vertex1.y) + vertex1.x);
}

static CowColor priv_cow_color_from_int(unsigned int code) {
    CowColor color;
    color.r = (code >> 24) & 0xFF;
    color.g = (code >> 16) & 0xFF;
    color.b = (code >> 8) & 0xFF;
    color.a = code & 0xFF;
    return color;
}

static void priv_cow_add_line_to_path(CowPath *path, float x, float y) {
    path->point_counts += 2;
    path->points = realloc(path->points, sizeof(float) * path->point_counts);
    path->points[path->point_counts - 2] = x;
    path->points[path->point_counts - 1] = y;
}

static CowPoint priv_cow_compute_bezier_point(CowPoint p0, CowPoint p1, CowPoint p2, CowPoint p3, float t) {
    double u = 1 - t;
    double tt = t * t;
    double uu = u * u;
    double uuu = uu * u;
    double ttt = tt * t;

    CowPoint point;
    point.x = uuu * p0.x + 3 * uu * t * p1.x + 3 * u * tt * p2.x + ttt * p3.x;
    point.y = uuu * p0.y + 3 * uu * t * p1.y + 3 * u * tt * p2.y + ttt * p3.y;
    return point;
}

static void priv_cow_add_bezier_to_path(CowPath *path, CowPoint p0, CowPoint pc, CowPoint pc1, CowPoint p1) {
    const int segments = 15;
    for (int i = 0; i < segments; ++i) {
        float t = (float) (i + 1) / (float) (segments + 1);
        CowPoint p = priv_cow_compute_bezier_point(p0, pc, pc1, p1, t);
        priv_cow_add_line_to_path(path, p.x, p.y);
    }

    priv_cow_add_line_to_path(path, p1.x, p1.y);
}

static void priv_cow_add_path_to_shape(CowShape *shape, CowPath *path) {
    ++shape->path_counts;
    shape->paths = realloc(shape->paths, sizeof(CowPath) * shape->path_counts);
    shape->paths[shape->path_counts - 1] = *path;
}

static CowShape *priv_cow_make_lines_from_svg(CowSvg svg) {
    CowShape *current_shape = NULL;
    CowShape *first_shape = NULL;
    CowPoint lastPos;
    for (NSVGshape *shape = svg.svg->shapes; shape != NULL; shape = shape->next) {
        if (current_shape == NULL) {
            current_shape = calloc(sizeof(CowShape), 1);
            current_shape->color = priv_cow_color_from_int(shape->fill.color);
            first_shape = current_shape;
        } else {
            current_shape->next = calloc(sizeof(CowShape), 1);
            current_shape = (CowShape *) current_shape->next;
        }

        for (NSVGpath *path = shape->paths; path != NULL; path = path->next) {
            CowPath current_path = {0};
            lastPos = (CowPoint) {path->pts[0], path->pts[1]};
            priv_cow_add_line_to_path(&current_path, lastPos.x, lastPos.y);
            for (int i = 0; i < path->npts - 1; i += 3) {
                float *p = &path->pts[i * 2];
                priv_cow_add_bezier_to_path(&current_path, lastPos, (CowPoint) {p[2], p[3]},
                                            (CowPoint) {p[4], p[5]}, (CowPoint) {p[6], p[7]});
                lastPos = (CowPoint) {p[6], p[7]};
            }
            priv_cow_add_path_to_shape(current_shape, &current_path);
        }
    }

    return first_shape;
}

CowSvg cowvg_load_svg(char *svg_str) {
    CowSvg svg;
    svg.svg = nsvgParse(svg_str, "px", 96);
    return svg;
}

void cowvg_free_svg(CowSvg svg) {
    nsvgDelete(svg.svg);
}

static CowPoint priv_cow_get_smallest_point(CowShape *shape) {
    CowPoint smallest = {FLT_MAX, FLT_MAX};
    CowShape *current = shape;
    while (current != NULL) {

        for (int i = 0; i < current->path_counts; ++i) {
            CowPath *path = &current->paths[i];
            for (int j = 0; j < path->point_counts - 1; j += 2) {
                if (smallest.x > path->points[j]) {
                    smallest.x = path->points[j];
                }

                if (smallest.y > path->points[j + 1]) {
                    smallest.y = path->points[j + 1];
                }
            }
        }

        current = (CowShape *) current->next;
    }

    return smallest;
}

static CowPoint priv_cow_get_largest_point(CowShape *shape) {
    CowPoint largest = {FLT_MIN, FLT_MIN};
    CowShape *current = shape;
    while (current != NULL) {

        for (int i = 0; i < current->path_counts; ++i) {
            CowPath *path = &current->paths[i];
            for (int j = 0; j < path->point_counts - 1; j += 2) {
                if (largest.x < path->points[j]) {
                    largest.x = path->points[j];
                }

                if (largest.y < path->points[j + 1]) {
                    largest.y = path->points[j + 1];
                }
            }
        }

        current = (CowShape *) current->next;
    }

    return largest;
}

static void priv_cow_apply_transform(CowShape *shape, float translate_x, float translate_y, float scale) {
    CowShape *current = shape;
    while (current != NULL) {
        for (int i = 0; i < current->path_counts; ++i) {
            CowPath *path = &current->paths[i];
            for (int j = 0; j < path->point_counts - 1; j += 2) {
                path->points[j] += translate_x;
                path->points[j] *= scale;
                path->points[j + 1] += translate_y;
                path->points[j + 1] *= scale;
            }
        }
        current = (CowShape *) current->next;
    }
}

static int priv_cow_scanline_sub_pixel(CowShape *current_shape, CowPoint pixel) {
    int sub_hit = 0;
    //for every path
    for (int i = 0; i < current_shape->path_counts; ++i) {
        CowPath *current_path = &current_shape->paths[i];
        CowPoint *points = (CowPoint *) current_path->points;
        const int points_count = current_path->point_counts / 2;
        for (int j = 1; j < points_count; j += 1) {
            if (is_point_inside_edge(pixel, points[j - 1], points[j])) {
                ++sub_hit;
            }
        }

        if (is_point_inside_edge(pixel, points[points_count - 1], points[0])) {
            ++sub_hit;
        }
    }

    return sub_hit;
}

static float priv_cow_scanline_pixel(CowShape *current_shape, int x, int y) {
    const int subpixels = 5;
    int sub_hit = 0;
    for (int sub_x = 0; sub_x < subpixels; ++sub_x) {
        for (int sub_y = 0; sub_y < subpixels; ++sub_y) {
            CowPoint pixel;
            pixel.x = (float) x;
            pixel.x += (float) (sub_x + 1) / (float) (subpixels + 1);
            pixel.y = (float) y;
            pixel.y += (float) (sub_y + 1) / (float) (subpixels + 1);
            int hit = priv_cow_scanline_sub_pixel(current_shape, pixel);

            if (hit % 2 != 0) {
                ++sub_hit;
            }
        }
    }

    return (float) sub_hit / (float) (subpixels * subpixels);
}

static void cow_traverse(CowSurface *surface, CowShape *current_shape) {
    CowColor *pixels = (CowColor *) surface->pixels;
    const int w = surface->width;
    const int h = surface->height;

    //for every pixel
    for (int x = 0; x < w; ++x) {
        for (int y = 0; y < h; ++y) {
            float subpixel = priv_cow_scanline_pixel(current_shape, x, y);
            CowColor color = current_shape->color;
            color.a = (unsigned char) ((float) color.a * subpixel);
            pixels[y * w + x] = color;
        }
    }
}

static CowShape *priv_cow_make_outlines(CowShape *shape_to_outline, float thickness) {
    CowShape *new_shapes = NULL;
    CowShape *current = shape_to_outline;
    while (current != NULL) {
        for (int i = 0; i < current->path_counts; ++i) {
            CowPath *path = &current->paths[i];
            for (int j = 0; j < path->point_counts - 1; j += 2) {
            }
        }
        current = (CowShape *) current->next;
    }
}

CowSurface *cowvg_rasterize(CowSvg svg, int surface_width, int surface_height) {
    CowSurface *surface = malloc(sizeof(CowSurface));
    surface->width = surface_width;
    surface->height = surface_height;
    surface->pixels = calloc(1, surface_width * surface_height * 4);

    CowShape *shape = priv_cow_make_lines_from_svg(svg);

    CowPoint translation = priv_cow_get_smallest_point(shape);
    CowPoint largest = priv_cow_get_largest_point(shape);

    float scale_to_fill_x = (float) surface_width / (largest.x - translation.x);
    float scale_to_fill_y = (float) surface_height / (largest.y - translation.y);

    priv_cow_apply_transform(shape, -translation.x, -translation.y,
                             scale_to_fill_x < scale_to_fill_y ? scale_to_fill_x : scale_to_fill_y);

    CowShape *current_shape = shape;
    while (current_shape != NULL) {
        cow_traverse(surface, current_shape);
        current_shape = (CowShape *) current_shape->next;
    }

    return surface;
}