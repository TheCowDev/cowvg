#ifndef COWVG_COWVG_H
#define COWVG_COWVG_H

typedef struct {
    unsigned char *pixels;
    int width;
    int height;
} CowSurface;

typedef struct {
    struct NSVGimage *svg;
} CowSvg;


CowSvg cowvg_load_svg(char *svg_str);

void cowvg_free_svg(CowSvg svg);

CowSurface *cowvg_rasterize(CowSvg svg, int surface_width, int surface_height);

#endif //COWVG_COWVG_H