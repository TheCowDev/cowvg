#include <stdio.h>
#include <stdlib.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION

#include "stb_image_write.h"
#include "cowvg.h"

char *loadFileAsString(const char *filename) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        printf("Failed to open the file: %s\n", filename);
        return NULL;
    }

    fseek(file, 0, SEEK_END);
    long fileSize = ftell(file);
    rewind(file);

    char *buffer = (char *) malloc(fileSize + 1);
    if (buffer == NULL) {
        printf("Failed to allocate memory for file content.\n");
        fclose(file);
        return NULL;
    }

    size_t bytesRead = fread(buffer, sizeof(char), fileSize, file);
    if (bytesRead != fileSize) {
        printf("Error reading file: %s\n", filename);
        fclose(file);
        free(buffer);
        return NULL;
    }

    buffer[fileSize] = '\0';

    fclose(file);
    return buffer;
}


int main() {
    CowSvg svg = cowvg_load_svg(loadFileAsString("test.svg"));
    CowSurface *surface = cowvg_rasterize(svg, 400, 400);
    stbi_write_png("test.png", surface->width, surface->height, 4, surface->pixels, surface->width * 4);
    return 0;
}