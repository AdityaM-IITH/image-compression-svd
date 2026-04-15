#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "image.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

int img_read_rgb(const char *path, int *h, int *w, double **R, double **G, double **B) {
    int width, height, orig_comp;
    
    unsigned char *data = stbi_load(path, &width, &height, &orig_comp, 3);
    
    if(!data){
        fprintf(stderr, "STB Error: %s\n", stbi_failure_reason());
        return -1;
    }

    *w = width;
    *h = height;
    size_t n = (size_t)width * height;

    double *r_data = malloc(n * sizeof(double));
    double *g_data = malloc(n * sizeof(double));
    double *b_data = malloc(n * sizeof(double));

    if(!r_data || !g_data || !b_data){
        if(r_data) free(r_data);
        if(g_data) free(g_data);
        if(b_data) free(b_data);
        stbi_image_free(data);
        return -3;
    }

    for(size_t i = 0; i < n; i++){
        r_data[i] = (double)data[i * 3 + 0];
        g_data[i] = (double)data[i * 3 + 1];
        b_data[i] = (double)data[i * 3 + 2];
    }

    stbi_image_free(data);
    
    *R = r_data;
    *G = g_data;
    *B = b_data;
    return 0;
}

int jpg_read_rgb(const char *path, int *h, int *w, double **R, double **G, double **B){
    return img_read_rgb(path, h, w, R, G, B);
}

int png_read_rgb(const char *path, int *h, int *w, double **R, double **G, double **B){
    return img_read_rgb(path, h, w, R, G, B);
}

static unsigned char* pack_rgb_data(int h, int w, double* R, double* G, double* B) {
    size_t n = (size_t)h * w;
    unsigned char* data = malloc(n * 3 * sizeof(unsigned char)); 
    if (!data) return NULL;

    for (size_t i = 0; i < n; i++) {
        double r = R[i], g = G[i], b = B[i];
        
        if (r < 0) r = 0; if (r > 255) r = 255;
        if (g < 0) g = 0; if (g > 255) g = 255;
        if (b < 0) b = 0; if (b > 255) b = 255;
        
        data[i * 3 + 0] = (unsigned char)(r + 0.5);
        data[i * 3 + 1] = (unsigned char)(g + 0.5);
        data[i * 3 + 2] = (unsigned char)(b + 0.5);
    }
    return data;
}

int jpg_write_rgb(const char *path, int h, int w, int quality, double *R, double *G, double *B){
    unsigned char* data = pack_rgb_data(h, w, R, G, B); 
    if (!data) return -1;
    
    int result = stbi_write_jpg(path, w, h, 3, data, quality); 
    
    free(data);
    return (result == 0) ? -2 : 0;
}

int png_write_rgb(const char *path, int h, int w, double *R, double *G, double *B){
    unsigned char* data = pack_rgb_data(h, w, R, G, B); 
    if (!data) return -1;
    
    int stride_in_bytes = w * 3;
    int result = stbi_write_png(path, w, h, 3, data, stride_in_bytes);
    
    free(data);
    return (result == 0) ? -2 : 0;
}
int ppm_write(const char *path, int h, int w, int mx, double *R, double *G, double *B){
    FILE *f = fopen(path, "wb");
    if(!f) return -1;
    if(mx<=0) mx = 255;    
    fprintf(f, "P6\n%d %d\n%d\n", w, h, mx);
    size_t n = (size_t)h * w;
    double r, g, b;
    for(size_t i=0; i<n; i++){
        r = R[i]; g = G[i]; b = B[i];
        if(r < 0) r = 0; if(r > mx) r = mx;
        if(g < 0) g = 0; if(g > mx) g = mx;
        if(b < 0) b = 0; if(b > mx) b = mx;
        fputc((unsigned char)(r + 0.5), f);
        fputc((unsigned char)(g + 0.5), f);
        fputc((unsigned char)(b + 0.5), f);
    }
    fclose(f);
    return 0;
}