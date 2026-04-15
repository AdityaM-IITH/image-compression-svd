#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "image.h"
#include "svd.h"

int gray(double *R, double *G, double *B, size_t pix, int m, int n) {
    int err = 0;
    for (int i = 0; i < m; i++) {
        if (R[i] != G[i] || R[i] != B[i]) {
            err++;
        }
    }
    if (err < 0.01 * m * n) return 1;
    return 0;
}

int main() {
    int m, n, mxg = 255;
    double *Rin, *Gin, *Bin;
    char ifname[100], ofname[100];
    int k;
    printf("Enter the value of K: ");
    scanf("%d", &k);
    printf("Enter Image File name with format: ");
    scanf("%s", ifname);
    printf("Enter where to output the image: ");
    scanf("%s", ofname);

    printf("Opening Image: %s\n", ifname);
    int rstat = img_read_rgb(ifname, &m, &n, &Rin, &Gin, &Bin);
    if (rstat != 0) {
        printf("Error! Couldn't open image file.\n");
        return 1;
    }

    size_t pix = (size_t)m * n;
    if (gray(Rin, Gin, Bin, pix, m, n)) {
        printf("Image is Grayscale. Compressing once...\n");
        double *AhR = channel(Rin, m, n, k, mxg);
        double *AhG = calloc(pix, sizeof(double));
        double *AhB = calloc(pix, sizeof(double));
        cpyv(AhR, AhG, pix);
        cpyv(AhR, AhB, pix);
        printf("Grayscale Compression Complete. Output file %s produced\n", ofname);
        ppm_write(ofname, m, n, mxg, AhR, AhG, AhB);
        free(AhR); free(AhG); free(AhB);
    } else {
        printf("Image is RGB coloured. Compressing channel-wise...\n");
        double *AhR = channel(Rin, m, n, k, mxg);
        double *AhG = channel(Gin, m, n, k, mxg);
        double *AhB = channel(Bin, m, n, k, mxg);
        printf("RGB Compression Complete. Output file %s produced\n", ofname);
        ppm_write(ofname, m, n, mxg, AhR, AhG, AhB);
        free(AhR); free(AhG); free(AhB);
    }

    free(Rin); free(Gin); free(Bin);
    return 0;
}
