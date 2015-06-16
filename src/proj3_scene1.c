//
//  proj3_scene1.c
//  
//
//  Created by Will Kearney on 4/4/15.
//
//

#include <stdio.h>
#include <time.h>
#include "graphics.h"

#define NUM_WINDOWS 22

#define SIZE_color 5
#define SIGMA_color 1.4

#define SIZE_edge 5
#define SIGMA_edge 2.0 // use 0.4 for IMAGE_3 (stained glass)

#define IMAGE_1 "../images/rita_hayworth1.ppm"
#define IMAGE_2 "../images/rita_hayworth2.ppm"
#define IMAGE_3 "../images/glass.ppm"


int main(void) {
    
    clock_t start, end;
    double cpu_time_used;
    
    start = clock();
    
    Image *rita;
    Image *edges;
    
    rita = image_read(IMAGE_1);
    edges = image_create(rita->rows, rita->cols);
    // copy image
    for (int row=0;row<rita->rows;row++) {
        for (int col=0;col<rita->cols;col++) {
            fpixel_copy(&edges->data[row][col], &rita->data[row][col]);
        }
    }
    
    printf("here0\n");
    

    // perform smoothing
    remove_noise(rita, SIZE_color, SIGMA_color);
    remove_noise(edges, SIZE_edge, SIGMA_edge);
    
    printf("here1\n");
    
    // color segmentation
    mean_shift(rita, NUM_WINDOWS);
    
    printf("here2\n");
    
    // find gradiant strengths
    canny_edge_detector(edges);
    
    printf("here3\n");
    
    
    // combine!
    for (int row=0;row<rita->rows;row++) {
        for (int col=0;col<rita->cols;col++) {
            if (edges->data[row][col].rgb[0] == 0.0 || edges->data[row][col].rgb[0] == 0.0 || edges->data[row][col].rgb[0] == 0.0 ) {
                fpixel_copy(&rita->data[row][col], &edges->data[row][col]);
            }
        }
    }
    

    // write out image and convert
    image_write(rita, "../images/rita.ppm");
    system("convert ../images/rita.ppm ../images/rita.png");
    system("open ../images/rita.png");

    
    image_free(edges);
    image_free(rita);
    
    end = clock();
    
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    
    printf("Time: %f\n", cpu_time_used);

}





