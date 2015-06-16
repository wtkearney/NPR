//
//  npr.c
//  
//
//  Created by Will Kearney on 4/4/15.
//
//
//
// Contains functions used to perform Non-Photorealistic (NPR) techniques on an image including noise removal, edge detection, and image segmentation
//
//
//
// Canny Edge detection guidance from http://dasl.mem.drexel.edu/alumni/bGreen/www.pages.drexel.edu/_weg22/can_tut.html
// and http://en.wikipedia.org/wiki/Canny_edge_detector
//
// Mean shift based on http://stackoverflow.com/questions/4831813/image-segmentation-using-mean-shift-explained
//

#include <stdio.h>
#include <math.h>
#include "graphics.h"


#define T1 0.5
#define T2 0.2

#define EPSILON 0.001




/* perform mean shift algorithm in 3D color space */
void mean_shift(Image *src, int numWindows) {
    // Seed (only once)
    srand48(arc4random());
    
    float distance, tempDistance, r, g, b;
    int convergence = 0;
    
    // localize variables
    int rows = src->rows;
    int cols = src->cols;
    
    // mode array will hold latest window number of each pixel
    int *mode;
    mode = malloc( sizeof(int)*src->imagesize );
    
    // store the feature space coordinates for each window
    float windowCoord[numWindows*3];
    
    // store the previous coord for each window (to test for convergence)
    float tempCoord[numWindows*3];
    
    // counter for each window (used when calculating means)
    float counter[numWindows];
    
    // init window coordinates in feature space
    // TODO: what should this distribution look like?
    for (int winIdx=0;winIdx<numWindows;winIdx++) {
        windowCoord[winIdx*3] = drand48();
        windowCoord[winIdx*3+1] = drand48();
        windowCoord[winIdx*3+2] = drand48();
    }
    
    while (convergence == 0) {
        for (int row=0;row<rows;row++) {
            for (int col=0;col<cols;col++) {
                // localize color info
                r = src->data[row][col].rgb[0];
                g = src->data[row][col].rgb[1];
                b = src->data[row][col].rgb[2];
                distance = HUGE_VALF; // init distance really big (huge, even)
                
                for (int winIdx=0;winIdx<numWindows;winIdx++) {
                    
                    tempDistance = sqrt((r - windowCoord[winIdx*3])*(r - windowCoord[winIdx*3])
                                    + (g - windowCoord[winIdx*3+1])*(g - windowCoord[winIdx*3+1])
                                    + (b - windowCoord[winIdx*3+2])*(b - windowCoord[winIdx*3+2]));
                    
                    if (tempDistance < distance) {
                        distance = tempDistance;
                        mode[row*cols + col] = winIdx;
                    }
                    
                } // end window loop
                
                
            }
        } // end image loop
        
        // re-init counter array and temp means
        for (int winIdx=0;winIdx<numWindows;winIdx++) {
            counter[winIdx] = 0;
            tempCoord[winIdx*3] = 0;
            tempCoord[winIdx*3+1] = 0;
            tempCoord[winIdx*3+2] = 0;
        }
        
        // calculate means for each window
        for (int row=0;row<rows;row++) {
            for (int col=0;col<cols;col++) {
                tempCoord[ mode[row*cols + col]*3 ] += src->data[row][col].rgb[0];
                tempCoord[ mode[row*cols + col]*3+1 ] += src->data[row][col].rgb[1];
                tempCoord[ mode[row*cols + col]*3+2 ] += src->data[row][col].rgb[2];
                
                counter[ mode[row*cols + col] ] += 1;
            }
        }
        
        convergence = 1;
        for (int winIdx=0;winIdx<numWindows;winIdx++) {
            if (counter[winIdx] != 0) {
                tempCoord[winIdx*3] /= counter[winIdx];
                tempCoord[winIdx*3+1] /= counter[winIdx];
                tempCoord[winIdx*3+2] /= counter[winIdx];
            }
            
            printf("Convergence: %f, %f\n", tempCoord[winIdx*3], windowCoord[winIdx*3]);
            printf("Convergence1: %f, %f\n", tempCoord[winIdx*3+1], windowCoord[winIdx*3+1]);
            printf("Convergence2: %f, %f\n", tempCoord[winIdx*3+2], windowCoord[winIdx*3+2]);
            
            // check if all the windows have converged
            // TODO: This does not technically work. Something strange happening with comparing the two coordinates...
            if (abs(tempCoord[winIdx*3] - windowCoord[winIdx*3]) > EPSILON
                || abs(tempCoord[winIdx*3+1] - windowCoord[winIdx*3+1]) > EPSILON
                || abs(tempCoord[winIdx*3+2] - windowCoord[winIdx*3+2]) > EPSILON) {
                convergence = 0;
            } else {
                // update new coordinates with mean
                windowCoord[winIdx*3] = tempCoord[winIdx*3];
                windowCoord[winIdx*3+1] = tempCoord[winIdx*3+1];
                windowCoord[winIdx*3+2] = tempCoord[winIdx*3+2];
            }
        }
        
        
    } // end while loop
    
    
    for (int row=0;row<rows;row++) {
        for (int col=0;col<cols;col++) {
            
            src->data[row][col].rgb[0] = windowCoord[ mode[row*cols + col]*3 ];
            src->data[row][col].rgb[1] = windowCoord[ mode[row*cols + col]*3+1 ];
            src->data[row][col].rgb[2] = windowCoord[ mode[row*cols + col]*3+2 ];

        }
    }
    
    free(mode);
}


/* Finds the edge strength by taking the gradiant of the image. Uses the Sobel operator. */
void canny_edge_detector(Image *src) {
    float Gx, Gy, theta;
    
    // localize variables
    int rows = src->rows;
    int cols = src->cols;
    
    // interesting, have to put this on the heap because larger images files overflow the stack...
    float *edgeStrength;
    edgeStrength = malloc( sizeof(float)*src->imagesize );
    
    int *edgeDirection;
    edgeDirection = malloc( sizeof(int)*src->imagesize );
    
    int *edgeHysteresis;
    edgeHysteresis = malloc( sizeof(int)*src->imagesize );
    
    for (int row=0;row<rows;row++) {
        for (int col=0;col<cols;col++) {
            
            if (row < 1 || row >= rows-1 || col < 1 || col >= cols-1) {
                //printf("EXCLUDE -- row: %d, col: %d\n", row, col);
                continue;
            }
            
            /* Edge Strength */
            
            Gx = 0.0;
            Gy = 0.0;
            
//              Gx              Gy
//            -1 0 +1        +1 +2 +1
//            -2 0 +2         0  0  0
//            -1 0 +1        -1 -2 -1
            
            // red
            Gx += (-src->data[row-1][col-1].rgb[0]) + (src->data[row-1][col+1].rgb[0])
            + (-2*src->data[row][col-1].rgb[0]) + (2*src->data[row][col+1].rgb[0])
            + (-src->data[row+1][col-1].rgb[0]) + (src->data[row+1][col+1].rgb[0]);
            // green
            Gx += (-src->data[row-1][col-1].rgb[1]) + (src->data[row-1][col+1].rgb[1])
            + (-2*src->data[row][col-1].rgb[1]) + (2*src->data[row][col+1].rgb[1])
            + (-src->data[row+1][col-1].rgb[1]) + (src->data[row+1][col+1].rgb[1]);
            // blue
            Gx += (-src->data[row-1][col-1].rgb[2]) + (src->data[row-1][col+1].rgb[2])
            + (-2*src->data[row][col-1].rgb[2]) + (2*src->data[row][col+1].rgb[2])
            + (-src->data[row+1][col-1].rgb[2]) + (src->data[row+1][col+1].rgb[2]);
            Gx = Gx / 3;
            if (Gx < 0) Gx = -Gx;
            
            // red
            Gy = (src->data[row-1][col-1].rgb[0]) + (2*src->data[row-1][col].rgb[0]) + (src->data[row-1][col+1].rgb[0])
            + (-src->data[row+1][col-1].rgb[0]) + (-2*src->data[row+1][col].rgb[0]) + (-src->data[row+1][col+1].rgb[0]);
            // green
            Gy = (src->data[row-1][col-1].rgb[1]) + (2*src->data[row-1][col].rgb[1]) + (src->data[row-1][col+1].rgb[1])
            + (-src->data[row+1][col-1].rgb[1]) + (-2*src->data[row+1][col].rgb[1]) + (-src->data[row+1][col+1].rgb[1]);
            // blue
            Gy = (src->data[row-1][col-1].rgb[2]) + (2*src->data[row-1][col].rgb[2]) + (src->data[row-1][col+1].rgb[2])
            + (-src->data[row+1][col-1].rgb[2]) + (-2*src->data[row+1][col].rgb[2]) + (-src->data[row+1][col+1].rgb[2]);
            Gy = Gy / 3;
            if (Gy < 0) Gy = -Gy;
            
            edgeStrength[row*cols + col] = Gx + Gy;
            
            /* Edge Direction */
            
            if (Gx == 0) {
                if (Gy == 0) {
                    theta = 0.0;
                } else {
                    theta = 90.0;
                }
                
            } else {
                theta = atan(Gy / Gx);
            }
            
            // round edge detection
            if ((0 <= theta && theta < 22.5) || (157.5 <= theta && theta < 180)) {
                edgeDirection[row*cols + col] = 0;
            } else if (22.5 <= theta < 67.5) {
                edgeDirection[row*cols + col] = 45;
            } else if (67.5 <= theta < 112.5) {
                edgeDirection[row*cols + col] = 90;
            } else { // 112.5 <= theta < 157.5
                edgeDirection[row*cols + col] = 135;
            }
            
        }
    }
    
    /* Non-maximum Suppression */
    
    for (int row=0;row<rows;row++) {
        for (int col=0;col<cols;col++) {
            
            if (row < 1 || row >= rows-1 || col < 1 || col >= cols-1) {
                continue;
            }
            
            // Compare the edge strength of the current pixel with the edge strength of the pixel in the positive and negative gradient directions. If the edge strength of the current pixel is the largest compared to the other pixels in the mask with the same direction, the value will be preserved. Otherwise, the value will be suppressed.
            if (edgeDirection[row*cols + col] == 0
                && edgeStrength[row*cols + col] >= T2
                && edgeStrength[row*cols + col] > edgeStrength[row*cols + col-1]
                && edgeStrength[row*cols + col] > edgeStrength[row*cols + col+1]) {
                if (edgeStrength[row*cols + col] >= T1) {
                    edgeHysteresis[row*cols + col] = 2; // strong edge
                } else {
                    edgeHysteresis[row*cols + col] = 1; // weak edge
                }
            } else if (edgeDirection[row*cols + col] == 45
                       && edgeStrength[row*cols + col] >= T2
                       && edgeStrength[row*cols + col] > edgeStrength[(row-1)*cols + col+1]
                       && edgeStrength[row*cols + col] > edgeStrength[(row+1)*cols + col-1]) {
                if (edgeStrength[row*cols + col] >= T1) {
                    edgeHysteresis[row*cols + col] = 2; // strong edge
                } else {
                    edgeHysteresis[row*cols + col] = 1; // weak edge
                }
            } else if (edgeDirection[row*cols + col] == 90
                       && edgeStrength[row*cols + col] >= T2
                       && edgeStrength[row*cols + col] > edgeStrength[(row-1)*cols + col]
                       && edgeStrength[row*cols + col] > edgeStrength[(row+1)*cols + col]) {
                if (edgeStrength[row*cols + col] >= T1) {
                    edgeHysteresis[row*cols + col] = 2; // strong edge
                } else {
                    edgeHysteresis[row*cols + col] = 1; // weak edge
                }
            } else if (edgeDirection[row*cols + col] == 135
                       && edgeStrength[row*cols + col] >= T2
                       && edgeStrength[row*cols + col] > edgeStrength[(row-1)*cols + col-1]
                       && edgeStrength[row*cols + col] > edgeStrength[(row+1)*cols + col+1]) {
                if (edgeStrength[row*cols + col] >= T1) {
                    edgeHysteresis[row*cols + col] = 2; // strong edge
                } else {
                    edgeHysteresis[row*cols + col] = 1; // weak edge
                }
            } else {
                edgeHysteresis[row*cols + col] = 0; // not an edge
            }
            
            
        }
    } // end non-maximum suppression
    
    
    /* Hysteresis */
    for (int row=0;row<rows;row++) {
        for (int col=0;col<cols;col++) {
            
            if (row < 1 || row >= rows-1 || col < 1 || col >= cols-1) {
                src->data[row][col].rgb[0] = 1.0;
                src->data[row][col].rgb[1] = 1.0;
                src->data[row][col].rgb[2] = 1.0;
                continue;
            }

            
            if (edgeHysteresis[row*cols + col] == 2) {
                src->data[row][col].rgb[0] = 0.0;
                src->data[row][col].rgb[1] = 0.0;
                src->data[row][col].rgb[2] = 0.0;
                src->data[row][col].a = 1.0;
                src->data[row][col].z = 1.0;
            } else if (edgeHysteresis[row*cols + col] == 1) {
                // check surrounding pixels to see if any are "strong"
                if (edgeHysteresis[(row-1)*cols + col-1] == 2       // NW
                    || edgeHysteresis[(row-1)*cols + col] == 2      // N
                    || edgeHysteresis[(row-1)*cols + col+1] == 2    // NE
                    || edgeHysteresis[row*cols + col-1] == 2        // W
                    || edgeHysteresis[row*cols + col+1] == 2        // E
                    || edgeHysteresis[(row+1)*cols + col-1] == 2    // SW
                    || edgeHysteresis[(row+1)*cols + col] == 2      // S
                    || edgeHysteresis[(row+1)*cols + col+1] == 2) { // SE
                    src->data[row][col].rgb[0] = 0.0;
                    src->data[row][col].rgb[1] = 0.0;
                    src->data[row][col].rgb[2] = 0.0;
                    src->data[row][col].a = 1.0;
                    src->data[row][col].z = 1.0;
                } else {
                    src->data[row][col].rgb[0] = 1.0;
                    src->data[row][col].rgb[1] = 1.0;
                    src->data[row][col].rgb[2] = 1.0;
                    src->data[row][col].a = 1.0;
                    src->data[row][col].z = 1.0;
                }
            } else {
                src->data[row][col].rgb[0] = 1.0;
                src->data[row][col].rgb[1] = 1.0;
                src->data[row][col].rgb[2] = 1.0;
                src->data[row][col].a = 1.0;
                src->data[row][col].z = 1.0;
            }
            
        }
    }
    
    free(edgeStrength);
    free(edgeDirection);
    free(edgeHysteresis);
    
}

/* use convolution mask to smooth specified image */
void remove_noise(Image *src, int size, float sigma) {
    
    // get kernel
    float** kernel = gaussian_kernel(size, sigma);
    
    for (int i=0;i<size;i++) {
        for (int j=0;j<size;j++) {
            printf("%f ", kernel[i][j]);
            
        }
        printf("\n");
    }

    // localize variables
    int rows = src->rows;
    int cols = src->cols;
    float accumulation_r, accumulation_g, accumulation_b = 0.0;
    
    //printf("Total rows: %d, Total cols: %d\n", rows, cols);
    
    // create temp image to write to
    Image *temp = image_create(rows, cols);
    
    for (int row=0;row<rows;row++) {
        for (int col=0;col<cols;col++) {
            accumulation_r = 0.0;
            accumulation_g = 0.0;
            accumulation_b = 0.0;
            
            if (row < size/2 || row >= rows-(size/2) || col < size/2 || col >= cols-(size/2)) {
                //printf("EXCLUDE -- row: %d, col: %d\n", row, col);
                continue;
            }
            
            //printf("row: %d, col: %d\n", row, col);
            
            for (int k_row=0;k_row<size;k_row++) {
                for (int k_col=0;k_col<size;k_col++) {
                    //if (row >= 1023) printf("kernel row: %d, kernel col: %d\n", k_row, k_col);
                    
                    // red band
                    accumulation_r += src->data[row-(size/2)+k_row][col-(size/2)+k_col].rgb[0] * kernel[k_row][k_col];
                    // green band
                    accumulation_g += src->data[row-(size/2)+k_row][col-(size/2)+k_col].rgb[1] * kernel[k_row][k_col];
                    // blue band
                    accumulation_b += src->data[row-(size/2)+k_row][col-(size/2)+k_col].rgb[2] * kernel[k_row][k_col];

                }
            }
            
            //printf("(%f, %f, %f)\n", accumulation_r, accumulation_g, accumulation_b);
            
            temp->data[row][col].rgb[0] = 1-accumulation_r/256;
            temp->data[row][col].rgb[1] = 1-accumulation_g/256;
            temp->data[row][col].rgb[2] = 1-accumulation_b/256;
            
            
        }
    } // end image loop
    
    for (int row=0;row<rows;row++) {
        for (int col=0;col<cols;col++) {
            
            fpixel_copy(&src->data[row][col], &temp->data[row][col]);
            
        }
    }
    
    image_free(temp);
    destroy_kernel(kernel);
    
}


/* Create a Gaussian filter of the specified size */
float** gaussian_kernel(int size, float sigma) {
    int ub = size/2;
    int lb = -ub;
    float r, sum = 0.0;
    
//    printf("UB: %d\n", ub);
//    printf("LB: %d\n", lb);
    
    float** kernel = malloc( size*sizeof(float*) );
    
    for (int i=lb;i<=ub;i++) {
        kernel[i+ub] = (float *) malloc( size*sizeof(float) );
        for (int j=lb;j<=ub;j++) {
            // calculate value
            r = sqrt(i*i + j*j);
            kernel[i+ub][j+ub] = (exp( -(r*r) / 2*sigma*sigma )); // (2*M_PI*sigma*sigma) don't need because of normalization
            
            // sum used for normalizing
            sum += kernel[i+ub][j+ub];
        }
    }
    
    // normalize distribution
    for (int i=0;i<size;i++) {
        for (int j=0;j<size;j++) {
            kernel[i][j] /= sum;
        }
    }
    
//    // testing
//    for (int i=0;i<size;i++) {
//        for (int j=0;j<size;j++) {
//            if (i == size/2 && j == size/2) {
//                kernel[i][j] = 1;
//            } else {
//                kernel[i][j] = 0;
//            }
//        }
//    }
    
    return kernel;
}

/* Deallocate mem used for gaussian kernel */
void destroy_kernel(float** kernel) {
    free(*kernel);
    free(kernel);
}



