//
//  npr.h
//  
//
//  Created by Will Kearney on 4/4/15.
//
//

#ifndef _npr_h
#define _npr_h

void mean_shift(Image *src, int numWindows);
void canny_edge_detector(Image *src);
float** gaussian_kernel(int size, float sigma);
void destroy_kernel(float** kernel);
void remove_noise(Image *src, int size, float sigma);


#endif
