#include <cutil_inline.h>

#include "medianFilter_kernel.h"
#include "medianFilter_gold.h"

void loadImage(
               unsigned char **img,
               unsigned int *width,
               unsigned int *height,
               const char* image_filename,
               const char* exec_path)
{
    const char* image_path = cutFindFilePath(image_filename, exec_path);
    if (image_path == 0) {
        fprintf(stderr, "Error finding image file '%s'\n", image_filename);
        exit(-1);
    }

    cutilCheckError(cutLoadPGMub(image_path, img, width, height));

    if (!*img) {
        fprintf(stderr, "Error opening file '%s'\n", image_path);
        exit(-1);
    }

    printf("Loaded '%s', %d x %d pixels\n", image_path, *width, *height);
}

void saveImage(
               unsigned char* img,
               const unsigned int w,
               const unsigned int h,
               const char* file_path)
{
    cutSavePGMub(file_path, img, w, h);
}

bool test(
          const unsigned char* result1,
          const unsigned char* result2,
          const unsigned int w,
          const unsigned int h,
          const unsigned int r)
{
    for (unsigned int x = r; x < w-r; ++x) {
        for (unsigned int y = r; y < h-r; ++y) {
            const unsigned int idx = y*w+x;
            if (result1[idx] != result2[idx]) 
                return false;
        }
    }
    return true;
}

int main(int, char* argv[])
{
    // In assignment 2.1, we will fix filter_radius as 1,
    // and therefore it's filter size is 3x3.
    const unsigned int filter_radius = 1;

    unsigned int width, height;
    unsigned char* h_img = NULL;
    unsigned char* your_result = NULL;
    unsigned char* gold_result = NULL;

    // This is the filename of the image file in PGM format
    // and the file should be put under the subdirectory, "data".
    const char* image_filename = "test.pgm";
    
    loadImage(&h_img, &width, &height, image_filename, argv[0]);

    if (width <= 2*filter_radius || height <= 2*filter_radius) {
        fprintf(stderr, "Filter radius is too large.\n");
        exit(-1);
    }

    your_result = (unsigned char*) malloc(width*height*sizeof(unsigned char));
    gold_result = (unsigned char*) malloc(width*height*sizeof(unsigned char));


    // Run your median filter
    {
        unsigned int timer = 0;
        cutilCheckError(cutCreateTimer(&timer));
        cutilCheckError(cutStartTimer(timer));
        
        // You should implemnt medianFilter() in medianFilter_kernel.cu
        medianFilter(h_img, your_result, width, height, filter_radius);

        cutilCheckError(cutStopTimer(timer));
        printf("[Yours] Processing time: %f (ms) \n", cutGetTimerValue(timer));
        cutilCheckError(cutDeleteTimer(timer));
    }

    // Run the reference median filter
    {
        unsigned int timer = 0;
        cutilCheckError(cutCreateTimer(&timer));
        cutilCheckError(cutStartTimer(timer));
        
        medianFilter_gold(h_img, gold_result, width, height, filter_radius);

        cutilCheckError(cutStopTimer(timer));
        printf("[Gold]  Processing time: %f (ms) \n", cutGetTimerValue(timer));
        cutilCheckError(cutDeleteTimer(timer));
    }

    // You can use saveImage() to save the result.
    // Under Windows, you can use IrfanView (http://www.irfanview.com/) to view PGM files.
    // Or you can convert PGM into JPEG, PNG, and etc. 
    
    // saveImage(your_result, width, height, "output.pgm");
    // saveImage(gold_result, width, height, "gold.pgm");

    // Compare your result and the reference result.
    {
        if (test(your_result, gold_result, width, height, filter_radius)) {
            printf("PASSED\n");
        } else {
            printf("FAILED\n");
        }
    }

    free(gold_result);
    free(your_result);
    return 0;
}
