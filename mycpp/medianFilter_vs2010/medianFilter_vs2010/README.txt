--------------------------------------------------------------------------------
Assignment #2.I sample codes, GPU Programming 2012
--------------------------------------------------------------------------------

===============================================================================
FILES
===============================================================================

  medianFilter.cpp        : the main program which you don't need to modify
  medianFilter_kernel.cu  : you should implement your median filter with CUDA
                            here
  medianFilter_gold.cpp   : a simple implementation of median filter in CPU for
                            your reference
  medianFilter_kernel.h   : header of your median filter
  medianFilter_gold.h     : header of the reference median filter
  data/	                  : testing image files should be put under this
                            directory
    test.pgm              : one of our testing image
  Makefile                : the Makefile for gcc under Linux and Mac OS X
  medianFilter_vs2010.sln : the solution file for Visual Studio 2010 
  medianFilter_vs2010.vcproj : the project file for Visual Studio 2010 

===============================================================================
WHAT YOU NEED TO DO IS
===============================================================================

  What you need to do is to write a median filter in medianFilter_kernel.cu with
  CUDA. The function of your medianFilter() should be the same as
  medianFilter_gold() in medianFilter_gold.cpp. However, you can always take the
  parameter "r" as 1 in assignment 2.I. Finally you will only (need to) submit
  medianFilter_kenel.cu to us.

===============================================================================
HOW TO USE THE SAMPLE CODES
===============================================================================
  
  * Unzip this zip file and put the directory under your SDK sample code
    directory (For example, NVIDIA GPU Computing SDK/C/src/medianFilter)
  * For Windows users, you can use Visual Studio 2010 with medianFilter_vs2010.sln
  * For Linux or Mac OS X users, you can use gcc with Makefile
  
If you have any problem, you are always welcome to contact me.

TA,
Ken-Yi Lee (feis.tw@gmail.com)
