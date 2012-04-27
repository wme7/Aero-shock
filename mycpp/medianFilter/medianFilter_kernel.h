#ifndef _MEDIAN_FILTER_KERNEL_H_
#define _MEDIAN_FILTER_KERNEL_H_

extern "C" void medianFilter(
		const unsigned char *src, unsigned char *dst,
		unsigned int w, unsigned int h,
		unsigned int r);

#endif
