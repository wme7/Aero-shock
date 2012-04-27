#ifndef _MEDIAN_FILTER_GOLD_H_
#define _MEDIAN_FILTER_GOLD_H_

extern "C" void medianFilter_gold(
		const unsigned char *src, unsigned char *dst,
		unsigned int width, unsigned int height,
		unsigned int filter_radius);

#endif
