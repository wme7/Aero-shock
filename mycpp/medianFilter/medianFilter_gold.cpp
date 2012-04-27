#include <vector>
#include <algorithm>

#include "medianFilter_gold.h"

extern "C" void medianFilter_gold(
    const unsigned char *src, unsigned char *dst,
    unsigned int w, unsigned int h, unsigned int r)
{
    const int sr = (signed int)r;
    for (int px = sr; px < (signed int)w-sr; ++px) {
        for (int py = sr; py < (signed int)h-sr; ++py) {
            std::vector<unsigned char> list;
            for (int dx = -sr; dx <= sr; ++dx) {
                for (int dy = -sr; dy <= sr; ++dy) {
                    const int qx = px+dx;
                    const int qy = py+dy;
                    list.push_back(src[qy*w+qx]);
                }
            }
            sort(list.begin(), list.end());
            dst[py*w+px] = list[list.size()/2];
        }
    }
}
