/*
        Copyright (C) 2014-2020 Igor van den Hoven ivdhoven@gmail.com
*/

/*
        This program is free software: you can redistribute it and/or modify
        it under the terms of the GNU General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.

        This program is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        GNU General Public License for more details.

        You should have received a copy of the GNU General Public License
        along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/*
        quadsort 1.1.1.1
*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <algorithm>

typedef int CMPFUNC(const void* a, const void* b);

template <CMPFUNC cmp>
void quad_swap32(void* array, void* swap, size_t nmemb)
{
    size_t offset;

    int* pta = (int*)array;
    int* pts = (int*)swap;

    for (offset = 0; offset + 4 <= nmemb; offset += 4) {
        if (cmp(&pta[0], &pta[1]) > 0) {
            pts[0] = pta[1];
            pts[1] = pta[0];
        } else {
            pts[0] = pta[0];
            pts[1] = pta[1];
        }

        if (cmp(&pta[2], &pta[3]) > 0) {
            pts[2] = pta[3];
            pts[3] = pta[2];
        } else {
            pts[2] = pta[2];
            pts[3] = pta[3];
        }

        if (cmp(&pts[1], &pts[2]) <= 0) {
            *pta++ = pts[0];
            *pta++ = pts[1];
            *pta++ = pts[2];
            *pta++ = pts[3];
        } else if (cmp(&pts[0], &pts[3]) > 0) {
            *pta++ = pts[2];
            *pta++ = pts[3];
            *pta++ = pts[0];
            *pta++ = pts[1];
        } else if (cmp(&pts[0], &pts[2]) > 0) {
            *pta++ = pts[2];
            *pta++ = pts[0];

            if (cmp(&pts[1], &pts[3]) > 0) {
                *pta++ = pts[3];
                *pta++ = pts[1];
            } else {
                *pta++ = pts[1];
                *pta++ = pts[3];
            }
        } else {
            *pta++ = pts[0];
            *pta++ = pts[2];

            if (cmp(&pts[1], &pts[3]) > 0) {
                *pta++ = pts[3];
                *pta++ = pts[1];
            } else {
                *pta++ = pts[1];
                *pta++ = pts[3];
            }
        }
    }

    switch (nmemb - offset) {
        case 0:
        case 1: return;
        case 2:
            if (cmp(&pta[0], &pta[1]) > 0) {
                pts[0] = pta[0];
                pta[0] = pta[1];
                pta[1] = pts[0];
            }
            return;
        case 3:
            if (cmp(&pta[0], &pta[1]) > 0) {
                pts[0] = pta[0];
                pta[0] = pta[1];
                pta[1] = pts[0];
            }
            if (cmp(&pta[1], &pta[2]) > 0) {
                pts[0] = pta[1];
                pta[1] = pta[2];
                pta[2] = pts[0];
            }
            if (cmp(&pta[0], &pta[1]) > 0) {
                pts[0] = pta[0];
                pta[0] = pta[1];
                pta[1] = pts[0];
            }
            return;
        default: assert(nmemb < 1 && nmemb > 3);
    }
}

template <CMPFUNC cmp>
void quad_swap64(void* array, void* swap, size_t nmemb)
{
    size_t     offset;
    long long *pts, *pta;

    pta = (long long*)array;
    pts = (long long*)swap;

    for (offset = 0; offset + 4 <= nmemb; offset += 4) {
        if (cmp(&pta[0], &pta[1]) > 0) {
            pts[0] = pta[1];
            pts[1] = pta[0];
        } else {
            pts[0] = pta[0];
            pts[1] = pta[1];
        }

        if (cmp(&pta[2], &pta[3]) > 0) {
            pts[2] = pta[3];
            pts[3] = pta[2];
        } else {
            pts[2] = pta[2];
            pts[3] = pta[3];
        }

        if (cmp(&pts[1], &pts[2]) <= 0) {
            *pta++ = pts[0];
            *pta++ = pts[1];
            *pta++ = pts[2];
            *pta++ = pts[3];
        } else if (cmp(&pts[0], &pts[3]) > 0) {
            *pta++ = pts[2];
            *pta++ = pts[3];
            *pta++ = pts[0];
            *pta++ = pts[1];
        } else if (cmp(&pts[0], &pts[2]) > 0) {
            *pta++ = pts[2];
            *pta++ = pts[0];

            if (cmp(&pts[1], &pts[3]) > 0) {
                *pta++ = pts[3];
                *pta++ = pts[1];
            } else {
                *pta++ = pts[1];
                *pta++ = pts[3];
            }
        } else {
            *pta++ = pts[0];
            *pta++ = pts[2];

            if (cmp(&pts[1], &pts[3]) > 0) {
                *pta++ = pts[3];
                *pta++ = pts[1];
            } else {
                *pta++ = pts[1];
                *pta++ = pts[3];
            }
        }
    }

    switch (nmemb - offset) {
        case 0:
        case 1: return;
        case 2:
            if (cmp(&pta[0], &pta[1]) > 0) {
                pts[0] = pta[0];
                pta[0] = pta[1];
                pta[1] = pts[0];
            }
            return;
        case 3:
            if (cmp(&pta[0], &pta[1]) > 0) {
                pts[0] = pta[0];
                pta[0] = pta[1];
                pta[1] = pts[0];
            }
            if (cmp(&pta[1], &pta[2]) > 0) {
                pts[0] = pta[1];
                pta[1] = pta[2];
                pta[2] = pts[0];
            }
            if (cmp(&pta[0], &pta[1]) > 0) {
                pts[0] = pta[0];
                pta[0] = pta[1];
                pta[1] = pts[0];
            }
            return;
        default: assert(nmemb < 1 && nmemb > 3);
    }
}

template <CMPFUNC cmp>
void quad_sort32(void* array, void* swap, size_t nmemb)
{
    size_t offset, block = 4;
    int *  pta, *pts, *c, *c_max, *d, *d_max, *end;

    end = (int*)array;
    end += nmemb;

    while (block < nmemb) {
        offset = 0;

        while (offset + block < nmemb) {
            pta = (decltype(pta))array;
            pta += offset;

            d_max = pta + block;

            if (cmp(d_max - 1, d_max) <= 0) {
                if (offset + block * 3 < nmemb) {
                    d_max = pta + block * 3;

                    if (cmp(d_max - 1, d_max) <= 0) {
                        d_max = pta + block * 2;

                        if (cmp(d_max - 1, d_max) <= 0) {
                            offset += block * 4;
                            continue;
                        }
                        pts = (decltype(pts))swap;

                        c     = pta;
                        c_max = pta + block * 2;
                        d     = c_max;
                        d_max = offset + block * 4 <= nmemb ? d + block * 2 : end;

                        while (c < c_max) *pts++ = *c++;

                        while (d < d_max) *pts++ = *d++;

                        goto step3;
                    }
                    pts = (int*)swap;

                    c     = pta;
                    c_max = pta + block * 2;

                    while (c < c_max) *pts++ = *c++;

                    goto step2;
                } else if (offset + block * 2 < nmemb) {
                    d_max = pta + block * 2;

                    if (cmp(d_max - 1, d_max) <= 0) {
                        offset += block * 4;
                        continue;
                    }
                    pts = (int*)swap;

                    c     = pta;
                    c_max = pta + block * 2;

                    while (c < c_max) *pts++ = *c++;

                    goto step2;
                } else {
                    offset += block * 4;
                    continue;
                }
            }

            // step1:

            pts = (int*)swap;

            c     = pta;
            c_max = pta + block;

            d     = c_max;
            d_max = offset + block * 2 <= nmemb ? d + block : end;

            if (cmp(c_max - 1, d_max - 1) <= 0) {
                while (c < c_max) {
                    while (cmp(c, d) > 0) { *pts++ = *d++; }
                    *pts++ = *c++;
                }
                while (d < d_max) *pts++ = *d++;
            } else if (cmp(c, d_max - 1) > 0) {
                while (d < d_max) *pts++ = *d++;

                while (c < c_max) *pts++ = *c++;
            } else {
                while (d < d_max) {
                    while (cmp(c, d) <= 0) { *pts++ = *c++; }
                    *pts++ = *d++;
                }

                while (c < c_max) { *pts++ = *c++; }
            }

        step2:

            if (offset + block * 2 < nmemb) {
                c = pta + block * 2;

                if (offset + block * 3 < nmemb) {
                    c_max = c + block;
                    d     = c_max;
                    d_max = offset + block * 4 <= nmemb ? d + block : end;

                    if (cmp(c_max - 1, d_max - 1) <= 0) {
                        while (c < c_max) {
                            while (cmp(c, d) > 0) { *pts++ = *d++; }
                            *pts++ = *c++;
                        }
                        while (d < d_max) *pts++ = *d++;
                    } else if (cmp(c, d_max - 1) > 0) {
                        while (d < d_max) *pts++ = *d++;
                        while (c < c_max) *pts++ = *c++;
                    } else {
                        while (d < d_max) {
                            while (cmp(c, d) <= 0) { *pts++ = *c++; }
                            *pts++ = *d++;
                        }
                        while (c < c_max) *pts++ = *c++;
                    }
                } else {
                    while (c < end) *pts++ = *c++;
                }
            }

        step3:

            pts = (int*)swap;

            c = pts;

            if (offset + block * 2 < nmemb) {
                c_max = c + block * 2;

                d     = c_max;
                d_max = offset + block * 4 <= nmemb ? d + block * 2 : pts + nmemb - offset;

                if (cmp(c_max - 1, d_max - 1) <= 0) {
                    while (c < c_max) {
                        while (cmp(c, d) > 0) { *pta++ = *d++; }
                        *pta++ = *c++;
                    }
                    while (d < d_max) *pta++ = *d++;
                } else if (cmp(c, d_max - 1) > 0) {
                    while (d < d_max) *pta++ = *d++;
                    while (c < c_max) *pta++ = *c++;
                } else {
                    while (d < d_max) {
                        while (cmp(d, c) > 0) { *pta++ = *c++; }
                        *pta++ = *d++;
                    }
                    while (c < c_max) *pta++ = *c++;
                }
            } else {
                d_max = pts + nmemb - offset;

                while (c < d_max) *pta++ = *c++;
            }
            offset += block * 4;
        }
        block *= 4;
    }
}

template <CMPFUNC cmp>
void quad_sort64(void* array, void* swap, size_t nmemb)
{
    size_t     offset, block = 4;
    long long *pta, *pts, *c, *c_max, *d, *d_max, *end;

    end = (decltype(end))array;
    end += nmemb;

    while (block < nmemb) {
        offset = 0;

        while (offset + block < nmemb) {
            pta = (decltype(pta))array;
            pta += offset;

            d_max = pta + block;

            if (cmp(d_max - 1, d_max) <= 0) {
                if (offset + block * 3 < nmemb) {
                    d_max = pta + block * 3;

                    if (cmp(d_max - 1, d_max) <= 0) {
                        d_max = pta + block * 2;

                        if (cmp(d_max - 1, d_max) <= 0) {
                            offset += block * 4;
                            continue;
                        }
                        pts = (decltype(pts))swap;

                        c     = pta;
                        c_max = pta + block * 2;
                        d     = c_max;
                        d_max = offset + block * 4 <= nmemb ? d + block * 2 : end;

                        while (c < c_max) *pts++ = *c++;

                        while (d < d_max) *pts++ = *d++;

                        goto step3;
                    }
                    pts = (decltype(pts))swap;

                    c     = pta;
                    c_max = pta + block * 2;

                    while (c < c_max) *pts++ = *c++;

                    goto step2;
                } else if (offset + block * 2 < nmemb) {
                    d_max = pta + block * 2;

                    if (cmp(d_max - 1, d_max) <= 0) {
                        offset += block * 4;
                        continue;
                    }
                    pts = (decltype(pts))swap;

                    c     = pta;
                    c_max = pta + block * 2;

                    while (c < c_max) *pts++ = *c++;

                    goto step2;
                } else {
                    offset += block * 4;
                    continue;
                }
            }

            // step1:

            pts = (decltype(pts))swap;

            c     = pta;
            c_max = pta + block;

            d     = c_max;
            d_max = offset + block * 2 <= nmemb ? d + block : end;

            if (cmp(c_max - 1, d_max - 1) <= 0) {
                while (c < c_max) {
                    while (cmp(c, d) > 0) { *pts++ = *d++; }
                    *pts++ = *c++;
                }
                while (d < d_max) *pts++ = *d++;
            } else if (cmp(c, d_max - 1) > 0) {
                while (d < d_max) *pts++ = *d++;

                while (c < c_max) *pts++ = *c++;
            } else {
                while (d < d_max) {
                    while (cmp(c, d) <= 0) { *pts++ = *c++; }
                    *pts++ = *d++;
                }

                while (c < c_max) { *pts++ = *c++; }
            }

        step2:

            if (offset + block * 2 < nmemb) {
                c = pta + block * 2;

                if (offset + block * 3 < nmemb) {
                    c_max = c + block;
                    d     = c_max;
                    d_max = offset + block * 4 <= nmemb ? d + block : end;

                    if (cmp(c_max - 1, d_max - 1) <= 0) {
                        while (c < c_max) {
                            while (cmp(c, d) > 0) { *pts++ = *d++; }
                            *pts++ = *c++;
                        }
                        while (d < d_max) *pts++ = *d++;
                    } else if (cmp(c, d_max - 1) > 0) {
                        while (d < d_max) *pts++ = *d++;
                        while (c < c_max) *pts++ = *c++;
                    } else {
                        while (d < d_max) {
                            while (cmp(c, d) <= 0) { *pts++ = *c++; }
                            *pts++ = *d++;
                        }
                        while (c < c_max) *pts++ = *c++;
                    }
                } else {
                    while (c < end) *pts++ = *c++;
                }
            }

        step3:

            pts = (long long*)swap;

            c = pts;

            if (offset + block * 2 < nmemb) {
                c_max = c + block * 2;

                d     = c_max;
                d_max = offset + block * 4 <= nmemb ? d + block * 2 : pts + nmemb - offset;

                if (cmp(c_max - 1, d_max - 1) <= 0) {
                    while (c < c_max) {
                        while (cmp(c, d) > 0) { *pta++ = *d++; }
                        *pta++ = *c++;
                    }
                    while (d < d_max) *pta++ = *d++;
                } else if (cmp(c, d_max - 1) > 0) {
                    while (d < d_max) *pta++ = *d++;
                    while (c < c_max) *pta++ = *c++;
                } else {
                    while (d < d_max) {
                        while (cmp(d, c) > 0) { *pta++ = *c++; }
                        *pta++ = *d++;
                    }
                    while (c < c_max) *pta++ = *c++;
                }
            } else {
                d_max = pts + nmemb - offset;

                while (c < d_max) *pta++ = *c++;
            }
            offset += block * 4;
        }
        block *= 4;
    }
}

template <CMPFUNC cmp>
void quadsort(void* array, size_t nmemb, size_t size)
{
    void* swap;

    swap = malloc(nmemb * size);

    if (size == sizeof(int)) {
        quad_swap32<cmp>(array, swap, nmemb);

        quad_sort32<cmp>(array, swap, nmemb);
    } else if (size == sizeof(long long)) {
        quad_swap64<cmp>(array, swap, nmemb);

        quad_sort64<cmp>(array, swap, nmemb);
    } else {
        assert(size == 4 || size == 8);
    }

    free(swap);
}

static inline int cmp_int(const void* a, const void* b) { return *(int*)a - *(int*)b; }

static inline int cmp_str(const void* a, const void* b)
{
    return strcmp(*(const char**)a, *(const char**)b);
}

static inline int cmp_float(const void* a, const void* b) { return *(float*)a - *(float*)b; }



