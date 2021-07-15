#include <stdio.h>

#define PROF_USER_EVENTS_ONLY
#define PROF_EVENT_LIST \
    PROF_EVENT_CACHE(L1D, READ, MISS) \
    PROF_EVENT_CACHE(L1D, WRITE, MISS)
#include "prof.h"

int main()
{
    uint64_t faults[2] = { 0 };

    PROF_START();
    // slow code goes here...
    PROF_DO(faults[index] += counter);

    // fast or uninteresting code goes here...

    PROF_START();
    // slow code goes here...
    PROF_DO(faults[index] += counter);

    printf("Total L1 faults: R = %lu; W = %lu\n", faults[0], faults[1]);
}