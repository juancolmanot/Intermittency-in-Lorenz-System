#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

bool reinjection_1d(
    long double xn,
    long double xn1,
    long double xfp,
    long double c
){
    bool reinj = 0;
    if (xn < xfp - c || xn > xfp + c) {
        if (xn1 >= xfp - c && xn1 <= xfp + c) {
            reinj = 1;

        }
    }
    return reinj;
}