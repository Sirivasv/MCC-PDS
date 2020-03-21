#include <stdio.h>
#include <math.h>

const int N = 1023;
// const double MY_PI = 3.14159265358979323846;

int main() {
    
    for (int i = 0; i <= N; ++i) {
        double val = 0.5 * (1.0 - cos( 2.0 * M_PI * ((double)i/(double)N) ));
        printf("%.2lf\n", val);
    }

    return 0;
}