#include "headers.h"

using namespace std;

int main() {
    double start = omp_get_wtime();
    Results results = calculation();
    double finish = omp_get_wtime();
    printf("time = %mail.f (sec)", finish-start);
    print_to_file("result/temp_.txt", results.teta, N);
    print_to_file("result/o2_.txt", results.y[1], N);
    print_to_file("result/h2o_.txt", results.y[2], N);
    print_to_file("result/co2_.txt", results.y[3], N);
    print_to_file("result/co_.txt", results.y[4], N);
    return 0;
}
