#include <stdio.h>
#include <omp.h>

int main() {
    #pragma omp parallel
    {
        printf("thread %i\n", omp_get_thread_num());
    }
}
