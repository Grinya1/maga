#include "headers_mpi.h"
#include "mpi.h"

using namespace std;

int main(int argc, char **argv) {
    int numproc, rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    double start = MPI_Wtime();
    Results results = calculation(numproc, rank);
    double finish = MPI_Wtime();

    printf("rank %i finished\ttime = %.f (sec)", rank, finish-start);
    
    MPI_Finalize();

    printf("%i", rank);
    if (rank == 0) {
        print_to_file("c:\\Projects\\result\\temp_.txt", results.teta, N);
        print_to_file("result/o2_.txt", results.y[1], N);
        print_to_file("result/h2o_.txt", results.y[2], N);
        print_to_file("result/co2_.txt", results.y[3], N);
        print_to_file("result/co_.txt", results.y[4], N);
    }
    
    return 0;
}