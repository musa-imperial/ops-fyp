#include <math.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <time.h>

#include <mpi.h>

#define IDX(I,J) ((J)*(Nx) + (I))

void printMatrix(const double* array, int rows, int cols, const std::string& filename, bool p2terminal) {
    std::ofstream outFile(filename);
    if (!outFile) {
        std::cerr << "Error: Couldn't open the file " << filename << std::endl;
        return;
    }

    for (int i = rows-1; i > -1; --i) {
        for (int j = 0; j < cols; ++j) {
            outFile << array[i * cols + j] << " ";
            if (p2terminal) std::cout << array[i * cols + j] << " ";
        }
        outFile << std::endl;
        if (p2terminal) std::cout << std::endl;
    }
    outFile.close();
}

void populateArray(double* arr, int rows, int cols, int rank) {
    // Seed the random number generator with current time
    srand(time(NULL)+rank);
    
    // Populate the array

    for (int i = rows-1; i > -1; --i) {
        for (int j = 0; j < cols; ++j) {
            arr[i * cols + j] = 100*rank+10*j+i;//int(rand() % 100);
        }
    }


    // for (int i = 0; i < rows * cols; ++i) {
    //     arr[i] = rand() % 100; // Generate random numbers between 0 and 99
    // }
}

int main(int argc, char **argv)
{

  MPI_Init(&argc, &argv);

  int rank, procs;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &procs);

  //assume
  int blocks_x = sqrt(procs);
  int blocks_y = sqrt(procs);

  if (blocks_x * blocks_y != procs) {
    throw std::runtime_error(
        "the total number of process should equal to blocks_x*blocks_y");
  }

  int dims[2] = {blocks_x, blocks_y};

  int periodical[2] = {0, 0};
  int reorder = 0;

  int coord[2] = {};
  MPI_Comm cart_comm;
  MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periodical, reorder, &cart_comm);
  MPI_Cart_coords(cart_comm, rank, 2, coord);

  std::cout << "x" << ", " << "y" << ", " << "rank" << std::endl;
  std::cout << coord[0] << ", " << coord[1] << ", " << rank << std::endl;
 
  int leftId, rightId, downId, upId;
  MPI_Cart_shift(cart_comm, 0, 1, &leftId, &rightId);
  MPI_Cart_shift(cart_comm, 1, 1, &downId, &upId);

  int x = coord[0];
  int y = coord[1];

  int Nx, Ny;

  Nx = 10;
  Ny = 10;

  int    Npts = (Nx)*(Ny);

  double *A = nullptr;
  double *Anew = nullptr;

  A   = new double[Npts];
  Anew   = new double[Npts];
  
  populateArray(A, Nx, Ny, rank);

  std::cout << A[0] << std::endl;

  std::string filename = "output" + std::to_string(rank) + ".txt";

  MPI_Status status;

  int tag = 1;

  for (int k = 0; k < 10; k++) {
    for(int  j = 1; j < Ny-1; j++ )
        {

        if (leftId != MPI_PROC_NULL)
        {
            //send to proccess above the array pointing to x = 1, y=j, w/ ghosts
            //send comp domain as single double
            //send to left
            
            MPI_Send(&A[IDX(1, j)], 1, MPI_DOUBLE, leftId, tag, MPI_COMM_WORLD);
            
        }
        if (rightId != MPI_PROC_NULL)
        {
            //store in the array pointing to x = Nx-1, y = j
            //store is ghost points
            //rcv from right
            //put into ghost
            
            MPI_Recv(&A[IDX(Nx-1, j)], 1, MPI_DOUBLE, rightId, tag, MPI_COMM_WORLD, &status);
        }

        // send to right + receive from left
        tag = 2;
        //get right value
        if (rightId != MPI_PROC_NULL)//(x != blocks_x - 1)
        {
            //send to proccess above the array pointing to x = Nx-2, y=j, w/ ghosts
            //send comp domain as single double
            //send to right
            //send the comp domain
            
            MPI_Send(&A[IDX(Nx-2, j)], 1, MPI_DOUBLE, rightId, tag, MPI_COMM_WORLD);
        }

        if (leftId != MPI_PROC_NULL)
        {
            //store in the array pointing to x = 0, y = j
            //store is ghost points
            //rcv from left
            //put into ghost
            
            MPI_Recv(&A[IDX(0, j)], 1, MPI_DOUBLE, leftId, tag, MPI_COMM_WORLD, &status);
        }
        
        }


    //send up + receive from below
        tag = 3;
        if (upId != MPI_PROC_NULL)
        {
            //send to proccess above the array pointing to x = 0, y=Ny-1, w/o ghost
            //send to proccess above the array pointing to x = 1, y=Ny-2, w/ ghosts
            //send comp domain as vector
            
            MPI_Send(&A[IDX(1, Ny-2)], Nx-2, MPI_DOUBLE, upId, tag, MPI_COMM_WORLD); //Fix
        }

        //rcv from down value
        if (downId != MPI_PROC_NULL)
        {
            //store in the array pointing to x = 1, y = 0
            //store in ghost points
        
            MPI_Recv(&A[IDX(1, 0)], Nx-2, MPI_DOUBLE, downId, tag, MPI_COMM_WORLD, &status);
        }

        // send down + receive from above
        tag = 4;
        if (downId != MPI_PROC_NULL)
        { 
            //send to proccess below the array pointing to x = 0, y= 0, w/o ghost
            //send to proccess below the array pointing to x = 1, y= 1, w/ ghosts
            // send to down
            
            //std::cout << interior_range_idx[0] << ", "<< interior_range_idx[2] << std::endl;
            MPI_Send(&A[IDX(1, 1)], Nx-2, MPI_DOUBLE, downId, tag, MPI_COMM_WORLD);
        }
        if (upId != MPI_PROC_NULL)
        {
            //store in the ghost array pointing to x = 1, y = Nx-1 w ghosts
            //rcv from up
            //std::cout << "rank: " << rank << "Entered 8 " << std::endl;
            MPI_Recv(&A[IDX(1, Nx-1)], Nx-2, MPI_DOUBLE, upId, tag, MPI_COMM_WORLD, &status);
            
        }

  }

  printMatrix(A, Nx, Ny, filename,false);

  delete[] A;
  delete[] Anew;
  MPI_Finalize();

  return 0;
  
}

