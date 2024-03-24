#include <math.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <iomanip>
#include <mpi.h>


#define IDX(I,J) ((J)*(Nx+2) + (I))

void printMatrix(const double* array, int rows, int cols, int interior_range_idx[], const std::string& filename) {
    std::ofstream outFile(filename);
    if (!outFile) {
        std::cerr << "Error: Couldn't open the file " << filename << std::endl;
        return;
    }

    outFile << std::setprecision(5) << std::fixed;

    for(int j = interior_range_idx[3]; j > interior_range_idx[2]-1; j--)
    {
      for(int i = interior_range_idx[0]; i <interior_range_idx[1]+1; i++)
      {
        outFile << array[j * cols + i] << " ";
      }
      outFile << std::endl;
    }

    // for (int i = rows-1; i > -1; --i) {
    //     for (int j = 0; j < cols; ++j) {
    //         outFile << array[i * cols + j] << " ";
    //     }
    //     outFile << std::endl;
    // }
    outFile.close();
}

void set_local_N(int &Nx, int &Global_Nx, int &blocks_x, int &rank) {

  int    rx = Global_Nx % blocks_x; //remainder
  
  Nx = (Global_Nx-rx) / blocks_x; //minimum number of points in a block

  if (rank < (Nx % blocks_x)) { // for ranks < r, Nx_local is k + 1
    Nx++;
  }
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

  // std::cout << "x" << ", " << "y" << ", " << "rank" << std::endl;
  // std::cout << coord[0] << ", " << coord[1] << ", " << rank << std::endl;
 
  int leftId, rightId, downId, upId;
  MPI_Cart_shift(cart_comm, 0, 1, &leftId, &rightId);
  MPI_Cart_shift(cart_comm, 1, 1, &downId, &upId);

  int x = coord[0];
  int y = coord[1];

  // if (leftId < 0) {
  //   leftId = MPI_PROC_NULL;
  // }
  // if (rightId < 0) {
  //   leftId = MPI_PROC_NULL;
  // }
  // if (upId < 0) {
  //   leftId = MPI_PROC_NULL;
  // }
  // if (downId < 0) {
  //   leftId = MPI_PROC_NULL;
  // }

  //std::cout << "leftId" << ", " << "rightId" << ", " << "downId" << ", " << "upId" << std::endl;
  //std::cout << leftId << ", " << rightId << ", " << downId << ", " << upId << std::endl;


  double t;
  double dt   = 0.01;
  double T    = 1;
  double dx;
  double dy;
  int    Global_Nx   = 100;
  int    Global_Ny   = 100;

  int Nx, Ny;

  Nx = Global_Nx/blocks_x;
  Ny = Global_Ny/blocks_y;
  //set_local_N(Nx, Global_Nx, blocks_x, rank);
  //set_local_N(Ny, Global_Ny, blocks_y, rank);

  std::cout << "Rank: " << rank << " Nx: " << Nx <<" Ny: "<< Ny << std::endl;

  int    Npts = (Nx+2)*(Ny+2);
  double Lx   = 1.0;
  double Ly   = 1.0;
  double nu   = 0.1;

  //calculate dx, dy
  dx = Lx / (Global_Nx-1);
  dy = Ly / (Global_Ny-1);

  double hnudt = nu*dt/dx;


  int i, j;

  double pi  = 2.0 * asin(1.0);

  double *A = nullptr;
  double *Anew = nullptr;
  

  A   = new double[Npts];
  Anew   = new double[Npts];
  
  // = 2;
  //std::cout << A[-1] << std::endl;
  //std::cout << A[Npts] << std::endl;
  // set boundary conditions
  clock_t begin = clock();

  
  double x_corner = Global_Nx/dims[0]*coord[0]*dx; // Fix this for non equal distribution of grid points
  double y_corner = Global_Nx/dims[1]*coord[1]*dy;

  //std::cout << "x_corner" << ", " << "y_corner" << std::endl;
  //std::cout << x_corner << ", " << y_corner << std::endl;

  // {left, right, bottom, top}
  int interior_range_idx[4] = {1, Nx, 1, Ny}; //for non boundary process

  //only for boundary processes 
  if (leftId == MPI_PROC_NULL) {
    interior_range_idx[0]+=2; 
    //std::cout << "True check" << std::endl;
  }

  if (rightId == MPI_PROC_NULL) {
    interior_range_idx[1]-=2; 
  }

  if (downId == MPI_PROC_NULL) {
    interior_range_idx[2]+=2; 
  }

  if (upId == MPI_PROC_NULL) {
    interior_range_idx[3]-=2; 
  }


  if (coord[1] == 0) {
  //bottom
  for (i = 1; i < Nx; i++)
    A[IDX(i,1)]   = 0; //replace 1 with boundary interior ranges
    Anew[IDX(1,j)] = 0;
  }

  if (coord[1] == blocks_y-1) {
  //top
  for (i = 1; i < Nx; i++)
    A[IDX(i,Ny-1)] = 0;
    Anew[IDX(i,Ny-1)] = 0;
  }
  
  if (coord[0] == 0) {
  //left
  for (j = 1; j < Ny; j++)
    A[IDX(1,j)] = 0;//sin(pi * j / (Ny-1));
    Anew[IDX(1,j)] = 0;
}

  if (coord[0] == blocks_x-1) {
  //right
  for (j = 1; j < Ny; j++)
    A[IDX(Nx-1,j)] = 0;//sin(pi * j / (Ny-1))*exp(-pi);
    Anew[IDX(Nx-1,j)] = 0;
}



  //initial condition

  for( j = interior_range_idx[2]; j < interior_range_idx[3]+1; j++ )
    {
      for(i = interior_range_idx[0]; i <interior_range_idx[1]+1; i++)
      {
        A[IDX(i,j)] = 5*sin(pi*(x_corner+dx*i)/Lx)*sin(pi*(y_corner+dy*j)/Ly);
        
        Anew[IDX(i,j)] = 5*sin(pi*(x_corner+dx*i)/Lx)*sin(pi*(y_corner+dy*j)/Ly);
        
      }
    }

  printMatrix(A, Ny+2, Nx+2, interior_range_idx, "initial_condition.txt");

  

  // for (j = 0; j < Nx; j++)
  // {
  //   Anew[IDX(Nx-1,j)] = sin(pi * j / (Ny-1))*exp(-pi);
  // }
  

  for (t = 0; t < T; t+=dt)
  {
    if (rank==0) std::cout << "Current time step: " <<t << std::endl;
    MPI_Status status;

    //use differnet tag for different exchange

    // send to left + receive from right
    int tag = 1;
    //for left edge, there is no element need to send to left
    for( j = interior_range_idx[2]; j < interior_range_idx[3]+1; j++ )
    {
      //if (rank==0) std::cout << "Entered for loop " << std::endl;
      if (leftId != MPI_PROC_NULL)
      {
          //send to left
          //send the comp domain
          
          MPI_Send(&A[IDX(interior_range_idx[0], j)], 1, MPI_DOUBLE, leftId, tag, MPI_COMM_WORLD);
          
      }
      
      if (rightId != MPI_PROC_NULL)
      {
          //rcv from right
          //put into ghost
          
          MPI_Recv(&A[IDX(interior_range_idx[1]+1, j)], 1, MPI_DOUBLE, rightId, tag, MPI_COMM_WORLD, &status);
      }

      // send to right + receive from left
      tag = 2;
      //get right value
      if (rightId != MPI_PROC_NULL)//(x != blocks_x - 1)
      {
          //send to right
          //send the comp domain
          
          MPI_Send(&A[IDX(interior_range_idx[1], j)], 1, MPI_DOUBLE, rightId, tag, MPI_COMM_WORLD);
      }

      if (leftId != MPI_PROC_NULL)
      {
          //rcv from left
          //put into ghost
          
          MPI_Recv(&A[IDX(interior_range_idx[0]-1, j)], 1, MPI_DOUBLE, leftId, tag, MPI_COMM_WORLD, &status);
      }
    
    }
    
    //std::cout << "rank: " << rank << "middle " << std::endl;
   
    //send up + receive from below
    tag = 3;
    if (upId != MPI_PROC_NULL)
    {
        //send to up
        //send comp domain as vector
        
        MPI_Send(&A[IDX(interior_range_idx[0], interior_range_idx[3])], interior_range_idx[1]-interior_range_idx[0]+1, MPI_DOUBLE, upId, tag, MPI_COMM_WORLD); //Fix
    }

    //rcv down value
    if (downId != MPI_PROC_NULL)
    {
        MPI_Recv(&A[IDX(interior_range_idx[0], interior_range_idx[2]-1)], interior_range_idx[1]-interior_range_idx[0]+1, MPI_DOUBLE, downId, tag, MPI_COMM_WORLD, &status);
    }

    // send down + receive from above
    tag = 4;
    if (downId != MPI_PROC_NULL)
    {
        // send to down
      
        MPI_Send(&A[IDX(interior_range_idx[0], interior_range_idx[2])], interior_range_idx[1]-interior_range_idx[0]+1, MPI_DOUBLE, downId, tag, MPI_COMM_WORLD);
    }
    if (upId != MPI_PROC_NULL)
    {
        //rcv from up

        std::cout << "rank: " <<rank << " " << interior_range_idx[0] << ", "<< interior_range_idx[3] << std::endl;
        MPI_Recv(&A[IDX(interior_range_idx[0], interior_range_idx[3]+1)], interior_range_idx[1]-interior_range_idx[0]+1, MPI_DOUBLE, upId, tag, MPI_COMM_WORLD, &status);
        
    }

     //MPI_Barrier(MPI_COMM_WORLD);

     for( j = interior_range_idx[2]; j < interior_range_idx[3]+1; j++ )
    {
      for(i = interior_range_idx[0]; i <interior_range_idx[1]+1; i++)
      {
        Anew[IDX(i,j)] = A[IDX(i,j)]+hnudt*(A[IDX(i+1,j)]+A[IDX(i-1,j)]+A[IDX(i,j+1)]+A[IDX(i,j-1)]-4*A[IDX(i,j)]); 
        
      }
    }
    
    if (rank==0) std::cout << "Entered fd loop " << std::endl;
    for( j = interior_range_idx[2]; j < interior_range_idx[3]+1; j++ )
    {
      for(i = interior_range_idx[0]; i <interior_range_idx[1]+1; i++)
      {
        A[IDX(i,j)] = Anew[IDX(i,j)];    
      }
    }

    //MPI_Barrier(MPI_COMM_WORLD);    
    
  }

  clock_t end = clock();

  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

  std::string filename = "output" + std::to_string(rank) + ".txt";
  printMatrix(A, Ny+2, Nx+2, interior_range_idx ,filename);

	std::cout << "-- Run-time: " << 
			time_spent << " ms --\n";


  delete[] A;
  delete[] Anew;
  MPI_Finalize();

  return 0;
  
}

