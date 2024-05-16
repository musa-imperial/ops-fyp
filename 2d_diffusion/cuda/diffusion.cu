#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define IDX(I, J) ((J)*Nx+ (I))

cudaError_t cudaCheck(cudaError_t result)
{
#if defined(DEBUG) || defined(_DEBUG)
  if (result != cudaSuccess) {
    fprintf(stderr, "CUDA Runtime Error: %s\n", 
            cudaGetErrorString(result));
    assert(result == cudaSuccess);
  }
#endif
  return result;
}

// CUDA diffusion kernel
__global__ void interior_kernel(double *Anew, double *A, int Npts, int Nx, double hnudt) 
{
    //Get our global thread ID 1D block and grid
    int id = blockIdx.x*blockDim.x+threadIdx.x;

    //at id we are at (i, j)
    //id+1 -> (i+1, j)
    //id+blockIdx

    // Add for n number of elements
    if (id < Npts) 
    {
        if (id%Nx!=0 && (id+1)%Nx!=0 && id<Npts-Nx && id>Nx-1)
        Anew[id] = A[id]+hnudt*(A[id+1]+A[id-1]+A[id+Nx]+A[id-Nx]-4*A[id]);
    }
}

__global__ void copy(double *Anew, double *A, int Npts) {

    int id = blockIdx.x*blockDim.x+threadIdx.x;
    if (id < Npts) {
        A[id] = Anew[id];
    }
}

int main()
{
    double t;
    double dt = 0.001;
    double T = 1.0;

    int Nx = 40;
    int Ny = 40;
    int Npts = Nx*Ny;

    double dx = 1;
    double dy = 1;
    double Lx = dx*(Nx-1);
    double Ly = dy*(Ny-1);

    double nu = 0.1;

    double hnudt = nu*dt/dx/dx;

    int i, j;

    double pi = 2.0 * asin(1.0);

    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    double *h_A, *h_Anew; // CPU (host) vectors
    double *d_A, *d_Anew; // GPU(host) (device vectors)
    size_t size = Npts*sizeof(double); //size of each vector

    h_Anew = (double*)malloc(size); //Allocate CPU vectors
    h_A = (double*)malloc(size);

    cudaCheck( cudaMalloc(&d_Anew, size));
    cudaCheck( cudaMalloc(&d_A, size));

    for ( i = 0; i < Nx; i++) {
        h_A[IDX(i, 0)] = 0.0;
        h_Anew[IDX(i, 0)] = 0.0;
        h_A[IDX(i, Ny-1)] = 0.0;
        h_Anew[IDX(i, Ny-1)] = 0.0;
    }

    for ( j = 0; j < Ny; j++) {
        h_A[IDX(0, j)] = 0.0;
        h_Anew[IDX(0, j)] = 0.0;
        h_A[IDX(Nx-1, j)] = 0.0;
        h_Anew[IDX(Nx-1, j)] = 0.0;
    }

    for ( j = 1; j < Ny-1; j++) {
        for (i = 1;i < Nx-1;i++) {
            h_A[IDX(i,j)] = 5*sin(pi*(dx*i)/Lx)*sin(pi*(dy*j)/Ly);
            h_Anew[IDX(i,j)] = 5*sin(pi*(dx*i)/Lx)*sin(pi*(dy*j)/Ly);
        }
    }

    //Copy host vectors to device vector
    cudaCheck( cudaMemcpy( d_Anew, h_Anew, size, cudaMemcpyHostToDevice));
    cudaCheck( cudaMemcpy( d_A, h_A, size, cudaMemcpyHostToDevice));

    int blockSize, gridSize;
    blockSize = 1024; // One-dimensional block- number of threads
    gridSize = (int)ceil((double)Npts/blockSize); //One-d grid
                                                // number of grids

    cudaEventRecord(start);
    for (t = 0; t < T; t+=dt) {
        // Execute the kernel
        interior_kernel<<<gridSize, blockSize>>>(d_Anew, d_A, Npts, Nx, hnudt);

        copy<<<gridSize, blockSize>>>(d_Anew, d_A, Npts);

        // Copy array back output vector from device to host
        cudaCheck( cudaMemcpy( h_A, d_A, size, cudaMemcpyDeviceToHost));
    }
    cudaEventRecord(stop);

    cudaEventSynchronize(stop);
    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);
    // Sum up vector c and print result diveded by N: 1- no error
    double u;
    double error = 0.0;
    double max_error = 0.0;
    for( j = 1; j < Ny-1; j++ )
    {
      for(i = 1; i < Nx-1; i++)
      {
        u = 5*exp(-nu*pi*pi*(1/Lx/Lx+1/Ly/Ly)*T)*sin(pi/Lx*(dx*(i)))*sin(pi/Ly*(dy*(j)));

        error = error + sqrt(abs(h_A[IDX(i,j)]*h_A[IDX(i,j)]-u*u));
        max_error = fmax(max_error, fabs((h_A[IDX(i,j)]-u)/u));

        //std::cout << dx*(i) << ", "<< dy*(j) << ", " << (A[IDX(i,j)]) << ", " << u << ", " << fabs((A[IDX(i,j)]-u))/u << std::endl;
      }
    }

    printf("Max percentage error: %le\n", 100*max_error);
    printf("Runtime %fs\n ", milliseconds/10e2);
    cudaCheck( cudaFree(d_Anew)); // Release device memory
    cudaCheck( cudaFree(d_A));
    cudaCheck( cudaFree(h_Anew)); // Release host memory
    cudaCheck( cudaFree(h_A));
   

    return 0;



}