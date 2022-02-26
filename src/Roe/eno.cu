#include <cstdio>
#include <vector>
using std::vector;

#include "eno.cuh"
#include "util.cuh"

dim3 grid, block;
dim3 grid2, block2;
dim3 grid3, block3;
const int block_size = 32;
timer func_timer;

int offset1 = (Ny + 7) * 4;
int offset2 = 4;
double *U_cuda, *F_cuda, *Fp_cuda, *Fd_cuda, *F_p_cuda, *F_d_cuda, *F__cuda, *q3p_cuda, *q3d_cuda;
void ENO_Solver(double ***U, double ***U1, double ***U2, double ***F, double ***Fp, double ***Fd, double ***F_p,
                double ***F_d, double ***F_, double ***G, double ***Gp, double ***Gd, double ***G_p, double ***G_d,
                double ***G_, double ****LAMDA_, double ****q3p, double ****q3d, double dx, double dy, double &dt) {
  bound(U, dx, dy);
  LF_x(U, LAMDA_, F, Fp, Fd);
  func_timer.start("ENO_x");
  ENO_x(U, F, Fp, Fd, F_p, F_d, F_, q3p, q3d, dx, dy, dt);
  func_timer.stop("ENO_x");

  func_timer.start("ENO_x_cuda");
  eno_x_cuda<<<grid, block>>>(U_cuda, F_cuda, Fp_cuda, Fd_cuda, F_p_cuda, F_d_cuda, F__cuda, q3p_cuda, q3d_cuda, dx, dy,
                              dt, offset1, offset2);

  // __global__ void eno_x_cuda_2(double* U, double* F_, double dy, double r, int offset1, int offset2) {
  grid2.x = (Nx + 1 + 31) / 32;
  grid2.y = (int(0.5 / dy) + 1 + 31) / 32;
  eno_x_cuda_2<<<grid, block>>>(U_cuda, F_cuda, dy, 1.f, offset1, offset2);

  // __global__ void eno_x_cuda_3(double* U, double* F_, double dx, double dy, double r, int offset1, int offset2) {
  grid3.x = ((2.0 / dx - 1.0 / dx) + 1 + 31) / 32;
  grid3.y = (Ny - int(0.5 / dy) + 1 + 31) / 32;
  eno_x_cuda_3<<<grid, block>>>(U_cuda, F_cuda, dx, dy, 1.f, offset1, offset2);
  cudaDeviceSynchronize();
  func_timer.stop("ENO_x_cuda");

  bound(U, dx, dy);
  LF_y(U, LAMDA_, G, Gp, Gd);
  ENO_y(U, G, Gp, Gd, G_p, G_d, G_, q3p, q3d, dx, dy, dt);
  bound(U, dx, dy);
  LF_y(U, LAMDA_, G, Gp, Gd);
  ENO_y(U, G, Gp, Gd, G_p, G_d, G_, q3p, q3d, dx, dy, dt);
  bound(U, dx, dy);
  LF_x(U, LAMDA_, F, Fp, Fd);
  func_timer.start("ENO_x");
  ENO_x(U, F, Fp, Fd, F_p, F_d, F_, q3p, q3d, dx, dy, dt);
  func_timer.stop("ENO_x");
  bound(U, dx, dy);
}

void alloc(double **ptr, const vector<int> &dims, int depth) {
  if (depth + 1 >= dims.size()) return;
  for (int i = 0; i < dims[depth]; ++i) {
    cudaMallocManaged(&(ptr[i]), dims[depth + 1] * sizeof(void *));
    alloc((double **)ptr[i], dims, depth + 1);
  }
}

double *alloc_nd(const vector<int> &dims) {
  double *output = nullptr;
  cudaMallocManaged(&output, dims[0] * sizeof(void *));
  alloc((double **)output, dims, 0);
  return output;
}

int main(int argc, char **argv) {
  double dx, dy, dt = 0, T = 0;

  vector<int> shape1 = {Nx + 7, Ny + 7, 4};
  int shape1_ = (Nx + 7) * (Ny + 7) * 4;
  vector<int> shape2 = {Nx + 7, Ny + 7, 4, 4};
  int shape2_ = (Nx + 7) * (Ny + 7) * 4 * 4;
  vector<int> shape3 = {Nx + 7, Ny + 7, 4, 3};
  int shape3_ = (Nx + 7) * (Ny + 7) * 4 * 3;
  vector<int> shape4 = {Nx + 7, Ny + 7, 1};
  int shape4_ = (Nx + 7) * (Ny + 7) * 1;

  double ***U = (double ***)alloc_nd(shape1);
  double ***U_ = (double ***)alloc_nd(shape1);
  double ***Ut = (double ***)alloc_nd(shape1);
  double ***U1 = (double ***)alloc_nd(shape1);
  double ***U2 = (double ***)alloc_nd(shape1);
  double ***G = (double ***)alloc_nd(shape1);
  double ***G_ = (double ***)alloc_nd(shape1);
  double ***F = (double ***)alloc_nd(shape1);
  double ***F_ = (double ***)alloc_nd(shape1);
  double ***G_p = (double ***)alloc_nd(shape1);
  double ***G_d = (double ***)alloc_nd(shape1);
  double ***Gp = (double ***)alloc_nd(shape1);
  double ***Gd = (double ***)alloc_nd(shape1);
  double ***F_p = (double ***)alloc_nd(shape1);
  double ***F_d = (double ***)alloc_nd(shape1);
  double ***Fp = (double ***)alloc_nd(shape1);
  double ***Fd = (double ***)alloc_nd(shape1);

  double ****L_ = (double ****)alloc_nd(shape2);
  double ****R_ = (double ****)alloc_nd(shape2);
  double ****LAMDA_ = (double ****)alloc_nd(shape2);

  double ****q3p = (double ****)alloc_nd(shape3);
  double ****q3d = (double ****)alloc_nd(shape3);

  double ***a_ = (double ***)alloc_nd(shape4);

  cudaMalloc((void **)&U_cuda, shape1_ * sizeof(double));
  cudaMalloc((void **)&F_cuda, shape1_ * sizeof(double));
  cudaMalloc((void **)&Fp_cuda, shape1_ * sizeof(double));
  cudaMalloc((void **)&Fd_cuda, shape1_ * sizeof(double));
  cudaMalloc((void **)&F_p_cuda, shape1_ * sizeof(double));
  cudaMalloc((void **)&F_d_cuda, shape1_ * sizeof(double));
  cudaMalloc((void **)&F__cuda, shape1_ * sizeof(double));
  cudaMalloc((void **)&q3p_cuda, shape3_ * sizeof(double));
  cudaMalloc((void **)&q3d_cuda, shape3_ * sizeof(double));

  const int blocksize = 32;
  grid.x = (Nx + 2 + 31) / blocksize;
  grid.y = (Ny + 2 + 31) / blocksize;
  block.x = blocksize;
  block.y = blocksize;
  block2.x = blocksize;
  block2.y = blocksize;
  block3.x = blocksize;
  block3.y = blocksize;

  // int shape1_ = (Nx + 7) * (Ny + 7) * 4;

  for (int iter = 0; iter < 10; ++iter) {
    initial(U, dx, dy);
    int n = 0;
    timeval s, e;
    gettimeofday(&s, NULL);
    while (T <= TT) {
      dt = CFL(U, dx, dy, ENOCFL);
      ENO_Solver(U, U1, U2, F, Fp, Fd, F_p, F_d, F_, G, Gp, Gd, F_p, F_d, F_, LAMDA_, q3p, q3d, dx, dy, dt);
      T += dt;
      n++;
      virtual_clear(U, dx, dy);
    }
    gettimeofday(&e, NULL);
    double ms = get_elapsed_time_ms(s, e);
    printf("ENO total time: %lf ms\n", ms);
    printf("ENO iter: %d\n", n);
  }

  func_timer.show_all();

  if (argc == 2) {
    // result check
    printf("check dat file: %s\n", argv[1]);
    check(argv[1], U);
  }
  return 0;
}