#include <cstdio>
#include <vector>
using std::vector;

#include "eno.cuh"
#include "util.cuh"
// #define CUDAV2

dim3 grid, block;
dim3 grid2, block2;
dim3 grid3, block3;
const int blocksize = 32;
#ifdef CUDAV2
const int blocksize_Y = 8;
#else
const int blocksize_Y = 32;
#endif

timer func_timer;

int offset1 = (Ny + 7) * 4;
int offset2 = 4;
double *U_cuda, *F_cuda, *Fp_cuda, *Fd_cuda, *F_p_cuda, *F_d_cuda, *F__cuda, *q3p_cuda, *q3d_cuda;
void ENO_Solver(Tensor &U, Tensor &U1, Tensor &U2, Tensor &F, Tensor &Fp, Tensor &Fd, Tensor &F_p, Tensor &F_d,
                Tensor &F_, Tensor &G, Tensor &Gp, Tensor &Gd, Tensor &G_p, Tensor &G_d, Tensor &G_, Tensor &LAMDA_,
                Tensor &q3p, Tensor &q3d, double dx, double dy, double &dt) {
  bound_Tensor(U, dx, dy);
  LF_x_Tensor(U, LAMDA_, F, Fp, Fd);
  func_timer.start("ENO_x");
  ENO_x_Tensor(U, F, Fp, Fd, F_p, F_d, F_, q3p, q3d, dx, dy, dt);
  func_timer.stop("ENO_x");

  func_timer.start("ENO_x_cuda");
#ifdef CUDAV2
  eno_x_cuda_v2<<<grid, block>>>(U_cuda, F_cuda, Fp_cuda, Fd_cuda, F_p_cuda, F_d_cuda, F__cuda, q3p_cuda, q3d_cuda, dx,
                                 dy, dt, offset1, offset2);
#else
  eno_x_cuda<<<grid, block>>>(U_cuda, F_cuda, Fp_cuda, Fd_cuda, F_p_cuda, F_d_cuda, F__cuda, q3p_cuda, q3d_cuda, dx, dy,
                              dt, offset1, offset2);
#endif

  // __global__ void eno_x_cuda_2(double* U, double* F_, double dy, double r, int offset1, int offset2) {
  grid2.x = (Nx + 1 + blocksize - 1) / blocksize;
  grid2.y = (int(0.5 / dy) + 1 + blocksize_Y - 1) / blocksize_Y;
#ifdef CUDAV2
  eno_x_cuda_2_v2<<<grid2, block2>>>(U_cuda, F_cuda, dy, 1.f, offset1, offset2);
#else
  eno_x_cuda_2<<<grid2, block2>>>(U_cuda, F_cuda, dy, 1.f, offset1, offset2);
#endif

  // __global__ void eno_x_cuda_3(double* U, double* F_, double dx, double dy, double r, int offset1, int offset2) {
  grid3.x = ((2.0 / dx - 1.0 / dx) + 1 + blocksize - 1) / blocksize;
  grid3.y = (Ny - int(0.5 / dy) + 1 + blocksize_Y - 1) / blocksize_Y;
#ifdef CUDAV2
  eno_x_cuda_3_v2<<<grid3, block3>>>(U_cuda, F_cuda, dx, dy, 1.f, offset1, offset2);
#else
  eno_x_cuda_3<<<grid3, block3>>>(U_cuda, F_cuda, dx, dy, 1.f, offset1, offset2);
#endif
  cudaDeviceSynchronize();
  func_timer.stop("ENO_x_cuda");

  bound_Tensor(U, dx, dy);
  LF_y_Tensor(U, LAMDA_, G, Gp, Gd);
  ENO_y_Tensor(U, G, Gp, Gd, G_p, G_d, G_, q3p, q3d, dx, dy, dt);
  bound_Tensor(U, dx, dy);
  LF_y_Tensor(U, LAMDA_, G, Gp, Gd);
  ENO_y_Tensor(U, G, Gp, Gd, G_p, G_d, G_, q3p, q3d, dx, dy, dt);
  bound_Tensor(U, dx, dy);
  LF_x_Tensor(U, LAMDA_, F, Fp, Fd);
  // func_timer.start("ENO_x");
  ENO_x_Tensor(U, F, Fp, Fd, F_p, F_d, F_, q3p, q3d, dx, dy, dt);
  // func_timer.stop("ENO_x");
  bound_Tensor(U, dx, dy);
}

void alloc(double **ptr, const vector<int> &dims, int depth) {
  if (depth + 1 >= dims.size()) return;
  for (int i = 0; i < dims[depth]; ++i) {
    ptr[i] = (double *)malloc(dims[depth + 1] * sizeof(void *));
    // cudaMallocManaged(&(ptr[i]), dims[depth + 1] * sizeof(void *));
    alloc((double **)ptr[i], dims, depth + 1);
  }
}

double *alloc_nd(const vector<int> &dims) {
  double *output = nullptr;
  output = (double *)malloc(dims[0] * sizeof(void *));
  // cudaMallocManaged(&output, dims[0] * sizeof(void *));
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

  Tensor U(shape1, shape1_);
  Tensor U_(shape1, shape1_);
  Tensor Ut(shape1, shape1_);
  Tensor U1(shape1, shape1_);
  Tensor U2(shape1, shape1_);
  Tensor G(shape1, shape1_);
  Tensor G_(shape1, shape1_);
  Tensor F(shape1, shape1_);
  Tensor F_(shape1, shape1_);
  Tensor G_p(shape1, shape1_);
  Tensor G_d(shape1, shape1_);
  Tensor Gp(shape1, shape1_);
  Tensor Gd(shape1, shape1_);
  Tensor F_p(shape1, shape1_);
  Tensor F_d(shape1, shape1_);
  Tensor Fp(shape1, shape1_);
  Tensor Fd(shape1, shape1_);
  // double ***U = (double ***)alloc_nd(shape1);
  // double ***U_ = (double ***)alloc_nd(shape1);
  // double ***Ut = (double ***)alloc_nd(shape1);
  // double ***U1 = (double ***)alloc_nd(shape1);
  // double ***U2 = (double ***)alloc_nd(shape1);
  // double ***G = (double ***)alloc_nd(shape1);
  // double ***G_ = (double ***)alloc_nd(shape1);
  // double ***F = (double ***)alloc_nd(shape1);
  // double ***F_ = (double ***)alloc_nd(shape1);
  // double ***G_p = (double ***)alloc_nd(shape1);
  // double ***G_d = (double ***)alloc_nd(shape1);
  // double ***Gp = (double ***)alloc_nd(shape1);
  // double ***Gd = (double ***)alloc_nd(shape1);
  // double ***F_p = (double ***)alloc_nd(shape1);
  // double ***F_d = (double ***)alloc_nd(shape1);
  // double ***Fp = (double ***)alloc_nd(shape1);
  // double ***Fd = (double ***)alloc_nd(shape1);

  Tensor L_(shape2, shape2_);
  Tensor R_(shape2, shape2_);
  Tensor LAMDA_(shape2, shape2_);
  // double ****L_ = (double ****)alloc_nd(shape2);
  // double ****R_ = (double ****)alloc_nd(shape2);
  // double ****LAMDA_ = (double ****)alloc_nd(shape2);

  Tensor q3p(shape3, shape3_);
  Tensor q3d(shape3, shape3_);
  // double ****q3p = (double ****)alloc_nd(shape3);
  // double ****q3d = (double ****)alloc_nd(shape3);

  Tensor a_(shape4, shape4_);
  // double ***a_ = (double ***)alloc_nd(shape4);

  cudaMalloc((void **)&U_cuda, shape1_ * sizeof(double));
  cudaMalloc((void **)&F_cuda, shape1_ * sizeof(double));
  cudaMalloc((void **)&Fp_cuda, shape1_ * sizeof(double));
  cudaMalloc((void **)&Fd_cuda, shape1_ * sizeof(double));
  cudaMalloc((void **)&F_p_cuda, shape1_ * sizeof(double));
  cudaMalloc((void **)&F_d_cuda, shape1_ * sizeof(double));
  cudaMalloc((void **)&F__cuda, shape1_ * sizeof(double));
  cudaMalloc((void **)&q3p_cuda, shape3_ * sizeof(double));
  cudaMalloc((void **)&q3d_cuda, shape3_ * sizeof(double));

  grid.x = (Nx + 2 + blocksize - 1) / blocksize;
  grid.y = (Ny + 2 + blocksize_Y - 1) / blocksize_Y;
  block.x = blocksize;
  block.y = blocksize_Y;
  block2.x = blocksize;
  block2.y = blocksize_Y;
  block3.x = blocksize;
  block3.y = blocksize_Y;
#ifdef CUDAV2
  block.z = 3;
  block2.z = 3;
  block3.z = 3;
#endif

  // int shape1_ = (Nx + 7) * (Ny + 7) * 4;

  initial_Tensor(U, dx, dy);
  int n = 0;
  timeval s, e;
  gettimeofday(&s, NULL);
  while (T <= TT) {
    dt = CFL_Tensor(U, dx, dy, ENOCFL);
    ENO_Solver(U, U1, U2, F, Fp, Fd, F_p, F_d, F_, G, Gp, Gd, F_p, F_d, F_, LAMDA_, q3p, q3d, dx, dy, dt);
    T += dt;
    n++;
    virtual_clear_Tensor(U, dx, dy);
  }
  gettimeofday(&e, NULL);
  double ms = get_elapsed_time_ms(s, e);
  printf("ENO total time: %lf ms\n", ms);
  printf("ENO iter: %d\n", n);

  func_timer.show_all();

  if (argc == 2) {
    // result check
    printf("check dat file: %s\n", argv[1]);
    check_Tensor(argv[1], U);
  }
  return 0;
}