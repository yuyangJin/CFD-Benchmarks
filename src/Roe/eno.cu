#include <cstdio>
#include <vector>
using std::vector;

#include "eno.h"
#include "util.cuh"

timer func_timer;
void ENO_Solver(double ***U, double ***U1, double ***U2, double ***F, double ***Fp, double ***Fd, double ***F_p,
                double ***F_d, double ***F_, double ***G, double ***Gp, double ***Gd, double ***G_p, double ***G_d,
                double ***G_, double ****LAMDA_, double ****q3p, double ****q3d, double dx, double dy, double &dt) {
  bound(U, dx, dy);
  LF_x(U, LAMDA_, F, Fp, Fd);
  func_timer.start("ENO_x");
  ENO_x(U, F, Fp, Fd, F_p, F_d, F_, q3p, q3d, dx, dy, dt);
  func_timer.stop("ENO_x");
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
  vector<int> shape2 = {Nx + 7, Ny + 7, 4, 4};
  vector<int> shape3 = {Nx + 7, Ny + 7, 4, 3};
  vector<int> shape4 = {Nx + 7, Ny + 7, 1};

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

  func_timer.show_all();

  if (argc == 2) {
    // result check
    printf("check dat file: %s\n", argv[1]);
    check(argv[1], U);
  }
  return 0;
}