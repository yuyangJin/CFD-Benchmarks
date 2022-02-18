#include "nnd.h"
#include <cstdio>
#include "util.h"

double U[Nx + 7][Ny + 7][4], U_[Nx + 7][Ny + 7][4];
double Ut[Nx + 7][Ny + 7][4], U1[Nx + 7][Ny + 7][4], U2[Nx + 7][Ny + 7][4];
double L_[Nx + 7][Ny + 7][4][4], R_[Nx + 7][Ny + 7][4][4];
double a_[Nx + 7][Ny + 7][1], LAMDA_[Nx + 7][Ny + 7][4][4];
double G[Nx + 7][Ny + 7][4], G_[Nx + 7][Ny + 7][4], F[Nx + 7][Ny + 7][4], F_[Nx + 7][Ny + 7][4];

// NND
double Fp[Nx + 7][Ny + 7][4], Fd[Nx + 7][Ny + 7][4];
double Gp[Nx + 7][Ny + 7][4], Gd[Nx + 7][Ny + 7][4];

timer func_timer;

void NND_Solver(double U[Nx + 7][Ny + 7][4], double Fp[Nx + 7][Ny + 7][4], double Fd[Nx + 7][Ny + 7][4],
                double Gp[Nx + 7][Ny + 7][4], double Gd[Nx + 7][Ny + 7][4], double dx, double dy, double& dt) {
  bound(U, dx, dy);
  func_timer.start("NND_x");
  NND_x(U, Fp, Fd, dx, dy, dt);
  func_timer.stop("NND_x");
  bound(U, dx, dy);
  NND_y(U, Gp, Gd, dx, dy, dt);
  bound(U, dx, dy);
  NND_y(U, Gp, Gd, dx, dy, dt);
  bound(U, dx, dy);
  func_timer.start("NND_x");
  NND_x(U, Fp, Fd, dx, dy, dt);
  func_timer.stop("NND_x");
}

int main(int argc, char** argv) {
  double dx, dy, dt = 0, T = 0;
  initial(U, dx, dy);

  int n = 0;
  timeval s, e;
  gettimeofday(&s, NULL);
  while (T <= TT) {
    dt = CFL(U, dx, dy, NNDCFL);
    NND_Solver(U, Fp, Fd, Gp, Gd, dx, dy, dt);
    T += dt;
    n++;
    virtual_clear(U, dx, dy);
  }
  gettimeofday(&e, NULL);
  double ms = get_elapsed_time_ms(s, e);
  printf("NND total time: %lf ms\n", ms);
  printf("NND iter: %d\n", n);

  func_timer.show_all();

  if (argc == 2) {
    // result check
    printf("check dat file: %s\n", argv[1]);
    check(argv[1], U);
  }
  return 0;
}