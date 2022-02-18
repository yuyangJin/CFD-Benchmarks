#include "muscl.h"
#include <cstdio>
#include "util.h"

double U[Nx + 7][Ny + 7][4], U_[Nx + 7][Ny + 7][4];
double Ut[Nx + 7][Ny + 7][4], U1[Nx + 7][Ny + 7][4], U2[Nx + 7][Ny + 7][4];
double L_[Nx + 7][Ny + 7][4][4], R_[Nx + 7][Ny + 7][4][4];
double a_[Nx + 7][Ny + 7][1], LAMDA_[Nx + 7][Ny + 7][4][4];
double G[Nx + 7][Ny + 7][4], G_[Nx + 7][Ny + 7][4], F[Nx + 7][Ny + 7][4], F_[Nx + 7][Ny + 7][4];

// MUSCL
double U_L[Nx + 7][Ny + 7][4], U_R[Nx + 7][Ny + 7][4];
double Fp[Nx + 7][Ny + 7][4], Fd[Nx + 7][Ny + 7][4];
double Gp[Nx + 7][Ny + 7][4], Gd[Nx + 7][Ny + 7][4];

timer func_timer;

void MUSCL_Solver(double U[Nx + 7][Ny + 7][4], double U_L[Nx + 7][Ny + 7][4], double U_R[Nx + 7][Ny + 7][4],
                  double Fp[Nx + 7][Ny + 7][4], double Fd[Nx + 7][Ny + 7][4], double F_[Nx + 7][Ny + 7][4],
                  double Gp[Nx + 7][Ny + 7][4], double Gd[Nx + 7][Ny + 7][4], double G_[Nx + 7][Ny + 7][4], double dx,
                  double dy, double& dt) {
  bound(U, dx, dy);
  dt = CFL(U, dx, dy, MUSCLCFL);
  func_timer.start("MUSCL_x");
  MUSCL_x(U, U_L, U_R, F_, Fp, Fd, dx, dy, dt / 2.0);
  func_timer.stop("MUSCL_x");
  bound(U, dx, dy);
  dt = CFL(U, dx, dy, MUSCLCFL);
  MUSCL_y(U, U_L, U_R, G_, Gp, Gd, dx, dy, dt / 2.0);
  bound(U, dx, dy);
  dt = CFL(U, dx, dy, MUSCLCFL);
  MUSCL_y(U, U_L, U_R, G_, Gp, Gd, dx, dy, dt / 2.0);
  bound(U, dx, dy);
  dt = CFL(U, dx, dy, MUSCLCFL);
  func_timer.start("MUSCL_x");
  MUSCL_x(U, U_L, U_R, F_, Fp, Fd, dx, dy, dt / 2.0);
  func_timer.stop("MUSCL_x");
  bound(U, dx, dy);
}

int main(int argc, char** argv) {
  double dx, dy, dt = 0, T = 0;
  initial(U, dx, dy);

  int n = 0;
  timeval s, e;
  gettimeofday(&s, NULL);
  while (T <= TT) {
    dt = CFL(U, dx, dy, MUSCLCFL);
    MUSCL_Solver(U, U_L, U_R, Fp, Fd, F_, Gp, Gd, G_, dx, dy, dt);
    T += dt;
    n++;
    virtual_clear(U, dx, dy);
  }
  gettimeofday(&e, NULL);
  double ms = get_elapsed_time_ms(s, e);
  printf("MUSCL total time: %lf ms\n", ms);
  printf("MUSCL iter: %d\n", n);

  func_timer.show_all();

  if (argc == 2) {
    // result check
    printf("check dat file: %s\n", argv[1]);
    check(argv[1], U);
  }
  return 0;
}