#include "eno.h"
#include <cstdio>
#include "util.h"

double U[Nx + 7][Ny + 7][4], U_[Nx + 7][Ny + 7][4];
double Ut[Nx + 7][Ny + 7][4], U1[Nx + 7][Ny + 7][4], U2[Nx + 7][Ny + 7][4];
double L_[Nx + 7][Ny + 7][4][4], R_[Nx + 7][Ny + 7][4][4];
double a_[Nx + 7][Ny + 7][1], LAMDA_[Nx + 7][Ny + 7][4][4];
double G[Nx + 7][Ny + 7][4], G_[Nx + 7][Ny + 7][4], F[Nx + 7][Ny + 7][4], F_[Nx + 7][Ny + 7][4];

// ENO
double q3p[Nx + 7][Ny + 7][4][3], q3d[Nx + 7][Ny + 7][4][3];
double F_p[Nx + 7][Ny + 7][4], F_d[Nx + 7][Ny + 7][4];
double G_p[Nx + 7][Ny + 7][4], G_d[Nx + 7][Ny + 7][4];
double Fp[Nx + 7][Ny + 7][4], Fd[Nx + 7][Ny + 7][4];
double Gp[Nx + 7][Ny + 7][4], Gd[Nx + 7][Ny + 7][4];

timer func_timer;

void ENO_Solver(double U[Nx + 7][Ny + 7][4], double U1[Nx + 7][Ny + 7][4], double U2[Nx + 7][Ny + 7][4],
                double F[Nx + 7][Ny + 7][4], double Fp[Nx + 7][Ny + 7][4], double Fd[Nx + 7][Ny + 7][4],
                double F_p[Nx + 7][Ny + 7][4], double F_d[Nx + 7][Ny + 7][4], double F_[Nx + 7][Ny + 7][4],
                double G[Nx + 7][Ny + 7][4], double Gp[Nx + 7][Ny + 7][4], double Gd[Nx + 7][Ny + 7][4],
                double G_p[Nx + 7][Ny + 7][4], double G_d[Nx + 7][Ny + 7][4], double G_[Nx + 7][Ny + 7][4],
                double LAMDA_[Nx + 7][Ny + 7][4][4], double q3p[Nx + 7][Ny + 7][4][3], double q3d[Nx + 7][Ny + 7][4][3],
                double dx, double dy, double& dt) {
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

int main(int argc, char** argv) {
  double dx, dy, dt = 0, T = 0;
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

  double enox_ms;
  int enox_n;
  func_timer.info("ENO_x", enox_ms, enox_n);
  printf("ENO_x: %lf ms / %d count, avg %lf ms\n", enox_ms, enox_n, enox_ms / enox_n);

  if (argc == 2) {
    // result check
    printf("check dat file: %s\n", argv[1]);
    check(argv[1], U);
  }
  return 0;
}