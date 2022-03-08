#include "compact.h"
#include <cstdio>
#include "util.h"

double U[Nx + 7][Ny + 7][4], U_[Nx + 7][Ny + 7][4];
double Ut[Nx + 7][Ny + 7][4], U1[Nx + 7][Ny + 7][4], U2[Nx + 7][Ny + 7][4];
double L_[Nx + 7][Ny + 7][4][4], R_[Nx + 7][Ny + 7][4][4];
double a_[Nx + 7][Ny + 7][1], LAMDA_[Nx + 7][Ny + 7][4][4];
double G[Nx + 7][Ny + 7][4], G_[Nx + 7][Ny + 7][4], F[Nx + 7][Ny + 7][4], F_[Nx + 7][Ny + 7][4];

// COMPACT
double fp[Nx + 7][Ny + 7][4], fd[Nx + 7][Ny + 7][4];                     //紧致格式中用小写的f表示F
double gp[Nx + 7][Ny + 7][4], gd[Nx + 7][Ny + 7][4], p[Nx + 7][Ny + 7];  // p是压力
double Fp[Nx + 7][Ny + 7][4],
    Fd[Nx + 7][Ny + 7][4];  // Fp和Fd表示对F作Sterger-Warming分裂以后的正负F 其中Fp=Fplus Fd=Fdecrease
double Gp[Nx + 7][Ny + 7][4], Gd[Nx + 7][Ny + 7][4];  // G表示y方向的量

timer func_timer;
int case_id;

void Compact_Solver(double U[Nx + 7][Ny + 7][4], double Ut[Nx + 7][Ny + 7][4], double U1[Nx + 7][Ny + 7][4],
                    double U2[Nx + 7][Ny + 7][4], double Fp[Nx + 7][Ny + 7][4], double Fd[Nx + 7][Ny + 7][4],
                    double fp[Nx + 7][Ny + 7][4], double fd[Nx + 7][Ny + 7][4], double Gp[Nx + 7][Ny + 7][4],
                    double Gd[Nx + 7][Ny + 7][4], double gp[Nx + 7][Ny + 7][4], double gd[Nx + 7][Ny + 7][4], double dx,
                    double dy, double dt, double p[Nx + 7][Ny + 7]) {
  bound(U, dx, dy, case_id);
  Compact_x(U, Ut, Fp, Fd, fp, fd, dx, dy, dt, p);
  func_timer.start("Compact_RK_x");
  Compact_RK_x(U, Ut, U1, U2, Fp, Fd, fp, fd, dx, dy, dt, p, LAMDA_);
  func_timer.stop("Compact_RK_x");
  bound(U, dx, dy, case_id);
  Compact_y(U, Ut, Gp, Gd, gp, gd, dx, dy, dt, p);
  Compact_RK_y(U, Ut, U1, U2, Gp, Gd, gp, gd, dx, dy, dt, p, LAMDA_);
  bound(U, dx, dy, case_id);
  Compact_y(U, Ut, Gp, Gd, gp, gd, dx, dy, dt, p);
  Compact_RK_y(U, Ut, U1, U2, Gp, Gd, gp, gd, dx, dy, dt, p, LAMDA_);
  bound(U, dx, dy, case_id);
  Compact_x(U, Ut, Fp, Fd, fp, fd, dx, dy, dt, p);
  func_timer.start("Compact_RK_x");
  Compact_RK_x(U, Ut, U1, U2, Fp, Fd, fp, fd, dx, dy, dt, p, LAMDA_);
  func_timer.stop("Compact_RK_x");
}

int main(int argc, char** argv) {
  if (argc != 2) {
    printf("Please enter: \"./eno [case id]\" ");
    return 0;
  }

  case_id = atoi(argv[1]);

  double dx, dy, dt = 0, T = 0;
  initial(U, dx, dy);

  int n = 0;
  timeval s, e;
  gettimeofday(&s, NULL);
  while (T <= TT) {
    dt = CFL(U, dx, dy, COMPACTCFL);

    Compact_Solver(U, Ut, U1, U2, Fp, Fd, fp, fd, Gp, Gd, gp, gd, dx, dy, dt, p);
    T += dt;
    n++;
    virtual_clear(U, dx, dy);
  }
  gettimeofday(&e, NULL);
  double ms = get_elapsed_time_ms(s, e);
  printf("COMPACT total time: %lf ms\n", ms);
  printf("COMPACT iter: %d\n", n);

  func_timer.show_all();

  if (argc == 2) {
    // result check
    printf("check dat file: %s\n", argv[1]);
    check(argv[1], U);
  }
  return 0;
}