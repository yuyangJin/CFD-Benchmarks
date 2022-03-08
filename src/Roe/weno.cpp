#include "weno.h"
#include <cstdio>
#include "util.h"

double U[Nx + 7][Ny + 7][4], U_[Nx + 7][Ny + 7][4];
double Ut[Nx + 7][Ny + 7][4], U1[Nx + 7][Ny + 7][4], U2[Nx + 7][Ny + 7][4];
double L_[Nx + 7][Ny + 7][4][4], R_[Nx + 7][Ny + 7][4][4];
double a_[Nx + 7][Ny + 7][1], LAMDA_[Nx + 7][Ny + 7][4][4];
double G[Nx + 7][Ny + 7][4], G_[Nx + 7][Ny + 7][4], F[Nx + 7][Ny + 7][4], F_[Nx + 7][Ny + 7][4];

// WENO
double q3p[Nx + 7][Ny + 7][4][3], q3d[Nx + 7][Ny + 7][4][3];  //这里q3p和q3d 表示模板的正负差值f
double ISd[Nx + 7][Ny + 7][4][3], ISp[Nx + 7][Ny + 7][4][3];  //这里ISd和omegap以及alphap的定义见教材
double Fp[Nx + 7][Ny + 7][4],
    Fd[Nx + 7][Ny + 7][4];  // Fp和Fd表示对F作Sterger-Warming分裂以后的正负F 其中Fp=Fplus Fd=Fdecrease
double Gp[Nx + 7][Ny + 7][4], Gd[Nx + 7][Ny + 7][4];    // G表示y方向的量
double F_p[Nx + 7][Ny + 7][4], F_d[Nx + 7][Ny + 7][4];  //这里F_p表示正的流通量F[i+1/2]
double G_p[Nx + 7][Ny + 7][4], G_d[Nx + 7][Ny + 7][4];
double omegap[Nx + 7][Ny + 7][4][3], omegad[Nx + 7][Ny + 7][4][3];
double alphap[Nx + 7][Ny + 7][4][3], alphad[Nx + 7][Ny + 7][4][3];

timer func_timer;
int case_id;

void WENO_Solver(double U[Nx + 7][Ny + 7][4], double U1[Nx + 7][Ny + 7][4], double U2[Nx + 7][Ny + 7][4],
                 double ISp[Nx + 7][Ny + 7][4][3], double ISd[Nx + 7][Ny + 7][4][3],
                 double omegap[Nx + 7][Ny + 7][4][3], double omegad[Nx + 7][Ny + 7][4][3],
                 double alphap[Nx + 7][Ny + 7][4][3], double alphad[Nx + 7][Ny + 7][4][3],
                 double q3p[Nx + 7][Ny + 7][4][3], double q3d[Nx + 7][Ny + 7][4][3],
                 double LAMDA_[Nx + 7][Ny + 7][4][4], double Fp[Nx + 7][Ny + 7][4], double Fd[Nx + 7][Ny + 7][4],
                 double F_p[Nx + 7][Ny + 7][4], double F_d[Nx + 7][Ny + 7][4], double F[Nx + 7][Ny + 7][4],
                 double F_[Nx + 7][Ny + 7][4], double G[Nx + 7][Ny + 7][4], double G_[Nx + 7][Ny + 7][4],
                 double Gp[Nx + 7][Ny + 7][4], double Gd[Nx + 7][Ny + 7][4], double G_p[Nx + 7][Ny + 7][4],
                 double G_d[Nx + 7][Ny + 7][4], double dx, double dy, double& dt) {
  bound(U, dx, dy, case_id);
  LF_x(U, LAMDA_, F, Fp, Fd);
  func_timer.start("WENO_x");
  WENO_x(U, ISp, ISd, omegap, omegad, alphap, alphad, q3p, q3d, Fp, Fd, F_p, F_d, F_, dx, dy, dt);
  func_timer.stop("WENO_x");
  bound(U, dx, dy, case_id);
  LF_y(U, LAMDA_, G, Gp, Gd);
  WENO_y(U, ISp, ISd, omegap, omegad, alphap, alphad, q3p, q3d, Gp, Gd, G_p, G_d, G_, dx, dy, dt);
  bound(U, dx, dy, case_id);
  LF_y(U, LAMDA_, G, Gp, Gd);
  WENO_y(U, ISp, ISd, omegap, omegad, alphap, alphad, q3p, q3d, Gp, Gd, G_p, G_d, G_, dx, dy, dt);
  bound(U, dx, dy, case_id);
  LF_x(U, LAMDA_, F, Fp, Fd);
  func_timer.start("WENO_x");
  WENO_x(U, ISp, ISd, omegap, omegad, alphap, alphad, q3p, q3d, Fp, Fd, F_p, F_d, F_, dx, dy, dt);
  func_timer.stop("WENO_x");
  bound(U, dx, dy, case_id);
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
    dt = CFL(U, dx, dy, WENOCFL);
    WENO_Solver(U, U1, U2, ISp, ISd, omegap, omegad, alphap, alphad, q3p, q3d, LAMDA_, Fp, Fd, F_p, F_d, F, F_, G, G_,
                Gp, Gd, G_p, G_d, dx, dy, dt);
    T += dt;
    n++;
    virtual_clear(U, dx, dy);
  }
  gettimeofday(&e, NULL);
  double ms = get_elapsed_time_ms(s, e);
  printf("WENO total time: %lf ms\n", ms);
  printf("WENO iter: %d\n", n);

  func_timer.show_all();

  if (argc == 2) {
    // result check
    printf("check dat file: %s\n", argv[1]);
    check(argv[1], U);
  }
  return 0;
}