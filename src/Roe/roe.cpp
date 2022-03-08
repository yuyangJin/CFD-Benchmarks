#include "roe.h"
#include <cstdio>
#include "util.h"

double U[Nx + 7][Ny + 7][4], U_[Nx + 7][Ny + 7][4];
double Ut[Nx + 7][Ny + 7][4], U1[Nx + 7][Ny + 7][4], U2[Nx + 7][Ny + 7][4];
double L_[Nx + 7][Ny + 7][4][4], R_[Nx + 7][Ny + 7][4][4];
double a_[Nx + 7][Ny + 7][1], LAMDA_[Nx + 7][Ny + 7][4][4];
double G[Nx + 7][Ny + 7][4], G_[Nx + 7][Ny + 7][4], F[Nx + 7][Ny + 7][4], F_[Nx + 7][Ny + 7][4];

// ROE
// alpha_表示Roe方法中的alpha=L*(U[i+1]-U[i]).theta表示特征向量叠加的和
double alpha_[Nx + 7][Ny + 7][4], theta[Nx + 7][Ny + 7][4];

timer func_timer;
int case_id;


void Roe_Solver(double U[Nx + 7][Ny + 7][4], double U_[Nx + 7][Ny + 7][4], double LAMDA_[Nx + 7][Ny + 7][4][4],
                double L_[Nx + 7][Ny + 7][4][4], double R_[Nx + 7][Ny + 7][4][4], double alpha_[Nx + 7][Ny + 7][4],
                double theta[Nx + 7][Ny + 7][4], double F[Nx + 7][Ny + 7][4], double F_[Nx + 7][Ny + 7][4],
                double G[Nx + 7][Ny + 7][4], double G_[Nx + 7][Ny + 7][4], double dx, double dy, double& dt,
                double a_[Nx + 7][Ny + 7][1])  // Roe求解器
{
  bound(U, dx, dy, case_id);
  func_timer.start("Roe_x");
  Roe_x(U, U_, L_, R_, F_, F, a_, LAMDA_, alpha_, theta, dx, dy, dt);
  func_timer.stop("Roe_x");
  bound(U, dx, dy, case_id);
  Roe_y(U, U_, L_, R_, G_, G, a_, LAMDA_, alpha_, theta, dx, dy, dt);
  bound(U, dx, dy, case_id);
  Roe_y(U, U_, L_, R_, G_, G, a_, LAMDA_, alpha_, theta, dx, dy, dt);
  bound(U, dx, dy, case_id);
  func_timer.start("Roe_x");
  Roe_x(U, U_, L_, R_, F_, F, a_, LAMDA_, alpha_, theta, dx, dy, dt);
  func_timer.stop("Roe_x");
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
    dt = CFL(U, dx, dy, ROECFL);
    Roe_Solver(U, U_, LAMDA_, L_, R_, alpha_, theta, F, F_, G, G_, dx, dy, dt, a_);
    T += dt;
    n++;
    virtual_clear(U, dx, dy);
  }
  gettimeofday(&e, NULL);
  double ms = get_elapsed_time_ms(s, e);
  printf("Roe total time: %lf ms\n", ms);
  printf("Roe iter: %d\n", n);

  func_timer.show_all();

  if (argc == 2) {
    // result check
    printf("check dat file: %s\n", argv[1]);
    check(argv[1], U);
  }
  return 0;
}