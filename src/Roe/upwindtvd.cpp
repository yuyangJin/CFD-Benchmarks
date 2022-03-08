#include "upwindtvd.h"
#include <cstdio>
#include "util.h"

double U[Nx + 7][Ny + 7][4], U_[Nx + 7][Ny + 7][4];
double Ut[Nx + 7][Ny + 7][4], U1[Nx + 7][Ny + 7][4], U2[Nx + 7][Ny + 7][4];
double L_[Nx + 7][Ny + 7][4][4], R_[Nx + 7][Ny + 7][4][4];
double a_[Nx + 7][Ny + 7][1], LAMDA_[Nx + 7][Ny + 7][4][4];
double G[Nx + 7][Ny + 7][4], G_[Nx + 7][Ny + 7][4], F[Nx + 7][Ny + 7][4], F_[Nx + 7][Ny + 7][4];

// UpWindTVD
double g_[Nx + 7][Ny + 7][4];  //这里alpha_,g_,g,gama_,Q_,theta分别表示迎风TVD算法中用到的中间变量,其意义见教材
double g[Nx + 7][Ny + 7][4], gama_[Nx + 7][Ny + 7][4], Q_[Nx + 7][Ny + 7][4];  //其中theta定义于Roe方法中的意义不同
double alpha_[Nx + 7][Ny + 7][4],
    theta[Nx + 7][Ny + 7][4];  // alpha_表示Roe方法中的alpha=L*(U[i+1]-U[i]).theta表示特征向量叠加的和

timer func_timer;
int case_id;


void UpWindTVD_Solver(double U[Nx + 7][Ny + 7][4], double U_[Nx + 7][Ny + 7][4], double LAMDA_[Nx + 7][Ny + 7][4][4],
                      double L_[Nx + 7][Ny + 7][4][4], double R_[Nx + 7][Ny + 7][4][4],
                      double alpha_[Nx + 7][Ny + 7][4], double g_[Nx + 7][Ny + 7][4], double g[Nx + 7][Ny + 7][4],
                      double gama_[Nx + 7][Ny + 7][4], double Q_[Nx + 7][Ny + 7][4], double theta[Nx + 7][Ny + 7][4],
                      double F[Nx + 7][Ny + 7][4], double F_[Nx + 7][Ny + 7][4], double G[Nx + 7][Ny + 7][4],
                      double G_[Nx + 7][Ny + 7][4], double dx, double dy, double& dt,
                      double a_[Nx + 7][Ny + 7][1])  // TVD求解器
{
  bound(U, dx, dy, case_id);  //先对边界以及虚拟节点赋值
  func_timer.start("UpWindTVD_x");
  UpWindTVD_x(U, U_, a_, LAMDA_, L_, R_, alpha_, g_, g, gama_, Q_, theta, F, F_, dx, dy, dt);  // TVD算法求解U
  func_timer.stop("UpWindTVD_x");
  bound(U, dx, dy, case_id);
  UpWindTVD_y(U, U_, a_, LAMDA_, L_, R_, alpha_, g_, g, gama_, Q_, theta, G, G_, dx, dy, dt);
  bound(U, dx, dy, case_id);
  UpWindTVD_y(U, U_, a_, LAMDA_, L_, R_, alpha_, g_, g, gama_, Q_, theta, G, G_, dx, dy, dt);
  bound(U, dx, dy, case_id);
  func_timer.start("UpWindTVD_x");
  UpWindTVD_x(U, U_, a_, LAMDA_, L_, R_, alpha_, g_, g, gama_, Q_, theta, F, F_, dx, dy, dt);
  func_timer.stop("UpWindTVD_x");
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
    dt = CFL(U, dx, dy, UPWINDTVDCFL);
    UpWindTVD_Solver(U, U_, LAMDA_, L_, R_, alpha_, g_, g, gama_, Q_, theta, F, F_, G, G_, dx, dy, dt, a_);
    T += dt;
    n++;
    virtual_clear(U, dx, dy);
  }
  gettimeofday(&e, NULL);
  double ms = get_elapsed_time_ms(s, e);
  printf("UpWindTVD total time: %lf ms\n", ms);
  printf("UpWindTVD iter: %d\n", n);

  func_timer.show_all();

  if (argc == 2) {
    // result check
    printf("check dat file: %s\n", argv[1]);
    check(argv[1], U);
  }
  return 0;
}