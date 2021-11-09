#ifndef RIEMANN2D_H
#define RIEMANN2D_H

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

using namespace std;

#define GAMA 1.4 //气体常数
#define PI 3.1415926

#define MIN(x,y) (((x)<(y))?(x):(y))
#define MAX(x,y) (((x)>(y))?(x):(y))


#define Lx 2.0 //计算区域
#define Ly 2.0

#define TT 0.6 //总时间
#define Sf 0.2 //时间步长因子
#define epsilon 0.12 //人工粘性项小参数

#define Jx 100 
#define Jy 100 //网格数

//全局变量
double U[Jx+2][Jy+2][4],U_old[Jx+2][Jy+2][4], F[Jx+2][Jy+2][4],G[Jx+2][Jy+2][4];
double UF[Jx+2][Jy+2][4],UG[Jx+2][Jy+2][4];
double rou[Jx+2][Jy+2],u[Jx+2][Jy+2],v[Jx+2][Jy+2],h[Jx+2][Jy+2],p[Jx+2][Jy+2];

void Init(double U[Jx+2][Jy+2][4],double& dx,double& dy);
void bound(double U[Jx+2][Jy+2][4]);
double CFL(double U[Jx+2][Jy+2][4],double dx,double dy);
double Q_z(double x);
double delta(double x,double r);
double minmod(double x,double y,double z);
void U2F(double UF[4],double F[4]);
void U2G(double UG[4],double G[4]);
void Calculate_field_p(double U[Jx+2][Jy+2][4],double rou[Jx+2][Jy+2],double u[Jx+2][Jy+2],double v[Jx+2][Jy+2],double h[Jx+2][Jy+2],double p[Jx+2][Jy+2]);
void lamda_x(double lamda_x[4],double u,double v,double a);
void lamda_y(double lamda_y[4],double u,double v,double a);
void LR_x(double R[4][4],double L[4][4],double u,double v,double a );
void LR_y(double R[4][4],double L[4][4],double u,double v,double a);
void save_U_last();
void LLx(double U[Jx+2][Jy+2][4],double dx,double dt);
void LLy(double U[Jx+2][Jy+2][4],double dy,double dt);
void symmetry_TVD(double U[Jx+2][Jy+2][4],double dx,double dy,double dt);
double error(double U1[Jx+2][Jy+2][4],double U2[Jx+2][Jy+2][4]);
void Output(double U[Jx+2][Jy+2][4],double dx,double dy,double T);

#endif