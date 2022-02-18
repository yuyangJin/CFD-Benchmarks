#ifndef __MUSCL_H
#define __MUSCL_H

#include "util.h"

// MUSCL算法中用到的两个常数
#define K -1
#define BETA 1

#ifndef MUSCL_FUSION
//作用:求Fp和Fd
void U2FpFd_AVG(double U_L[4], double U_R[4], double Fp[4], double Fd[4]) {
  double H, rou, u, v, p, a, lambda1, lambda2, lambda3, lambda4, lambda1p, lambda2p, lambda3p, lambda4p, lambda1d,
      lambda2d, lambda3d, lambda4d;
  //先算Fp即F+
  rou = U_L[0];
  u = U_L[1] / U_L[0];
  v = U_L[2] / U_L[0];
  p = (GAMA - 1) * (U_L[3] - 0.5 * rou * (u * u + v * v));
  a = sqrt(GAMA * p / rou);
  H = a * a / (GAMA - 1) + 0.5 * (u * u + v * v);
  lambda1 = u;
  lambda2 = u;
  lambda3 = u - a;
  lambda4 = u + a;
  //计算正特征值
  lambda1p = 0.5 * (fabs(lambda1) + lambda1);
  lambda2p = 0.5 * (fabs(lambda2) + lambda2);
  lambda3p = 0.5 * (fabs(lambda3) + lambda3);
  lambda4p = 0.5 * (fabs(lambda4) + lambda4);
  //特征值修正
  lambda1p = 0.5 * (lambda1p + sqrt(lambda1p * lambda1p + 1e-16));
  lambda2p = 0.5 * (lambda2p + sqrt(lambda2p * lambda2p + 1e-16));
  lambda3p = 0.5 * (lambda3p + sqrt(lambda3p * lambda3p + 1e-16));
  lambda4p = 0.5 * (lambda4p + sqrt(lambda4p * lambda4p + 1e-16));

  Fp[0] = rou / 2.0 / GAMA * (2 * (GAMA - 1) * lambda1p + lambda3p + lambda4p);
  Fp[1] = rou / 2.0 / GAMA * (2 * u * (GAMA - 1) * lambda1p + (u - a) * lambda3p + (u + a) * lambda4p);
  Fp[2] = rou / 2.0 / GAMA * (v * 2 * (GAMA - 1) * lambda1p + v * lambda3p + v * lambda4p);
  Fp[3] =
      rou / 2.0 / GAMA * (2 * ((GAMA - 1) * H - a * a) * lambda1p + (H - a * u) * lambda3p + (H + a * u) * lambda4p);
  //现在就算Fd即F-
  rou = U_R[0];
  u = U_R[1] / U_R[0];
  v = U_R[2] / U_R[0];
  p = (GAMA - 1) * (U_R[3] - 0.5 * rou * (u * u + v * v));
  a = sqrt(GAMA * p / rou);
  H = a * a / (GAMA - 1) + 0.5 * (u * u + v * v);
  lambda1 = u;
  lambda2 = u;
  lambda3 = u - a;
  lambda4 = u + a;
  lambda1d = -0.5 * (fabs(lambda1) - lambda1);
  lambda2d = -0.5 * (fabs(lambda2) - lambda2);
  lambda3d = -0.5 * (fabs(lambda3) - lambda3);
  lambda4d = -0.5 * (fabs(lambda4) - lambda4);
  lambda1d = 0.5 * (lambda1d - sqrt(lambda1d * lambda1d + 1e-16));
  lambda2d = 0.5 * (lambda2d - sqrt(lambda2d * lambda2d + 1e-16));
  lambda3d = 0.5 * (lambda3d - sqrt(lambda3d * lambda3d + 1e-16));
  lambda4d = 0.5 * (lambda4d - sqrt(lambda4d * lambda4d + 1e-16));

  Fd[0] = rou / 2.0 / GAMA * (2 * (GAMA - 1) * lambda1d + lambda3d + lambda4d);
  Fd[1] = rou / 2.0 / GAMA * (2 * u * (GAMA - 1) * lambda1d + (u - a) * lambda3d + (u + a) * lambda4d);
  Fd[2] = rou / 2.0 / GAMA * (v * 2 * (GAMA - 1) * lambda1d + v * lambda3d + v * lambda4d);
  Fd[3] =
      rou / 2.0 / GAMA * (2 * ((GAMA - 1) * H - a * a) * lambda1d + (H - a * u) * lambda3d + (H + a * u) * lambda4d);
}
void MUSCL_x(double U[Nx + 7][Ny + 7][4], double U_L[Nx + 7][Ny + 7][4], double U_R[Nx + 7][Ny + 7][4],
             double F_[Nx + 7][Ny + 7][4], double Fp[Nx + 7][Ny + 7][4], double Fd[Nx + 7][Ny + 7][4], double dx,
             double dy, double dt) {
  int i, j, k;
  double r = dt / dx;

  //计算U_L和U_R
  for (i = 2; i <= Nx + 3; i++)
    for (j = 2; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) {
          U_L[i][j][k] = U[i][j][k] +
                         0.25 * (1 - K) * minmod(U[i][j][k] - U[i - 1][j][k], BETA * (U[i + 1][j][k] - U[i][j][k])) +
                         0.25 * (1 + K) * minmod(U[i + 1][j][k] - U[i][j][k], BETA * (U[i][j][k] - U[i - 1][j][k]));
          U_R[i][j][k] =
              U[i + 1][j][k] -
              0.25 * (1 - K) * minmod(U[i + 2][j][k] - U[i + 1][j][k], BETA * (U[i + 1][j][k] - U[i][j][k])) -
              0.25 * (1 + K) * minmod(U[i + 1][j][k] - U[i][j][k], BETA * (U[i + 2][j][k] - U[i + 1][j][k]));
        }
  //计算Fp和Fd
  for (i = 2; i <= Nx + 3; i++)
    for (j = 2; j <= Ny + 3; j++)
      if (U[i][j][0] != 0) U2FpFd_AVG(U_L[i][j], U_R[i][j], Fp[i][j], Fd[i][j]);
  //计算F_
  for (i = 2; i <= Nx + 3; i++)
    for (j = 2; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) F_[i][j][k] = Fp[i][j][k] + Fd[i][j][k];
  //计算U
  for (i = 3; i <= Nx + 3; i++)
    for (j = 3; j <= int(0.5 / dy) + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) U[i][j][k] = U[i][j][k] - r * (F_[i][j][k] - F_[i - 1][j][k]);

  for (i = int(1.0 / dx) + 3; i <= int(2.0 / dx) + 3; i++)
    for (j = int(0.5 / dy) + 4; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) U[i][j][k] = U[i][j][k] - r * (F_[i][j][k] - F_[i - 1][j][k]);
}
//算Gp和Gd
void U2GpGd_AVG(double U_L[4], double U_R[4], double Gp[4], double Gd[4]) {
  double H, rou, u, v, p, a, lambda1, lambda2, lambda3, lambda4, lambda1p, lambda2p, lambda3p, lambda4p;
  double lambda1d, lambda2d, lambda3d, lambda4d;
  //计算Gp即G+的函数
  rou = U_L[0];
  u = U_L[1] / U_L[0];
  v = U_L[2] / U_L[0];
  p = (GAMA - 1) * (U_L[3] - 0.5 * rou * (u * u + v * v));
  a = sqrt(GAMA * p / rou);
  H = a * a / (GAMA - 1) + 0.5 * (u * u + v * v);
  lambda1 = v;
  lambda2 = v;
  lambda3 = v - a;
  lambda4 = v + a;
  lambda1p = 0.5 * (fabs(lambda1) + lambda1);
  lambda2p = 0.5 * (fabs(lambda2) + lambda2);
  lambda3p = 0.5 * (fabs(lambda3) + lambda3);
  lambda4p = 0.5 * (fabs(lambda4) + lambda4);
  lambda1p = 0.5 * (lambda1p + sqrt(lambda1p * lambda1p + 1e-16));
  lambda2p = 0.5 * (lambda2p + sqrt(lambda2p * lambda2p + 1e-16));
  lambda3p = 0.5 * (lambda3p + sqrt(lambda3p * lambda3p + 1e-16));
  lambda4p = 0.5 * (lambda4p + sqrt(lambda4p * lambda4p + 1e-16));

  Gp[0] = rou / 2.0 / GAMA * (2 * (GAMA - 1) * lambda1p + lambda3p + lambda4p);
  Gp[1] = rou / 2.0 / GAMA * (u * 2 * (GAMA - 1) * lambda1p + u * lambda3p + u * lambda4p);
  Gp[2] = rou / 2.0 / GAMA * (2 * v * (GAMA - 1) * lambda1p + (v - a) * lambda3p + (v + a) * lambda4p);
  Gp[3] =
      rou / 2.0 / GAMA * (2 * ((GAMA - 1) * H - a * a) * lambda1p + (H - a * v) * lambda3p + (H + a * v) * lambda4p);
  //计算Gd即G-的函数
  rou = U_R[0];
  u = U_R[1] / U_R[0];
  v = U_R[2] / U_R[0];
  p = (GAMA - 1) * (U_R[3] - 0.5 * rou * (u * u + v * v));
  a = sqrt(GAMA * p / rou);
  H = a * a / (GAMA - 1) + 0.5 * (u * u + v * v);
  lambda1 = v;
  lambda2 = v;
  lambda3 = v - a;
  lambda4 = v + a;
  lambda1d = -0.5 * (fabs(lambda1) - lambda1);
  lambda2d = -0.5 * (fabs(lambda2) - lambda2);
  lambda3d = -0.5 * (fabs(lambda3) - lambda3);
  lambda4d = -0.5 * (fabs(lambda4) - lambda4);
  lambda1d = 0.5 * (lambda1d - sqrt(lambda1d * lambda1d + 1e-16));
  lambda2d = 0.5 * (lambda2d - sqrt(lambda2d * lambda2d + 1e-16));
  lambda3d = 0.5 * (lambda3d - sqrt(lambda3d * lambda3d + 1e-16));
  lambda4d = 0.5 * (lambda4d - sqrt(lambda4d * lambda4d + 1e-16));

  Gd[0] = rou / 2.0 / GAMA * (2 * (GAMA - 1) * lambda1d + lambda3d + lambda4d);
  Gd[1] = rou / 2.0 / GAMA * (u * 2 * (GAMA - 1) * lambda1d + u * lambda3d + u * lambda4d);
  Gd[2] = rou / 2.0 / GAMA * (2 * v * (GAMA - 1) * lambda1d + (v - a) * lambda3d + (v + a) * lambda4d);
  Gd[3] =
      rou / 2.0 / GAMA * (2 * ((GAMA - 1) * H - a * a) * lambda1d + (H - a * v) * lambda3d + (H + a * v) * lambda4d);
}
void MUSCL_y(double U[Nx + 7][Ny + 7][4], double U_L[Nx + 7][Ny + 7][4], double U_R[Nx + 7][Ny + 7][4],
             double G_[Nx + 7][Ny + 7][4], double Gp[Nx + 7][Ny + 7][4], double Gd[Nx + 7][Ny + 7][4], double dx,
             double dy, double dt) {
  int i, j, k;
  double r = dt / dy;
  //计算U_L和U_R
  for (i = 2; i <= Nx + 3; i++)
    for (j = 2; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) {
          U_L[i][j][k] = U[i][j][k] +
                         0.25 * (1 - K) * minmod(U[i][j][k] - U[i][j - 1][k], BETA * (U[i][j + 1][k] - U[i][j][k])) +
                         0.25 * (1 + K) * minmod(U[i][j + 1][k] - U[i][j][k], BETA * (U[i][j][k] - U[i][j - 1][k]));
          U_R[i][j][k] =
              U[i][j + 1][k] -
              0.25 * (1 - K) * minmod(U[i][j + 2][k] - U[i][j + 1][k], BETA * (U[i][j + 1][k] - U[i][j][k])) -
              0.25 * (1 + K) * minmod(U[i][j + 1][k] - U[i][j][k], BETA * (U[i][j + 2][k] - U[i][j + 1][k]));
        }
  //计算Gp和Gd
  for (i = 2; i <= Nx + 3; i++)
    for (j = 2; j <= Ny + 3; j++)
      if (U[i][j][0] != 0) U2GpGd_AVG(U_L[i][j], U_R[i][j], Gp[i][j], Gd[i][j]);
  //计算G_
  for (i = 2; i <= Nx + 3; i++)
    for (j = 2; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) G_[i][j][k] = Gp[i][j][k] + Gd[i][j][k];
  //计算U
  for (i = 3; i <= Nx + 3; i++)
    for (j = 3; j <= int(0.5 / dy) + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) U[i][j][k] = U[i][j][k] - r * (G_[i][j][k] - G_[i][j - 1][k]);

  for (i = int(1.0 / dx) + 3; i <= int(2.0 / dx) + 3; i++)
    for (j = int(0.5 / dy) + 4; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) U[i][j][k] = U[i][j][k] - r * (G_[i][j][k] - G_[i][j - 1][k]);
}

#else

//作用:求Fp和Fd
void U2FpFd_AVG(double U_L[4], double U_R[4], double F_[4]) {
  double H, rou, u, v, p, a, lambda1, lambda2, lambda3, lambda4, lambda1p, lambda2p, lambda3p, lambda4p, lambda1d,
      lambda2d, lambda3d, lambda4d;
  double f0,f1,f2,f3;
  //先算Fp即F+
  rou = U_L[0];
  u = U_L[1] / U_L[0];
  v = U_L[2] / U_L[0];
  p = (GAMA - 1) * (U_L[3] - 0.5 * rou * (u * u + v * v));
  a = sqrt(GAMA * p / rou);
  H = a * a / (GAMA - 1) + 0.5 * (u * u + v * v);
  lambda1 = u;
  //lambda2 = u;
  lambda3 = u - a;
  lambda4 = u + a;
  //计算正特征值
  /*
  lambda1p = 0.5 * (fabs(lambda1) + lambda1);
  //lambda2p = 0.5 * (fabs(lambda2) + lambda2);
  lambda3p = 0.5 * (fabs(lambda3) + lambda3);
  lambda4p = 0.5 * (fabs(lambda4) + lambda4);
  //特征值修正
  lambda1p = 0.5 * (lambda1p + sqrt(lambda1p * lambda1p + 1e-16));
  //lambda2p = 0.5 * (lambda2p + sqrt(lambda2p * lambda2p + 1e-16));
  lambda3p = 0.5 * (lambda3p + sqrt(lambda3p * lambda3p + 1e-16));
  lambda4p = 0.5 * (lambda4p + sqrt(lambda4p * lambda4p + 1e-16));
  */
  //lambda1p = 0.5 * (0.5 * (fabs(lambda1) + lambda1) + sqrt(0.25 * (fabs(lambda1) + lambda1) * (fabs(lambda1) + lambda1) + 1e-16));
  lambda1p = 0.5 * (0.5 * (fabs(lambda1) + lambda1) + sqrt(0.5 * lambda1 * (fabs(lambda1) + lambda1) + 1e-16));
  lambda3p = 0.5 * (0.5 * (fabs(lambda3) + lambda3) + sqrt(0.5 * lambda3 * (fabs(lambda3) + lambda3) + 1e-16));
  lambda4p = 0.5 * (0.5 * (fabs(lambda4) + lambda4) + sqrt(0.5 * lambda4 * (fabs(lambda4) + lambda4) + 1e-16));

  f0 = rou / 2.0 / GAMA * (2 * (GAMA - 1) * lambda1p + lambda3p + lambda4p);
  f1 = rou / 2.0 / GAMA * (2 * u * (GAMA - 1) * lambda1p + (u - a) * lambda3p + (u + a) * lambda4p);
  f2 = rou / 2.0 / GAMA * (v * 2 * (GAMA - 1) * lambda1p + v * lambda3p + v * lambda4p);
  f3 =
      rou / 2.0 / GAMA * (2 * ((GAMA - 1) * H - a * a) * lambda1p + (H - a * u) * lambda3p + (H + a * u) * lambda4p);
  //现在就算Fd即F-
  rou = U_R[0];
  u = U_R[1] / U_R[0];
  v = U_R[2] / U_R[0];
  p = (GAMA - 1) * (U_R[3] - 0.5 * rou * (u * u + v * v));
  a = sqrt(GAMA * p / rou);
  H = a * a / (GAMA - 1) + 0.5 * (u * u + v * v);
  lambda1 = u;
  //lambda2 = u;
  lambda3 = u - a;
  lambda4 = u + a;
  /*
  lambda1d = -0.5 * (fabs(lambda1) - lambda1);
  //lambda2d = -0.5 * (fabs(lambda2) - lambda2);
  lambda3d = -0.5 * (fabs(lambda3) - lambda3);
  lambda4d = -0.5 * (fabs(lambda4) - lambda4);
  lambda1d = 0.5 * (lambda1d - sqrt(lambda1d * lambda1d + 1e-16));
  //lambda2d = 0.5 * (lambda2d - sqrt(lambda2d * lambda2d + 1e-16));
  lambda3d = 0.5 * (lambda3d - sqrt(lambda3d * lambda3d + 1e-16));
  lambda4d = 0.5 * (lambda4d - sqrt(lambda4d * lambda4d + 1e-16));
  */
  lambda1d = 0.5 * (-0.5 * (fabs(lambda1) - lambda1) - sqrt(-0.5 * lambda1 * (fabs(lambda1) - lambda1) + 1e-16));
  lambda3d = 0.5 * (-0.5 * (fabs(lambda3) - lambda3) - sqrt(-0.5 * lambda3 * (fabs(lambda3) - lambda3) + 1e-16));
  lambda4d = 0.5 * (-0.5 * (fabs(lambda4) - lambda4) - sqrt(-0.5 * lambda4 * (fabs(lambda4) - lambda4) + 1e-16));

  F_[0] = f0 + rou / 2.0 / GAMA * (2 * (GAMA - 1) * lambda1d + lambda3d + lambda4d);
  F_[1] = f1 + rou / 2.0 / GAMA * (2 * u * (GAMA - 1) * lambda1d + (u - a) * lambda3d + (u + a) * lambda4d);
  F_[2] = f2 + rou / 2.0 / GAMA * (v * 2 * (GAMA - 1) * lambda1d + v * lambda3d + v * lambda4d);
  F_[3] = f3 +
      rou / 2.0 / GAMA * (2 * ((GAMA - 1) * H - a * a) * lambda1d + (H - a * u) * lambda3d + (H + a * u) * lambda4d);
}
void MUSCL_x(double U[Nx + 7][Ny + 7][4], double U_L[Nx + 7][Ny + 7][4], double U_R[Nx + 7][Ny + 7][4],
             double F_[Nx + 7][Ny + 7][4], double Fp[Nx + 7][Ny + 7][4], double Fd[Nx + 7][Ny + 7][4], double dx,
             double dy, double dt) {
  int i, j, k;
  double r = dt / dx;

  for (i = 2; i <= Nx + 3; i++)
    for (j = 2; j <= Ny + 3; j++)
      if (U[i][j][0] != 0) {
        for (k = 0; k <= 3; k++) {
          //计算U_L和U_R
          U_L[i][j][k] = U[i][j][k] +
                         0.25 * (1 - K) * minmod(U[i][j][k] - U[i - 1][j][k], BETA * (U[i + 1][j][k] - U[i][j][k])) +
                         0.25 * (1 + K) * minmod(U[i + 1][j][k] - U[i][j][k], BETA * (U[i][j][k] - U[i - 1][j][k]));
          U_R[i][j][k] =
              U[i + 1][j][k] -
              0.25 * (1 - K) * minmod(U[i + 2][j][k] - U[i + 1][j][k], BETA * (U[i + 1][j][k] - U[i][j][k])) -
              0.25 * (1 + K) * minmod(U[i + 1][j][k] - U[i][j][k], BETA * (U[i + 2][j][k] - U[i + 1][j][k]));
        }

        //计算Fp和Fd
        //计算F_
        //U2FpFd_AVG(U_L[i][j], U_R[i][j], Fp[i][j], Fd[i][j]);
        U2FpFd_AVG(U_L[i][j], U_R[i][j], F_[i][j]);
      }

  //计算U
  for (i = 3; i <= Nx + 3; i++)
    for (j = 3; j <= int(0.5 / dy) + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) U[i][j][k] = U[i][j][k] - r * (F_[i][j][k] - F_[i - 1][j][k]);

  for (i = int(1.0 / dx) + 3; i <= int(2.0 / dx) + 3; i++)
    for (j = int(0.5 / dy) + 4; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) U[i][j][k] = U[i][j][k] - r * (F_[i][j][k] - F_[i - 1][j][k]);
}
//算Gp和Gd
void U2GpGd_AVG(double U_L[4], double U_R[4], double G_[4]) {
  double H, rou, u, v, p, a, lambda1, lambda2, lambda3, lambda4, lambda1p, lambda2p, lambda3p, lambda4p;
  double lambda1d, lambda2d, lambda3d, lambda4d;
  double g0,g1,g2,g3;
  //计算Gp即G+的函数
  rou = U_L[0];
  u = U_L[1] / U_L[0];
  v = U_L[2] / U_L[0];
  p = (GAMA - 1) * (U_L[3] - 0.5 * rou * (u * u + v * v));
  a = sqrt(GAMA * p / rou);
  H = a * a / (GAMA - 1) + 0.5 * (u * u + v * v);
  lambda1 = v;
  //lambda2 = v;
  lambda3 = v - a;
  lambda4 = v + a;
  /*
  lambda1p = 0.5 * (fabs(lambda1) + lambda1);
  //lambda2p = 0.5 * (fabs(lambda2) + lambda2);
  lambda3p = 0.5 * (fabs(lambda3) + lambda3);
  lambda4p = 0.5 * (fabs(lambda4) + lambda4);
  lambda1p = 0.5 * (lambda1p + sqrt(lambda1p * lambda1p + 1e-16));
  //lambda2p = 0.5 * (lambda2p + sqrt(lambda2p * lambda2p + 1e-16));
  lambda3p = 0.5 * (lambda3p + sqrt(lambda3p * lambda3p + 1e-16));
  lambda4p = 0.5 * (lambda4p + sqrt(lambda4p * lambda4p + 1e-16));
  */
  lambda1p = 0.5 * (0.5 * (fabs(lambda1) + lambda1) + sqrt(0.5 * lambda1 * (fabs(lambda1) + lambda1) + 1e-16));
  lambda3p = 0.5 * (0.5 * (fabs(lambda3) + lambda3) + sqrt(0.5 * lambda3 * (fabs(lambda3) + lambda3) + 1e-16));
  lambda4p = 0.5 * (0.5 * (fabs(lambda4) + lambda4) + sqrt(0.5 * lambda4 * (fabs(lambda4) + lambda4) + 1e-16));


  g0 = rou / 2.0 / GAMA * (2 * (GAMA - 1) * lambda1p + lambda3p + lambda4p);
  g1 = rou / 2.0 / GAMA * (u * 2 * (GAMA - 1) * lambda1p + u * lambda3p + u * lambda4p);
  g2 = rou / 2.0 / GAMA * (2 * v * (GAMA - 1) * lambda1p + (v - a) * lambda3p + (v + a) * lambda4p);
  g3 =
      rou / 2.0 / GAMA * (2 * ((GAMA - 1) * H - a * a) * lambda1p + (H - a * v) * lambda3p + (H + a * v) * lambda4p);
  //计算Gd即G-的函数
  rou = U_R[0];
  u = U_R[1] / U_R[0];
  v = U_R[2] / U_R[0];
  p = (GAMA - 1) * (U_R[3] - 0.5 * rou * (u * u + v * v));
  a = sqrt(GAMA * p / rou);
  H = a * a / (GAMA - 1) + 0.5 * (u * u + v * v);
  lambda1 = v;
  //lambda2 = v;
  lambda3 = v - a;
  lambda4 = v + a;
  /*
  lambda1d = -0.5 * (fabs(lambda1) - lambda1);
  //lambda2d = -0.5 * (fabs(lambda2) - lambda2);
  lambda3d = -0.5 * (fabs(lambda3) - lambda3);
  lambda4d = -0.5 * (fabs(lambda4) - lambda4);
  lambda1d = 0.5 * (lambda1d - sqrt(lambda1d * lambda1d + 1e-16));
  //lambda2d = 0.5 * (lambda2d - sqrt(lambda2d * lambda2d + 1e-16));
  lambda3d = 0.5 * (lambda3d - sqrt(lambda3d * lambda3d + 1e-16));
  lambda4d = 0.5 * (lambda4d - sqrt(lambda4d * lambda4d + 1e-16));
  */
  lambda1d = 0.5 * (-0.5 * (fabs(lambda1) - lambda1) - sqrt(-0.5 * lambda1 * (fabs(lambda1) - lambda1) + 1e-16));
  lambda3d = 0.5 * (-0.5 * (fabs(lambda3) - lambda3) - sqrt(-0.5 * lambda3 * (fabs(lambda3) - lambda3) + 1e-16));
  lambda4d = 0.5 * (-0.5 * (fabs(lambda4) - lambda4) - sqrt(-0.5 * lambda4 * (fabs(lambda4) - lambda4) + 1e-16));


  G_[0] = g0 + rou / 2.0 / GAMA * (2 * (GAMA - 1) * lambda1d + lambda3d + lambda4d);
  G_[1] = g1 + rou / 2.0 / GAMA * (u * 2 * (GAMA - 1) * lambda1d + u * lambda3d + u * lambda4d);
  G_[2] = g2 + rou / 2.0 / GAMA * (2 * v * (GAMA - 1) * lambda1d + (v - a) * lambda3d + (v + a) * lambda4d);
  G_[3] = g3 +
      rou / 2.0 / GAMA * (2 * ((GAMA - 1) * H - a * a) * lambda1d + (H - a * v) * lambda3d + (H + a * v) * lambda4d);
}
void MUSCL_y(double U[Nx + 7][Ny + 7][4], double U_L[Nx + 7][Ny + 7][4], double U_R[Nx + 7][Ny + 7][4],
             double G_[Nx + 7][Ny + 7][4], double Gp[Nx + 7][Ny + 7][4], double Gd[Nx + 7][Ny + 7][4], double dx,
             double dy, double dt) {
  int i, j, k;
  double r = dt / dy;
  for (i = 2; i <= Nx + 3; i++)
    for (j = 2; j <= Ny + 3; j++)
      if (U[i][j][0] != 0) {
        for (k = 0; k <= 3; k++) {
          //计算U_L和U_R
          U_L[i][j][k] = U[i][j][k] +
                         0.25 * (1 - K) * minmod(U[i][j][k] - U[i][j - 1][k], BETA * (U[i][j + 1][k] - U[i][j][k])) +
                         0.25 * (1 + K) * minmod(U[i][j + 1][k] - U[i][j][k], BETA * (U[i][j][k] - U[i][j - 1][k]));
          U_R[i][j][k] =
              U[i][j + 1][k] -
              0.25 * (1 - K) * minmod(U[i][j + 2][k] - U[i][j + 1][k], BETA * (U[i][j + 1][k] - U[i][j][k])) -
              0.25 * (1 + K) * minmod(U[i][j + 1][k] - U[i][j][k], BETA * (U[i][j + 2][k] - U[i][j + 1][k]));
        }

        //计算Gp和Gd
        //计算G_
        U2GpGd_AVG(U_L[i][j], U_R[i][j], G_[i][j]);
      }

  //计算U
  for (i = 3; i <= Nx + 3; i++)
    for (j = 3; j <= int(0.5 / dy) + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) U[i][j][k] = U[i][j][k] - r * (G_[i][j][k] - G_[i][j - 1][k]);

  for (i = int(1.0 / dx) + 3; i <= int(2.0 / dx) + 3; i++)
    for (j = int(0.5 / dy) + 4; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) U[i][j][k] = U[i][j][k] - r * (G_[i][j][k] - G_[i][j - 1][k]);
}


#endif //MUSCL_FUSION

#endif