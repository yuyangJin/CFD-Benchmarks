#ifndef _NND_H
#define _NND_H
#include "util.h"

#ifndef NND_FUSION
void U2FpFd(double U[4], double Fp[4], double Fd[4]) {
  double H, rou, u, v, p, a, lambda1, lambda2, lambda3, lambda4, lambda1p, lambda2p, lambda3p, lambda4p, lambda1d,
      lambda2d, lambda3d, lambda4d;
  //先算Fp即F+
  rou = U[0];
  u = U[1] / U[0];
  v = U[2] / U[0];
  p = (GAMA - 1) * (U[3] - 0.5 * rou * (u * u + v * v));
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
  Fp[0] = rou / 2 / GAMA * (2 * (GAMA - 1) * lambda1p + lambda3p + lambda4p);
  Fp[1] = rou / 2 / GAMA * (2 * u * (GAMA - 1) * lambda1p + (u - a) * lambda3p + (u + a) * lambda4p);
  Fp[2] = rou / 2 / GAMA * (v * 2 * (GAMA - 1) * lambda1p + v * lambda3p + v * lambda4p);
  Fp[3] = rou / 2 / GAMA * (2 * ((GAMA - 1) * H - a * a) * lambda1p + (H - a * u) * lambda3p + (H + a * u) * lambda4p);
  //现在就算Fd即F-
  lambda1d = -0.5 * (fabs(lambda1) - lambda1);
  lambda2d = -0.5 * (fabs(lambda2) - lambda2);
  lambda3d = -0.5 * (fabs(lambda3) - lambda3);
  lambda4d = -0.5 * (fabs(lambda4) - lambda4);
  lambda1d = 0.5 * (lambda1d - sqrt(lambda1d * lambda1d + 1e-16));
  lambda2d = 0.5 * (lambda2d - sqrt(lambda2d * lambda2d + 1e-16));
  lambda3d = 0.5 * (lambda3d - sqrt(lambda3d * lambda3d + 1e-16));
  lambda4d = 0.5 * (lambda4d - sqrt(lambda4d * lambda4d + 1e-16));
  Fd[0] = rou / 2 / GAMA * (2 * (GAMA - 1) * lambda1d + lambda3d + lambda4d);
  Fd[1] = rou / 2 / GAMA * (2 * u * (GAMA - 1) * lambda1d + (u - a) * lambda3d + (u + a) * lambda4d);
  Fd[2] = rou / 2 / GAMA * (v * 2 * (GAMA - 1) * lambda1d + v * lambda3d + v * lambda4d);
  Fd[3] = rou / 2 / GAMA * (2 * ((GAMA - 1) * H - a * a) * lambda1d + (H - a * u) * lambda3d + (H + a * u) * lambda4d);
}
void U2GpGd(double U[4], double Gp[4], double Gd[4]) {
  double H, rou, u, v, p, a, lambda1, lambda2, lambda3, lambda4, lambda1p, lambda2p, lambda3p, lambda4p, lambda1d,
      lambda2d, lambda3d, lambda4d;
  //计算Gp即G+的函数
  rou = U[0];
  u = U[1] / U[0];
  v = U[2] / U[0];
  p = (GAMA - 1) * (U[3] - 0.5 * rou * (u * u + v * v));
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
  Gp[0] = rou / 2 / GAMA * (2 * (GAMA - 1) * lambda1p + lambda3p + lambda4p);
  Gp[1] = rou / 2 / GAMA * (u * 2 * (GAMA - 1) * lambda1p + u * lambda3p + u * lambda4p);
  Gp[2] = rou / 2 / GAMA * (2 * v * (GAMA - 1) * lambda1p + (v - a) * lambda3p + (v + a) * lambda4p);
  Gp[3] = rou / 2 / GAMA * (2 * ((GAMA - 1) * H - a * a) * lambda1p + (H - a * v) * lambda3p + (H + a * v) * lambda4p);
  //计算Gd即G-的函数
  lambda1d = -0.5 * (fabs(lambda1) - lambda1);
  lambda2d = -0.5 * (fabs(lambda2) - lambda2);
  lambda3d = -0.5 * (fabs(lambda3) - lambda3);
  lambda4d = -0.5 * (fabs(lambda4) - lambda4);
  lambda1d = 0.5 * (lambda1d - sqrt(lambda1d * lambda1d + 1e-16));
  lambda2d = 0.5 * (lambda2d - sqrt(lambda2d * lambda2d + 1e-16));
  lambda3d = 0.5 * (lambda3d - sqrt(lambda3d * lambda3d + 1e-16));
  lambda4d = 0.5 * (lambda4d - sqrt(lambda4d * lambda4d + 1e-16));
  Gd[0] = rou / 2 / GAMA * (2 * (GAMA - 1) * lambda1d + lambda3d + lambda4d);
  Gd[1] = rou / 2 / GAMA * (u * 2 * (GAMA - 1) * lambda1d + u * lambda3d + u * lambda4d);
  Gd[2] = rou / 2 / GAMA * (2 * v * (GAMA - 1) * lambda1d + (v - a) * lambda3d + (v + a) * lambda4d);
  Gd[3] = rou / 2 / GAMA * (2 * ((GAMA - 1) * H - a * a) * lambda1d + (H - a * v) * lambda3d + (H + a * v) * lambda4d);
}

void NND_x(double U[Nx + 7][Ny + 7][4], double Fp[Nx + 7][Ny + 7][4], double Fd[Nx + 7][Ny + 7][4], double dx,
           double dy, double &dt) {
  int i, j, k;
  double r;
  dt = CFL(U, dx, dy, NNDCFL);
  r = dt / dx;
  for (i = 0; i <= Nx + 6; i++)
    for (j = 0; j <= Ny + 6; j++)
      if (U[i][j][0] != 0) U2FpFd(U[i][j], Fp[i][j], Fd[i][j]);
  //分区计算U
  for (i = 3; i <= Nx + 3; i++)
    for (j = 3; j <= int(0.5 / dy) + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) {
          U[i][j][k] =
              U[i][j][k] -
              r * (Fp[i][j][k] + 0.5 * minmod(Fp[i][j][k] - Fp[i - 1][j][k], Fp[i + 1][j][k] - Fp[i][j][k]) +
                   Fd[i + 1][j][k] - 0.5 * minmod(Fd[i + 1][j][k] - Fd[i][j][k], Fd[i + 2][j][k] - Fd[i + 1][j][k]) -
                   Fp[i - 1][j][k] - 0.5 * minmod(Fp[i - 1][j][k] - Fp[i - 2][j][k], Fp[i][j][k] - Fp[i - 1][j][k]) -
                   Fd[i][j][k] + 0.5 * minmod(Fd[i][j][k] - Fd[i - 1][j][k], Fd[i + 1][j][k] - Fd[i][j][k]));
        }

  for (i = int(1.0 / dx) + 3; i <= int(2.0 / dx) + 3; i++)
    for (j = int(0.5 / dy) + 4; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) {
          U[i][j][k] =
              U[i][j][k] -
              r * (Fp[i][j][k] + 0.5 * minmod(Fp[i][j][k] - Fp[i - 1][j][k], Fp[i + 1][j][k] - Fp[i][j][k]) +
                   Fd[i + 1][j][k] - 0.5 * minmod(Fd[i + 1][j][k] - Fd[i][j][k], Fd[i + 2][j][k] - Fd[i + 1][j][k]) -
                   Fp[i - 1][j][k] - 0.5 * minmod(Fp[i - 1][j][k] - Fp[i - 2][j][k], Fp[i][j][k] - Fp[i - 1][j][k]) -
                   Fd[i][j][k] + 0.5 * minmod(Fd[i][j][k] - Fd[i - 1][j][k], Fd[i + 1][j][k] - Fd[i][j][k]));
        }
}
void NND_y(double U[Nx + 7][Ny + 7][4], double Gp[Nx + 7][Ny + 7][4], double Gd[Nx + 7][Ny + 7][4], double dx,
           double dy, double &dt) {
  int i, j, k;
  double r;
  dt = CFL(U, dx, dy, NNDCFL);
  r = dt / dy;
  for (i = 0; i <= Nx + 6; i++)
    for (j = 0; j <= Ny + 6; j++)
      if (U[i][j][0] != 0) U2GpGd(U[i][j], Gp[i][j], Gd[i][j]);
  //计算U
  for (i = 3; i <= Nx + 3; i++)
    for (j = 3; j <= int(0.5 / dy) + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) {
          U[i][j][k] =
              U[i][j][k] -
              r * (Gp[i][j][k] + 0.5 * minmod(Gp[i][j][k] - Gp[i][j - 1][k], Gp[i][j + 1][k] - Gp[i][j][k]) +
                   Gd[i][j + 1][k] - 0.5 * minmod(Gd[i][j + 1][k] - Gd[i][j][k], Gd[i][j + 2][k] - Gd[i][j + 1][k]) -
                   Gp[i][j - 1][k] - 0.5 * minmod(Gp[i][j - 1][k] - Gp[i][j - 2][k], Gp[i][j][k] - Gp[i][j - 1][k]) -
                   Gd[i][j][k] + 0.5 * minmod(Gd[i][j][k] - Gd[i][j - 1][k], Gd[i][j + 1][k] - Gd[i][j][k]));
        }

  for (i = int(1.0 / dx) + 3; i <= int(2.0 / dx) + 3; i++)
    for (j = int(0.5 / dy) + 4; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) {
          U[i][j][k] =
              U[i][j][k] -
              r * (Gp[i][j][k] + 0.5 * minmod(Gp[i][j][k] - Gp[i][j - 1][k], Gp[i][j + 1][k] - Gp[i][j][k]) +
                   Gd[i][j + 1][k] - 0.5 * minmod(Gd[i][j + 1][k] - Gd[i][j][k], Gd[i][j + 2][k] - Gd[i][j + 1][k]) -
                   Gp[i][j - 1][k] - 0.5 * minmod(Gp[i][j - 1][k] - Gp[i][j - 2][k], Gp[i][j][k] - Gp[i][j - 1][k]) -
                   Gd[i][j][k] + 0.5 * minmod(Gd[i][j][k] - Gd[i][j - 1][k], Gd[i][j + 1][k] - Gd[i][j][k]));
        }
}

#else

void U2FpFd(double U[4], double Fp[4], double Fd[4]) {
  double H, rou, u, v, p, a, lambda1, lambda2, lambda3, lambda4, lambda1p, lambda2p, lambda3p, lambda4p, lambda1d,
      lambda2d, lambda3d, lambda4d;
  //先算Fp即F+
  rou = U[0];
  u = U[1] / U[0];
  v = U[2] / U[0];
  p = (GAMA - 1) * (U[3] - 0.5 * rou * (u * u + v * v));
  a = sqrt(GAMA * p / rou);
  H = a * a / (GAMA - 1) + 0.5 * (u * u + v * v);
  lambda1 = u;
  //lambda2 = u;
  lambda3 = u - a;
  lambda4 = u + a;
  //计算正特征值
  lambda1p = 0.5 * (fabs(lambda1) + lambda1);
  //lambda2p = 0.5 * (fabs(lambda2) + lambda2);
  lambda3p = 0.5 * (fabs(lambda3) + lambda3);
  lambda4p = 0.5 * (fabs(lambda4) + lambda4);
  //特征值修正
  lambda1p = 0.5 * (lambda1p + sqrt(lambda1p * lambda1p + 1e-16));
  //lambda2p = 0.5 * (lambda2p + sqrt(lambda2p * lambda2p + 1e-16));
  lambda3p = 0.5 * (lambda3p + sqrt(lambda3p * lambda3p + 1e-16));
  lambda4p = 0.5 * (lambda4p + sqrt(lambda4p * lambda4p + 1e-16));
  Fp[0] = rou / 2 / GAMA * (2 * (GAMA - 1) * lambda1p + lambda3p + lambda4p);
  Fp[1] = rou / 2 / GAMA * (2 * u * (GAMA - 1) * lambda1p + (u - a) * lambda3p + (u + a) * lambda4p);
  Fp[2] = rou / 2 / GAMA * (v * 2 * (GAMA - 1) * lambda1p + v * lambda3p + v * lambda4p);
  Fp[3] = rou / 2 / GAMA * (2 * ((GAMA - 1) * H - a * a) * lambda1p + (H - a * u) * lambda3p + (H + a * u) * lambda4p);
  //现在就算Fd即F-
  lambda1d = -0.5 * (fabs(lambda1) - lambda1);
  //lambda2d = -0.5 * (fabs(lambda2) - lambda2);
  lambda3d = -0.5 * (fabs(lambda3) - lambda3);
  lambda4d = -0.5 * (fabs(lambda4) - lambda4);
  lambda1d = 0.5 * (lambda1d - sqrt(lambda1d * lambda1d + 1e-16));
  //lambda2d = 0.5 * (lambda2d - sqrt(lambda2d * lambda2d + 1e-16));
  lambda3d = 0.5 * (lambda3d - sqrt(lambda3d * lambda3d + 1e-16));
  lambda4d = 0.5 * (lambda4d - sqrt(lambda4d * lambda4d + 1e-16));
  Fd[0] = rou / 2 / GAMA * (2 * (GAMA - 1) * lambda1d + lambda3d + lambda4d);
  Fd[1] = rou / 2 / GAMA * (2 * u * (GAMA - 1) * lambda1d + (u - a) * lambda3d + (u + a) * lambda4d);
  Fd[2] = rou / 2 / GAMA * (v * 2 * (GAMA - 1) * lambda1d + v * lambda3d + v * lambda4d);
  Fd[3] = rou / 2 / GAMA * (2 * ((GAMA - 1) * H - a * a) * lambda1d + (H - a * u) * lambda3d + (H + a * u) * lambda4d);
}
void U2GpGd(double U[4], double Gp[4], double Gd[4]) {
  double H, rou, u, v, p, a, lambda1, lambda2, lambda3, lambda4, lambda1p, lambda2p, lambda3p, lambda4p, lambda1d,
      lambda2d, lambda3d, lambda4d;
  //计算Gp即G+的函数
  rou = U[0];
  u = U[1] / U[0];
  v = U[2] / U[0];
  p = (GAMA - 1) * (U[3] - 0.5 * rou * (u * u + v * v));
  a = sqrt(GAMA * p / rou);
  H = a * a / (GAMA - 1) + 0.5 * (u * u + v * v);
  lambda1 = v;
  //lambda2 = v;
  lambda3 = v - a;
  lambda4 = v + a;
  lambda1p = 0.5 * (fabs(lambda1) + lambda1);
  //lambda2p = 0.5 * (fabs(lambda2) + lambda2);
  lambda3p = 0.5 * (fabs(lambda3) + lambda3);
  lambda4p = 0.5 * (fabs(lambda4) + lambda4);
  lambda1p = 0.5 * (lambda1p + sqrt(lambda1p * lambda1p + 1e-16));
  //lambda2p = 0.5 * (lambda2p + sqrt(lambda2p * lambda2p + 1e-16));
  lambda3p = 0.5 * (lambda3p + sqrt(lambda3p * lambda3p + 1e-16));
  lambda4p = 0.5 * (lambda4p + sqrt(lambda4p * lambda4p + 1e-16));
  Gp[0] = rou / 2 / GAMA * (2 * (GAMA - 1) * lambda1p + lambda3p + lambda4p);
  Gp[1] = rou / 2 / GAMA * (u * 2 * (GAMA - 1) * lambda1p + u * lambda3p + u * lambda4p);
  Gp[2] = rou / 2 / GAMA * (2 * v * (GAMA - 1) * lambda1p + (v - a) * lambda3p + (v + a) * lambda4p);
  Gp[3] = rou / 2 / GAMA * (2 * ((GAMA - 1) * H - a * a) * lambda1p + (H - a * v) * lambda3p + (H + a * v) * lambda4p);
  //计算Gd即G-的函数
  lambda1d = -0.5 * (fabs(lambda1) - lambda1);
  //lambda2d = -0.5 * (fabs(lambda2) - lambda2);
  lambda3d = -0.5 * (fabs(lambda3) - lambda3);
  lambda4d = -0.5 * (fabs(lambda4) - lambda4);
  lambda1d = 0.5 * (lambda1d - sqrt(lambda1d * lambda1d + 1e-16));
  //lambda2d = 0.5 * (lambda2d - sqrt(lambda2d * lambda2d + 1e-16));
  lambda3d = 0.5 * (lambda3d - sqrt(lambda3d * lambda3d + 1e-16));
  lambda4d = 0.5 * (lambda4d - sqrt(lambda4d * lambda4d + 1e-16));
  Gd[0] = rou / 2 / GAMA * (2 * (GAMA - 1) * lambda1d + lambda3d + lambda4d);
  Gd[1] = rou / 2 / GAMA * (u * 2 * (GAMA - 1) * lambda1d + u * lambda3d + u * lambda4d);
  Gd[2] = rou / 2 / GAMA * (2 * v * (GAMA - 1) * lambda1d + (v - a) * lambda3d + (v + a) * lambda4d);
  Gd[3] = rou / 2 / GAMA * (2 * ((GAMA - 1) * H - a * a) * lambda1d + (H - a * v) * lambda3d + (H + a * v) * lambda4d);
}

void NND_x(double U[Nx + 7][Ny + 7][4], double Fp[Nx + 7][Ny + 7][4], double Fd[Nx + 7][Ny + 7][4], double dx,
           double dy, double &dt) {
  int i, j, k;
  double r;
  dt = CFL(U, dx, dy, NNDCFL);
  r = dt / dx;
  for (i = 0; i <= Nx + 6; i++)
    for (j = 0; j <= Ny + 6; j++)
      if (U[i][j][0] != 0) U2FpFd(U[i][j], Fp[i][j], Fd[i][j]);
  //分区计算U
  for (i = 3; i <= Nx + 3; i++)
    for (j = 3; j <= int(0.5 / dy) + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) {
          U[i][j][k] =
              U[i][j][k] -
              r * (Fp[i][j][k] + 0.5 * minmod(Fp[i][j][k] - Fp[i - 1][j][k], Fp[i + 1][j][k] - Fp[i][j][k]) +
                   Fd[i + 1][j][k] - 0.5 * minmod(Fd[i + 1][j][k] - Fd[i][j][k], Fd[i + 2][j][k] - Fd[i + 1][j][k]) -
                   Fp[i - 1][j][k] - 0.5 * minmod(Fp[i - 1][j][k] - Fp[i - 2][j][k], Fp[i][j][k] - Fp[i - 1][j][k]) -
                   Fd[i][j][k] + 0.5 * minmod(Fd[i][j][k] - Fd[i - 1][j][k], Fd[i + 1][j][k] - Fd[i][j][k]));
        }

  for (i = int(1.0 / dx) + 3; i <= int(2.0 / dx) + 3; i++)
    for (j = int(0.5 / dy) + 4; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) {
          U[i][j][k] =
              U[i][j][k] -
              r * (Fp[i][j][k] + 0.5 * minmod(Fp[i][j][k] - Fp[i - 1][j][k], Fp[i + 1][j][k] - Fp[i][j][k]) +
                   Fd[i + 1][j][k] - 0.5 * minmod(Fd[i + 1][j][k] - Fd[i][j][k], Fd[i + 2][j][k] - Fd[i + 1][j][k]) -
                   Fp[i - 1][j][k] - 0.5 * minmod(Fp[i - 1][j][k] - Fp[i - 2][j][k], Fp[i][j][k] - Fp[i - 1][j][k]) -
                   Fd[i][j][k] + 0.5 * minmod(Fd[i][j][k] - Fd[i - 1][j][k], Fd[i + 1][j][k] - Fd[i][j][k]));
        }
}
void NND_y(double U[Nx + 7][Ny + 7][4], double Gp[Nx + 7][Ny + 7][4], double Gd[Nx + 7][Ny + 7][4], double dx,
           double dy, double &dt) {
  int i, j, k;
  double r;
  dt = CFL(U, dx, dy, NNDCFL);
  r = dt / dy;
  for (i = 0; i <= Nx + 6; i++)
    for (j = 0; j <= Ny + 6; j++)
      if (U[i][j][0] != 0) U2GpGd(U[i][j], Gp[i][j], Gd[i][j]);
  //计算U
  for (i = 3; i <= Nx + 3; i++)
    for (j = 3; j <= int(0.5 / dy) + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) {
          U[i][j][k] =
              U[i][j][k] -
              r * (Gp[i][j][k] + 0.5 * minmod(Gp[i][j][k] - Gp[i][j - 1][k], Gp[i][j + 1][k] - Gp[i][j][k]) +
                   Gd[i][j + 1][k] - 0.5 * minmod(Gd[i][j + 1][k] - Gd[i][j][k], Gd[i][j + 2][k] - Gd[i][j + 1][k]) -
                   Gp[i][j - 1][k] - 0.5 * minmod(Gp[i][j - 1][k] - Gp[i][j - 2][k], Gp[i][j][k] - Gp[i][j - 1][k]) -
                   Gd[i][j][k] + 0.5 * minmod(Gd[i][j][k] - Gd[i][j - 1][k], Gd[i][j + 1][k] - Gd[i][j][k]));
        }

  for (i = int(1.0 / dx) + 3; i <= int(2.0 / dx) + 3; i++)
    for (j = int(0.5 / dy) + 4; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) {
          U[i][j][k] =
              U[i][j][k] -
              r * (Gp[i][j][k] + 0.5 * minmod(Gp[i][j][k] - Gp[i][j - 1][k], Gp[i][j + 1][k] - Gp[i][j][k]) +
                   Gd[i][j + 1][k] - 0.5 * minmod(Gd[i][j + 1][k] - Gd[i][j][k], Gd[i][j + 2][k] - Gd[i][j + 1][k]) -
                   Gp[i][j - 1][k] - 0.5 * minmod(Gp[i][j - 1][k] - Gp[i][j - 2][k], Gp[i][j][k] - Gp[i][j - 1][k]) -
                   Gd[i][j][k] + 0.5 * minmod(Gd[i][j][k] - Gd[i][j - 1][k], Gd[i][j + 1][k] - Gd[i][j][k]));
        }
}


#endif // NND_FUSION

#endif