#ifndef __COMPACT_H
#define __COMPACT_h

#include "util.h"

#ifndef COMPACT_FUSION

void CPT_fpfd(double U[4], double fp[4], double fd[4]) {
  double H, rou, u, v, p, a, lambda1, lambda2, lambda3, lambda4, lambda1p, lambda2p, lambda3p, lambda4p, lambda1d,
      lambda2d, lambda3d, lambda4d;
  //先算fp即f+
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
  fp[0] = rou / 2.0 / GAMA * (2 * (GAMA - 1) * lambda1p + lambda3p + lambda4p);
  fp[1] = rou / 2.0 / GAMA * (2 * u * (GAMA - 1) * lambda1p + (u - a) * lambda3p + (u + a) * lambda4p);
  fp[2] = rou / 2.0 / GAMA * (v * 2 * (GAMA - 1) * lambda1p + v * lambda3p + v * lambda4p);
  fp[3] =
      rou / 2.0 / GAMA * (2 * ((GAMA - 1) * H - a * a) * lambda1p + (H - a * u) * lambda3p + (H + a * u) * lambda4p);
  //现在就算fd即f-
  lambda1d = -0.5 * (fabs(lambda1) - lambda1);
  lambda2d = -0.5 * (fabs(lambda2) - lambda2);
  lambda3d = -0.5 * (fabs(lambda3) - lambda3);
  lambda4d = -0.5 * (fabs(lambda4) - lambda4);
  lambda1d = 0.5 * (lambda1d - sqrt(lambda1d * lambda1d + 1e-16));
  lambda2d = 0.5 * (lambda2d - sqrt(lambda2d * lambda2d + 1e-16));
  lambda3d = 0.5 * (lambda3d - sqrt(lambda3d * lambda3d + 1e-16));
  lambda4d = 0.5 * (lambda4d - sqrt(lambda4d * lambda4d + 1e-16));
  fd[0] = rou / 2.0 / GAMA * (2 * (GAMA - 1) * lambda1d + lambda3d + lambda4d);
  fd[1] = rou / 2.0 / GAMA * (2 * u * (GAMA - 1) * lambda1d + (u - a) * lambda3d + (u + a) * lambda4d);
  fd[2] = rou / 2.0 / GAMA * (v * 2 * (GAMA - 1) * lambda1d + v * lambda3d + v * lambda4d);
  fd[3] =
      rou / 2.0 / GAMA * (2 * ((GAMA - 1) * H - a * a) * lambda1d + (H - a * u) * lambda3d + (H + a * u) * lambda4d);
}

void Compact_x(double U[Nx + 7][Ny + 7][4], double Ut[Nx + 7][Ny + 7][4], double Fp[Nx + 7][Ny + 7][4],
               double Fd[Nx + 7][Ny + 7][4], double fp[Nx + 7][Ny + 7][4], double fd[Nx + 7][Ny + 7][4], double dx,
               double dy, double dt, double p[Nx + 7][Ny + 7]) {
  int i, j, k;
  double r = dt / dx;
  double q, rou, u, v;
  double nu;
  for (i = 0; i <= Nx + 6; i++)
    for (j = 0; j <= Ny + 6; j++)
      if (U[i][j][0] != 0) {
        rou = U[i][j][0];
        u = U[i][j][1] / U[i][j][0];
        v = U[i][j][2] / U[i][j][0];
        p[i][j] = (GAMA - 1) * (U[i][j][3] - 0.5 * rou * (u * u + v * v));
      }

  nu = 5 * r * (1 - 5 * r);
  //人工黏性
  for (i = 3; i <= Nx + 3; i++)
    for (j = 3; j <= int(0.5 / dy) + 3; j++) {
      q = fabs(fabs(U[i + 1][j][0] - U[i][j][0]) - fabs(U[i][j][0] - U[i - 1][j][0])) /
          (fabs(U[i + 1][j][0] - U[i][j][0]) + fabs(U[i][j][0] - U[i - 1][j][0]) + 1e-100);
      for (k = 0; k <= 3; k++) {
        Ut[i][j][k] = U[i][j][k] + 0.5 * nu * q * (U[i + 1][j][k] - 2 * U[i][j][k] + U[i - 1][j][k]);
      }
    }

  for (i = int(1.0 / dx) + 3; i <= int(2.0 / dx) + 3; i++)
    for (j = int(0.5 / dy) + 4; j <= Ny + 3; j++) {
      q = fabs(fabs(U[i + 1][j][0] - U[i][j][0]) - fabs(U[i][j][0] - U[i - 1][j][0])) /
          (fabs(U[i + 1][j][0] - U[i][j][0]) + fabs(U[i][j][0] - U[i - 1][j][0]) + 1e-100);
      for (k = 0; k <= 3; k++) {
        Ut[i][j][k] = U[i][j][k] + 0.5 * nu * q * (U[i + 1][j][k] - 2 * U[i][j][k] + U[i - 1][j][k]);
      }
    }

  for (i = 3; i <= Nx + 3; i++)
    for (j = 3; j <= int(0.5 / dy) + 3; j++) {
      for (k = 0; k <= 3; k++) {
        U[i][j][k] = Ut[i][j][k];
      }
    }

  for (i = int(1.0 / dx) + 3; i <= int(2.0 / dx) + 3; i++)
    for (j = int(0.5 / dy) + 4; j <= Ny + 3; j++)
      for (k = 0; k <= 3; k++) U[i][j][k] = Ut[i][j][k];

  //计算fp和fd
  for (i = 0; i <= Nx + 6; i++)
    for (j = 0; j <= Ny + 6; j++)
      if (U[i][j][0] != 0) CPT_fpfd(U[i][j], fp[i][j], fd[i][j]);

  //分区计算F+的左值
  for (j = 3; j <= int(0.5 / dy) + 3; j++)
    for (k = 0; k <= 3; k++) Fp[2][j][k] = 0.5 * (3 * (fp[3][j][k] - fp[2][j][k]) - (fp[4][j][k] - fp[3][j][k]));

  for (j = int(0.5 / dy) + 4; j <= Ny + 3; j++)
    for (k = 0; k <= 3; k++)
      Fp[int(1.0 / dx) + 2][j][k] = 0.5 * (3 * (fp[int(1.0 / dx) + 3][j][k] - fp[int(1.0 / dx) + 2][j][k]) -
                                           (fp[int(1.0 / dx) + 4][j][k] - fp[int(1.0 / dx) + 3][j][k]));

  //分区计算全场的F+
  for (i = 3; i <= Nx + 3; i++)
    for (j = 3; j <= int(0.5 / dy) + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++)
          Fp[i][j][k] = -0.5 * Fp[i - 1][j][k] + 5.0 / 4.0 * (fp[i][j][k] - fp[i - 1][j][k]) +
                        0.25 * (fp[i + 1][j][k] - fp[i][j][k]);

  for (i = int(1.0 / dx) + 3; i <= int(2.0 / dx) + 3; i++)
    for (j = int(0.5 / dy) + 4; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++)
          Fp[i][j][k] = -0.5 * Fp[i - 1][j][k] + 5.0 / 4.0 * (fp[i][j][k] - fp[i - 1][j][k]) +
                        0.25 * (fp[i + 1][j][k] - fp[i][j][k]);

  //分区计算F-的右值
  for (j = 3; j <= int(0.5 / dy) + 3; j++)
    for (k = 0; k <= 3; k++)
      Fd[Nx + 4][j][k] = 0.5 * (3 * (fd[Nx + 4][j][k] - fd[Nx + 3][j][k]) - (fd[Nx + 3][j][k] - fd[Nx + 2][j][k]));

  for (j = int(0.5 / dy) + 4; j <= Ny + 3; j++)
    for (k = 0; k <= 3; k++)
      Fd[int(2.0 / dx) + 4][j][k] = 0.5 * (3 * (fd[int(2.0 / dx) + 4][j][k] - fd[int(2.0 / dx) + 3][j][k]) -
                                           (fd[int(2.0 / dx) + 3][j][k] - fd[int(2.0 / dx) + 2][j][k]));

  //分区计算全场的F-值
  for (i = Nx + 3; i >= 3; i--)
    for (j = 3; j <= int(0.5 / dy) + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++)
          Fd[i][j][k] = -0.5 * Fd[i + 1][j][k] + 5.0 / 4.0 * (fd[i + 1][j][k] - fd[i][j][k]) +
                        0.25 * (fd[i][j][k] - fd[i - 1][j][k]);

  for (i = int(2.0 / dx) + 3; i >= int(1.0 / dx) + 3; i--)
    for (j = int(0.5 / dy) + 4; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++)
          Fd[i][j][k] = -0.5 * Fd[i + 1][j][k] + 5.0 / 4.0 * (fd[i + 1][j][k] - fd[i][j][k]) +
                        0.25 * (fd[i][j][k] - fd[i - 1][j][k]);
}

void Compact_RK_x(double U[Nx + 7][Ny + 7][4], double Ut[Nx + 7][Ny + 7][4], double U1[Nx + 7][Ny + 7][4],
                  double U2[Nx + 7][Ny + 7][4], double Fp[Nx + 7][Ny + 7][4], double Fd[Nx + 7][Ny + 7][4],
                  double fp[Nx + 7][Ny + 7][4], double fd[Nx + 7][Ny + 7][4], double dx, double dy, double dt,
                  double p[Nx + 7][Ny + 7], double LAMDA_[Nx + 7][Ny + 7][4][4]) {
  int i, j, k;
  dt = CFL(U, dx, dy, COMPACTCFL);
  double r = dt / dx;
  //分区计算U1
  for (i = 3; i <= Nx + 3; i++)
    for (j = 3; j <= int(0.5 / dy) + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) U1[i][j][k] = U[i][j][k] - r * (Fp[i][j][k] + Fd[i][j][k]);

  for (i = int(1.0 / dx) + 3; i <= int(2.0 / dx) + 3; i++)
    for (j = int(0.5 / dy) + 4; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) U1[i][j][k] = U[i][j][k] - r * (Fp[i][j][k] + Fd[i][j][k]);

  // U1求边界条件
  bound(U1, dx, dy);
  //计算dt
  dt = CFL(U1, dx, dy, COMPACTCFL);
  r = dt / dx;
  for (i = 0; i <= Nx + 6; i++)
    for (j = 0; j <= Ny + 6; j++)
      if (U[i][j][0] != 0) CPT_fpfd(U1[i][j], fp[i][j], fd[i][j]);
  //用U1算Fpd
  Compact_x(U1, Ut, Fp, Fd, fp, fd, dx, dy, dt, p);
  // RK算U2
  for (i = 3; i <= Nx + 3; i++)
    for (j = 3; j <= int(0.5 / dy) + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) U2[i][j][k] = 0.75 * U[i][j][k] + 0.25 * U1[i][j][k] - r * (Fp[i][j][k] + Fd[i][j][k]);

  for (i = int(1.0 / dx) + 3; i <= int(2.0 / dx) + 3; i++)
    for (j = int(0.5 / dy) + 4; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) U2[i][j][k] = 0.75 * U[i][j][k] + 0.25 * U1[i][j][k] - r * (Fp[i][j][k] + Fd[i][j][k]);

  bound(U2, dx, dy);
  dt = CFL(U2, dx, dy, COMPACTCFL);
  r = dt / dx;
  for (i = 0; i <= Nx + 6; i++)
    for (j = 0; j <= Ny + 6; j++)
      if (U[i][j][0] != 0) CPT_fpfd(U2[i][j], fp[i][j], fd[i][j]);
  Compact_x(U2, Ut, Fp, Fd, fp, fd, dx, dy, dt, p);
  // RK算U3-即U
  for (i = 3; i <= Nx + 3; i++)
    for (j = 3; j <= int(0.5 / dy) + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++)
          U[i][j][k] = 1.0 / 3.0 * U[i][j][k] + 2.0 / 3.0 * U2[i][j][k] - 2.0 / 3.0 * r * (Fp[i][j][k] + Fd[i][j][k]);

  for (i = int(1.0 / dx) + 3; i <= int(2.0 / dx) + 3; i++)
    for (j = int(0.5 / dy) + 4; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++)
          U[i][j][k] = 1.0 / 3.0 * U[i][j][k] + 2.0 / 3.0 * U2[i][j][k] - 2.0 / 3.0 * r * (Fp[i][j][k] + Fd[i][j][k]);
}
///////////-------------------------------y
void CPT_gpgd(double U[4], double gp[4], double gd[4]) {
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
  gp[0] = rou / 2.0 / GAMA * (2 * (GAMA - 1) * lambda1p + lambda3p + lambda4p);
  gp[1] = rou / 2.0 / GAMA * (u * 2 * (GAMA - 1) * lambda1p + u * lambda3p + u * lambda4p);
  gp[2] = rou / 2.0 / GAMA * (2 * v * (GAMA - 1) * lambda1p + (v - a) * lambda3p + (v + a) * lambda4p);
  gp[3] =
      rou / 2.0 / GAMA * (2 * ((GAMA - 1) * H - a * a) * lambda1p + (H - a * v) * lambda3p + (H + a * v) * lambda4p);
  //计算Gd即G-的函数
  lambda1d = -0.5 * (fabs(lambda1) - lambda1);
  lambda2d = -0.5 * (fabs(lambda2) - lambda2);
  lambda3d = -0.5 * (fabs(lambda3) - lambda3);
  lambda4d = -0.5 * (fabs(lambda4) - lambda4);
  lambda1d = 0.5 * (lambda1d - sqrt(lambda1d * lambda1d + 1e-16));
  lambda2d = 0.5 * (lambda2d - sqrt(lambda2d * lambda2d + 1e-16));
  lambda3d = 0.5 * (lambda3d - sqrt(lambda3d * lambda3d + 1e-16));
  lambda4d = 0.5 * (lambda4d - sqrt(lambda4d * lambda4d + 1e-16));
  gd[0] = rou / 2.0 / GAMA * (2 * (GAMA - 1) * lambda1d + lambda3d + lambda4d);
  gd[1] = rou / 2.0 / GAMA * (u * 2 * (GAMA - 1) * lambda1d + u * lambda3d + u * lambda4d);
  gd[2] = rou / 2.0 / GAMA * (2 * v * (GAMA - 1) * lambda1d + (v - a) * lambda3d + (v + a) * lambda4d);
  gd[3] =
      rou / 2.0 / GAMA * (2 * ((GAMA - 1) * H - a * a) * lambda1d + (H - a * v) * lambda3d + (H + a * v) * lambda4d);
}

void Compact_y(double U[Nx + 7][Ny + 7][4], double Ut[Nx + 7][Ny + 7][4], double Gp[Nx + 7][Ny + 7][4],
               double Gd[Nx + 7][Ny + 7][4], double gp[Nx + 7][Ny + 7][4], double gd[Nx + 7][Ny + 7][4], double dx,
               double dy, double dt, double p[Nx + 7][Ny + 7]) {
  int i, j, k;
  double r = dt / dy;
  double q, rou, u, v;
  double nu;
  for (i = 0; i <= Nx + 6; i++)
    for (j = 0; j <= Ny + 6; j++)
      if (U[i][j][0] != 0) {
        rou = U[i][j][0];
        u = U[i][j][1] / U[i][j][0];
        v = U[i][j][2] / U[i][j][0];
        p[i][j] = (GAMA - 1) * (U[i][j][3] - 0.5 * rou * (u * u + v * v));
      }

  nu = 5 * r * (1 - 5 * r);
  //人工黏性
  for (i = 3; i <= Nx + 3; i++)
    for (j = 3; j <= int(0.5 / dy) + 3; j++) {
      q = fabs(fabs(U[i][j + 1][0] - U[i][j][0]) - fabs(U[i][j][0] - U[i][j - 1][0])) /
          (fabs(U[i][j + 1][0] - U[i][j][0]) + fabs(U[i][j][0] - U[i][j - 1][0]) + 1e-100);
      for (k = 0; k <= 3; k++) {
        Ut[i][j][k] = U[i][j][k] + 0.5 * nu * q * (U[i][j + 1][k] - 2 * U[i][j][k] + U[i][j - 1][k]);
      }
    }
  for (i = int(1.0 / dx) + 3; i <= int(2.0 / dx) + 3; i++)
    for (j = int(0.5 / dy) + 4; j <= Ny + 3; j++) {
      q = fabs(fabs(U[i][j + 1][0] - U[i][j][0]) - fabs(U[i][j][0] - U[i][j - 1][0])) /
          (fabs(U[i][j + 1][0] - U[i][j][0]) + fabs(U[i][j][0] - U[i][j - 1][0]) + 1e-100);
      for (k = 0; k <= 3; k++) {
        Ut[i][j][k] = U[i][j][k] + 0.5 * nu * q * (U[i][j + 1][k] - 2 * U[i][j][k] + U[i][j - 1][k]);
      }
    }

  for (i = 3; i <= Nx + 3; i++)
    for (j = 3; j <= int(0.5 / dy) + 3; j++)
      for (k = 0; k <= 3; k++) U[i][j][k] = Ut[i][j][k];

  for (i = int(1.0 / dx) + 3; i <= int(2.0 / dx) + 3; i++)
    for (j = int(0.5 / dy) + 4; j <= Ny + 3; j++)
      for (k = 0; k <= 3; k++) U[i][j][k] = Ut[i][j][k];

  //计算gp和gd
  for (i = 0; i <= Nx + 6; i++)
    for (j = 0; j <= Ny + 6; j++)
      if (U[i][j][0] != 0) CPT_gpgd(U[i][j], gp[i][j], gd[i][j]);

  //分区计算G+的下值
  for (i = 3; i <= Nx + 3; i++)
    for (k = 0; k <= 3; k++) Gp[i][2][k] = 0.5 * (3 * (gp[i][3][k] - gp[i][2][k]) - (gp[i][4][k] - gp[i][3][k]));

  //分区计算全场的G+
  for (j = 3; j <= Ny + 3; j++)
    for (i = 3; i <= Nx + 3; i++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++)
          Gp[i][j][k] = -0.5 * Gp[i][j - 1][k] + 5.0 / 4.0 * (gp[i][j][k] - gp[i][j - 1][k]) +
                        0.25 * (gp[i][j + 1][k] - gp[i][j][k]);

  //分区计算G-的上值
  for (i = 3; i <= int(1.0 / dx) + 2; i++)
    for (k = 0; k <= 3; k++)
      Gd[i][int(0.5 / dy) + 4][k] = 0.5 * (3 * (gd[i][int(0.5 / dy) + 4][k] - gd[i][int(0.5 / dy) + 3][k]) -
                                           (gd[i][int(0.5 / dy) + 3][k] - gd[i][int(0.5 / dy) + 2][k]));

  for (i = int(1.0 / dx) + 3; i <= int(2.0 / dx) + 3; i++)
    for (k = 0; k <= 3; k++)
      Gd[i][int(1.0 / dy) + 4][k] = 0.5 * (3 * (gd[i][int(1.0 / dy) + 4][k] - gd[i][int(1.0 / dy) + 3][k]) -
                                           (gd[i][int(1.0 / dy) + 3][k] - gd[i][int(1.0 / dy) + 2][k]));

  for (i = int(2.0 / dx) + 4; i <= Nx + 3; i++)
    for (k = 0; k <= 3; k++)
      Gd[i][int(0.5 / dy) + 4][k] = 0.5 * (3 * (gd[i][int(0.5 / dy) + 4][k] - gd[i][int(0.5 / dy) + 3][k]) -
                                           (gd[i][int(0.5 / dy) + 3][k] - gd[i][int(0.5 / dy) + 2][k]));

  //分区计算全场的F-值
  for (j = int(0.5 / dy) + 3; j >= 3; j--)
    for (i = 3; i <= int(1.0 / dx) + 2; i++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++)
          Gd[i][j][k] = -0.5 * Gd[i][j + 1][k] + 5.0 / 4.0 * (gd[i][j + 1][k] - gd[i][j][k]) +
                        0.25 * (gd[i][j][k] - gd[i][j - 1][k]);

  for (j = int(0.5 / dy) + 3; j >= 3; j--)
    for (i = int(2.0 / dx) + 4; i <= Nx + 3; i++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++)
          Gd[i][j][k] = -0.5 * Gd[i][j + 1][k] + 5.0 / 4.0 * (gd[i][j + 1][k] - gd[i][j][k]) +
                        0.25 * (gd[i][j][k] - gd[i][j - 1][k]);

  for (j = int(1.0 / dy) + 3; j >= 3; j--)
    for (i = int(1.0 / dx) + 3; i <= int(2.0 / dx) + 3; i++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++)
          Gd[i][j][k] = -0.5 * Gd[i][j + 1][k] + 5.0 / 4.0 * (gd[i][j + 1][k] - gd[i][j][k]) +
                        0.25 * (gd[i][j][k] - gd[i][j - 1][k]);
}
void Compact_RK_y(double U[Nx + 7][Ny + 7][4], double Ut[Nx + 7][Ny + 7][4], double U1[Nx + 7][Ny + 7][4],
                  double U2[Nx + 7][Ny + 7][4], double Gp[Nx + 7][Ny + 7][4], double Gd[Nx + 7][Ny + 7][4],
                  double gp[Nx + 7][Ny + 7][4], double gd[Nx + 7][Ny + 7][4], double dx, double dy, double dt,
                  double p[Nx + 7][Ny + 7], double LAMDA_[Nx + 7][Ny + 7][4][4]) {
  int i, j, k;
  dt = CFL(U, dx, dy, COMPACTCFL);
  double r = dt / dy;
  //分区计算U
  for (i = 3; i <= Nx + 3; i++)
    for (j = 3; j <= int(0.5 / dy) + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) U1[i][j][k] = U[i][j][k] - r * (Gp[i][j][k] + Gd[i][j][k]);

  for (i = int(1.0 / dx) + 3; i <= int(2.0 / dx) + 3; i++)
    for (j = int(0.5 / dy) + 4; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) U1[i][j][k] = U[i][j][k] - r * (Gp[i][j][k] + Gd[i][j][k]);

  bound(U1, dx, dy);
  dt = CFL(U1, dx, dy, COMPACTCFL);
  r = dt / dy;
  for (i = 0; i <= Nx + 6; i++)
    for (j = 0; j <= Ny + 6; j++)
      if (U[i][j][0] != 0) CPT_gpgd(U1[i][j], gp[i][j], gd[i][j]);
  Compact_y(U1, Ut, Gp, Gd, gp, gd, dx, dy, dt, p);
  for (i = 3; i <= Nx + 3; i++)
    for (j = 3; j <= int(0.5 / dy) + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) U2[i][j][k] = 0.75 * U[i][j][k] + 0.25 * U1[i][j][k] - r * (Gp[i][j][k] + Gd[i][j][k]);

  for (i = int(1.0 / dx) + 3; i <= int(2.0 / dx) + 3; i++)
    for (j = int(0.5 / dy) + 4; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) U2[i][j][k] = 0.75 * U[i][j][k] + 0.25 * U1[i][j][k] - r * (Gp[i][j][k] + Gd[i][j][k]);

  bound(U2, dx, dy);
  dt = CFL(U2, dx, dy, COMPACTCFL);
  r = dt / dy;
  for (i = 0; i <= Nx + 6; i++)
    for (j = 0; j <= Ny + 6; j++)
      if (U[i][j][0] != 0) CPT_gpgd(U2[i][j], gp[i][j], gd[i][j]);
  Compact_y(U2, Ut, Gp, Gd, gp, gd, dx, dy, dt, p);
  for (i = 3; i <= Nx + 3; i++)
    for (j = 3; j <= int(0.5 / dy) + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++)
          U[i][j][k] = 1.0 / 3.0 * U[i][j][k] + 2.0 / 3.0 * U2[i][j][k] - 2.0 / 3.0 * r * (Gp[i][j][k] + Gd[i][j][k]);

  for (i = int(1.0 / dx) + 3; i <= int(2.0 / dx) + 3; i++)
    for (j = int(0.5 / dy) + 4; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++)
          U[i][j][k] = 1.0 / 3.0 * U[i][j][k] + 2.0 / 3.0 * U2[i][j][k] - 2.0 / 3.0 * r * (Gp[i][j][k] + Gd[i][j][k]);
}

#else

void CPT_fpfd(double U[4], double fp[4], double fd[4]) {
  double H, rou, u, v, p, a, lambda1, lambda2, lambda3, lambda4, lambda1p, lambda2p, lambda3p, lambda4p, lambda1d,
      lambda2d, lambda3d, lambda4d;
  //先算fp即f+
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
  fp[0] = rou / 2.0 / GAMA * (2 * (GAMA - 1) * lambda1p + lambda3p + lambda4p);
  fp[1] = rou / 2.0 / GAMA * (2 * u * (GAMA - 1) * lambda1p + (u - a) * lambda3p + (u + a) * lambda4p);
  fp[2] = rou / 2.0 / GAMA * (v * 2 * (GAMA - 1) * lambda1p + v * lambda3p + v * lambda4p);
  fp[3] =
      rou / 2.0 / GAMA * (2 * ((GAMA - 1) * H - a * a) * lambda1p + (H - a * u) * lambda3p + (H + a * u) * lambda4p);
  //现在就算fd即f-
  lambda1d = -0.5 * (fabs(lambda1) - lambda1);
  lambda2d = -0.5 * (fabs(lambda2) - lambda2);
  lambda3d = -0.5 * (fabs(lambda3) - lambda3);
  lambda4d = -0.5 * (fabs(lambda4) - lambda4);
  lambda1d = 0.5 * (lambda1d - sqrt(lambda1d * lambda1d + 1e-16));
  lambda2d = 0.5 * (lambda2d - sqrt(lambda2d * lambda2d + 1e-16));
  lambda3d = 0.5 * (lambda3d - sqrt(lambda3d * lambda3d + 1e-16));
  lambda4d = 0.5 * (lambda4d - sqrt(lambda4d * lambda4d + 1e-16));
  fd[0] = rou / 2.0 / GAMA * (2 * (GAMA - 1) * lambda1d + lambda3d + lambda4d);
  fd[1] = rou / 2.0 / GAMA * (2 * u * (GAMA - 1) * lambda1d + (u - a) * lambda3d + (u + a) * lambda4d);
  fd[2] = rou / 2.0 / GAMA * (v * 2 * (GAMA - 1) * lambda1d + v * lambda3d + v * lambda4d);
  fd[3] =
      rou / 2.0 / GAMA * (2 * ((GAMA - 1) * H - a * a) * lambda1d + (H - a * u) * lambda3d + (H + a * u) * lambda4d);
}

void Compact_x(double U[Nx + 7][Ny + 7][4], double Ut[Nx + 7][Ny + 7][4], double Fp[Nx + 7][Ny + 7][4],
               double Fd[Nx + 7][Ny + 7][4], double fp[Nx + 7][Ny + 7][4], double fd[Nx + 7][Ny + 7][4], double dx,
               double dy, double dt, double p[Nx + 7][Ny + 7]) {
  int i, j, k;
  double r = dt / dx;
  double q, rou, u, v;
  double nu;
  for (i = 0; i <= Nx + 6; i++)
    for (j = 0; j <= Ny + 6; j++)
      if (U[i][j][0] != 0) {
        rou = U[i][j][0];
        u = U[i][j][1] / U[i][j][0];
        v = U[i][j][2] / U[i][j][0];
        p[i][j] = (GAMA - 1) * (U[i][j][3] - 0.5 * rou * (u * u + v * v));
      }

  nu = 5 * r * (1 - 5 * r);
  //人工黏性
  for (i = 3; i <= Nx + 3; i++)
    for (j = 3; j <= int(0.5 / dy) + 3; j++) {
      q = fabs(fabs(U[i + 1][j][0] - U[i][j][0]) - fabs(U[i][j][0] - U[i - 1][j][0])) /
          (fabs(U[i + 1][j][0] - U[i][j][0]) + fabs(U[i][j][0] - U[i - 1][j][0]) + 1e-100);
      for (k = 0; k <= 3; k++) {
        Ut[i][j][k] = U[i][j][k] + 0.5 * nu * q * (U[i + 1][j][k] - 2 * U[i][j][k] + U[i - 1][j][k]);
      }
    }

  for (i = int(1.0 / dx) + 3; i <= int(2.0 / dx) + 3; i++)
    for (j = int(0.5 / dy) + 4; j <= Ny + 3; j++) {
      q = fabs(fabs(U[i + 1][j][0] - U[i][j][0]) - fabs(U[i][j][0] - U[i - 1][j][0])) /
          (fabs(U[i + 1][j][0] - U[i][j][0]) + fabs(U[i][j][0] - U[i - 1][j][0]) + 1e-100);
      for (k = 0; k <= 3; k++) {
        Ut[i][j][k] = U[i][j][k] + 0.5 * nu * q * (U[i + 1][j][k] - 2 * U[i][j][k] + U[i - 1][j][k]);
      }
    }

  for (i = 3; i <= Nx + 3; i++)
    for (j = 3; j <= int(0.5 / dy) + 3; j++) {
      for (k = 0; k <= 3; k++) {
        U[i][j][k] = Ut[i][j][k];
      }
    }

  for (i = int(1.0 / dx) + 3; i <= int(2.0 / dx) + 3; i++)
    for (j = int(0.5 / dy) + 4; j <= Ny + 3; j++)
      for (k = 0; k <= 3; k++) U[i][j][k] = Ut[i][j][k];

  //计算fp和fd
  for (i = 0; i <= Nx + 6; i++)
    for (j = 0; j <= Ny + 6; j++)
      if (U[i][j][0] != 0) CPT_fpfd(U[i][j], fp[i][j], fd[i][j]);

  //分区计算F+的左值
  for (j = 3; j <= int(0.5 / dy) + 3; j++)
    for (k = 0; k <= 3; k++) Fp[2][j][k] = 0.5 * (3 * (fp[3][j][k] - fp[2][j][k]) - (fp[4][j][k] - fp[3][j][k]));

  for (j = int(0.5 / dy) + 4; j <= Ny + 3; j++)
    for (k = 0; k <= 3; k++)
      Fp[int(1.0 / dx) + 2][j][k] = 0.5 * (3 * (fp[int(1.0 / dx) + 3][j][k] - fp[int(1.0 / dx) + 2][j][k]) -
                                           (fp[int(1.0 / dx) + 4][j][k] - fp[int(1.0 / dx) + 3][j][k]));

  //分区计算全场的F+
  for (i = 3; i <= Nx + 3; i++)
    for (j = 3; j <= int(0.5 / dy) + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++)
          Fp[i][j][k] = -0.5 * Fp[i - 1][j][k] + 5.0 / 4.0 * (fp[i][j][k] - fp[i - 1][j][k]) +
                        0.25 * (fp[i + 1][j][k] - fp[i][j][k]);

  for (i = int(1.0 / dx) + 3; i <= int(2.0 / dx) + 3; i++)
    for (j = int(0.5 / dy) + 4; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++)
          Fp[i][j][k] = -0.5 * Fp[i - 1][j][k] + 5.0 / 4.0 * (fp[i][j][k] - fp[i - 1][j][k]) +
                        0.25 * (fp[i + 1][j][k] - fp[i][j][k]);

  //分区计算F-的右值
  for (j = 3; j <= int(0.5 / dy) + 3; j++)
    for (k = 0; k <= 3; k++)
      Fd[Nx + 4][j][k] = 0.5 * (3 * (fd[Nx + 4][j][k] - fd[Nx + 3][j][k]) - (fd[Nx + 3][j][k] - fd[Nx + 2][j][k]));

  for (j = int(0.5 / dy) + 4; j <= Ny + 3; j++)
    for (k = 0; k <= 3; k++)
      Fd[int(2.0 / dx) + 4][j][k] = 0.5 * (3 * (fd[int(2.0 / dx) + 4][j][k] - fd[int(2.0 / dx) + 3][j][k]) -
                                           (fd[int(2.0 / dx) + 3][j][k] - fd[int(2.0 / dx) + 2][j][k]));

  //分区计算全场的F-值
  for (i = Nx + 3; i >= 3; i--)
    for (j = 3; j <= int(0.5 / dy) + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++)
          Fd[i][j][k] = -0.5 * Fd[i + 1][j][k] + 5.0 / 4.0 * (fd[i + 1][j][k] - fd[i][j][k]) +
                        0.25 * (fd[i][j][k] - fd[i - 1][j][k]);

  for (i = int(2.0 / dx) + 3; i >= int(1.0 / dx) + 3; i--)
    for (j = int(0.5 / dy) + 4; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++)
          Fd[i][j][k] = -0.5 * Fd[i + 1][j][k] + 5.0 / 4.0 * (fd[i + 1][j][k] - fd[i][j][k]) +
                        0.25 * (fd[i][j][k] - fd[i - 1][j][k]);
}

void Compact_RK_x(double U[Nx + 7][Ny + 7][4], double Ut[Nx + 7][Ny + 7][4], double U1[Nx + 7][Ny + 7][4],
                  double U2[Nx + 7][Ny + 7][4], double Fp[Nx + 7][Ny + 7][4], double Fd[Nx + 7][Ny + 7][4],
                  double fp[Nx + 7][Ny + 7][4], double fd[Nx + 7][Ny + 7][4], double dx, double dy, double dt,
                  double p[Nx + 7][Ny + 7], double LAMDA_[Nx + 7][Ny + 7][4][4]) {
  int i, j, k;
  dt = CFL(U, dx, dy, COMPACTCFL);
  double r = dt / dx;
  //分区计算U1
  for (i = 3; i <= Nx + 3; i++)
    for (j = 3; j <= int(0.5 / dy) + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) U1[i][j][k] = U[i][j][k] - r * (Fp[i][j][k] + Fd[i][j][k]);

  for (i = int(1.0 / dx) + 3; i <= int(2.0 / dx) + 3; i++)
    for (j = int(0.5 / dy) + 4; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) U1[i][j][k] = U[i][j][k] - r * (Fp[i][j][k] + Fd[i][j][k]);

  // U1求边界条件
  bound(U1, dx, dy);
  //计算dt
  dt = CFL(U1, dx, dy, COMPACTCFL);
  r = dt / dx;
  for (i = 0; i <= Nx + 6; i++)
    for (j = 0; j <= Ny + 6; j++)
      if (U[i][j][0] != 0) CPT_fpfd(U1[i][j], fp[i][j], fd[i][j]);
  //用U1算Fpd
  Compact_x(U1, Ut, Fp, Fd, fp, fd, dx, dy, dt, p);
  // RK算U2
  for (i = 3; i <= Nx + 3; i++)
    for (j = 3; j <= int(0.5 / dy) + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) U2[i][j][k] = 0.75 * U[i][j][k] + 0.25 * U1[i][j][k] - r * (Fp[i][j][k] + Fd[i][j][k]);

  for (i = int(1.0 / dx) + 3; i <= int(2.0 / dx) + 3; i++)
    for (j = int(0.5 / dy) + 4; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) U2[i][j][k] = 0.75 * U[i][j][k] + 0.25 * U1[i][j][k] - r * (Fp[i][j][k] + Fd[i][j][k]);

  bound(U2, dx, dy);
  dt = CFL(U2, dx, dy, COMPACTCFL);
  r = dt / dx;
  for (i = 0; i <= Nx + 6; i++)
    for (j = 0; j <= Ny + 6; j++)
      if (U[i][j][0] != 0) CPT_fpfd(U2[i][j], fp[i][j], fd[i][j]);
  Compact_x(U2, Ut, Fp, Fd, fp, fd, dx, dy, dt, p);
  // RK算U3-即U
  for (i = 3; i <= Nx + 3; i++)
    for (j = 3; j <= int(0.5 / dy) + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++)
          U[i][j][k] = 1.0 / 3.0 * U[i][j][k] + 2.0 / 3.0 * U2[i][j][k] - 2.0 / 3.0 * r * (Fp[i][j][k] + Fd[i][j][k]);

  for (i = int(1.0 / dx) + 3; i <= int(2.0 / dx) + 3; i++)
    for (j = int(0.5 / dy) + 4; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++)
          U[i][j][k] = 1.0 / 3.0 * U[i][j][k] + 2.0 / 3.0 * U2[i][j][k] - 2.0 / 3.0 * r * (Fp[i][j][k] + Fd[i][j][k]);
}
///////////-------------------------------y
void CPT_gpgd(double U[4], double gp[4], double gd[4]) {
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
  gp[0] = rou / 2.0 / GAMA * (2 * (GAMA - 1) * lambda1p + lambda3p + lambda4p);
  gp[1] = rou / 2.0 / GAMA * (u * 2 * (GAMA - 1) * lambda1p + u * lambda3p + u * lambda4p);
  gp[2] = rou / 2.0 / GAMA * (2 * v * (GAMA - 1) * lambda1p + (v - a) * lambda3p + (v + a) * lambda4p);
  gp[3] =
      rou / 2.0 / GAMA * (2 * ((GAMA - 1) * H - a * a) * lambda1p + (H - a * v) * lambda3p + (H + a * v) * lambda4p);
  //计算Gd即G-的函数
  lambda1d = -0.5 * (fabs(lambda1) - lambda1);
  lambda2d = -0.5 * (fabs(lambda2) - lambda2);
  lambda3d = -0.5 * (fabs(lambda3) - lambda3);
  lambda4d = -0.5 * (fabs(lambda4) - lambda4);
  lambda1d = 0.5 * (lambda1d - sqrt(lambda1d * lambda1d + 1e-16));
  lambda2d = 0.5 * (lambda2d - sqrt(lambda2d * lambda2d + 1e-16));
  lambda3d = 0.5 * (lambda3d - sqrt(lambda3d * lambda3d + 1e-16));
  lambda4d = 0.5 * (lambda4d - sqrt(lambda4d * lambda4d + 1e-16));
  gd[0] = rou / 2.0 / GAMA * (2 * (GAMA - 1) * lambda1d + lambda3d + lambda4d);
  gd[1] = rou / 2.0 / GAMA * (u * 2 * (GAMA - 1) * lambda1d + u * lambda3d + u * lambda4d);
  gd[2] = rou / 2.0 / GAMA * (2 * v * (GAMA - 1) * lambda1d + (v - a) * lambda3d + (v + a) * lambda4d);
  gd[3] =
      rou / 2.0 / GAMA * (2 * ((GAMA - 1) * H - a * a) * lambda1d + (H - a * v) * lambda3d + (H + a * v) * lambda4d);
}

void Compact_y(double U[Nx + 7][Ny + 7][4], double Ut[Nx + 7][Ny + 7][4], double Gp[Nx + 7][Ny + 7][4],
               double Gd[Nx + 7][Ny + 7][4], double gp[Nx + 7][Ny + 7][4], double gd[Nx + 7][Ny + 7][4], double dx,
               double dy, double dt, double p[Nx + 7][Ny + 7]) {
  int i, j, k;
  double r = dt / dy;
  double q, rou, u, v;
  double nu;
  for (i = 0; i <= Nx + 6; i++)
    for (j = 0; j <= Ny + 6; j++)
      if (U[i][j][0] != 0) {
        rou = U[i][j][0];
        u = U[i][j][1] / U[i][j][0];
        v = U[i][j][2] / U[i][j][0];
        p[i][j] = (GAMA - 1) * (U[i][j][3] - 0.5 * rou * (u * u + v * v));
      }

  nu = 5 * r * (1 - 5 * r);
  //人工黏性
  for (i = 3; i <= Nx + 3; i++)
    for (j = 3; j <= int(0.5 / dy) + 3; j++) {
      q = fabs(fabs(U[i][j + 1][0] - U[i][j][0]) - fabs(U[i][j][0] - U[i][j - 1][0])) /
          (fabs(U[i][j + 1][0] - U[i][j][0]) + fabs(U[i][j][0] - U[i][j - 1][0]) + 1e-100);
      for (k = 0; k <= 3; k++) {
        Ut[i][j][k] = U[i][j][k] + 0.5 * nu * q * (U[i][j + 1][k] - 2 * U[i][j][k] + U[i][j - 1][k]);
      }
    }
  for (i = int(1.0 / dx) + 3; i <= int(2.0 / dx) + 3; i++)
    for (j = int(0.5 / dy) + 4; j <= Ny + 3; j++) {
      q = fabs(fabs(U[i][j + 1][0] - U[i][j][0]) - fabs(U[i][j][0] - U[i][j - 1][0])) /
          (fabs(U[i][j + 1][0] - U[i][j][0]) + fabs(U[i][j][0] - U[i][j - 1][0]) + 1e-100);
      for (k = 0; k <= 3; k++) {
        Ut[i][j][k] = U[i][j][k] + 0.5 * nu * q * (U[i][j + 1][k] - 2 * U[i][j][k] + U[i][j - 1][k]);
      }
    }

  for (i = 3; i <= Nx + 3; i++)
    for (j = 3; j <= int(0.5 / dy) + 3; j++)
      for (k = 0; k <= 3; k++) U[i][j][k] = Ut[i][j][k];

  for (i = int(1.0 / dx) + 3; i <= int(2.0 / dx) + 3; i++)
    for (j = int(0.5 / dy) + 4; j <= Ny + 3; j++)
      for (k = 0; k <= 3; k++) U[i][j][k] = Ut[i][j][k];

  //计算gp和gd
  for (i = 0; i <= Nx + 6; i++)
    for (j = 0; j <= Ny + 6; j++)
      if (U[i][j][0] != 0) CPT_gpgd(U[i][j], gp[i][j], gd[i][j]);

  //分区计算G+的下值
  for (i = 3; i <= Nx + 3; i++)
    for (k = 0; k <= 3; k++) Gp[i][2][k] = 0.5 * (3 * (gp[i][3][k] - gp[i][2][k]) - (gp[i][4][k] - gp[i][3][k]));

  //分区计算全场的G+
  for (j = 3; j <= Ny + 3; j++)
    for (i = 3; i <= Nx + 3; i++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++)
          Gp[i][j][k] = -0.5 * Gp[i][j - 1][k] + 5.0 / 4.0 * (gp[i][j][k] - gp[i][j - 1][k]) +
                        0.25 * (gp[i][j + 1][k] - gp[i][j][k]);

  //分区计算G-的上值
  for (i = 3; i <= int(1.0 / dx) + 2; i++)
    for (k = 0; k <= 3; k++)
      Gd[i][int(0.5 / dy) + 4][k] = 0.5 * (3 * (gd[i][int(0.5 / dy) + 4][k] - gd[i][int(0.5 / dy) + 3][k]) -
                                           (gd[i][int(0.5 / dy) + 3][k] - gd[i][int(0.5 / dy) + 2][k]));

  for (i = int(1.0 / dx) + 3; i <= int(2.0 / dx) + 3; i++)
    for (k = 0; k <= 3; k++)
      Gd[i][int(1.0 / dy) + 4][k] = 0.5 * (3 * (gd[i][int(1.0 / dy) + 4][k] - gd[i][int(1.0 / dy) + 3][k]) -
                                           (gd[i][int(1.0 / dy) + 3][k] - gd[i][int(1.0 / dy) + 2][k]));

  for (i = int(2.0 / dx) + 4; i <= Nx + 3; i++)
    for (k = 0; k <= 3; k++)
      Gd[i][int(0.5 / dy) + 4][k] = 0.5 * (3 * (gd[i][int(0.5 / dy) + 4][k] - gd[i][int(0.5 / dy) + 3][k]) -
                                           (gd[i][int(0.5 / dy) + 3][k] - gd[i][int(0.5 / dy) + 2][k]));

  //分区计算全场的F-值
  for (j = int(0.5 / dy) + 3; j >= 3; j--)
    for (i = 3; i <= int(1.0 / dx) + 2; i++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++)
          Gd[i][j][k] = -0.5 * Gd[i][j + 1][k] + 5.0 / 4.0 * (gd[i][j + 1][k] - gd[i][j][k]) +
                        0.25 * (gd[i][j][k] - gd[i][j - 1][k]);

  for (j = int(0.5 / dy) + 3; j >= 3; j--)
    for (i = int(2.0 / dx) + 4; i <= Nx + 3; i++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++)
          Gd[i][j][k] = -0.5 * Gd[i][j + 1][k] + 5.0 / 4.0 * (gd[i][j + 1][k] - gd[i][j][k]) +
                        0.25 * (gd[i][j][k] - gd[i][j - 1][k]);

  for (j = int(1.0 / dy) + 3; j >= 3; j--)
    for (i = int(1.0 / dx) + 3; i <= int(2.0 / dx) + 3; i++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++)
          Gd[i][j][k] = -0.5 * Gd[i][j + 1][k] + 5.0 / 4.0 * (gd[i][j + 1][k] - gd[i][j][k]) +
                        0.25 * (gd[i][j][k] - gd[i][j - 1][k]);
}
void Compact_RK_y(double U[Nx + 7][Ny + 7][4], double Ut[Nx + 7][Ny + 7][4], double U1[Nx + 7][Ny + 7][4],
                  double U2[Nx + 7][Ny + 7][4], double Gp[Nx + 7][Ny + 7][4], double Gd[Nx + 7][Ny + 7][4],
                  double gp[Nx + 7][Ny + 7][4], double gd[Nx + 7][Ny + 7][4], double dx, double dy, double dt,
                  double p[Nx + 7][Ny + 7], double LAMDA_[Nx + 7][Ny + 7][4][4]) {
  int i, j, k;
  dt = CFL(U, dx, dy, COMPACTCFL);
  double r = dt / dy;
  //分区计算U
  for (i = 3; i <= Nx + 3; i++)
    for (j = 3; j <= int(0.5 / dy) + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) U1[i][j][k] = U[i][j][k] - r * (Gp[i][j][k] + Gd[i][j][k]);

  for (i = int(1.0 / dx) + 3; i <= int(2.0 / dx) + 3; i++)
    for (j = int(0.5 / dy) + 4; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) U1[i][j][k] = U[i][j][k] - r * (Gp[i][j][k] + Gd[i][j][k]);

  bound(U1, dx, dy);
  dt = CFL(U1, dx, dy, COMPACTCFL);
  r = dt / dy;
  for (i = 0; i <= Nx + 6; i++)
    for (j = 0; j <= Ny + 6; j++)
      if (U[i][j][0] != 0) CPT_gpgd(U1[i][j], gp[i][j], gd[i][j]);
  Compact_y(U1, Ut, Gp, Gd, gp, gd, dx, dy, dt, p);
  for (i = 3; i <= Nx + 3; i++)
    for (j = 3; j <= int(0.5 / dy) + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) U2[i][j][k] = 0.75 * U[i][j][k] + 0.25 * U1[i][j][k] - r * (Gp[i][j][k] + Gd[i][j][k]);

  for (i = int(1.0 / dx) + 3; i <= int(2.0 / dx) + 3; i++)
    for (j = int(0.5 / dy) + 4; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) U2[i][j][k] = 0.75 * U[i][j][k] + 0.25 * U1[i][j][k] - r * (Gp[i][j][k] + Gd[i][j][k]);

  bound(U2, dx, dy);
  dt = CFL(U2, dx, dy, COMPACTCFL);
  r = dt / dy;
  for (i = 0; i <= Nx + 6; i++)
    for (j = 0; j <= Ny + 6; j++)
      if (U[i][j][0] != 0) CPT_gpgd(U2[i][j], gp[i][j], gd[i][j]);
  Compact_y(U2, Ut, Gp, Gd, gp, gd, dx, dy, dt, p);
  for (i = 3; i <= Nx + 3; i++)
    for (j = 3; j <= int(0.5 / dy) + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++)
          U[i][j][k] = 1.0 / 3.0 * U[i][j][k] + 2.0 / 3.0 * U2[i][j][k] - 2.0 / 3.0 * r * (Gp[i][j][k] + Gd[i][j][k]);

  for (i = int(1.0 / dx) + 3; i <= int(2.0 / dx) + 3; i++)
    for (j = int(0.5 / dy) + 4; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++)
          U[i][j][k] = 1.0 / 3.0 * U[i][j][k] + 2.0 / 3.0 * U2[i][j][k] - 2.0 / 3.0 * r * (Gp[i][j][k] + Gd[i][j][k]);
}


#endif // COMPACT_FUSION

#endif