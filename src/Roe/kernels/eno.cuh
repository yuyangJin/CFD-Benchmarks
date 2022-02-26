#ifndef __ENO_H
#define __ENO_H

#include "util.cuh"

#ifndef ENO_FUSION
//利用Roe平均以后的值计算x方向的特征值与特征向量
void LF_x(double U[Nx + 7][Ny + 7][4], double LAMDA_[Nx + 7][Ny + 7][4][4], double F[Nx + 7][Ny + 7][4],
          double Fp[Nx + 7][Ny + 7][4], double Fd[Nx + 7][Ny + 7][4]) {
  int i, j, k;
  double rou, u, v, p, a, maxlamda = 1e-100;
  //对LAMDA_赋值
  for (i = 0; i <= Nx + 6; i++)
    for (j = 0; j <= Ny + 6; j++)
      if (U[i][j][0] != 0) {
        rou = U[i][j][0];
        u = U[i][j][1] / U[i][j][0];
        v = U[i][j][2] / U[i][j][0];
        p = (GAMA - 1) * (U[i][j][3] - 0.5 * rou * (u * u + v * v));
        a = sqrt(p * GAMA / rou);
        LAMDA_[i][j][0][0] = u - a;
        LAMDA_[i][j][1][1] = u;
        LAMDA_[i][j][2][2] = u + a;
        LAMDA_[i][j][3][3] = u;
        for (k = 0; k <= 3; k++) {
          if (fabs(LAMDA_[i][j][k][k]) >= maxlamda) maxlamda = fabs(LAMDA_[i][j][k][k]);
        }
      }
  // U2F
  for (i = 0; i <= Nx + 6; i++)
    for (j = 0; j <= Ny + 6; j++)
      if (U[i][j][0] != 0) U2F(U[i][j], F[i][j]);

  //计算Fpd
  for (i = 0; i <= Nx + 6; i++)
    for (j = 0; j <= Ny + 6; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) {
          Fp[i][j][k] = 0.5 * (F[i][j][k] + maxlamda * U[i][j][k]);
          Fd[i][j][k] = 0.5 * (F[i][j][k] - maxlamda * U[i][j][k]);
        }
}

//此函数将之前计算的特征之余特征向量通过TVD算法计算U
void ENO_x(double U[Nx + 7][Ny + 7][4], double F[Nx + 7][Ny + 7][4], double Fp[Nx + 7][Ny + 7][4],
           double Fd[Nx + 7][Ny + 7][4], double F_p[Nx + 7][Ny + 7][4], double F_d[Nx + 7][Ny + 7][4],
           double F_[Nx + 7][Ny + 7][4], double q3p[Nx + 7][Ny + 7][4][3], double q3d[Nx + 7][Ny + 7][4][3], double dx,
           double dy, double dt) {
  int i, j, k;
  double r;
  dt = CFL(U, dx, dy, ENOCFL);
  r = dt / dx;

  //计算q3pd q5pd
  for (i = 2; i <= Nx + 3; i++)
    for (j = 2; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) {
          q3p[i][j][k][0] = 1.0 / 3.0 * Fp[i - 2][j][k] - 7.0 / 6.0 * Fp[i - 1][j][k] + 11.0 / 6.0 * Fp[i][j][k];
          q3p[i][j][k][1] = -1.0 / 6.0 * Fp[i - 1][j][k] + 5.0 / 6.0 * Fp[i][j][k] + 1.0 / 3.0 * Fp[i + 1][j][k];
          q3p[i][j][k][2] = 1.0 / 3.0 * Fp[i][j][k] + 5.0 / 6.0 * Fp[i + 1][j][k] - 1.0 / 6.0 * Fp[i + 2][j][k];

          q3d[i][j][k][0] = -1.0 / 6.0 * Fd[i - 1][j][k] + 5.0 / 6.0 * Fd[i][j][k] + 1.0 / 3.0 * Fd[i + 1][j][k];
          q3d[i][j][k][1] = 1.0 / 3.0 * Fd[i][j][k] + 5.0 / 6.0 * Fd[i + 1][j][k] - 1.0 / 6.0 * Fd[i + 2][j][k];
          q3d[i][j][k][2] = 11.0 / 6.0 * Fd[i + 1][j][k] - 7.0 / 6.0 * Fd[i + 2][j][k] + 1.0 / 3.0 * Fd[i + 3][j][k];
        }

  //判断差商大小并且赋值
  for (i = 2; i <= Nx + 3; i++)
    for (j = 2; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) {
          if (fabs(U[i + 1][j][k] - U[i][j][k]) > fabs(U[i][j][k] - U[i - 1][j][k]) &&
              fabs(U[i + 1][j][k] - U[i][j][k]) > fabs(U[i - 1][j][k] - U[i - 2][j][k]))
            F_p[i][j][k] = q3p[i][j][k][0];
          else {
            if (fabs(U[i][j][k] - U[i - 1][j][k]) > fabs(U[i + 1][j][k] - U[i][j][k]) &&
                fabs(U[i][j][k] - U[i - 1][j][k]) > fabs(U[i + 2][j][k] - U[i + 1][j][k]))
              F_p[i][j][k] = q3p[i][j][k][2];
            else
              F_p[i][j][k] = q3p[i][j][k][1];
          }
        }

  for (i = 2; i <= Nx + 3; i++)
    for (j = 2; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) {
          if (fabs(U[i + 2][j][k] - U[i + 1][j][k]) > fabs(U[i + 1][j][k] - U[i][j][k]) &&
              fabs(U[i + 2][j][k] - U[i + 1][j][k]) > fabs(U[i][j][k] - U[i - 1][j][k]))
            F_d[i][j][k] = q3d[i][j][k][0];
          else {
            if (fabs(U[i + 1][j][k] - U[i][j][k]) > fabs(U[i + 2][j][k] - U[i + 1][j][k]) &&
                fabs(U[i + 1][j][k] - U[i][j][k]) > fabs(U[i + 3][j][k] - U[i + 2][j][k]))
              F_d[i][j][k] = q3d[i][j][k][2];
            else
              F_d[i][j][k] = q3d[i][j][k][1];
          }
        }

  //计算F_
  for (i = 2; i <= Nx + 3; i++)
    for (j = 2; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) F_[i][j][k] = F_p[i][j][k] + F_d[i][j][k];
  //分区计算U
  for (i = 3; i <= Nx + 3; i++)
    for (j = 3; j <= int(0.5 / dy) + 3; j++)
      for (k = 0; k <= 3; k++) U[i][j][k] = U[i][j][k] - r * (F_[i][j][k] - F_[i - 1][j][k]);

  for (i = int(1.0 / dx) + 3; i <= int(2.0 / dx) + 3; i++)
    for (j = int(0.5 / dy) + 4; j <= Ny + 3; j++)
      for (k = 0; k <= 3; k++) U[i][j][k] = U[i][j][k] - r * (F_[i][j][k] - F_[i - 1][j][k]);
}

//利用Roe平均以后的值计算x方向的特征值与特征向量
void LF_y(double U[Nx + 7][Ny + 7][4], double LAMDA_[Nx + 7][Ny + 7][4][4], double G[Nx + 7][Ny + 7][4],
          double Gp[Nx + 7][Ny + 7][4], double Gd[Nx + 7][Ny + 7][4]) {
  int i, j, k;
  double rou, u, v, a, p, maxlamda = 1e-100;
  ;

  //对LAMDA_赋值
  for (i = 0; i <= Nx + 5; i++)
    for (j = 0; j <= Ny + 5; j++)
      if (U[i][j][0] != 0) {
        rou = U[i][j][0];
        u = U[i][j][1] / U[i][j][0];
        v = U[i][j][2] / U[i][j][0];
        p = (GAMA - 1) * (U[i][j][3] - 0.5 * rou * (u * u + v * v));
        // a=pow(p*GAMA/rou,0.5);
        a = sqrt(p * GAMA / rou);
        LAMDA_[i][j][0][0] = v - a;
        LAMDA_[i][j][1][1] = v;
        LAMDA_[i][j][2][2] = v + a;
        LAMDA_[i][j][3][3] = v;
        for (k = 0; k <= 3; k++) {
          if (fabs(LAMDA_[i][j][k][k]) >= maxlamda) maxlamda = fabs(LAMDA_[i][j][k][k]);
        }
      }

  // U2G
  for (i = 0; i <= Nx + 6; i++)
    for (j = 0; j <= Ny + 6; j++)
      if (U[i][j][0] != 0) U2G(U[i][j], G[i][j]);

  //计算Gpd
  for (i = 0; i <= Nx + 6; i++)
    for (j = 0; j <= Ny + 6; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) {
          Gp[i][j][k] = 0.5 * (G[i][j][k] + maxlamda * U[i][j][k]);
          Gd[i][j][k] = 0.5 * (G[i][j][k] - maxlamda * U[i][j][k]);
        }
}

//此函数将之前计算的特征之余特征向量通过TVD算法计算U
void ENO_y(double U[Nx + 7][Ny + 7][4], double G[Nx + 7][Ny + 7][4], double Gp[Nx + 7][Ny + 7][4],
           double Gd[Nx + 7][Ny + 7][4], double G_p[Nx + 7][Ny + 7][4], double G_d[Nx + 7][Ny + 7][4],
           double G_[Nx + 7][Ny + 7][4], double q3p[Nx + 7][Ny + 7][4][3], double q3d[Nx + 7][Ny + 7][4][3], double dx,
           double dy, double dt) {
  int i, j, k;
  double r;
  dt = CFL(U, dx, dy, ENOCFL);
  r = dt / dy;

  //计算q3pd q5pd
  for (i = 2; i <= Nx + 3; i++)
    for (j = 2; j <= Ny + 3; j++)
      if (U[i][j][0] != 0) {
        for (k = 0; k <= 3; k++) {
          q3p[i][j][k][0] = 1.0 / 3.0 * Gp[i][j - 2][k] - 7.0 / 6.0 * Gp[i][j - 1][k] + 11.0 / 6.0 * Gp[i][j][k];
          q3p[i][j][k][1] = -1.0 / 6.0 * Gp[i][j - 1][k] + 5.0 / 6.0 * Gp[i][j][k] + 1.0 / 3.0 * Gp[i][j + 1][k];
          q3p[i][j][k][2] = 1.0 / 3.0 * Gp[i][j][k] + 5.0 / 6.0 * Gp[i][j + 1][k] - 1.0 / 6.0 * Gp[i][j + 2][k];

          q3d[i][j][k][0] = -1.0 / 6.0 * Gd[i][j - 1][k] + 5.0 / 6.0 * Gd[i][j][k] + 1.0 / 3.0 * Gd[i][j + 1][k];
          q3d[i][j][k][1] = 1.0 / 3.0 * Gd[i][j][k] + 5.0 / 6.0 * Gd[i][j + 1][k] - 1.0 / 6.0 * Gd[i][j + 2][k];
          q3d[i][j][k][2] = 11.0 / 6.0 * Gd[i][j + 1][k] - 7.0 / 6.0 * Gd[i][j + 2][k] + 1.0 / 3.0 * Gd[i][j + 3][k];
        }
      }

  //判断差商大小并且赋值
  for (i = 2; i <= Nx + 3; i++)
    for (j = 2; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) {
          if (fabs(U[i][j + 1][k] - U[i][j][k]) > fabs(U[i][j][k] - U[i][j - 1][k]) &&
              fabs(U[i][j + 1][k] - U[i][j][k]) > fabs(U[i][j - 1][k] - U[i][j - 2][k]))
            G_p[i][j][k] = q3p[i][j][k][0];
          else {
            if (fabs(U[i][j][k] - U[i][j - 1][k]) > fabs(U[i][j + 1][k] - U[i][j][k]) &&
                fabs(U[i][j][k] - U[i][j - 1][k]) > fabs(U[i][j + 2][k] - U[i][j + 1][k]))
              G_p[i][j][k] = q3p[i][j][k][2];
            else
              G_p[i][j][k] = q3p[i][j][k][1];
          }
        }

  for (i = 2; i <= Nx + 3; i++)
    for (j = 2; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) {
          if (fabs(U[i][j + 2][k] - U[i][j + 1][k]) > fabs(U[i][j + 1][k] - U[i][j][k]) &&
              fabs(U[i][j + 2][k] - U[i][j + 1][k]) > fabs(U[i][j][k] - U[i][j - 1][k]))
            G_d[i][j][k] = q3d[i][j][k][0];
          else {
            if (fabs(U[i][j + 1][k] - U[i][j][k]) > fabs(U[i][j + 2][k] - U[i][j + 1][k]) &&
                fabs(U[i][j + 1][k] - U[i][j][k]) > fabs(U[i][j + 3][k] - U[i][j + 2][k]))
              G_d[i][j][k] = q3d[i][j][k][2];
            else
              G_d[i][j][k] = q3d[i][j][k][1];
          }
        }

  //计算G_
  for (i = 2; i <= Nx + 3; i++)
    for (j = 2; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) G_[i][j][k] = G_p[i][j][k] + G_d[i][j][k];

  //分区计算U
  for (i = 3; i <= Nx + 3; i++)
    for (j = 3; j <= int(0.5 / dy) + 3; j++)
      for (k = 0; k <= 3; k++) U[i][j][k] = U[i][j][k] - r * (G_[i][j][k] - G_[i][j - 1][k]);

  for (i = int(1.0 / dx) + 3; i <= int(2.0 / dx) + 3; i++)
    for (j = int(0.5 / dy) + 4; j <= Ny + 3; j++)
      for (k = 0; k <= 3; k++) U[i][j][k] = U[i][j][k] - r * (G_[i][j][k] - G_[i][j - 1][k]);
}

#else

//利用Roe平均以后的值计算x方向的特征值与特征向量
void LF_x(double*** U, double**** LAMDA_, double*** F, double*** Fp, double*** Fd) {
  int i, j, k;
  double rou, u, v, p, a, maxlamda = 1e-100;
  //对LAMDA_赋值
  for (i = 0; i <= Nx + 6; i++)
    for (j = 0; j <= Ny + 6; j++)
      if (U[i][j][0] != 0) {
        rou = U[i][j][0];
        u = U[i][j][1] / U[i][j][0];
        v = U[i][j][2] / U[i][j][0];
        p = (GAMA - 1) * (U[i][j][3] - 0.5 * rou * (u * u + v * v));
        a = sqrt(p * GAMA / rou);
        LAMDA_[i][j][0][0] = u - a;
        LAMDA_[i][j][1][1] = u;
        LAMDA_[i][j][2][2] = u + a;
        LAMDA_[i][j][3][3] = u;
        for (k = 0; k <= 3; k++) {
          if (fabs(LAMDA_[i][j][k][k]) >= maxlamda) maxlamda = fabs(LAMDA_[i][j][k][k]);
        }
      }
  for (i = 0; i <= Nx + 6; i++)
    for (j = 0; j <= Ny + 6; j++)
      if (U[i][j][0] != 0) {
        // U2F
        U2F(U[i][j], F[i][j]);
        //计算Fpd
        for (k = 0; k <= 3; k++) {
          Fp[i][j][k] = 0.5 * (F[i][j][k] + maxlamda * U[i][j][k]);
          Fd[i][j][k] = 0.5 * (F[i][j][k] - maxlamda * U[i][j][k]);
        }
      }
}

__global__ void eno_x_cuda(double* U, double* F, double* Fp, double* Fd, double* F_p, double* F_d, double* F_,
                           double* q3p, double* q3d, double dx, double dy, double dt, int offset1, int offset2) {
  int i = blockIdx.x * blockDim.x + threadIdx.x + 2;
  int j = blockIdx.y * blockDim.y + threadIdx.y + 2;

  // for (i = 2; i <= Nx + 3; i++)
  //   for (j = 2; j <= Ny + 3; j++)
  if (i >= Nx + 2 || j >= Ny + 2) return;
  if (1)
    // if (U[i * offset1 + j * offset2] != 0)
    for (int k = 0; k <= 3; k++) {
      //判断差商大小并且赋值
      if (fabs(U[(i + 1) * offset1 + j * offset2 + k] - U[i * offset1 + j * offset2 + k]) >
              fabs(U[i * offset1 + j * offset2 + k] - U[(i - 1) * offset1 + j * offset2 + k]) &&
          fabs(U[(i + 1) * offset1 + j * offset2 + k] - U[i * offset1 + j * offset2 + k]) >
              fabs(U[(i - 1) * offset1 + j * offset2 + k] - U[(i - 2) * offset1 + j * offset2 + k]))
        // F_p[i*offset1 + j*offset2 + k] = q3p[i*offset1*3 + j*offset2*3 + k*3+0];
        F_p[i * offset1 + j * offset2 + k] = 1.0 / 3.0 * Fp[(i - 2) * offset1 + j * offset2 + k] -
                                             7.0 / 6.0 * Fp[(i - 1) * offset1 + j * offset2 + k] +
                                             11.0 / 6.0 * Fp[i * offset1 + j * offset2 + k];
      else {
        if (fabs(U[i * offset1 + j * offset2 + k] - U[(i - 1) * offset1 + j * offset2 + k]) >
                fabs(U[(i + 1) * offset1 + j * offset2 + k] - U[i * offset1 + j * offset2 + k]) &&
            fabs(U[i * offset1 + j * offset2 + k] - U[(i - 1) * offset1 + j * offset2 + k]) >
                fabs(U[(i + 2) * offset1 + j * offset2 + k] - U[(i + 1) * offset1 + j * offset2 + k]))
          // F_p[i*offset1 + j*offset2 + k] = q3p[i*offset1*3 + j*offset2*3 + k*3+2];
          F_p[i * offset1 + j * offset2 + k] = 1.0 / 3.0 * Fp[i * offset1 + j * offset2 + k] +
                                               5.0 / 6.0 * Fp[(i + 1) * offset1 + j * offset2 + k] -
                                               1.0 / 6.0 * Fp[(i + 2) * offset1 + j * offset2 + k];
        else
          // F_p[i*offset1 + j*offset2 + k] = q3p[i*offset1*3 + j*offset2*3 + k*3+1];
          F_p[i * offset1 + j * offset2 + k] = -1.0 / 6.0 * Fp[(i - 1) * offset1 + j * offset2 + k] +
                                               5.0 / 6.0 * Fp[i * offset1 + j * offset2 + k] +
                                               1.0 / 3.0 * Fp[(i + 1) * offset1 + j * offset2 + k];
      }
      if (fabs(U[(i + 2) * offset1 + j * offset2 + k] - U[(i + 1) * offset1 + j * offset2 + k]) >
              fabs(U[(i + 1) * offset1 + j * offset2 + k] - U[i * offset1 + j * offset2 + k]) &&
          fabs(U[(i + 2) * offset1 + j * offset2 + k] - U[(i + 1) * offset1 + j * offset2 + k]) >
              fabs(U[i * offset1 + j * offset2 + k] - U[(i - 1) * offset1 + j * offset2 + k]))
        // F_d[i*offset1 + j*offset2 + k] = q3d[i*offset1*3 + j*offset2*3 + k*3+0];
        F_d[i * offset1 + j * offset2 + k] = -1.0 / 6.0 * Fd[(i - 1) * offset1 + j * offset2 + k] +
                                             5.0 / 6.0 * Fd[i * offset1 + j * offset2 + k] +
                                             1.0 / 3.0 * Fd[(i + 1) * offset1 + j * offset2 + k];
      else {
        if (fabs(U[(i + 1) * offset1 + j * offset2 + k] - U[i * offset1 + j * offset2 + k]) >
                fabs(U[(i + 2) * offset1 + j * offset2 + k] - U[(i + 1) * offset1 + j * offset2 + k]) &&
            fabs(U[(i + 1) * offset1 + j * offset2 + k] - U[i * offset1 + j * offset2 + k]) >
                fabs(U[(i + 3) * offset1 + j * offset2 + k] - U[(i + 2) * offset1 + j * offset2 + k]))
          // F_d[i*offset1 + j*offset2 + k] = q3d[i*offset1*3 + j*offset2*3 + k*3+2];
          F_d[i * offset1 + j * offset2 + k] = 11.0 / 6.0 * Fd[(i + 1) * offset1 + j * offset2 + k] -
                                               7.0 / 6.0 * Fd[(i + 2) * offset1 + j * offset2 + k] +
                                               1.0 / 3.0 * Fd[(i + 3) * offset1 + j * offset2 + k];
        else
          // F_d[i*offset1 + j*offset2 + k] = q3d[i*offset1*3 + j*offset2*3 + k*3+1];
          F_d[i * offset1 + j * offset2 + k] = 1.0 / 3.0 * Fd[i * offset1 + j * offset2 + k] +
                                               5.0 / 6.0 * Fd[(i + 1) * offset1 + j * offset2 + k] -
                                               1.0 / 6.0 * Fd[(i + 2) * offset1 + j * offset2 + k];
      }

      // 计算F_
      F_[i * offset1 + j * offset2 + k] = F_p[i * offset1 + j * offset2 + k] + F_d[i * offset1 + j * offset2 + k];
    }
}

// for (i = 3; i <= Nx + 3; i++)
//   for (j = 3; j <= int(0.5 / dy) + 3; j++)
//     for (k = 0; k <= 3; k++) U[i][j][k] = U[i][j][k] - r * (F_[i][j][k] - F_[i - 1][j][k]);

// for (i = int(1.0 / dx) + 3; i <= int(2.0 / dx) + 3; i++)
//   for (j = int(0.5 / dy) + 4; j <= Ny + 3; j++)
//     for (k = 0; k <= 3; k++) U[i][j][k] = U[i][j][k] - r * (F_[i][j][k] - F_[i - 1][j][k]);

__global__ void eno_x_cuda_2(double* U, double* F_, double dy, double r, int offset1, int offset2) {
  // for (i = 3; i <= Nx + 3; i++)
  //   for (j = 3; j <= int(0.5 / dy) + 3; j++)
  int i = blockIdx.x * blockDim.x + threadIdx.x + 3;
  int j = blockIdx.y * blockDim.y + threadIdx.y + 3;
  if (i > Nx + 3 || j > int(0.5 / dy) + 3) return;
  for (int k = 0; k <= 3; k++)
    U[i * offset1 + j * offset2 + k] = U[i * offset1 + j * offset2 + k] - r * (F_[i * offset1 + j * offset2 + k] -
                                                                               F_[(i - 1) * offset1 + j * offset2 + k]);
}

__global__ void eno_x_cuda_3(double* U, double* F_, double dx, double dy, double r, int offset1, int offset2) {
  // for (i = 3; i <= Nx + 3; i++)
  //   for (j = 3; j <= int(0.5 / dy) + 3; j++)
  int i = blockIdx.x * blockDim.x + threadIdx.x + 3 + int(1.f / dx);
  int j = blockIdx.y * blockDim.y + threadIdx.y + 4 + int(0.5 / dy);
  if (i > int(2.0 / dx) + 3 || j > Ny + 3) return;
  for (int k = 0; k <= 3; k++)
    U[i * offset1 + j * offset2 + k] = U[i * offset1 + j * offset2 + k] - r * (F_[i * offset1 + j * offset2 + k] -
                                                                               F_[(i - 1) * offset1 + j * offset2 + k]);
}
//此函数将之前计算的特征之余特征向量通过TVD算法计算U
void ENO_x(double*** U, double*** F, double*** Fp, double*** Fd, double*** F_p, double*** F_d, double*** F_,
           double**** q3p, double**** q3d, double dx, double dy, double dt) {
  int i, j, k;
  double r;
  dt = CFL(U, dx, dy, ENOCFL);
  r = dt / dx;

  for (i = 2; i <= Nx + 3; i++)
    for (j = 2; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) {
          //计算q3pd q5pd
          /*
          q3p[i][j][k][0] = 1.0 / 3.0 * Fp[i - 2][j][k] - 7.0 / 6.0 * Fp[i - 1][j][k] + 11.0 / 6.0 * Fp[i][j][k];
          q3p[i][j][k][1] = -1.0 / 6.0 * Fp[i - 1][j][k] + 5.0 / 6.0 * Fp[i][j][k] + 1.0 / 3.0 * Fp[i + 1][j][k];
          q3p[i][j][k][2] = 1.0 / 3.0 * Fp[i][j][k] + 5.0 / 6.0 * Fp[i + 1][j][k] - 1.0 / 6.0 * Fp[i + 2][j][k];

          q3d[i][j][k][0] = -1.0 / 6.0 * Fd[i - 1][j][k] + 5.0 / 6.0 * Fd[i][j][k] + 1.0 / 3.0 * Fd[i + 1][j][k];
          q3d[i][j][k][1] = 1.0 / 3.0 * Fd[i][j][k] + 5.0 / 6.0 * Fd[i + 1][j][k] - 1.0 / 6.0 * Fd[i + 2][j][k];
          q3d[i][j][k][2] = 11.0 / 6.0 * Fd[i + 1][j][k] - 7.0 / 6.0 * Fd[i + 2][j][k] + 1.0 / 3.0 * Fd[i + 3][j][k];
          */

          //判断差商大小并且赋值
          if (fabs(U[i + 1][j][k] - U[i][j][k]) > fabs(U[i][j][k] - U[i - 1][j][k]) &&
              fabs(U[i + 1][j][k] - U[i][j][k]) > fabs(U[i - 1][j][k] - U[i - 2][j][k]))
            // F_p[i][j][k] = q3p[i][j][k][0];
            F_p[i][j][k] = 1.0 / 3.0 * Fp[i - 2][j][k] - 7.0 / 6.0 * Fp[i - 1][j][k] + 11.0 / 6.0 * Fp[i][j][k];
          else {
            if (fabs(U[i][j][k] - U[i - 1][j][k]) > fabs(U[i + 1][j][k] - U[i][j][k]) &&
                fabs(U[i][j][k] - U[i - 1][j][k]) > fabs(U[i + 2][j][k] - U[i + 1][j][k]))
              // F_p[i][j][k] = q3p[i][j][k][2];
              F_p[i][j][k] = 1.0 / 3.0 * Fp[i][j][k] + 5.0 / 6.0 * Fp[i + 1][j][k] - 1.0 / 6.0 * Fp[i + 2][j][k];
            else
              // F_p[i][j][k] = q3p[i][j][k][1];
              F_p[i][j][k] = -1.0 / 6.0 * Fp[i - 1][j][k] + 5.0 / 6.0 * Fp[i][j][k] + 1.0 / 3.0 * Fp[i + 1][j][k];
          }
          if (fabs(U[i + 2][j][k] - U[i + 1][j][k]) > fabs(U[i + 1][j][k] - U[i][j][k]) &&
              fabs(U[i + 2][j][k] - U[i + 1][j][k]) > fabs(U[i][j][k] - U[i - 1][j][k]))
            // F_d[i][j][k] = q3d[i][j][k][0];
            F_d[i][j][k] = -1.0 / 6.0 * Fd[i - 1][j][k] + 5.0 / 6.0 * Fd[i][j][k] + 1.0 / 3.0 * Fd[i + 1][j][k];
          else {
            if (fabs(U[i + 1][j][k] - U[i][j][k]) > fabs(U[i + 2][j][k] - U[i + 1][j][k]) &&
                fabs(U[i + 1][j][k] - U[i][j][k]) > fabs(U[i + 3][j][k] - U[i + 2][j][k]))
              // F_d[i][j][k] = q3d[i][j][k][2];
              F_d[i][j][k] = 11.0 / 6.0 * Fd[i + 1][j][k] - 7.0 / 6.0 * Fd[i + 2][j][k] + 1.0 / 3.0 * Fd[i + 3][j][k];
            else
              // F_d[i][j][k] = q3d[i][j][k][1];
              F_d[i][j][k] = 1.0 / 3.0 * Fd[i][j][k] + 5.0 / 6.0 * Fd[i + 1][j][k] - 1.0 / 6.0 * Fd[i + 2][j][k];
          }

          // 计算F_
          F_[i][j][k] = F_p[i][j][k] + F_d[i][j][k];
        }

  //分区计算U
  for (i = 3; i <= Nx + 3; i++)
    for (j = 3; j <= int(0.5 / dy) + 3; j++)
      for (k = 0; k <= 3; k++) U[i][j][k] = U[i][j][k] - r * (F_[i][j][k] - F_[i - 1][j][k]);

  for (i = int(1.0 / dx) + 3; i <= int(2.0 / dx) + 3; i++)
    for (j = int(0.5 / dy) + 4; j <= Ny + 3; j++)
      for (k = 0; k <= 3; k++) U[i][j][k] = U[i][j][k] - r * (F_[i][j][k] - F_[i - 1][j][k]);
}

//利用Roe平均以后的值计算x方向的特征值与特征向量
void LF_y(double*** U, double**** LAMDA_, double*** G, double*** Gp, double*** Gd) {
  int i, j, k;
  double rou, u, v, a, p, maxlamda = 1e-100;
  ;

  //对LAMDA_赋值
  for (i = 0; i <= Nx + 5; i++)
    for (j = 0; j <= Ny + 5; j++)
      if (U[i][j][0] != 0) {
        rou = U[i][j][0];
        u = U[i][j][1] / U[i][j][0];
        v = U[i][j][2] / U[i][j][0];
        p = (GAMA - 1) * (U[i][j][3] - 0.5 * rou * (u * u + v * v));
        // a=pow(p*GAMA/rou,0.5);
        a = sqrt(p * GAMA / rou);
        LAMDA_[i][j][0][0] = v - a;
        LAMDA_[i][j][1][1] = v;
        LAMDA_[i][j][2][2] = v + a;
        LAMDA_[i][j][3][3] = v;
        for (k = 0; k <= 3; k++) {
          if (fabs(LAMDA_[i][j][k][k]) >= maxlamda) maxlamda = fabs(LAMDA_[i][j][k][k]);
        }
      }

  for (i = 0; i <= Nx + 6; i++)
    for (j = 0; j <= Ny + 6; j++)
      if (U[i][j][0] != 0) {
        // U2G
        U2G(U[i][j], G[i][j]);
        //计算Gpd
        for (k = 0; k <= 3; k++) {
          Gp[i][j][k] = 0.5 * (G[i][j][k] + maxlamda * U[i][j][k]);
          Gd[i][j][k] = 0.5 * (G[i][j][k] - maxlamda * U[i][j][k]);
        }
      }
}

//此函数将之前计算的特征之余特征向量通过TVD算法计算U
void ENO_y(double*** U, double*** G, double*** Gp, double*** Gd, double*** G_p, double*** G_d, double*** G_,
           double**** q3p, double**** q3d, double dx, double dy, double dt) {
  int i, j, k;
  double r;
  dt = CFL(U, dx, dy, ENOCFL);
  r = dt / dy;

  for (i = 2; i <= Nx + 3; i++)
    for (j = 2; j <= Ny + 3; j++)
      if (U[i][j][0] != 0) {
        for (k = 0; k <= 3; k++) {
          //计算q3pd q5pd
          /*
          q3p[i][j][k][0] = 1.0 / 3.0 * Gp[i][j - 2][k] - 7.0 / 6.0 * Gp[i][j - 1][k] + 11.0 / 6.0 * Gp[i][j][k];
          q3p[i][j][k][1] = -1.0 / 6.0 * Gp[i][j - 1][k] + 5.0 / 6.0 * Gp[i][j][k] + 1.0 / 3.0 * Gp[i][j + 1][k];
          q3p[i][j][k][2] = 1.0 / 3.0 * Gp[i][j][k] + 5.0 / 6.0 * Gp[i][j + 1][k] - 1.0 / 6.0 * Gp[i][j + 2][k];

          q3d[i][j][k][0] = -1.0 / 6.0 * Gd[i][j - 1][k] + 5.0 / 6.0 * Gd[i][j][k] + 1.0 / 3.0 * Gd[i][j + 1][k];
          q3d[i][j][k][1] = 1.0 / 3.0 * Gd[i][j][k] + 5.0 / 6.0 * Gd[i][j + 1][k] - 1.0 / 6.0 * Gd[i][j + 2][k];
          q3d[i][j][k][2] = 11.0 / 6.0 * Gd[i][j + 1][k] - 7.0 / 6.0 * Gd[i][j + 2][k] + 1.0 / 3.0 * Gd[i][j + 3][k];
          */

          //判断差商大小并且赋值
          if (fabs(U[i][j + 1][k] - U[i][j][k]) > fabs(U[i][j][k] - U[i][j - 1][k]) &&
              fabs(U[i][j + 1][k] - U[i][j][k]) > fabs(U[i][j - 1][k] - U[i][j - 2][k]))
            // G_p[i][j][k] = q3p[i][j][k][0];
            G_p[i][j][k] = 1.0 / 3.0 * Gp[i][j - 2][k] - 7.0 / 6.0 * Gp[i][j - 1][k] + 11.0 / 6.0 * Gp[i][j][k];
          else {
            if (fabs(U[i][j][k] - U[i][j - 1][k]) > fabs(U[i][j + 1][k] - U[i][j][k]) &&
                fabs(U[i][j][k] - U[i][j - 1][k]) > fabs(U[i][j + 2][k] - U[i][j + 1][k]))
              // G_p[i][j][k] = q3p[i][j][k][2];
              G_p[i][j][k] = 1.0 / 3.0 * Gp[i][j][k] + 5.0 / 6.0 * Gp[i][j + 1][k] - 1.0 / 6.0 * Gp[i][j + 2][k];
            else
              // G_p[i][j][k] = q3p[i][j][k][1];
              G_p[i][j][k] = -1.0 / 6.0 * Gp[i][j - 1][k] + 5.0 / 6.0 * Gp[i][j][k] + 1.0 / 3.0 * Gp[i][j + 1][k];
          }
          if (fabs(U[i][j + 2][k] - U[i][j + 1][k]) > fabs(U[i][j + 1][k] - U[i][j][k]) &&
              fabs(U[i][j + 2][k] - U[i][j + 1][k]) > fabs(U[i][j][k] - U[i][j - 1][k]))
            // G_d[i][j][k] = q3d[i][j][k][0];
            G_d[i][j][k] = -1.0 / 6.0 * Gd[i][j - 1][k] + 5.0 / 6.0 * Gd[i][j][k] + 1.0 / 3.0 * Gd[i][j + 1][k];
          else {
            if (fabs(U[i][j + 1][k] - U[i][j][k]) > fabs(U[i][j + 2][k] - U[i][j + 1][k]) &&
                fabs(U[i][j + 1][k] - U[i][j][k]) > fabs(U[i][j + 3][k] - U[i][j + 2][k]))
              // G_d[i][j][k] = q3d[i][j][k][2];
              G_d[i][j][k] = 11.0 / 6.0 * Gd[i][j + 1][k] - 7.0 / 6.0 * Gd[i][j + 2][k] + 1.0 / 3.0 * Gd[i][j + 3][k];
            else
              // G_d[i][j][k] = q3d[i][j][k][1];
              G_d[i][j][k] = 1.0 / 3.0 * Gd[i][j][k] + 5.0 / 6.0 * Gd[i][j + 1][k] - 1.0 / 6.0 * Gd[i][j + 2][k];
          }

          //计算G_
          G_[i][j][k] = G_p[i][j][k] + G_d[i][j][k];
        }
      }

  //分区计算U
  for (i = 3; i <= Nx + 3; i++)
    for (j = 3; j <= int(0.5 / dy) + 3; j++)
      for (k = 0; k <= 3; k++) U[i][j][k] = U[i][j][k] - r * (G_[i][j][k] - G_[i][j - 1][k]);

  for (i = int(1.0 / dx) + 3; i <= int(2.0 / dx) + 3; i++)
    for (j = int(0.5 / dy) + 4; j <= Ny + 3; j++)
      for (k = 0; k <= 3; k++) U[i][j][k] = U[i][j][k] - r * (G_[i][j][k] - G_[i][j - 1][k]);
}

#endif  // ENO_FUSION

#endif
