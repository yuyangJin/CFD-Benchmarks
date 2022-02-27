#ifndef __ENO_H
#define __ENO_H

#include "util.cuh"

#ifndef ENO_FUSION
//利用Roe平均以后的值计算x方向的特征值与特征向量
void LF_x(double U[Nx_new + 7][Ny_new + 7][4], double LAMDA_[Nx_new + 7][Ny_new + 7][4][4],
          double F[Nx_new + 7][Ny_new + 7][4], double Fp[Nx_new + 7][Ny_new + 7][4],
          double Fd[Nx_new + 7][Ny_new + 7][4]) {
  int i, j, k;
  double rou, u, v, p, a, maxlamda = 1e-100;
  //对LAMDA_赋值
  for (i = 0; i <= Nx_new + 6; i++)
    for (j = 0; j <= Ny_new + 6; j++)
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
  for (i = 0; i <= Nx_new + 6; i++)
    for (j = 0; j <= Ny_new + 6; j++)
      if (U[i][j][0] != 0) U2F(U[i][j], F[i][j]);

  //计算Fpd
  for (i = 0; i <= Nx_new + 6; i++)
    for (j = 0; j <= Ny_new + 6; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) {
          Fp[i][j][k] = 0.5 * (F[i][j][k] + maxlamda * U[i][j][k]);
          Fd[i][j][k] = 0.5 * (F[i][j][k] - maxlamda * U[i][j][k]);
        }
}

//此函数将之前计算的特征之余特征向量通过TVD算法计算U
void ENO_x(double U[Nx_new + 7][Ny_new + 7][4], double F[Nx_new + 7][Ny_new + 7][4],
           double Fp[Nx_new + 7][Ny_new + 7][4], double Fd[Nx_new + 7][Ny_new + 7][4],
           double F_p[Nx_new + 7][Ny_new + 7][4], double F_d[Nx_new + 7][Ny_new + 7][4],
           double F_[Nx_new + 7][Ny_new + 7][4], double q3p[Nx_new + 7][Ny_new + 7][4][3],
           double q3d[Nx_new + 7][Ny_new + 7][4][3], double dx, double dy, double dt) {
  int i, j, k;
  double r;
  dt = CFL(U, dx, dy, ENOCFL);
  r = dt / dx;

  //计算q3pd q5pd
  for (i = 2; i <= Nx_new + 3; i++)
    for (j = 2; j <= Ny_new + 3; j++)
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
  for (i = 2; i <= Nx_new + 3; i++)
    for (j = 2; j <= Ny_new + 3; j++)
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

  for (i = 2; i <= Nx_new + 3; i++)
    for (j = 2; j <= Ny_new + 3; j++)
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
  for (i = 2; i <= Nx_new + 3; i++)
    for (j = 2; j <= Ny_new + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) F_[i][j][k] = F_p[i][j][k] + F_d[i][j][k];
  //分区计算U
  for (i = 3; i <= Nx_new + 3; i++)
    for (j = 3; j <= int(0.5 / dy) + 3; j++)
      for (k = 0; k <= 3; k++) U[i][j][k] = U[i][j][k] - r * (F_[i][j][k] - F_[i - 1][j][k]);

  for (i = int(1.0 / dx) + 3; i <= int(2.0 / dx) + 3; i++)
    for (j = int(0.5 / dy) + 4; j <= Ny_new + 3; j++)
      for (k = 0; k <= 3; k++) U[i][j][k] = U[i][j][k] - r * (F_[i][j][k] - F_[i - 1][j][k]);
}

//利用Roe平均以后的值计算x方向的特征值与特征向量
void LF_y(double U[Nx_new + 7][Ny_new + 7][4], double LAMDA_[Nx_new + 7][Ny_new + 7][4][4],
          double G[Nx_new + 7][Ny_new + 7][4], double Gp[Nx_new + 7][Ny_new + 7][4],
          double Gd[Nx_new + 7][Ny_new + 7][4]) {
  int i, j, k;
  double rou, u, v, a, p, maxlamda = 1e-100;
  ;

  //对LAMDA_赋值
  for (i = 0; i <= Nx_new + 5; i++)
    for (j = 0; j <= Ny_new + 5; j++)
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
  for (i = 0; i <= Nx_new + 6; i++)
    for (j = 0; j <= Ny_new + 6; j++)
      if (U[i][j][0] != 0) U2G(U[i][j], G[i][j]);

  //计算Gpd
  for (i = 0; i <= Nx_new + 6; i++)
    for (j = 0; j <= Ny_new + 6; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) {
          Gp[i][j][k] = 0.5 * (G[i][j][k] + maxlamda * U[i][j][k]);
          Gd[i][j][k] = 0.5 * (G[i][j][k] - maxlamda * U[i][j][k]);
        }
}

//此函数将之前计算的特征之余特征向量通过TVD算法计算U
void ENO_y(double U[Nx_new + 7][Ny_new + 7][4], double G[Nx_new + 7][Ny_new + 7][4],
           double Gp[Nx_new + 7][Ny_new + 7][4], double Gd[Nx_new + 7][Ny_new + 7][4],
           double G_p[Nx_new + 7][Ny_new + 7][4], double G_d[Nx_new + 7][Ny_new + 7][4],
           double G_[Nx_new + 7][Ny_new + 7][4], double q3p[Nx_new + 7][Ny_new + 7][4][3],
           double q3d[Nx_new + 7][Ny_new + 7][4][3], double dx, double dy, double dt) {
  int i, j, k;
  double r;
  dt = CFL(U, dx, dy, ENOCFL);
  r = dt / dy;

  //计算q3pd q5pd
  for (i = 2; i <= Nx_new + 3; i++)
    for (j = 2; j <= Ny_new + 3; j++)
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
  for (i = 2; i <= Nx_new + 3; i++)
    for (j = 2; j <= Ny_new + 3; j++)
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

  for (i = 2; i <= Nx_new + 3; i++)
    for (j = 2; j <= Ny_new + 3; j++)
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
  for (i = 2; i <= Nx_new + 3; i++)
    for (j = 2; j <= Ny_new + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) G_[i][j][k] = G_p[i][j][k] + G_d[i][j][k];

  //分区计算U
  for (i = 3; i <= Nx_new + 3; i++)
    for (j = 3; j <= int(0.5 / dy) + 3; j++)
      for (k = 0; k <= 3; k++) U[i][j][k] = U[i][j][k] - r * (G_[i][j][k] - G_[i][j - 1][k]);

  for (i = int(1.0 / dx) + 3; i <= int(2.0 / dx) + 3; i++)
    for (j = int(0.5 / dy) + 4; j <= Ny_new + 3; j++)
      for (k = 0; k <= 3; k++) U[i][j][k] = U[i][j][k] - r * (G_[i][j][k] - G_[i][j - 1][k]);
}

#else

//利用Roe平均以后的值计算x方向的特征值与特征向量
void LF_x(double*** U, double**** LAMDA_, double*** F, double*** Fp, double*** Fd) {
  int i, j, k;
  double rou, u, v, p, a, maxlamda = 1e-100;
  //对LAMDA_赋值
  for (i = 0; i <= Nx_new + 6; i++)
    for (j = 0; j <= Ny_new + 6; j++)
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
  for (i = 0; i <= Nx_new + 6; i++)
    for (j = 0; j <= Ny_new + 6; j++)
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
void LF_x_Tensor(Tensor& U, Tensor& LAMDA_, Tensor& F, Tensor& Fp, Tensor& Fd) {
  int i, j, k;
  double rou, u, v, p, a, maxlamda = 1e-100;
  //对LAMDA_赋值
  for (i = 0; i <= Nx_new + 6; i++)
    for (j = 0; j <= Ny_new + 6; j++)
      if (U.indexing(i, j, 0) != 0) {
        rou = U.indexing(i, j, 0);
        u = U.indexing(i, j, 1) / U.indexing(i, j, 0);
        v = U.indexing(i, j, 2) / U.indexing(i, j, 0);
        p = (GAMA - 1) * (U.indexing(i, j, 3) - 0.5 * rou * (u * u + v * v));
        a = sqrt(p * GAMA / rou);
        LAMDA_.indexing(i, j, 0, 0) = u - a;
        LAMDA_.indexing(i, j, 1, 1) = u;
        LAMDA_.indexing(i, j, 2, 2) = u + a;
        LAMDA_.indexing(i, j, 3, 3) = u;
        for (k = 0; k <= 3; k++) {
          if (fabs(LAMDA_.indexing(i, j, k, k)) >= maxlamda) maxlamda = fabs(LAMDA_.indexing(i, j, k, k));
        }
      }
  for (i = 0; i <= Nx_new + 6; i++)
    for (j = 0; j <= Ny_new + 6; j++)
      if (U.indexing(i, j, 0) != 0) {
        // U2F
        U2F(U.indexing(i, j), F.indexing(i, j));
        //计算Fpd
        for (k = 0; k <= 3; k++) {
          Fp.indexing(i, j, k) = 0.5 * (F.indexing(i, j, k) + maxlamda * U.indexing(i, j, k));
          Fd.indexing(i, j, k) = 0.5 * (F.indexing(i, j, k) - maxlamda * U.indexing(i, j, k));
        }
      }
}

__global__ void eno_x_cuda(double* U, double* F, double* Fp, double* Fd, double* F_p, double* F_d, double* F_,
                           double* q3p, double* q3d, double dx, double dy, double dt, int offset1, int offset2,
                           int Nx_new, int Ny_new) {
  int i = blockIdx.x * blockDim.x + threadIdx.x + 2;
  int j = blockIdx.y * blockDim.y + threadIdx.y + 2;

  // for (i = 2; i <= Nx_new + 3; i++)
  //   for (j = 2; j <= Ny_new + 3; j++)
  if (i >= Nx_new + 2 || j >= Ny_new + 2) return;
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

// for (i = 3; i <= Nx_new + 3; i++)
//   for (j = 3; j <= int(0.5 / dy) + 3; j++)
//     for (k = 0; k <= 3; k++) U[i][j][k] = U[i][j][k] - r * (F_[i][j][k] - F_[i - 1][j][k]);

// for (i = int(1.0 / dx) + 3; i <= int(2.0 / dx) + 3; i++)
//   for (j = int(0.5 / dy) + 4; j <= Ny_new + 3; j++)
//     for (k = 0; k <= 3; k++) U[i][j][k] = U[i][j][k] - r * (F_[i][j][k] - F_[i - 1][j][k]);

__global__ void eno_x_cuda_2(double* U, double* F_, double dy, double r, int offset1, int offset2, int Nx_new,
                             int Ny_new) {
  // for (i = 3; i <= Nx_new + 3; i++)
  //   for (j = 3; j <= int(0.5 / dy) + 3; j++)
  int i = blockIdx.x * blockDim.x + threadIdx.x + 3;
  int j = blockIdx.y * blockDim.y + threadIdx.y + 3;
  if (i > Nx_new + 3 || j > int(0.5 / dy) + 3) return;
  for (int k = 0; k <= 3; k++)
    U[i * offset1 + j * offset2 + k] = U[i * offset1 + j * offset2 + k] - r * (F_[i * offset1 + j * offset2 + k] -
                                                                               F_[(i - 1) * offset1 + j * offset2 + k]);
}

__global__ void eno_x_cuda_3(double* U, double* F_, double dx, double dy, double r, int offset1, int offset2,
                             int Nx_new, int Ny_new) {
  // for (i = 3; i <= Nx_new + 3; i++)
  //   for (j = 3; j <= int(0.5 / dy) + 3; j++)
  int i = blockIdx.x * blockDim.x + threadIdx.x + 3 + int(1.f / dx);
  int j = blockIdx.y * blockDim.y + threadIdx.y + 4 + int(0.5 / dy);
  if (i > int(2.0 / dx) + 3 || j > Ny_new + 3) return;
  for (int k = 0; k <= 3; k++)
    U[i * offset1 + j * offset2 + k] = U[i * offset1 + j * offset2 + k] - r * (F_[i * offset1 + j * offset2 + k] -
                                                                               F_[(i - 1) * offset1 + j * offset2 + k]);
}

__global__ void eno_x_cuda_v2(double* U, double* F, double* Fp, double* Fd, double* F_p, double* F_d, double* F_,
                              double* q3p, double* q3d, double dx, double dy, double dt, int offset1, int offset2,
                              int Nx_new, int Ny_new) {
  int i = blockIdx.x * blockDim.x + threadIdx.x + 2;
  int j = blockIdx.y * blockDim.y + threadIdx.y + 2;
  int k = threadIdx.z;

  // for (i = 2; i <= Nx_new + 3; i++)
  //   for (j = 2; j <= Ny_new + 3; j++)
  if (i >= Nx_new + 2 || j >= Ny_new + 2) return;
  if (1)
  // if (U[i * offset1 + j * offset2] != 0)
  // for (int k = 0; k <= 3; k++)
  {
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

// for (i = 3; i <= Nx_new + 3; i++)
//   for (j = 3; j <= int(0.5 / dy) + 3; j++)
//     for (k = 0; k <= 3; k++) U[i][j][k] = U[i][j][k] - r * (F_[i][j][k] - F_[i - 1][j][k]);

// for (i = int(1.0 / dx) + 3; i <= int(2.0 / dx) + 3; i++)
//   for (j = int(0.5 / dy) + 4; j <= Ny_new + 3; j++)
//     for (k = 0; k <= 3; k++) U[i][j][k] = U[i][j][k] - r * (F_[i][j][k] - F_[i - 1][j][k]);

__global__ void eno_x_cuda_2_v2(double* U, double* F_, double dy, double r, int offset1, int offset2, int Nx_new,
                                int Ny_new) {
  // for (i = 3; i <= Nx_new + 3; i++)
  //   for (j = 3; j <= int(0.5 / dy) + 3; j++)
  int i = blockIdx.x * blockDim.x + threadIdx.x + 3;
  int j = blockIdx.y * blockDim.y + threadIdx.y + 3;
  int k = threadIdx.z;
  if (i > Nx_new + 3 || j > int(0.5 / dy) + 3) return;
  // for (int k = 0; k <= 3; k++)
  U[i * offset1 + j * offset2 + k] = U[i * offset1 + j * offset2 + k] -
                                     r * (F_[i * offset1 + j * offset2 + k] - F_[(i - 1) * offset1 + j * offset2 + k]);
}

__global__ void eno_x_cuda_3_v2(double* U, double* F_, double dx, double dy, double r, int offset1, int offset2,
                                int Nx_new, int Ny_new) {
  // for (i = 3; i <= Nx_new + 3; i++)
  //   for (j = 3; j <= int(0.5 / dy) + 3; j++)
  int i = blockIdx.x * blockDim.x + threadIdx.x + 3 + int(1.f / dx);
  int j = blockIdx.y * blockDim.y + threadIdx.y + 4 + int(0.5 / dy);
  int k = threadIdx.z;
  if (i > int(2.0 / dx) + 3 || j > Ny_new + 3) return;
  // for (int k = 0; k <= 3; k++)
  U[i * offset1 + j * offset2 + k] = U[i * offset1 + j * offset2 + k] -
                                     r * (F_[i * offset1 + j * offset2 + k] - F_[(i - 1) * offset1 + j * offset2 + k]);
}

//此函数将之前计算的特征之余特征向量通过TVD算法计算U
void ENO_x(double*** U, double*** F, double*** Fp, double*** Fd, double*** F_p, double*** F_d, double*** F_,
           double**** q3p, double**** q3d, double dx, double dy, double dt) {
  int i, j, k;
  double r;
  dt = CFL(U, dx, dy, ENOCFL);
  r = dt / dx;

  for (i = 2; i <= Nx_new + 3; i++)
    for (j = 2; j <= Ny_new + 3; j++)
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
  for (i = 3; i <= Nx_new + 3; i++)
    for (j = 3; j <= int(0.5 / dy) + 3; j++)
      for (k = 0; k <= 3; k++) U[i][j][k] = U[i][j][k] - r * (F_[i][j][k] - F_[i - 1][j][k]);

  for (i = int(1.0 / dx) + 3; i <= int(2.0 / dx) + 3; i++)
    for (j = int(0.5 / dy) + 4; j <= Ny_new + 3; j++)
      for (k = 0; k <= 3; k++) U[i][j][k] = U[i][j][k] - r * (F_[i][j][k] - F_[i - 1][j][k]);
}

void ENO_x_Tensor(Tensor& U, Tensor& F, Tensor& Fp, Tensor& Fd, Tensor& F_p, Tensor& F_d, Tensor& F_, Tensor& q3p,
                  Tensor& q3d, double dx, double dy, double dt) {
  int i, j, k;
  double r;
  dt = CFL_Tensor(U, dx, dy, ENOCFL);
  r = dt / dx;

  for (i = 2; i <= Nx_new + 3; i++)
    for (j = 2; j <= Ny_new + 3; j++)
      if (U.indexing(i, j, 0) != 0)
        for (k = 0; k <= 3; k++) {
          //计算q3pd q5pd
          /*
          q3p.indexing(i,j,k,0) = 1.0 / 3.0 * Fp.indexing(i - 2,j,k) - 7.0 / 6.0 * Fp.indexing(i - 1,j,k) + 11.0 / 6.0 *
          Fp.indexing(i,j,k); q3p.indexing(i,j,k,1) = -1.0 / 6.0 * Fp.indexing(i - 1,j,k) + 5.0 / 6.0 *
          Fp.indexing(i,j,k) + 1.0 / 3.0 * Fp.indexing(i + 1,j,k); q3p.indexing(i,j,k,2) = 1.0 / 3.0 *
          Fp.indexing(i,j,k) + 5.0 / 6.0 * Fp.indexing(i + 1,j,k) - 1.0 / 6.0 * Fp.indexing(i + 2,j,k);

          q3d.indexing(i,j,k,0) = -1.0 / 6.0 * Fd.indexing(i - 1,j,k) + 5.0 / 6.0 * Fd.indexing(i,j,k) + 1.0 / 3.0 *
          Fd.indexing(i + 1,j,k); q3d.indexing(i,j,k,1) = 1.0 / 3.0 * Fd.indexing(i,j,k) + 5.0 / 6.0 * Fd.indexing(i +
          1,j,k) - 1.0 / 6.0 * Fd.indexing(i + 2,j,k); q3d.indexing(i,j,k,2) = 11.0 / 6.0 * Fd.indexing(i + 1,j,k) - 7.0
          / 6.0 * Fd.indexing(i + 2,j,k) + 1.0 / 3.0 * Fd.indexing(i + 3,j,k);
          */

          //判断差商大小并且赋值
          if (fabs(U.indexing(i + 1, j, k) - U.indexing(i, j, k)) >
                  fabs(U.indexing(i, j, k) - U.indexing(i - 1, j, k)) &&
              fabs(U.indexing(i + 1, j, k) - U.indexing(i, j, k)) >
                  fabs(U.indexing(i - 1, j, k) - U.indexing(i - 2, j, k)))
            // F_p.indexing(i,j,k) = q3p.indexing(i,j,k,0);
            F_p.indexing(i, j, k) = 1.0 / 3.0 * Fp.indexing(i - 2, j, k) - 7.0 / 6.0 * Fp.indexing(i - 1, j, k) +
                                    11.0 / 6.0 * Fp.indexing(i, j, k);
          else {
            if (fabs(U.indexing(i, j, k) - U.indexing(i - 1, j, k)) >
                    fabs(U.indexing(i + 1, j, k) - U.indexing(i, j, k)) &&
                fabs(U.indexing(i, j, k) - U.indexing(i - 1, j, k)) >
                    fabs(U.indexing(i + 2, j, k) - U.indexing(i + 1, j, k)))
              // F_p.indexing(i,j,k) = q3p.indexing(i,j,k,2);
              F_p.indexing(i, j, k) = 1.0 / 3.0 * Fp.indexing(i, j, k) + 5.0 / 6.0 * Fp.indexing(i + 1, j, k) -
                                      1.0 / 6.0 * Fp.indexing(i + 2, j, k);
            else
              // F_p.indexing(i,j,k) = q3p.indexing(i,j,k,1);
              F_p.indexing(i, j, k) = -1.0 / 6.0 * Fp.indexing(i - 1, j, k) + 5.0 / 6.0 * Fp.indexing(i, j, k) +
                                      1.0 / 3.0 * Fp.indexing(i + 1, j, k);
          }
          if (fabs(U.indexing(i + 2, j, k) - U.indexing(i + 1, j, k)) >
                  fabs(U.indexing(i + 1, j, k) - U.indexing(i, j, k)) &&
              fabs(U.indexing(i + 2, j, k) - U.indexing(i + 1, j, k)) >
                  fabs(U.indexing(i, j, k) - U.indexing(i - 1, j, k)))
            // F_d.indexing(i,j,k) = q3d.indexing(i,j,k,0);
            F_d.indexing(i, j, k) = -1.0 / 6.0 * Fd.indexing(i - 1, j, k) + 5.0 / 6.0 * Fd.indexing(i, j, k) +
                                    1.0 / 3.0 * Fd.indexing(i + 1, j, k);
          else {
            if (fabs(U.indexing(i + 1, j, k) - U.indexing(i, j, k)) >
                    fabs(U.indexing(i + 2, j, k) - U.indexing(i + 1, j, k)) &&
                fabs(U.indexing(i + 1, j, k) - U.indexing(i, j, k)) >
                    fabs(U.indexing(i + 3, j, k) - U.indexing(i + 2, j, k)))
              // F_d.indexing(i,j,k) = q3d.indexing(i,j,k,2);
              F_d.indexing(i, j, k) = 11.0 / 6.0 * Fd.indexing(i + 1, j, k) - 7.0 / 6.0 * Fd.indexing(i + 2, j, k) +
                                      1.0 / 3.0 * Fd.indexing(i + 3, j, k);
            else
              // F_d.indexing(i,j,k) = q3d.indexing(i,j,k,1);
              F_d.indexing(i, j, k) = 1.0 / 3.0 * Fd.indexing(i, j, k) + 5.0 / 6.0 * Fd.indexing(i + 1, j, k) -
                                      1.0 / 6.0 * Fd.indexing(i + 2, j, k);
          }

          // 计算F_
          F_.indexing(i, j, k) = F_p.indexing(i, j, k) + F_d.indexing(i, j, k);
        }

  //分区计算U
  for (i = 3; i <= Nx_new + 3; i++)
    for (j = 3; j <= int(0.5 / dy) + 3; j++)
      for (k = 0; k <= 3; k++)
        U.indexing(i, j, k) = U.indexing(i, j, k) - r * (F_.indexing(i, j, k) - F_.indexing(i - 1, j, k));

  for (i = int(1.0 / dx) + 3; i <= int(2.0 / dx) + 3; i++)
    for (j = int(0.5 / dy) + 4; j <= Ny_new + 3; j++)
      for (k = 0; k <= 3; k++)
        U.indexing(i, j, k) = U.indexing(i, j, k) - r * (F_.indexing(i, j, k) - F_.indexing(i - 1, j, k));
}

//利用Roe平均以后的值计算x方向的特征值与特征向量
void LF_y(double*** U, double**** LAMDA_, double*** G, double*** Gp, double*** Gd) {
  int i, j, k;
  double rou, u, v, a, p, maxlamda = 1e-100;
  ;

  //对LAMDA_赋值
  for (i = 0; i <= Nx_new + 5; i++)
    for (j = 0; j <= Ny_new + 5; j++)
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

  for (i = 0; i <= Nx_new + 6; i++)
    for (j = 0; j <= Ny_new + 6; j++)
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

  for (i = 2; i <= Nx_new + 3; i++)
    for (j = 2; j <= Ny_new + 3; j++)
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
  for (i = 3; i <= Nx_new + 3; i++)
    for (j = 3; j <= int(0.5 / dy) + 3; j++)
      for (k = 0; k <= 3; k++) U[i][j][k] = U[i][j][k] - r * (G_[i][j][k] - G_[i][j - 1][k]);

  for (i = int(1.0 / dx) + 3; i <= int(2.0 / dx) + 3; i++)
    for (j = int(0.5 / dy) + 4; j <= Ny_new + 3; j++)
      for (k = 0; k <= 3; k++) U[i][j][k] = U[i][j][k] - r * (G_[i][j][k] - G_[i][j - 1][k]);
}

//利用Roe平均以后的值计算x方向的特征值与特征向量
void LF_y_Tensor(Tensor& U, Tensor& LAMDA_, Tensor& G, Tensor& Gp, Tensor& Gd) {
  int i, j, k;
  double rou, u, v, a, p, maxlamda = 1e-100;
  ;

  //对LAMDA_赋值
  for (i = 0; i <= Nx_new + 5; i++)
    for (j = 0; j <= Ny_new + 5; j++)
      if (U.indexing(i, j, 0) != 0) {
        rou = U.indexing(i, j, 0);
        u = U.indexing(i, j, 1) / U.indexing(i, j, 0);
        v = U.indexing(i, j, 2) / U.indexing(i, j, 0);
        p = (GAMA - 1) * (U.indexing(i, j, 3) - 0.5 * rou * (u * u + v * v));
        // a=pow(p*GAMA/rou,0.5);
        a = sqrt(p * GAMA / rou);
        LAMDA_.indexing(i, j, 0, 0) = v - a;
        LAMDA_.indexing(i, j, 1, 1) = v;
        LAMDA_.indexing(i, j, 2, 2) = v + a;
        LAMDA_.indexing(i, j, 3, 3) = v;
        for (k = 0; k <= 3; k++) {
          if (fabs(LAMDA_.indexing(i, j, k, k)) >= maxlamda) maxlamda = fabs(LAMDA_.indexing(i, j, k, k));
        }
      }

  for (i = 0; i <= Nx_new + 6; i++)
    for (j = 0; j <= Ny_new + 6; j++)
      if (U.indexing(i, j, 0) != 0) {
        // U2G
        U2G(U.indexing(i, j), G.indexing(i, j));
        //计算Gpd
        for (k = 0; k <= 3; k++) {
          Gp.indexing(i, j, k) = 0.5 * (G.indexing(i, j, k) + maxlamda * U.indexing(i, j, k));
          Gd.indexing(i, j, k) = 0.5 * (G.indexing(i, j, k) - maxlamda * U.indexing(i, j, k));
        }
      }
}

//此函数将之前计算的特征之余特征向量通过TVD算法计算U
void ENO_y_Tensor(Tensor& U, Tensor& G, Tensor& Gp, Tensor& Gd, Tensor& G_p, Tensor& G_d, Tensor& G_, Tensor& q3p,
                  Tensor& q3d, double dx, double dy, double dt) {
  int i, j, k;
  double r;
  dt = CFL_Tensor(U, dx, dy, ENOCFL);
  r = dt / dy;

  for (i = 2; i <= Nx_new + 3; i++)
    for (j = 2; j <= Ny_new + 3; j++)
      if (U.indexing(i, j, 0) != 0) {
        for (k = 0; k <= 3; k++) {
          //计算q3pd q5pd
          /*
          q3p.indexing(i,j,k,0) = 1.0 / 3.0 * Gp.indexing(i,j - 2,k) - 7.0 / 6.0 * Gp.indexing(i,j - 1,k) + 11.0 / 6.0 *
          Gp.indexing(i,j,k); q3p.indexing(i,j,k,1) = -1.0 / 6.0 * Gp.indexing(i,j - 1,k) + 5.0 / 6.0 *
          Gp.indexing(i,j,k) + 1.0 / 3.0 * Gp.indexing(i,j + 1,k); q3p.indexing(i,j,k,2) = 1.0 / 3.0 *
          Gp.indexing(i,j,k) + 5.0 / 6.0 * Gp.indexing(i,j + 1,k) - 1.0 / 6.0 * Gp.indexing(i,j + 2,k);

          q3d.indexing(i,j,k,0) = -1.0 / 6.0 * Gd.indexing(i,j - 1,k) + 5.0 / 6.0 * Gd.indexing(i,j,k) + 1.0 / 3.0 *
          Gd.indexing(i,j + 1,k); q3d.indexing(i,j,k,1) = 1.0 / 3.0 * Gd.indexing(i,j,k) + 5.0 / 6.0 * Gd.indexing(i,j +
          1,k) - 1.0 / 6.0 * Gd.indexing(i,j + 2,k); q3d.indexing(i,j,k,2) = 11.0 / 6.0 * Gd.indexing(i,j + 1,k) - 7.0
          / 6.0 * Gd.indexing(i,j + 2,k) + 1.0 / 3.0 * Gd.indexing(i,j + 3,k);
          */

          //判断差商大小并且赋值
          if (fabs(U.indexing(i, j + 1, k) - U.indexing(i, j, k)) >
                  fabs(U.indexing(i, j, k) - U.indexing(i, j - 1, k)) &&
              fabs(U.indexing(i, j + 1, k) - U.indexing(i, j, k)) >
                  fabs(U.indexing(i, j - 1, k) - U.indexing(i, j - 2, k)))
            // G_p.indexing(i,j,k) = q3p.indexing(i,j,k,0);
            G_p.indexing(i, j, k) = 1.0 / 3.0 * Gp.indexing(i, j - 2, k) - 7.0 / 6.0 * Gp.indexing(i, j - 1, k) +
                                    11.0 / 6.0 * Gp.indexing(i, j, k);
          else {
            if (fabs(U.indexing(i, j, k) - U.indexing(i, j - 1, k)) >
                    fabs(U.indexing(i, j + 1, k) - U.indexing(i, j, k)) &&
                fabs(U.indexing(i, j, k) - U.indexing(i, j - 1, k)) >
                    fabs(U.indexing(i, j + 2, k) - U.indexing(i, j + 1, k)))
              // G_p.indexing(i,j,k) = q3p.indexing(i,j,k,2);
              G_p.indexing(i, j, k) = 1.0 / 3.0 * Gp.indexing(i, j, k) + 5.0 / 6.0 * Gp.indexing(i, j + 1, k) -
                                      1.0 / 6.0 * Gp.indexing(i, j + 2, k);
            else
              // G_p.indexing(i,j,k) = q3p.indexing(i,j,k,1);
              G_p.indexing(i, j, k) = -1.0 / 6.0 * Gp.indexing(i, j - 1, k) + 5.0 / 6.0 * Gp.indexing(i, j, k) +
                                      1.0 / 3.0 * Gp.indexing(i, j + 1, k);
          }
          if (fabs(U.indexing(i, j + 2, k) - U.indexing(i, j + 1, k)) >
                  fabs(U.indexing(i, j + 1, k) - U.indexing(i, j, k)) &&
              fabs(U.indexing(i, j + 2, k) - U.indexing(i, j + 1, k)) >
                  fabs(U.indexing(i, j, k) - U.indexing(i, j - 1, k)))
            // G_d.indexing(i,j,k) = q3d.indexing(i,j,k,0);
            G_d.indexing(i, j, k) = -1.0 / 6.0 * Gd.indexing(i, j - 1, k) + 5.0 / 6.0 * Gd.indexing(i, j, k) +
                                    1.0 / 3.0 * Gd.indexing(i, j + 1, k);
          else {
            if (fabs(U.indexing(i, j + 1, k) - U.indexing(i, j, k)) >
                    fabs(U.indexing(i, j + 2, k) - U.indexing(i, j + 1, k)) &&
                fabs(U.indexing(i, j + 1, k) - U.indexing(i, j, k)) >
                    fabs(U.indexing(i, j + 3, k) - U.indexing(i, j + 2, k)))
              // G_d.indexing(i,j,k) = q3d.indexing(i,j,k,2);
              G_d.indexing(i, j, k) = 11.0 / 6.0 * Gd.indexing(i, j + 1, k) - 7.0 / 6.0 * Gd.indexing(i, j + 2, k) +
                                      1.0 / 3.0 * Gd.indexing(i, j + 3, k);
            else
              // G_d.indexing(i,j,k) = q3d.indexing(i,j,k,1);
              G_d.indexing(i, j, k) = 1.0 / 3.0 * Gd.indexing(i, j, k) + 5.0 / 6.0 * Gd.indexing(i, j + 1, k) -
                                      1.0 / 6.0 * Gd.indexing(i, j + 2, k);
          }

          //计算G_
          G_.indexing(i, j, k) = G_p.indexing(i, j, k) + G_d.indexing(i, j, k);
        }
      }

  //分区计算U
  for (i = 3; i <= Nx_new + 3; i++)
    for (j = 3; j <= int(0.5 / dy) + 3; j++)
      for (k = 0; k <= 3; k++)
        U.indexing(i, j, k) = U.indexing(i, j, k) - r * (G_.indexing(i, j, k) - G_.indexing(i, j - 1, k));

  for (i = int(1.0 / dx) + 3; i <= int(2.0 / dx) + 3; i++)
    for (j = int(0.5 / dy) + 4; j <= Ny_new + 3; j++)
      for (k = 0; k <= 3; k++)
        U.indexing(i, j, k) = U.indexing(i, j, k) - r * (G_.indexing(i, j, k) - G_.indexing(i, j - 1, k));
}

#endif  // ENO_FUSION

#endif
