#ifndef __WENO_H
#define __WENO_H

#include "util.h"

#define P 2  //模板不光滑性的放大系数
#define OPDT 4

#ifndef WENO_FUSION

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

//利用Roe平均以后的值计算x方向的特征值与特征向量
void LF_y(double U[Nx + 7][Ny + 7][4], double LAMDA_[Nx + 7][Ny + 7][4][4], double G[Nx + 7][Ny + 7][4],
          double Gp[Nx + 7][Ny + 7][4], double Gd[Nx + 7][Ny + 7][4]) {
  int i, j, k;
  double rou, u, v, a, p, maxlamda = 1e-100;

  //对LAMDA_赋值
  for (i = 0; i <= Nx + 5; i++)
    for (j = 0; j <= Ny + 5; j++)
      if (U[i][j][0] != 0) {
        rou = U[i][j][0];
        u = U[i][j][1] / U[i][j][0];
        v = U[i][j][2] / U[i][j][0];
        p = (GAMA - 1) * (U[i][j][3] - 0.5 * rou * (u * u + v * v));
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

void WENO_x(double U[Nx + 7][Ny + 7][4], double ISp[Nx + 7][Ny + 7][4][3], double ISd[Nx + 7][Ny + 7][4][3],
            double omegap[Nx + 7][Ny + 7][4][3], double omegad[Nx + 7][Ny + 7][4][3],
            double alphap[Nx + 7][Ny + 7][4][3], double alphad[Nx + 7][Ny + 7][4][3], double q3p[Nx + 7][Ny + 7][4][3],
            double q3d[Nx + 7][Ny + 7][4][3], double Fp[Nx + 7][Ny + 7][4], double Fd[Nx + 7][Ny + 7][4],
            double F_p[Nx + 7][Ny + 7][4], double F_d[Nx + 7][Ny + 7][4], double F_[Nx + 7][Ny + 7][4], double dx,
            double dy, double dt)  //此函数将之前计算的特征之余特征向量通过TVD算法计算U
{
  int i, j, k, l;
  dt = CFL(U, dx, dy, WENOCFL);
  double r = dt / dx;
  //计算ISpd
  for (i = 2; i <= Nx + 3; i++)
    for (j = 2; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) {
          double isp0_l = Fp[i - 2][j][k] - 2 * Fp[i - 1][j][k] + Fp[i][j][k];
          double isp0_r = Fp[i - 2][j][k] - 4 * Fp[i - 1][j][k] + 3 * Fp[i][j][k];
          double isp1_l = Fp[i - 1][j][k] - 2 * Fp[i][j][k] + Fp[i + 1][j][k];
          double isp1_r = Fp[i - 1][j][k] - Fp[i + 1][j][k];
          double isp2_l = Fp[i][j][k] - 2 * Fp[i + 1][j][k] + Fp[i + 2][j][k];
          double isp2_r = 3 * Fp[i][j][k] - 4 * Fp[i + 1][j][k] + Fp[i + 2][j][k];

          ISp[i][j][k][0] = 13.0 / 12.0 * isp0_l * isp0_l + 0.25 * isp0_r * isp0_r;
          ISp[i][j][k][1] = 13.0 / 12.0 * isp1_l * isp1_l + 0.25 * isp1_r * isp1_r;
          ISp[i][j][k][2] = 13.0 / 12.0 * isp2_l * isp2_l + 0.25 * isp2_r * isp2_r;

          double isd0_l = Fd[i + 1][j][k] - 2 * Fd[i][j][k] + Fd[i - 1][j][k];
          double isd0_r = 3 * Fd[i + 1][j][k] - 4 * Fd[i][j][k] + Fd[i - 1][j][k];
          double isd1_l = Fd[i + 2][j][k] - 2 * Fd[i + 1][j][k] + Fd[i][j][k];
          double isd1_r = Fd[i + 2][j][k] - Fd[i][j][k];
          double isd2_l = Fd[i + 3][j][k] - 2 * Fd[i + 2][j][k] + Fd[i + 1][j][k];
          double isd2_r = Fd[i + 3][j][k] - 4 * Fd[i + 2][j][k] + 3 * Fd[i + 1][j][k];

          ISd[i][j][k][0] = 13.0 / 12.0 * isd0_l * isd0_l + 0.25 * isd0_r * isd0_r;
          ISd[i][j][k][1] = 13.0 / 12.0 * isd1_l * isd1_l + 0.25 * isd1_r * isd1_r;
          ISd[i][j][k][2] = 13.0 / 12.0 * isd2_l * isd2_l + 0.25 * isd2_r * isd2_r;
        }
  //计算alphapd以及omegapd
  for (i = 2; i <= Nx + 3; i++)
    for (j = 2; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) {
          alphap[i][j][k][0] = 0.1 / pow(1e-16 + ISp[i][j][k][0], P);
          alphap[i][j][k][1] = 0.6 / pow(1e-16 + ISp[i][j][k][1], P);
          alphap[i][j][k][2] = 0.3 / pow(1e-16 + ISp[i][j][k][2], P);

          alphad[i][j][k][0] = 0.3 / pow(1e-16 + ISd[i][j][k][0], P);
          alphad[i][j][k][1] = 0.6 / pow(1e-16 + ISd[i][j][k][1], P);
          alphad[i][j][k][2] = 0.1 / pow(1e-16 + ISd[i][j][k][2], P);
          for (l = 0; l <= 2; l++) {
            omegap[i][j][k][l] = alphap[i][j][k][l] / (alphap[i][j][k][0] + alphap[i][j][k][1] + alphap[i][j][k][2]);
            omegad[i][j][k][l] = alphad[i][j][k][l] / (alphad[i][j][k][0] + alphad[i][j][k][1] + alphad[i][j][k][2]);
          }
        }
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

  //计算LF_pd
  for (i = 2; i <= Nx + 3; i++)
    for (j = 2; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) {
          F_p[i][j][k] = 0;
          F_d[i][j][k] = 0;
        }

  for (i = 2; i <= Nx + 3; i++)
    for (j = 2; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++)
          for (l = 0; l <= 2; l++) {
            F_p[i][j][k] += omegap[i][j][k][l] * q3p[i][j][k][l];
            F_d[i][j][k] += omegad[i][j][k][l] * q3d[i][j][k][l];
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

void WENO_y(double U[Nx + 7][Ny + 7][4], double ISp[Nx + 7][Ny + 7][4][3], double ISd[Nx + 7][Ny + 7][4][3],
            double omegap[Nx + 7][Ny + 7][4][3], double omegad[Nx + 7][Ny + 7][4][3],
            double alphap[Nx + 7][Ny + 7][4][3], double alphad[Nx + 7][Ny + 7][4][3], double q3p[Nx + 7][Ny + 7][4][3],
            double q3d[Nx + 7][Ny + 7][4][3], double Gp[Nx + 7][Ny + 7][4], double Gd[Nx + 7][Ny + 7][4],
            double G_p[Nx + 7][Ny + 7][4], double G_d[Nx + 7][Ny + 7][4], double G_[Nx + 7][Ny + 7][4], double dx,
            double dy, double dt)  //此函数将之前计算的特征之余特征向量通过TVD算法计算U
{
  int i, j, k, l;
  dt = CFL(U, dx, dy, WENOCFL);
  double r = dt / dy;
  //计算ISpd
  for (i = 2; i <= Nx + 3; i++)
    for (j = 2; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) {
          double isp0_l = Gp[i][j - 2][k] - 2 * Gp[i][j - 1][k] + Gp[i][j][k];
          double isp0_r = Gp[i][j - 2][k] - 4 * Gp[i][j - 1][k] + 3 * Gp[i][j][k];
          double isp1_l = Gp[i][j - 1][k] - 2 * Gp[i][j][k] + Gp[i][j + 1][k];
          double isp1_r = Gp[i][j - 1][k] - Gp[i][j + 1][k];
          double isp2_l = Gp[i][j][k] - 2 * Gp[i][j + 1][k] + Gp[i][j + 2][k];
          double isp2_r = 3 * Gp[i][j][k] - 4 * Gp[i][j + 1][k] + Gp[i][j + 2][k];

          ISp[i][j][k][0] = 13.0 / 12.0 * isp0_l * isp0_l + 0.25 * isp0_r * isp0_r;
          ISp[i][j][k][1] = 13.0 / 12.0 * isp1_l * isp1_l + 0.25 * isp1_r * isp1_r;
          ISp[i][j][k][2] = 13.0 / 12.0 * isp2_l * isp2_l + 0.25 * isp2_r * isp2_r;

          double isd0_l = Gd[i][j + 1][k] - 2 * Gd[i][j][k] + Gd[i][j - 1][k];
          double isd0_r = 3 * Gd[i][j + 1][k] - 4 * Gd[i][j][k] + Gd[i][j - 1][k];
          double isd1_l = Gd[i][j + 2][k] - 2 * Gd[i][j + 1][k] + Gd[i][j][k];
          double isd1_r = Gd[i][j + 2][k] - Gd[i][j][k];
          double isd2_l = Gd[i][j + 3][k] - 2 * Gd[i][j + 2][k] + Gd[i][j + 1][k];
          double isd2_r = Gd[i][j + 3][k] - 4 * Gd[i][j + 2][k] + 3 * Gd[i][j + 1][k];

          ISd[i][j][k][0] = 13.0 / 12.0 * isd0_l * isd0_l + 0.25 * isd0_r * isd0_r;
          ISd[i][j][k][1] = 13.0 / 12.0 * isd1_l * isd1_l + 0.25 * isd1_r * isd1_r;
          ISd[i][j][k][2] = 13.0 / 12.0 * isd2_l * isd2_l + 0.25 * isd2_r * isd2_r;
        }

  //计算alphapd以及omegapd
  for (i = 2; i <= Nx + 3; i++)
    for (j = 2; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) {
          alphap[i][j][k][0] = 0.1 / pow(1e-16 + ISp[i][j][k][0], P);
          alphap[i][j][k][1] = 0.6 / pow(1e-16 + ISp[i][j][k][1], P);
          alphap[i][j][k][2] = 0.3 / pow(1e-16 + ISp[i][j][k][2], P);

          alphad[i][j][k][0] = 0.3 / pow(1e-16 + ISd[i][j][k][0], P);
          alphad[i][j][k][1] = 0.6 / pow(1e-16 + ISd[i][j][k][1], P);
          alphad[i][j][k][2] = 0.1 / pow(1e-16 + ISd[i][j][k][2], P);
          for (l = 0; l <= 2; l++) {
            omegap[i][j][k][l] = alphap[i][j][k][l] / (alphap[i][j][k][0] + alphap[i][j][k][1] + alphap[i][j][k][2]);
            omegad[i][j][k][l] = alphad[i][j][k][l] / (alphad[i][j][k][0] + alphad[i][j][k][1] + alphad[i][j][k][2]);
          }
        }

  //计算q3pd q5pd
  for (i = 2; i <= Nx + 3; i++)
    for (j = 2; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) {
          q3p[i][j][k][0] = 1.0 / 3.0 * Gp[i][j - 2][k] - 7.0 / 6.0 * Gp[i][j - 1][k] + 11.0 / 6.0 * Gp[i][j][k];
          q3p[i][j][k][1] = -1.0 / 6.0 * Gp[i][j - 1][k] + 5.0 / 6.0 * Gp[i][j][k] + 1.0 / 3.0 * Gp[i][j + 1][k];
          q3p[i][j][k][2] = 1.0 / 3.0 * Gp[i][j][k] + 5.0 / 6.0 * Gp[i][j + 1][k] - 1.0 / 6.0 * Gp[i][j + 2][k];

          q3d[i][j][k][0] = -1.0 / 6.0 * Gd[i][j - 1][k] + 5.0 / 6.0 * Gd[i][j][k] + 1.0 / 3.0 * Gd[i][j + 1][k];
          q3d[i][j][k][1] = 1.0 / 3.0 * Gd[i][j][k] + 5.0 / 6.0 * Gd[i][j + 1][k] - 1.0 / 6.0 * Gd[i][j + 2][k];
          q3d[i][j][k][2] = 11.0 / 6.0 * Gd[i][j + 1][k] - 7.0 / 6.0 * Gd[i][j + 2][k] + 1.0 / 3.0 * Gd[i][j + 3][k];
        }

  //计算LG_pd
  for (i = 2; i <= Nx + 3; i++)
    for (j = 2; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) {
          G_p[i][j][k] = 0;
          G_d[i][j][k] = 0;
        }

  for (i = 2; i <= Nx + 3; i++)
    for (j = 2; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++)
          for (l = 0; l <= 2; l++) {
            G_p[i][j][k] += omegap[i][j][k][l] * q3p[i][j][k][l];
            G_d[i][j][k] += omegad[i][j][k][l] * q3d[i][j][k][l];
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

  for (i = 0; i <= Nx + 6; i++)
    for (j = 0; j <= Ny + 6; j++)
      if (U[i][j][0] != 0) {
        U2F(U[i][j], F[i][j]);
        for (k = 0; k <= 3; k++) {
          Fp[i][j][k] = 0.5 * (F[i][j][k] + maxlamda * U[i][j][k]);
          Fd[i][j][k] = 0.5 * (F[i][j][k] - maxlamda * U[i][j][k]);
        }
      }
}

//利用Roe平均以后的值计算x方向的特征值与特征向量
void LF_y(double U[Nx + 7][Ny + 7][4], double LAMDA_[Nx + 7][Ny + 7][4][4], double G[Nx + 7][Ny + 7][4],
          double Gp[Nx + 7][Ny + 7][4], double Gd[Nx + 7][Ny + 7][4]) {
  int i, j, k;
  double rou, u, v, a, p, maxlamda = 1e-100;

  //对LAMDA_赋值
  for (i = 0; i <= Nx + 5; i++)
    for (j = 0; j <= Ny + 5; j++)
      if (U[i][j][0] != 0) {
        rou = U[i][j][0];
        u = U[i][j][1] / U[i][j][0];
        v = U[i][j][2] / U[i][j][0];
        p = (GAMA - 1) * (U[i][j][3] - 0.5 * rou * (u * u + v * v));
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
        U2G(U[i][j], G[i][j]);
        for (k = 0; k <= 3; k++) {
          Gp[i][j][k] = 0.5 * (G[i][j][k] + maxlamda * U[i][j][k]);
          Gd[i][j][k] = 0.5 * (G[i][j][k] - maxlamda * U[i][j][k]);
        }
      }
}

void WENO_x(double U[Nx + 7][Ny + 7][4], double ISp[Nx + 7][Ny + 7][4][3], double ISd[Nx + 7][Ny + 7][4][3],
            double omegap[Nx + 7][Ny + 7][4][3], double omegad[Nx + 7][Ny + 7][4][3],
            double alphap[Nx + 7][Ny + 7][4][3], double alphad[Nx + 7][Ny + 7][4][3], double q3p[Nx + 7][Ny + 7][4][3],
            double q3d[Nx + 7][Ny + 7][4][3], double Fp[Nx + 7][Ny + 7][4], double Fd[Nx + 7][Ny + 7][4],
            double F_p[Nx + 7][Ny + 7][4], double F_d[Nx + 7][Ny + 7][4], double F_[Nx + 7][Ny + 7][4], double dx,
            double dy, double dt)  //此函数将之前计算的特征之余特征向量通过TVD算法计算U
{
  int i, j, k, l;
  dt = CFL(U, dx, dy, WENOCFL);
  double r = dt / dx;
  //计算ISpd
  for (i = 2; i <= Nx + 3; i++)
    for (j = 2; j <= Ny + 3; j++)
      if (U[i][j][0] != 0) {
        for (k = 0; k <= 3; k++) {
          double isp0_l = Fp[i - 2][j][k] - 2 * Fp[i - 1][j][k] + Fp[i][j][k];
          double isp0_r = Fp[i - 2][j][k] - 4 * Fp[i - 1][j][k] + 3 * Fp[i][j][k];
          double isp1_l = Fp[i - 1][j][k] - 2 * Fp[i][j][k] + Fp[i + 1][j][k];
          double isp1_r = Fp[i - 1][j][k] - Fp[i + 1][j][k];
          double isp2_l = Fp[i][j][k] - 2 * Fp[i + 1][j][k] + Fp[i + 2][j][k];
          double isp2_r = 3 * Fp[i][j][k] - 4 * Fp[i + 1][j][k] + Fp[i + 2][j][k];

          double ISp0 = 13.0 / 12.0 * isp0_l * isp0_l + 0.25 * isp0_r * isp0_r;
          double ISp1 = 13.0 / 12.0 * isp1_l * isp1_l + 0.25 * isp1_r * isp1_r;
          double ISp2 = 13.0 / 12.0 * isp2_l * isp2_l + 0.25 * isp2_r * isp2_r;

          double isd0_l = Fd[i + 1][j][k] - 2 * Fd[i][j][k] + Fd[i - 1][j][k];
          double isd0_r = 3 * Fd[i + 1][j][k] - 4 * Fd[i][j][k] + Fd[i - 1][j][k];
          double isd1_l = Fd[i + 2][j][k] - 2 * Fd[i + 1][j][k] + Fd[i][j][k];
          double isd1_r = Fd[i + 2][j][k] - Fd[i][j][k];
          double isd2_l = Fd[i + 3][j][k] - 2 * Fd[i + 2][j][k] + Fd[i + 1][j][k];
          double isd2_r = Fd[i + 3][j][k] - 4 * Fd[i + 2][j][k] + 3 * Fd[i + 1][j][k];

          double ISd0 = 13.0 / 12.0 * isd0_l * isd0_l + 0.25 * isd0_r * isd0_r;
          double ISd1 = 13.0 / 12.0 * isd1_l * isd1_l + 0.25 * isd1_r * isd1_r;
          double ISd2 = 13.0 / 12.0 * isd2_l * isd2_l + 0.25 * isd2_r * isd2_r;


          double local_alphap[3];
          double local_alphad[3];
          local_alphap[0] = 0.1 / pow(1e-16 + ISp0, P);
          local_alphap[1] = 0.6 / pow(1e-16 + ISp1, P);
          local_alphap[2] = 0.3 / pow(1e-16 + ISp2, P);

          local_alphad[0] = 0.3 / pow(1e-16 + ISd0, P);
          local_alphad[1] = 0.6 / pow(1e-16 + ISd1, P);
          local_alphad[2] = 0.1 / pow(1e-16 + ISd2, P);
          double local_omegap[3];
          double local_omegad[3];
          for (l = 0; l <= 2; l++) {
            local_omegap[l] = local_alphap[l] / (local_alphap[0] + local_alphap[1] + local_alphap[2]);
            local_omegad[l] = local_alphad[l] / (local_alphad[0] + local_alphad[1] + local_alphad[2]);
          }

          double local_q3p[3];
          double local_q3d[3];

          local_q3p[0] = 1.0 / 3.0 * Fp[i - 2][j][k] - 7.0 / 6.0 * Fp[i - 1][j][k] + 11.0 / 6.0 * Fp[i][j][k];
          local_q3p[1] = -1.0 / 6.0 * Fp[i - 1][j][k] + 5.0 / 6.0 * Fp[i][j][k] + 1.0 / 3.0 * Fp[i + 1][j][k];
          local_q3p[2] = 1.0 / 3.0 * Fp[i][j][k] + 5.0 / 6.0 * Fp[i + 1][j][k] - 1.0 / 6.0 * Fp[i + 2][j][k];

          local_q3d[0] = -1.0 / 6.0 * Fd[i - 1][j][k] + 5.0 / 6.0 * Fd[i][j][k] + 1.0 / 3.0 * Fd[i + 1][j][k];
          local_q3d[1] = 1.0 / 3.0 * Fd[i][j][k] + 5.0 / 6.0 * Fd[i + 1][j][k] - 1.0 / 6.0 * Fd[i + 2][j][k];
          local_q3d[2] = 11.0 / 6.0 * Fd[i + 1][j][k] - 7.0 / 6.0 * Fd[i + 2][j][k] + 1.0 / 3.0 * Fd[i + 3][j][k];


          double _F_p = 0;
          double _F_d = 0;


          for (l = 0; l <= 2; l++) {
            _F_p += local_omegap[l] * local_q3p[l];
            _F_d += local_omegad[l] * local_q3d[l];
          }

          F_[i][j][k] = _F_p + _F_d;
        }
      }

  //分区计算U
  for (i = 3; i <= Nx + 3; i++)
    for (j = 3; j <= int(0.5 / dy) + 3; j++)
      for (k = 0; k <= 3; k++) U[i][j][k] = U[i][j][k] - r * (F_[i][j][k] - F_[i - 1][j][k]);

  for (i = int(1.0 / dx) + 3; i <= int(2.0 / dx) + 3; i++)
    for (j = int(0.5 / dy) + 4; j <= Ny + 3; j++)
      for (k = 0; k <= 3; k++) U[i][j][k] = U[i][j][k] - r * (F_[i][j][k] - F_[i - 1][j][k]);
}

void WENO_y(double U[Nx + 7][Ny + 7][4], double ISp[Nx + 7][Ny + 7][4][3], double ISd[Nx + 7][Ny + 7][4][3],
            double omegap[Nx + 7][Ny + 7][4][3], double omegad[Nx + 7][Ny + 7][4][3],
            double alphap[Nx + 7][Ny + 7][4][3], double alphad[Nx + 7][Ny + 7][4][3], double q3p[Nx + 7][Ny + 7][4][3],
            double q3d[Nx + 7][Ny + 7][4][3], double Gp[Nx + 7][Ny + 7][4], double Gd[Nx + 7][Ny + 7][4],
            double G_p[Nx + 7][Ny + 7][4], double G_d[Nx + 7][Ny + 7][4], double G_[Nx + 7][Ny + 7][4], double dx,
            double dy, double dt)  //此函数将之前计算的特征之余特征向量通过TVD算法计算U
{
  int i, j, k, l;
  dt = CFL(U, dx, dy, WENOCFL);
  double r = dt / dy;
  //计算ISpd
  for (i = 2; i <= Nx + 3; i++)
    for (j = 2; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) {
          double isp0_l = Gp[i][j - 2][k] - 2 * Gp[i][j - 1][k] + Gp[i][j][k];
          double isp0_r = Gp[i][j - 2][k] - 4 * Gp[i][j - 1][k] + 3 * Gp[i][j][k];
          double isp1_l = Gp[i][j - 1][k] - 2 * Gp[i][j][k] + Gp[i][j + 1][k];
          double isp1_r = Gp[i][j - 1][k] - Gp[i][j + 1][k];
          double isp2_l = Gp[i][j][k] - 2 * Gp[i][j + 1][k] + Gp[i][j + 2][k];
          double isp2_r = 3 * Gp[i][j][k] - 4 * Gp[i][j + 1][k] + Gp[i][j + 2][k];

          double ISp0 = 13.0 / 12.0 * isp0_l * isp0_l + 0.25 * isp0_r * isp0_r;
          double ISp1 = 13.0 / 12.0 * isp1_l * isp1_l + 0.25 * isp1_r * isp1_r;
          double ISp2 = 13.0 / 12.0 * isp2_l * isp2_l + 0.25 * isp2_r * isp2_r;

          double isd0_l = Gd[i][j + 1][k] - 2 * Gd[i][j][k] + Gd[i][j - 1][k];
          double isd0_r = 3 * Gd[i][j + 1][k] - 4 * Gd[i][j][k] + Gd[i][j - 1][k];
          double isd1_l = Gd[i][j + 2][k] - 2 * Gd[i][j + 1][k] + Gd[i][j][k];
          double isd1_r = Gd[i][j + 2][k] - Gd[i][j][k];
          double isd2_l = Gd[i][j + 3][k] - 2 * Gd[i][j + 2][k] + Gd[i][j + 1][k];
          double isd2_r = Gd[i][j + 3][k] - 4 * Gd[i][j + 2][k] + 3 * Gd[i][j + 1][k];

          double ISd0 = 13.0 / 12.0 * isd0_l * isd0_l + 0.25 * isd0_r * isd0_r;
          double ISd1 = 13.0 / 12.0 * isd1_l * isd1_l + 0.25 * isd1_r * isd1_r;
          double ISd2 = 13.0 / 12.0 * isd2_l * isd2_l + 0.25 * isd2_r * isd2_r;

          double local_alphap[3];
          double local_alphad[3];
          local_alphap[0] = 0.1 / pow(1e-16 + ISp0, P);
          local_alphap[1] = 0.6 / pow(1e-16 + ISp1, P);
          local_alphap[2] = 0.3 / pow(1e-16 + ISp2, P);

          local_alphad[0] = 0.3 / pow(1e-16 + ISd0, P);
          local_alphad[1] = 0.6 / pow(1e-16 + ISd1, P);
          local_alphad[2] = 0.1 / pow(1e-16 + ISd2, P);
          double local_omegap[3];
          double local_omegad[3];
          for (l = 0; l <= 2; l++) {
            local_omegap[l] = local_alphap[l] / (local_alphap[0] + local_alphap[1] + local_alphap[2]);
            local_omegad[l] = local_alphad[l] / (local_alphad[0] + local_alphad[1] + local_alphad[2]);
          }

          double local_q3p[3];
          double local_q3d[3];
          local_q3p[0] = 1.0 / 3.0 * Gp[i][j - 2][k] - 7.0 / 6.0 * Gp[i][j - 1][k] + 11.0 / 6.0 * Gp[i][j][k];
          local_q3p[1] = -1.0 / 6.0 * Gp[i][j - 1][k] + 5.0 / 6.0 * Gp[i][j][k] + 1.0 / 3.0 * Gp[i][j + 1][k];
          local_q3p[2] = 1.0 / 3.0 * Gp[i][j][k] + 5.0 / 6.0 * Gp[i][j + 1][k] - 1.0 / 6.0 * Gp[i][j + 2][k];

          local_q3d[0] = -1.0 / 6.0 * Gd[i][j - 1][k] + 5.0 / 6.0 * Gd[i][j][k] + 1.0 / 3.0 * Gd[i][j + 1][k];
          local_q3d[1] = 1.0 / 3.0 * Gd[i][j][k] + 5.0 / 6.0 * Gd[i][j + 1][k] - 1.0 / 6.0 * Gd[i][j + 2][k];
          local_q3d[2] = 11.0 / 6.0 * Gd[i][j + 1][k] - 7.0 / 6.0 * Gd[i][j + 2][k] + 1.0 / 3.0 * Gd[i][j + 3][k];


          double _G_p = 0;
          double _G_d = 0;
          for (l = 0; l <= 2; l++) {
            _G_p += local_omegap[l] * local_q3p[l];
            _G_d += local_omegad[l] * local_q3d[l];
          }

          G_[i][j][k] = _G_p + _G_d;
        }


  //分区计算U
  for (i = 3; i <= Nx + 3; i++)
    for (j = 3; j <= int(0.5 / dy) + 3; j++)
      for (k = 0; k <= 3; k++) U[i][j][k] = U[i][j][k] - r * (G_[i][j][k] - G_[i][j - 1][k]);

  for (i = int(1.0 / dx) + 3; i <= int(2.0 / dx) + 3; i++)
    for (j = int(0.5 / dy) + 4; j <= Ny + 3; j++)
      for (k = 0; k <= 3; k++) U[i][j][k] = U[i][j][k] - r * (G_[i][j][k] - G_[i][j - 1][k]);
}

#endif  // WENO_FUSION

#endif