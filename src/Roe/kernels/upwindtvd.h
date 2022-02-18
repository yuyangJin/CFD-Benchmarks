#ifndef __UPWINDTVD_H
#define __UPWINDTVD_H
#include "util.h"

//迎风TVD
#define DELTA 0.5  // Q修正系数
#define OMEGA 0.3  // g修正系数

#ifndef UPWINDTVD_FUSION

// x方向的Roe平均
//作用:将U做Roe平均以后返回U_和声速a_
void RoeAVG_x(double U[Nx + 7][Ny + 7][4], double U_[Nx + 7][Ny + 7][4], double a_[Nx + 7][Ny + 7][1]) {
  int i, j;
  double uL, uR, vL, vR, rouL, rouR, pL, pR, D, HL, HR;  //这里uL和uR分别表示u(i+1/2)点左右的速度
  for (i = 0; i <= Nx + 5; i++)
    for (j = 0; j <= Ny + 5; j++)
      if (U[i][j][0] != 0 && U[i + 1][j][0] != 0) {
        rouL = U[i][j][0];
        rouR = U[i + 1][j][0];
        D = sqrt(rouR / rouL);
        uL = U[i][j][1] / U[i][j][0];
        uR = U[i + 1][j][1] / U[i + 1][j][0];
        vL = U[i][j][2] / U[i][j][0];
        vR = U[i + 1][j][2] / U[i + 1][j][0];
        pL = (GAMA - 1) * (U[i][j][3] - 0.5 * rouL * (uL * uL + vL * vL));
        pR = (GAMA - 1) * (U[i + 1][j][3] - 0.5 * rouR * (uR * uR + vR * vR));
        HL = pL * GAMA / (GAMA - 1) / rouL + (uL * uL + vL * vL) / 2;
        HR = pR * GAMA / (GAMA - 1) / rouR + (uR * uR + vR * vR) / 2;
        U_[i][j][0] = rouL * (1 + D) * (1 + D) / 4;  // U_表示roe平均以后的向量其4个分量分别是rou,u,v,H
        U_[i][j][1] = (uL + D * uR) / (1 + D);  //注意它和U分量不同,其中U_[i]表示U[i+0.5]处roe平均的结果
        U_[i][j][2] = (vL + D * vR) / (1 + D);
        U_[i][j][3] = (D * HR + HL) / (1 + D);
        a_[i][j][0] = sqrt((GAMA - 1) * (U_[i][j][3] - 0.5 * (U_[i][j][1] * U_[i][j][1] + U_[i][j][2] * U_[i][j][2])));
      }
}
void RoeAVG_y(double U[Nx + 7][Ny + 7][4], double U_[Nx + 7][Ny + 7][4],
              double a_[Nx + 7][Ny + 7][1])  // 这里U_[i][j][k]为在i+0.5处的Roe平均值
{
  int i, j;
  double uL, uR, vL, vR, rouL, rouR, pL, pR, D, HL, HR;
  for (i = 0; i <= Nx + 5; i++)
    for (j = 0; j <= Ny + 5; j++)
      if (U[i][j][0] != 0 && U[i][j + 1][0] != 0) {
        rouL = U[i][j][0];
        rouR = U[i][j + 1][0];
        D = sqrt(rouR / rouL);
        uL = U[i][j][1] / U[i][j][0];
        uR = U[i][j + 1][1] / U[i][j + 1][0];
        vL = U[i][j][2] / U[i][j][0];
        vR = U[i][j + 1][2] / U[i][j + 1][0];
        pL = (GAMA - 1) * (U[i][j][3] - 0.5 * rouL * (uL * uL + vL * vL));
        pR = (GAMA - 1) * (U[i][j + 1][3] - 0.5 * rouR * (uR * uR + vR * vR));
        HL = pL * GAMA / (GAMA - 1) / rouL + (uL * uL + vL * vL) / 2;
        HR = pR * GAMA / (GAMA - 1) / rouR + (uR * uR + vR * vR) / 2;
        U_[i][j][0] = rouL * (1 + D) * (1 + D) / 4;  // U_表示roe平均以后的向量其3个分量分别是rou,u,v,H注意它和U分量不同
        U_[i][j][1] = (uL + D * uR) / (1 + D);  //其中U[i]表示U[i+0.5]处roe平均的结果
        U_[i][j][2] = (vL + D * vR) / (1 + D);
        U_[i][j][3] = (D * HR + HL) / (1 + D);
        a_[i][j][0] = sqrt((GAMA - 1) * (U_[i][j][3] - 0.5 * (U_[i][j][1] * U_[i][j][1] + U_[i][j][2] * U_[i][j][2])));
      }
}
//作用:用Roe平均以后的U值求特征值与特征向量
void Roe_Eig_x(double U[Nx + 7][Ny + 7][4], double U_[Nx + 7][Ny + 7][4], double LAMDA_[Nx + 7][Ny + 7][4][4],
               double a_[Nx + 7][Ny + 7][1], double R_[Nx + 7][Ny + 7][4][4],
               double L_[Nx + 7][Ny + 7][4][4])  //利用Roe平均以后的值计算x方向的特征值与特征向量
{
  int i, j, k, m, n;
  for (i = 0; i <= Nx + 5; i++)  //先将LAMDA_全赋值0
    for (j = 0; j <= Ny + 5; j++)
      if (U[i][j][0] != 0)
        for (m = 0; m <= 3; m++)
          for (n = 0; n <= 3; n++) LAMDA_[i][j][m][n] = 0;

  for (i = 0; i <= Nx + 5; i++)  //对LAMDA_赋值
    for (j = 0; j <= Ny + 5; j++)
      if (U[i][j][0] != 0) {
        LAMDA_[i][j][0][0] = U_[i][j][1] - a_[i][j][0];
        LAMDA_[i][j][1][1] = U_[i][j][1];
        LAMDA_[i][j][2][2] = U_[i][j][1] + a_[i][j][0];
        LAMDA_[i][j][3][3] = U_[i][j][1];
        for (k = 0; k <= 3; k++)  //这里对lambda做一个修正,防止声速0点出现的问题
        {
          if (LAMDA_[i][j][k][k] >= 0) {
            LAMDA_[i][j][k][k] = 0.5 * (LAMDA_[i][j][k][k] + sqrt(LAMDA_[i][j][k][k] * LAMDA_[i][j][k][k] + 1e-8));
          } else
            LAMDA_[i][j][k][k] = 0.5 * (LAMDA_[i][j][k][k] - sqrt(LAMDA_[i][j][k][k] * LAMDA_[i][j][k][k] + 1e-8));
        }
      }
  //求右特征向量
  double u, v, a, rou, H;
  for (i = 0; i <= Nx + 5; i++)
    for (j = 0; j <= Ny + 5; j++)
      if (U[i][j][0] != 0) {
        rou = U_[i][j][0];
        u = U_[i][j][1];
        v = U_[i][j][2];
        a = a_[i][j][0];
        H = U_[i][j][3];
        R_[i][j][0][0] = 1;
        R_[i][j][0][1] = 1;
        R_[i][j][0][2] = 1;
        R_[i][j][0][3] = 0;
        R_[i][j][1][0] = -a + u;
        R_[i][j][1][1] = u;
        R_[i][j][1][2] = u + a;
        R_[i][j][1][3] = 0;
        R_[i][j][2][0] = v;
        R_[i][j][2][1] = v;
        R_[i][j][2][2] = v;
        R_[i][j][2][3] = 1;
        R_[i][j][3][0] = H - a * u;
        R_[i][j][3][1] = 0.5 * (u * u + v * v);
        R_[i][j][3][2] = H + u * a;
        R_[i][j][3][3] = v;
      }
  //求左特征向量
  double b1, b2;
  for (i = 0; i <= Nx + 5; i++)
    for (j = 0; j <= Ny + 5; j++)
      if (U[i][j][0] != 0) {
        rou = U_[i][j][0];
        u = U_[i][j][1];
        v = U_[i][j][2];
        a = a_[i][j][0];
        b2 = (GAMA - 1) / a / a;
        b1 = b2 * (u * u + v * v) / 2;
        L_[i][j][0][0] = (b1 + u / a) / 2;
        L_[i][j][0][1] = -(b2 * u + 1 / a) / 2;
        L_[i][j][0][2] = -b2 * v / 2;
        L_[i][j][0][3] = b2 / 2;
        L_[i][j][1][0] = 1 - b1;
        L_[i][j][1][1] = b2 * u;
        L_[i][j][1][2] = b2 * v;
        L_[i][j][1][3] = -b2;
        L_[i][j][2][0] = 0.5 * (b1 - u / a);
        L_[i][j][2][1] = 0.5 * (1 / a - b2 * u);
        L_[i][j][2][2] = -0.5 * b2 * v;
        L_[i][j][2][3] = b2 / 2;
        L_[i][j][3][0] = -v;
        L_[i][j][3][1] = 0;
        L_[i][j][3][2] = 1;
        L_[i][j][3][3] = 0;
      }
}
//用Roe平均以后的U值求特征值与特征向量
void Roe_Eig_y(double U[Nx + 7][Ny + 7][4], double U_[Nx + 7][Ny + 7][4], double LAMDA_[Nx + 7][Ny + 7][4][4],
               double a_[Nx + 7][Ny + 7][1], double R_[Nx + 7][Ny + 7][4][4], double L_[Nx + 7][Ny + 7][4][4]) {
  int i, j, k, m, n;
  for (i = 0; i <= Nx + 5; i++)  //先将LAMDA_全赋值0
    for (j = 0; j <= Ny + 5; j++)
      if (U[i][j][0] != 0)
        for (m = 0; m <= 3; m++)
          for (n = 0; n <= 3; n++) LAMDA_[i][j][m][n] = 0;

  for (i = 0; i <= Nx + 5; i++)  //对LAMDA_赋值
    for (j = 0; j <= Ny + 5; j++)
      if (U[i][j][0] != 0) {
        LAMDA_[i][j][0][0] = U_[i][j][2] - a_[i][j][0];
        LAMDA_[i][j][1][1] = U_[i][j][2];
        LAMDA_[i][j][2][2] = U_[i][j][2] + a_[i][j][0];
        LAMDA_[i][j][3][3] = U_[i][j][2];
        for (k = 0; k <= 3; k++) {
          if (LAMDA_[i][j][k][k] >= 0) {
            LAMDA_[i][j][k][k] = 0.5 * (LAMDA_[i][j][k][k] + sqrt(LAMDA_[i][j][k][k] * LAMDA_[i][j][k][k] + 1e-8));
          } else
            LAMDA_[i][j][k][k] = 0.5 * (LAMDA_[i][j][k][k] - sqrt(LAMDA_[i][j][k][k] * LAMDA_[i][j][k][k] + 1e-8));
        }
      }
  //求右特征向量
  double u, v, a, rou, H;
  for (i = 0; i <= Nx + 5; i++)
    for (j = 0; j <= Ny + 5; j++)
      if (U[i][j][0] != 0) {
        rou = U_[i][j][0];
        u = U_[i][j][1];
        v = U_[i][j][2];
        a = a_[i][j][0];
        H = U_[i][j][3];
        R_[i][j][0][0] = 1;
        R_[i][j][0][1] = 1;
        R_[i][j][0][2] = 1;
        R_[i][j][0][3] = 0;
        R_[i][j][1][0] = u;
        R_[i][j][1][1] = u;
        R_[i][j][1][2] = u;
        R_[i][j][1][3] = 1;
        R_[i][j][2][0] = v - a;
        R_[i][j][2][1] = v;
        R_[i][j][2][2] = v + a;
        R_[i][j][2][3] = 0;
        R_[i][j][3][0] = H - a * v;
        R_[i][j][3][1] = 0.5 * (u * u + v * v);
        R_[i][j][3][2] = H + v * a;
        R_[i][j][3][3] = u;
      }
  //求左特征向量
  double b1, b2;
  for (i = 0; i <= Nx + 5; i++)
    for (j = 0; j <= Ny + 5; j++)
      if (U[i][j][0] != 0) {
        rou = U_[i][j][0];
        u = U_[i][j][1];
        v = U_[i][j][2];
        a = a_[i][j][0];
        b2 = (GAMA - 1) / a / a;
        b1 = b2 * (u * u + v * v) / 2;
        L_[i][j][0][0] = (b1 + v / a) / 2;
        L_[i][j][0][1] = b2 * u / 2;
        L_[i][j][0][2] = -(b2 * v + 1 / a) / 2;
        L_[i][j][0][3] = b2 / 2;
        L_[i][j][1][0] = 1 - b1;
        L_[i][j][1][1] = b2 * u;
        L_[i][j][1][2] = b2 * v;
        L_[i][j][1][3] = -b2;
        L_[i][j][2][0] = 0.5 * (b1 - v / a);
        L_[i][j][2][1] = -0.5 * (b2 * u);
        L_[i][j][2][2] = 0.5 * (1 / a - b2 * v);
        L_[i][j][2][3] = b2 / 2;
        L_[i][j][3][0] = -u;
        L_[i][j][3][1] = 1;
        L_[i][j][3][2] = 0;
        L_[i][j][3][3] = 0;
      }
}

void UpWindTVD_x(double U[Nx + 7][Ny + 7][4], double U_[Nx + 7][Ny + 7][4], double a_[Nx + 7][Ny + 7][1],
                 double LAMDA_[Nx + 7][Ny + 7][4][4], double L_[Nx + 7][Ny + 7][4][4], double R_[Nx + 7][Ny + 7][4][4],
                 double alpha_[Nx + 7][Ny + 7][4], double g_[Nx + 7][Ny + 7][4], double g[Nx + 7][Ny + 7][4],
                 double gama_[Nx + 7][Ny + 7][4], double Q_[Nx + 7][Ny + 7][4], double theta[Nx + 7][Ny + 7][4],
                 double F[Nx + 7][Ny + 7][4], double F_[Nx + 7][Ny + 7][4], double dx, double dy,
                 double &dt)  //此函数将之前计算的特征之余特征向量通过TVD算法计算U
{
  int i, j, k, l;
  double r, z;
  RoeAVG_x(U, U_, a_);                      //求x方向的Roe平均
  dt = CFL_(U_, a_, dx, dy, UPWINDTVDCFL);  //用Roe平均以后的流场计算dt
  Roe_Eig_x(U, U_, LAMDA_, a_, R_, L_);
  r = dt / dx;
  //计算alpha_
  for (i = 0; i <= Nx + 5; i++)
    for (j = 0; j <= Ny + 5; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) alpha_[i][j][k] = 0;

  for (i = 0; i <= Nx + 5; i++)
    for (j = 0; j <= Ny + 5; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++)
          for (l = 0; l <= 3; l++) alpha_[i][j][k] += L_[i][j][k][l] * (U[i + 1][j][l] - U[i][j][l]);
  //计算g_
  for (i = 0; i <= Nx + 5; i++)
    for (j = 0; j <= Ny + 5; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) {
          z = LAMDA_[i][j][k][k];
          if (fabs(z) < DELTA) {
            g_[i][j][k] = alpha_[i][j][k] * 0.5 * (0.5 * (z * z / DELTA + DELTA) + r * z * z);
          } else
            g_[i][j][k] = alpha_[i][j][k] * 0.5 * (fabs(z) + r * z * z);
        }
  //用g_通过minmod来计算g
  for (i = 2; i <= Nx + 4; i++)
    for (j = 2; j <= Ny + 4; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) g[i][j][k] = minmod(g_[i][j][k], g_[i - 1][j][k]);
  //对g进行修正
  for (i = 2; i <= Nx + 4; i++)
    for (j = 2; j <= Ny + 4; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) {
          if (alpha_[i][j][k] != 0 || alpha_[i - 1][j][k] != 0) {
            g[i][j][k] = (1 +
                          OMEGA * fabs(alpha_[i][j][k] - alpha_[i - 1][j][k]) /
                              (fabs(alpha_[i][j][k]) + fabs(alpha_[i - 1][j][k]))) *
                         g[i][j][k];
          } else
            g[i][j][k] = 0;
        }
  //计算gama_
  for (i = 2; i <= Nx + 3; i++)
    for (j = 2; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) {
          if (alpha_[i][j][k] != 0) {
            gama_[i][j][k] = (g[i + 1][j][k] - g[i][j][k]) / alpha_[i][j][k];
          } else
            gama_[i][j][k] = 0;
        }
  //用gama_和lambda_来计算Q_
  for (i = 2; i <= Nx + 3; i++)
    for (j = 2; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) {
          z = LAMDA_[i][j][k][k] + gama_[i][j][k];
          if (fabs(z) < DELTA) {
            Q_[i][j][k] = 0.5 * (z * z / DELTA + DELTA);
          } else
            Q_[i][j][k] = fabs(z);
        }
  //计算反扩散粘性项theta=(g[i]+g[i+1]-Q(a+gama_)*alpha_)*R
  for (i = 2; i <= Nx + 3; i++)
    for (j = 2; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) theta[i][j][k] = 0;

  for (i = 2; i <= Nx + 3; i++)
    for (j = 2; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) {
          for (l = 0; l <= 3; l++) {
            theta[i][j][k] += R_[i][j][k][l] * (g[i][j][l] + g[i + 1][j][l] - Q_[i][j][l] * alpha_[i][j][l]);
          }
        }
  //计算F
  for (i = 0; i <= Nx + 6; i++)
    for (j = 0; j <= Ny + 6; j++)
      if (U[i][j][0] != 0) U2F(U[i][j], F[i][j]);
  //计算F_
  for (i = 2; i <= Nx + 3; i++)
    for (j = 2; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) F_[i][j][k] = 0.5 * (F[i][j][k] + F[i + 1][j][k]) + 0.5 * theta[i][j][k];
  //分区计算U
  for (i = 3; i <= Nx + 3; i++)
    for (j = 3; j <= int(0.5 / dy) + 3; j++)
      for (k = 0; k <= 3; k++) U[i][j][k] = U[i][j][k] - r * (F_[i][j][k] - F_[i - 1][j][k]);

  for (i = int(1.0 / dx) + 3; i <= int(2.0 / dx) + 3; i++)
    for (j = int(0.5 / dy) + 4; j <= Ny + 3; j++)
      for (k = 0; k <= 3; k++) U[i][j][k] = U[i][j][k] - r * (F_[i][j][k] - F_[i - 1][j][k]);
}
void UpWindTVD_y(double U[Nx + 7][Ny + 7][4], double U_[Nx + 7][Ny + 7][4], double a_[Nx + 7][Ny + 7][1],
                 double LAMDA_[Nx + 7][Ny + 7][4][4], double L_[Nx + 7][Ny + 7][4][4], double R_[Nx + 7][Ny + 7][4][4],
                 double alpha_[Nx + 7][Ny + 7][4], double g_[Nx + 7][Ny + 7][4], double g[Nx + 7][Ny + 7][4],
                 double gama_[Nx + 7][Ny + 7][4], double Q_[Nx + 7][Ny + 7][4], double theta[Nx + 7][Ny + 7][4],
                 double G[Nx + 7][Ny + 7][4], double G_[Nx + 7][Ny + 7][4], double dx, double dy, double &dt) {
  int i, j, k, l;
  double r, z;
  RoeAVG_y(U, U_, a_);
  dt = CFL_(U_, a_, dx, dy, UPWINDTVDCFL);
  Roe_Eig_y(U, U_, LAMDA_, a_, R_, L_);
  r = dt / dy;
  //计算alpha_
  for (i = 0; i <= Nx + 5; i++)
    for (j = 0; j <= Ny + 5; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) alpha_[i][j][k] = 0;

  for (i = 0; i <= Nx + 5; i++)
    for (j = 0; j <= Ny + 5; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++)
          for (l = 0; l <= 3; l++) alpha_[i][j][k] += L_[i][j][k][l] * (U[i][j + 1][l] - U[i][j][l]);
  //计算g_
  for (i = 0; i <= Nx + 5; i++)
    for (j = 0; j <= Ny + 5; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) {
          z = LAMDA_[i][j][k][k];
          if (fabs(z) < DELTA) {
            g_[i][j][k] = alpha_[i][j][k] * 0.5 * (0.5 * (z * z / DELTA + DELTA) + r * z * z);
          } else
            g_[i][j][k] = alpha_[i][j][k] * 0.5 * (fabs(z) + r * z * z);
        }
  //计算g
  for (i = 2; i <= Nx + 4; i++)
    for (j = 2; j <= Ny + 4; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) g[i][j][k] = minmod(g_[i][j][k], g_[i][j - 1][k]);
  //对g做修正
  for (i = 2; i <= Nx + 4; i++)
    for (j = 2; j <= Ny + 4; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) {
          if (alpha_[i][j][k] != 0 || alpha_[i][j - 1][k] != 0) {
            g[i][j][k] = (1 +
                          OMEGA * fabs(alpha_[i][j][k] - alpha_[i][j - 1][k]) /
                              (fabs(alpha_[i][j][k]) + fabs(alpha_[i][j - 1][k]))) *
                         g[i][j][k];
          } else
            g[i][j][k] = 0;
        }
  //计算gama_
  for (i = 2; i <= Nx + 3; i++)
    for (j = 2; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) {
          if (alpha_[i][j][k] != 0) {
            gama_[i][j][k] = (g[i][j + 1][k] - g[i][j][k]) / alpha_[i][j][k];
          } else
            gama_[i][j][k] = 0;
        }
  //计算Q_
  for (i = 2; i <= Nx + 3; i++)
    for (j = 2; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) {
          z = LAMDA_[i][j][k][k] + gama_[i][j][k];
          if (fabs(z) < DELTA) {
            Q_[i][j][k] = 0.5 * (z * z / DELTA + DELTA);
          } else
            Q_[i][j][k] = fabs(z);
        }
  //计算粘性项theta
  for (i = 2; i <= Nx + 3; i++)
    for (j = 2; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) theta[i][j][k] = 0;

  for (i = 2; i <= Nx + 3; i++)
    for (j = 2; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++)
          for (l = 0; l <= 3; l++)
            theta[i][j][k] += R_[i][j][k][l] * (g[i][j][l] + g[i][j + 1][l] - Q_[i][j][l] * alpha_[i][j][l]);
  //计算G
  for (i = 0; i <= Nx + 6; i++)
    for (j = 0; j <= Ny + 6; j++)
      if (U[i][j][0] != 0) U2G(U[i][j], G[i][j]);
  //计算G_
  for (i = 2; i <= Nx + 3; i++)
    for (j = 2; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) G_[i][j][k] = 0.5 * (G[i][j][k] + G[i][j + 1][k]) + 0.5 * theta[i][j][k];
  //分区计算U
  for (i = 3; i <= Nx + 3; i++)
    for (j = 3; j <= int(0.5 / dy) + 3; j++)
      for (k = 0; k <= 3; k++) U[i][j][k] = U[i][j][k] - r * (G_[i][j][k] - G_[i][j - 1][k]);

  for (i = int(1.0 / dx) + 3; i <= int(2.0 / dx) + 3; i++)
    for (j = int(0.5 / dy) + 4; j <= Ny + 3; j++)
      for (k = 0; k <= 3; k++) U[i][j][k] = U[i][j][k] - r * (G_[i][j][k] - G_[i][j - 1][k]);
}

#else

// x方向的Roe平均
//作用:将U做Roe平均以后返回U_和声速a_
void RoeAVG_x(double U[Nx + 7][Ny + 7][4], double U_[Nx + 7][Ny + 7][4], double a_[Nx + 7][Ny + 7][1]) {
  int i, j;
  double uL, uR, vL, vR, rouL, rouR, pL, pR, D, HL, HR;  //这里uL和uR分别表示u(i+1/2)点左右的速度
  for (i = 0; i <= Nx + 5; i++)
    for (j = 0; j <= Ny + 5; j++)
      if (U[i][j][0] != 0 && U[i + 1][j][0] != 0) {
        rouL = U[i][j][0];
        rouR = U[i + 1][j][0];
        D = sqrt(rouR / rouL);
        uL = U[i][j][1] / U[i][j][0];
        uR = U[i + 1][j][1] / U[i + 1][j][0];
        vL = U[i][j][2] / U[i][j][0];
        vR = U[i + 1][j][2] / U[i + 1][j][0];
        pL = (GAMA - 1) * (U[i][j][3] - 0.5 * rouL * (uL * uL + vL * vL));
        pR = (GAMA - 1) * (U[i + 1][j][3] - 0.5 * rouR * (uR * uR + vR * vR));
        HL = pL * GAMA / (GAMA - 1) / rouL + (uL * uL + vL * vL) / 2;
        HR = pR * GAMA / (GAMA - 1) / rouR + (uR * uR + vR * vR) / 2;
        U_[i][j][0] = rouL * (1 + D) * (1 + D) / 4;  // U_表示roe平均以后的向量其4个分量分别是rou,u,v,H
        U_[i][j][1] = (uL + D * uR) / (1 + D);  //注意它和U分量不同,其中U_[i]表示U[i+0.5]处roe平均的结果
        U_[i][j][2] = (vL + D * vR) / (1 + D);
        U_[i][j][3] = (D * HR + HL) / (1 + D);
        a_[i][j][0] = sqrt((GAMA - 1) * (U_[i][j][3] - 0.5 * (U_[i][j][1] * U_[i][j][1] + U_[i][j][2] * U_[i][j][2])));
      }
}
void RoeAVG_y(double U[Nx + 7][Ny + 7][4], double U_[Nx + 7][Ny + 7][4],
              double a_[Nx + 7][Ny + 7][1])  // 这里U_[i][j][k]为在i+0.5处的Roe平均值
{
  int i, j;
  double uL, uR, vL, vR, rouL, rouR, pL, pR, D, HL, HR;
  for (i = 0; i <= Nx + 5; i++)
    for (j = 0; j <= Ny + 5; j++)
      if (U[i][j][0] != 0 && U[i][j + 1][0] != 0) {
        rouL = U[i][j][0];
        rouR = U[i][j + 1][0];
        D = sqrt(rouR / rouL);
        uL = U[i][j][1] / U[i][j][0];
        uR = U[i][j + 1][1] / U[i][j + 1][0];
        vL = U[i][j][2] / U[i][j][0];
        vR = U[i][j + 1][2] / U[i][j + 1][0];
        pL = (GAMA - 1) * (U[i][j][3] - 0.5 * rouL * (uL * uL + vL * vL));
        pR = (GAMA - 1) * (U[i][j + 1][3] - 0.5 * rouR * (uR * uR + vR * vR));
        HL = pL * GAMA / (GAMA - 1) / rouL + (uL * uL + vL * vL) / 2;
        HR = pR * GAMA / (GAMA - 1) / rouR + (uR * uR + vR * vR) / 2;
        U_[i][j][0] = rouL * (1 + D) * (1 + D) / 4;  // U_表示roe平均以后的向量其3个分量分别是rou,u,v,H注意它和U分量不同
        U_[i][j][1] = (uL + D * uR) / (1 + D);  //其中U[i]表示U[i+0.5]处roe平均的结果
        U_[i][j][2] = (vL + D * vR) / (1 + D);
        U_[i][j][3] = (D * HR + HL) / (1 + D);
        a_[i][j][0] = sqrt((GAMA - 1) * (U_[i][j][3] - 0.5 * (U_[i][j][1] * U_[i][j][1] + U_[i][j][2] * U_[i][j][2])));
      }
}
//作用:用Roe平均以后的U值求特征值与特征向量
void Roe_Eig_x(double U[Nx + 7][Ny + 7][4], double U_[Nx + 7][Ny + 7][4], double LAMDA_[Nx + 7][Ny + 7][4][4],
               double a_[Nx + 7][Ny + 7][1], double R_[Nx + 7][Ny + 7][4][4],
               double L_[Nx + 7][Ny + 7][4][4])  //利用Roe平均以后的值计算x方向的特征值与特征向量
{
  int i, j, k, m, n;
  for (i = 0; i <= Nx + 5; i++)  //先将LAMDA_全赋值0
    for (j = 0; j <= Ny + 5; j++)
      if (U[i][j][0] != 0)
        for (m = 0; m <= 3; m++)
          for (n = 0; n <= 3; n++) LAMDA_[i][j][m][n] = 0;

  for (i = 0; i <= Nx + 5; i++)  //对LAMDA_赋值
    for (j = 0; j <= Ny + 5; j++)
      if (U[i][j][0] != 0) {
        LAMDA_[i][j][0][0] = U_[i][j][1] - a_[i][j][0];
        LAMDA_[i][j][1][1] = U_[i][j][1];
        LAMDA_[i][j][2][2] = U_[i][j][1] + a_[i][j][0];
        LAMDA_[i][j][3][3] = U_[i][j][1];
        for (k = 0; k <= 3; k++)  //这里对lambda做一个修正,防止声速0点出现的问题
        {
          if (LAMDA_[i][j][k][k] >= 0) {
            LAMDA_[i][j][k][k] = 0.5 * (LAMDA_[i][j][k][k] + sqrt(LAMDA_[i][j][k][k] * LAMDA_[i][j][k][k] + 1e-8));
          } else
            LAMDA_[i][j][k][k] = 0.5 * (LAMDA_[i][j][k][k] - sqrt(LAMDA_[i][j][k][k] * LAMDA_[i][j][k][k] + 1e-8));
        }
      }
  //求右特征向量
  for (i = 0; i <= Nx + 5; i++)
    for (j = 0; j <= Ny + 5; j++)
      if (U[i][j][0] != 0) {
        double u, v, a, rou, H;
        rou = U_[i][j][0];
        u = U_[i][j][1];
        v = U_[i][j][2];
        a = a_[i][j][0];
        H = U_[i][j][3];
        R_[i][j][0][0] = 1;
        R_[i][j][0][1] = 1;
        R_[i][j][0][2] = 1;
        R_[i][j][0][3] = 0;
        R_[i][j][1][0] = -a + u;
        R_[i][j][1][1] = u;
        R_[i][j][1][2] = u + a;
        R_[i][j][1][3] = 0;
        R_[i][j][2][0] = v;
        R_[i][j][2][1] = v;
        R_[i][j][2][2] = v;
        R_[i][j][2][3] = 1;
        R_[i][j][3][0] = H - a * u;
        R_[i][j][3][1] = 0.5 * (u * u + v * v);
        R_[i][j][3][2] = H + u * a;
        R_[i][j][3][3] = v;
      }
  //求左特征向量
  for (i = 0; i <= Nx + 5; i++)
    for (j = 0; j <= Ny + 5; j++)
      if (U[i][j][0] != 0) {
        double u, v, a, rou, H;
        double b1, b2;
        rou = U_[i][j][0];
        u = U_[i][j][1];
        v = U_[i][j][2];
        a = a_[i][j][0];
        b2 = (GAMA - 1) / a / a;
        b1 = b2 * (u * u + v * v) / 2;
        L_[i][j][0][0] = (b1 + u / a) / 2;
        L_[i][j][0][1] = -(b2 * u + 1 / a) / 2;
        L_[i][j][0][2] = -b2 * v / 2;
        L_[i][j][0][3] = b2 / 2;
        L_[i][j][1][0] = 1 - b1;
        L_[i][j][1][1] = b2 * u;
        L_[i][j][1][2] = b2 * v;
        L_[i][j][1][3] = -b2;
        L_[i][j][2][0] = 0.5 * (b1 - u / a);
        L_[i][j][2][1] = 0.5 * (1 / a - b2 * u);
        L_[i][j][2][2] = -0.5 * b2 * v;
        L_[i][j][2][3] = b2 / 2;
        L_[i][j][3][0] = -v;
        L_[i][j][3][1] = 0;
        L_[i][j][3][2] = 1;
        L_[i][j][3][3] = 0;
      }
}
//用Roe平均以后的U值求特征值与特征向量
void Roe_Eig_y(double U[Nx + 7][Ny + 7][4], double U_[Nx + 7][Ny + 7][4], double LAMDA_[Nx + 7][Ny + 7][4][4],
               double a_[Nx + 7][Ny + 7][1], double R_[Nx + 7][Ny + 7][4][4], double L_[Nx + 7][Ny + 7][4][4]) {
  int i, j, k, m, n;
  for (i = 0; i <= Nx + 5; i++)  //先将LAMDA_全赋值0
    for (j = 0; j <= Ny + 5; j++)
      if (U[i][j][0] != 0)
        for (m = 0; m <= 3; m++)
          for (n = 0; n <= 3; n++) LAMDA_[i][j][m][n] = 0;

  for (i = 0; i <= Nx + 5; i++)  //对LAMDA_赋值
    for (j = 0; j <= Ny + 5; j++)
      if (U[i][j][0] != 0) {
        LAMDA_[i][j][0][0] = U_[i][j][2] - a_[i][j][0];
        LAMDA_[i][j][1][1] = U_[i][j][2];
        LAMDA_[i][j][2][2] = U_[i][j][2] + a_[i][j][0];
        LAMDA_[i][j][3][3] = U_[i][j][2];
        for (k = 0; k <= 3; k++) {
          if (LAMDA_[i][j][k][k] >= 0) {
            LAMDA_[i][j][k][k] = 0.5 * (LAMDA_[i][j][k][k] + sqrt(LAMDA_[i][j][k][k] * LAMDA_[i][j][k][k] + 1e-8));
          } else
            LAMDA_[i][j][k][k] = 0.5 * (LAMDA_[i][j][k][k] - sqrt(LAMDA_[i][j][k][k] * LAMDA_[i][j][k][k] + 1e-8));
        }
      }
  //求右特征向量
  for (i = 0; i <= Nx + 5; i++)
    for (j = 0; j <= Ny + 5; j++)
      if (U[i][j][0] != 0) {
        double u, v, a, rou, H;
        rou = U_[i][j][0];
        u = U_[i][j][1];
        v = U_[i][j][2];
        a = a_[i][j][0];
        H = U_[i][j][3];
        R_[i][j][0][0] = 1;
        R_[i][j][0][1] = 1;
        R_[i][j][0][2] = 1;
        R_[i][j][0][3] = 0;
        R_[i][j][1][0] = u;
        R_[i][j][1][1] = u;
        R_[i][j][1][2] = u;
        R_[i][j][1][3] = 1;
        R_[i][j][2][0] = v - a;
        R_[i][j][2][1] = v;
        R_[i][j][2][2] = v + a;
        R_[i][j][2][3] = 0;
        R_[i][j][3][0] = H - a * v;
        R_[i][j][3][1] = 0.5 * (u * u + v * v);
        R_[i][j][3][2] = H + v * a;
        R_[i][j][3][3] = u;
      }
  //求左特征向量
  for (i = 0; i <= Nx + 5; i++)
    for (j = 0; j <= Ny + 5; j++)
      if (U[i][j][0] != 0) {
        double u, v, a, rou, H;
        double b1, b2;
        rou = U_[i][j][0];
        u = U_[i][j][1];
        v = U_[i][j][2];
        a = a_[i][j][0];
        b2 = (GAMA - 1) / a / a;
        b1 = b2 * (u * u + v * v) / 2;
        L_[i][j][0][0] = (b1 + v / a) / 2;
        L_[i][j][0][1] = b2 * u / 2;
        L_[i][j][0][2] = -(b2 * v + 1 / a) / 2;
        L_[i][j][0][3] = b2 / 2;
        L_[i][j][1][0] = 1 - b1;
        L_[i][j][1][1] = b2 * u;
        L_[i][j][1][2] = b2 * v;
        L_[i][j][1][3] = -b2;
        L_[i][j][2][0] = 0.5 * (b1 - v / a);
        L_[i][j][2][1] = -0.5 * (b2 * u);
        L_[i][j][2][2] = 0.5 * (1 / a - b2 * v);
        L_[i][j][2][3] = b2 / 2;
        L_[i][j][3][0] = -u;
        L_[i][j][3][1] = 1;
        L_[i][j][3][2] = 0;
        L_[i][j][3][3] = 0;
      }
}

void UpWindTVD_x(double U[Nx + 7][Ny + 7][4], double U_[Nx + 7][Ny + 7][4], double a_[Nx + 7][Ny + 7][1],
                 double LAMDA_[Nx + 7][Ny + 7][4][4], double L_[Nx + 7][Ny + 7][4][4], double R_[Nx + 7][Ny + 7][4][4],
                 double alpha_[Nx + 7][Ny + 7][4], double g_[Nx + 7][Ny + 7][4], double g[Nx + 7][Ny + 7][4],
                 double gama_[Nx + 7][Ny + 7][4], double Q_[Nx + 7][Ny + 7][4], double theta[Nx + 7][Ny + 7][4],
                 double F[Nx + 7][Ny + 7][4], double F_[Nx + 7][Ny + 7][4], double dx, double dy,
                 double &dt)  //此函数将之前计算的特征之余特征向量通过TVD算法计算U
{
  int i, j, k, l;
  double r, z;
  RoeAVG_x(U, U_, a_);                      //求x方向的Roe平均
  dt = CFL_(U_, a_, dx, dy, UPWINDTVDCFL);  //用Roe平均以后的流场计算dt
  r = dt / dx;

  static double local_LAMDA_[Nx + 7][Ny + 7][4];

  for (i = 2; i <= Nx + 3; i++)
    for (j = 2; j <= Ny + 3; j++)
      if (U[i][j][0] != 0) {
        local_LAMDA_[i][j][0] = U_[i][j][1] - a_[i][j][0];
        local_LAMDA_[i][j][1] = U_[i][j][1];
        local_LAMDA_[i][j][2] = U_[i][j][1] + a_[i][j][0];
        local_LAMDA_[i][j][3] = U_[i][j][1];

        double u, v, a, rou, H;
        double b1, b2;
        double local_L_[4][4];
        rou = U_[i][j][0];
        u = U_[i][j][1];
        v = U_[i][j][2];
        a = a_[i][j][0];
        b2 = (GAMA - 1) / a / a;
        b1 = b2 * (u * u + v * v) / 2;
        local_L_[0][0] = (b1 + u / a) / 2;
        local_L_[0][1] = -(b2 * u + 1 / a) / 2;
        local_L_[0][2] = -b2 * v / 2;
        local_L_[0][3] = b2 / 2;
        local_L_[1][0] = 1 - b1;
        local_L_[1][1] = b2 * u;
        local_L_[1][2] = b2 * v;
        local_L_[1][3] = -b2;
        local_L_[2][0] = 0.5 * (b1 - u / a);
        local_L_[2][1] = 0.5 * (1 / a - b2 * u);
        local_L_[2][2] = -0.5 * b2 * v;
        local_L_[2][3] = b2 / 2;
        local_L_[3][0] = -v;
        local_L_[3][1] = 0;
        local_L_[3][2] = 1;
        local_L_[3][3] = 0;

        for (k = 0; k <= 3; k++) {
          alpha_[i][j][k] = 0;
          for (l = 0; l <= 3; l++) alpha_[i][j][k] += local_L_[k][l] * (U[i + 1][j][l] - U[i][j][l]);

          if (local_LAMDA_[i][j][k] >= 0) {
            local_LAMDA_[i][j][k] =
                0.5 * (local_LAMDA_[i][j][k] + sqrt(local_LAMDA_[i][j][k] * local_LAMDA_[i][j][k] + 1e-8));
          } else
            local_LAMDA_[i][j][k] =
                0.5 * (local_LAMDA_[i][j][k] - sqrt(local_LAMDA_[i][j][k] * local_LAMDA_[i][j][k] + 1e-8));

          z = local_LAMDA_[i][j][k];

          if (fabs(z) < DELTA) {
            g_[i][j][k] = alpha_[i][j][k] * 0.5 * (0.5 * (z * z / DELTA + DELTA) + r * z * z);
          } else {
            g_[i][j][k] = alpha_[i][j][k] * 0.5 * (fabs(z) + r * z * z);
          }

          double t = minmod(g_[i][j][k], g_[i - 1][j][k]);

          if (alpha_[i][j][k] != 0 || alpha_[i - 1][j][k] != 0) {
            g[i][j][k] = (1 +
                          OMEGA * fabs(alpha_[i][j][k] - alpha_[i - 1][j][k]) /
                              (fabs(alpha_[i][j][k]) + fabs(alpha_[i - 1][j][k]))) *
                         t;
          } else
            g[i][j][k] = 0;
        }
      }
  //计算gama_
  for (i = 2; i <= Nx + 3; i++)
    for (j = 2; j <= Ny + 3; j++)
      if (U[i][j][0] != 0) {
        double u, v, a, H;
        u = U_[i][j][1];
        v = U_[i][j][2];
        a = a_[i][j][0];
        H = U_[i][j][3];
        double local_R_[4][4];
        local_R_[0][0] = 1;
        local_R_[0][1] = 1;
        local_R_[0][2] = 1;
        local_R_[0][3] = 0;
        local_R_[1][0] = -a + u;
        local_R_[1][1] = u;
        local_R_[1][2] = u + a;
        local_R_[1][3] = 0;
        local_R_[2][0] = v;
        local_R_[2][1] = v;
        local_R_[2][2] = v;
        local_R_[2][3] = 1;
        local_R_[3][0] = H - a * u;
        local_R_[3][1] = 0.5 * (u * u + v * v);
        local_R_[3][2] = H + u * a;
        local_R_[3][3] = v;

        double local_Q_[4];
        for (k = 0; k <= 3; k++) {
          double local_gama_ = 0;
          if (alpha_[i][j][k] != 0) {
            local_gama_ = (g[i + 1][j][k] - g[i][j][k]) / alpha_[i][j][k];
          }

          z = local_LAMDA_[i][j][k] + local_gama_;
          if (fabs(z) < DELTA) {
            local_Q_[k] = 0.5 * (z * z / DELTA + DELTA);
          } else
            local_Q_[k] = fabs(z);
        }
        for (k = 0; k <= 3; k++) {
          theta[i][j][k] = 0;
          for (l = 0; l <= 3; l++) {
            theta[i][j][k] += local_R_[k][l] * (g[i][j][l] + g[i + 1][j][l] - local_Q_[l] * alpha_[i][j][l]);
          }
        }
      }

  //计算F
  for (i = 2; i <= Nx + 4; i++)
    for (j = 2; j <= Ny + 4; j++)
      if (U[i][j][0] != 0) U2F(U[i][j], F[i][j]);
  //计算F_
  for (i = 2; i <= Nx + 3; i++)
    for (j = 2; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) F_[i][j][k] = 0.5 * (F[i][j][k] + F[i + 1][j][k]) + 0.5 * theta[i][j][k];
  //分区计算U
  for (i = 3; i <= Nx + 3; i++)
    for (j = 3; j <= int(0.5 / dy) + 3; j++)
      for (k = 0; k <= 3; k++) U[i][j][k] = U[i][j][k] - r * (F_[i][j][k] - F_[i - 1][j][k]);

  for (i = int(1.0 / dx) + 3; i <= int(2.0 / dx) + 3; i++)
    for (j = int(0.5 / dy) + 4; j <= Ny + 3; j++)
      for (k = 0; k <= 3; k++) U[i][j][k] = U[i][j][k] - r * (F_[i][j][k] - F_[i - 1][j][k]);
}
void UpWindTVD_y(double U[Nx + 7][Ny + 7][4], double U_[Nx + 7][Ny + 7][4], double a_[Nx + 7][Ny + 7][1],
                 double LAMDA_[Nx + 7][Ny + 7][4][4], double L_[Nx + 7][Ny + 7][4][4], double R_[Nx + 7][Ny + 7][4][4],
                 double alpha_[Nx + 7][Ny + 7][4], double g_[Nx + 7][Ny + 7][4], double g[Nx + 7][Ny + 7][4],
                 double gama_[Nx + 7][Ny + 7][4], double Q_[Nx + 7][Ny + 7][4], double theta[Nx + 7][Ny + 7][4],
                 double G[Nx + 7][Ny + 7][4], double G_[Nx + 7][Ny + 7][4], double dx, double dy, double &dt) {
  int i, j, k, l;
  double r, z;
  RoeAVG_y(U, U_, a_);
  dt = CFL_(U_, a_, dx, dy, UPWINDTVDCFL);
  r = dt / dy;

  static double local_LAMDA_[Nx + 7][Ny + 7][4];

  for (i = 2; i <= Nx + 3; i++)
    for (j = 2; j <= Ny + 3; j++)
      if (U[i][j][0] != 0) {
        local_LAMDA_[i][j][0] = U_[i][j][2] - a_[i][j][0];
        local_LAMDA_[i][j][1] = U_[i][j][2];
        local_LAMDA_[i][j][2] = U_[i][j][2] + a_[i][j][0];
        local_LAMDA_[i][j][3] = U_[i][j][2];

        double u, v, a, H;
        double b1, b2;
        double local_L_[4][4];
        u = U_[i][j][1];
        v = U_[i][j][2];
        a = a_[i][j][0];
        b2 = (GAMA - 1) / a / a;
        b1 = b2 * (u * u + v * v) / 2;
        local_L_[0][0] = (b1 + v / a) / 2;
        local_L_[0][1] = b2 * u / 2;
        local_L_[0][2] = -(b2 * v + 1 / a) / 2;
        local_L_[0][3] = b2 / 2;
        local_L_[1][0] = 1 - b1;
        local_L_[1][1] = b2 * u;
        local_L_[1][2] = b2 * v;
        local_L_[1][3] = -b2;
        local_L_[2][0] = 0.5 * (b1 - v / a);
        local_L_[2][1] = -0.5 * (b2 * u);
        local_L_[2][2] = 0.5 * (1 / a - b2 * v);
        local_L_[2][3] = b2 / 2;
        local_L_[3][0] = -u;
        local_L_[3][1] = 1;
        local_L_[3][2] = 0;
        local_L_[3][3] = 0;

        for (k = 0; k <= 3; k++) {
          alpha_[i][j][k] = 0;
          for (l = 0; l <= 3; l++) alpha_[i][j][k] += local_L_[k][l] * (U[i][j + 1][l] - U[i][j][l]);
          if (local_LAMDA_[i][j][k] >= 0) {
            local_LAMDA_[i][j][k] =
                0.5 * (local_LAMDA_[i][j][k] + sqrt(local_LAMDA_[i][j][k] * local_LAMDA_[i][j][k] + 1e-8));
          } else
            local_LAMDA_[i][j][k] =
                0.5 * (local_LAMDA_[i][j][k] - sqrt(local_LAMDA_[i][j][k] * local_LAMDA_[i][j][k] + 1e-8));

          z = local_LAMDA_[i][j][k];
          if (fabs(z) < DELTA) {
            g_[i][j][k] = alpha_[i][j][k] * 0.5 * (0.5 * (z * z / DELTA + DELTA) + r * z * z);
          } else
            g_[i][j][k] = alpha_[i][j][k] * 0.5 * (fabs(z) + r * z * z);

          g[i][j][k] = minmod(g_[i][j][k], g_[i][j - 1][k]);

          if (alpha_[i][j][k] != 0 || alpha_[i][j - 1][k] != 0) {
            g[i][j][k] = (1 +
                          OMEGA * fabs(alpha_[i][j][k] - alpha_[i][j - 1][k]) /
                              (fabs(alpha_[i][j][k]) + fabs(alpha_[i][j - 1][k]))) *
                         g[i][j][k];
          } else
            g[i][j][k] = 0;
        }
      }

  //计算gama_
  for (i = 2; i <= Nx + 3; i++)
    for (j = 2; j <= Ny + 3; j++)
      if (U[i][j][0] != 0) {
        double u, v, a, H;
        double local_R_[4][4];
        u = U_[i][j][1];
        v = U_[i][j][2];
        a = a_[i][j][0];
        H = U_[i][j][3];
        local_R_[0][0] = 1;
        local_R_[0][1] = 1;
        local_R_[0][2] = 1;
        local_R_[0][3] = 0;
        local_R_[1][0] = u;
        local_R_[1][1] = u;
        local_R_[1][2] = u;
        local_R_[1][3] = 1;
        local_R_[2][0] = v - a;
        local_R_[2][1] = v;
        local_R_[2][2] = v + a;
        local_R_[2][3] = 0;
        local_R_[3][0] = H - a * v;
        local_R_[3][1] = 0.5 * (u * u + v * v);
        local_R_[3][2] = H + v * a;
        local_R_[3][3] = u;

        double local_Q_[4];
        for (k = 0; k <= 3; k++) {
          double local_gama_ = 0;
          if (alpha_[i][j][k] != 0) {
            local_gama_ = (g[i][j + 1][k] - g[i][j][k]) / alpha_[i][j][k];
          }
          z = local_LAMDA_[i][j][k] + local_gama_;
          if (fabs(z) < DELTA) {
            local_Q_[k] = 0.5 * (z * z / DELTA + DELTA);
          } else
            local_Q_[k] = fabs(z);
        }
        for (k = 0; k <= 3; k++) theta[i][j][k] = 0;
        for (l = 0; l <= 3; l++)
          theta[i][j][k] += local_R_[k][l] * (g[i][j][l] + g[i][j + 1][l] - local_Q_[l] * alpha_[i][j][l]);
      }

  //计算G
  for (i = 2; i <= Nx + 4; i++)
    for (j = 2; j <= Ny + 4; j++)
      if (U[i][j][0] != 0) U2G(U[i][j], G[i][j]);
  //计算G_
  for (i = 2; i <= Nx + 3; i++)
    for (j = 2; j <= Ny + 3; j++)
      if (U[i][j][0] != 0)
        for (k = 0; k <= 3; k++) G_[i][j][k] = 0.5 * (G[i][j][k] + G[i][j + 1][k]) + 0.5 * theta[i][j][k];
  //分区计算U
  for (i = 3; i <= Nx + 3; i++)
    for (j = 3; j <= int(0.5 / dy) + 3; j++)
      for (k = 0; k <= 3; k++) U[i][j][k] = U[i][j][k] - r * (G_[i][j][k] - G_[i][j - 1][k]);

  for (i = int(1.0 / dx) + 3; i <= int(2.0 / dx) + 3; i++)
    for (j = int(0.5 / dy) + 4; j <= Ny + 3; j++)
      for (k = 0; k <= 3; k++) U[i][j][k] = U[i][j][k] - r * (G_[i][j][k] - G_[i][j - 1][k]);
}

#endif  // UPWINDTVD_FUSION
#endif