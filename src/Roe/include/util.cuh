#ifndef __UTIL_H
#define __UTIL_H

#include <sys/time.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <unordered_map>
#include <vector>
using std::vector;

#define ENO_FUSION
#define MUSCL_FUSION
#define ROE_FUSION
#define UPWINDTVD_FUSION
#define SYMTVD_FUSION
#define NND_FUSION
#define WENO_FUSION
#define COMPACT_FUSION

#define GAMA 1.4
#define PI 3.1415926
#define Lx 4.0
#define Ly 1.0
#define TT 0.1
#define Nx 400
#define Ny 100

#define ROECFL 0.4
#define MUSCLCFL 0.4
#define UPWINDTVDCFL 0.4
#define SYMTVDCFL 0.3
#define NNDCFL 0.4
#define ENOCFL 0.1
#define WENOCFL 0.2
#define COMPACTCFL 0.1

#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

#define get_elapsed_time_ms(_s, _e) (1000.0 * (_e.tv_sec - _s.tv_sec) + (_e.tv_usec - _s.tv_usec) / 1000.0)

class Tensor {
 public:
  vector<int> shape_;
  int shape_size_;
  vector<int> offset_;
  double* ptr_;
  double* cuda_ptr_;
  Tensor(const vector<int>& shape, int shape_size) : shape_(shape), shape_size_(shape_size) {
    for (int i = 1; i < shape_.size(); ++i) {
      int base = shape_[i];
      for (int j = i + 1; j < shape_.size(); ++j) base *= shape_[j];
      offset_.push_back(base);
    }
    for (auto& item : offset_) {
      printf("%d ", item);
    }
    printf("\n");
    ptr_ = new double[shape_size];
    cuda_ptr_ = new double[shape_size];
  }
  inline double* indexing(int i, int j) { return ptr_ + i * offset_[0] + j * offset_[1]; }
  inline double& indexing(int i, int j, int k) { return ptr_[i * offset_[0] + j * offset_[1]]; }
  inline double& indexing(int i, int j, int k, int f) { return ptr_[i * offset_[0] + j * offset_[1] + k * offset_[2]]; }
};

void save(const char* fn, double*** U) {
  std::ofstream out;
  out.open(fn, std::ofstream::binary);
  out.write(reinterpret_cast<const char*>(U), sizeof(double) * (Nx + 7) * (Ny + 7) * 4);
  out.close();
}

void check(const char* fn, double*** U) {
  double U_dat[Nx + 7][Ny + 7][4];
  std::ifstream in;
  in.open(fn, std::ifstream::binary);
  in.read(reinterpret_cast<char*>(U_dat), sizeof(double) * (Nx + 7) * (Ny + 7) * 4);
  in.close();
  for (int i = 0; i < Nx + 7; ++i) {
    for (int j = 0; j < Ny + 7; ++j) {
      for (int k = 0; k < 4; ++k) {
        if (fabs(U_dat[i][j][k] - U[i][j][k]) >= 1e-6) {
          printf("CHECK FAILED (%lf != %lf)\n", U_dat[i][j][k], U[i][j][k]);
          return;
        }
      }
    }
  }
  printf("CHECK PASS\n");
}

void check_Tensor(const char* fn, Tensor& U) {
  // double U_dat[Nx + 7][Ny + 7][4];
  Tensor U_dat({Nx + 7, Ny + 7, 4}, (Nx + 7) * (Ny + 7) * 4);
  std::ifstream in;
  in.open(fn, std::ifstream::binary);
  in.read((char*)U_dat.ptr_, sizeof(double) * (Nx + 7) * (Ny + 7) * 4);
  in.close();
  for (int i = 0; i < Nx + 7; ++i) {
    for (int j = 0; j < Ny + 7; ++j) {
      for (int k = 0; k < 4; ++k) {
        if (fabs(U_dat.indexing(i, j, k) - U.indexing(i, j, k)) >= 1e-6) {
          printf("CHECK FAILED (%lf != %lf)\n", U_dat.indexing(i, j, k) - U.indexing(i, j, k));
          return;
        }
      }
    }
  }
  printf("CHECK PASS\n");
}

//初始化函数
//作用:将全场赋初值
void initial(double*** U, double& dx, double& dy) {
  int i, j;
  dx = Lx / Nx;
  dy = Ly / Ny;
  double rou1 = 1, u1 = 0, v1 = 0, a1 = 1, p1 = 0.71429;
  double rou2 = 3.85714, p2 = 7.381, u2 = 2.22223, v2 = 0;
  for (i = 0; i <= Nx + 6; i++)  //初始条件  这里先把所有区域都赋值u1 v1 p1 rou1
    for (j = 0; j <= Ny + 6; j++) {
      U[i][j][0] = rou1;
      U[i][j][1] = rou1 * u1;
      U[i][j][2] = rou1 * v1;
      U[i][j][3] = p1 / (GAMA - 1) + 0.5 * rou1 * (u1 * u1 + v1 * v1);
    }
  //边界以外赋零
  for (i = 0; i <= int(1.0 / dx) - 1; i++)
    for (j = int(0.5 / dy) + 7; j <= Ny + 6; j++) {
      U[i][j][0] = 0;
      U[i][j][1] = 0;
      U[i][j][2] = 0;
      U[i][j][3] = 0;
    }
  for (i = int(2.0 / dx) + 7; i <= Nx + 6; i++)
    for (j = int(0.5 / dy) + 7; j <= Ny + 6; j++) {
      U[i][j][0] = 0;
      U[i][j][1] = 0;
      U[i][j][2] = 0;
      U[i][j][3] = 0;
    }
}

void initial_Tensor(Tensor& U, double& dx, double& dy) {
  int i, j;
  dx = Lx / Nx;
  dy = Ly / Ny;
  double rou1 = 1, u1 = 0, v1 = 0, a1 = 1, p1 = 0.71429;
  double rou2 = 3.85714, p2 = 7.381, u2 = 2.22223, v2 = 0;
  for (i = 0; i <= Nx + 6; i++)  //初始条件  这里先把所有区域都赋值u1 v1 p1 rou1
    for (j = 0; j <= Ny + 6; j++) {
      U.indexing(i, j, 0) = rou1;
      U.indexing(i, j, 1) = rou1 * u1;
      U.indexing(i, j, 2) = rou1 * v1;
      U.indexing(i, j, 3) = p1 / (GAMA - 1) + 0.5 * rou1 * (u1 * u1 + v1 * v1);
    }
  //边界以外赋零
  for (i = 0; i <= int(1.0 / dx) - 1; i++)
    for (j = int(0.5 / dy) + 7; j <= Ny + 6; j++) {
      U.indexing(i, j, 0) = 0;
      U.indexing(i, j, 1) = 0;
      U.indexing(i, j, 2) = 0;
      U.indexing(i, j, 3) = 0;
    }
  for (i = int(2.0 / dx) + 7; i <= Nx + 6; i++)
    for (j = int(0.5 / dy) + 7; j <= Ny + 6; j++) {
      U.indexing(i, j, 0) = 0;
      U.indexing(i, j, 1) = 0;
      U.indexing(i, j, 2) = 0;
      U.indexing(i, j, 3) = 0;
    }
}

// CFL稳定性条件
//入口:U(没有进行过Roe平均的U向量)
//出口:dt
double CFL(double*** U, double dx, double dy, double cfl) {
  int i, j;
  double u, v, rou, a, vel, maxvel, p;
  maxvel = 1e-100;
  for (i = 0; i <= Nx + 6; i++)
    for (j = 0; j <= Ny + 6; j++)
      if (U[i][j][0] != 0) {
        rou = U[i][j][0];
        u = U[i][j][1] / U[i][j][0];
        v = U[i][j][2] / U[i][j][0];
        p = (GAMA - 1) * (U[i][j][3] - 0.5 * rou * (u * u + v * v));
        a = sqrt(p * GAMA / rou);
        vel = a + fabs(u);
        if (vel >= maxvel) maxvel = vel;
        vel = a + fabs(v);
        if (vel >= maxvel) maxvel = vel;
      }
  return cfl * MIN(dx, dy) / maxvel;
}

double CFL_Tensor(Tensor& U, double dx, double dy, double cfl) {
  int i, j;
  double u, v, rou, a, vel, maxvel, p;
  maxvel = 1e-100;
  for (i = 0; i <= Nx + 6; i++)
    for (j = 0; j <= Ny + 6; j++)
      if (U.indexing(i, j, 0) != 0) {
        rou = U.indexing(i, j, 0);
        u = U.indexing(i, j, 1) / U.indexing(i, j, 0);
        v = U.indexing(i, j, 2) / U.indexing(i, j, 0);
        p = (GAMA - 1) * (U.indexing(i, j, 3) - 0.5 * rou * (u * u + v * v));
        a = sqrt(p * GAMA / rou);
        vel = a + fabs(u);
        if (vel >= maxvel) maxvel = vel;
        vel = a + fabs(v);
        if (vel >= maxvel) maxvel = vel;
      }
  return cfl * MIN(dx, dy) / maxvel;
}

// CFL_稳定性条件
//入口:U_是进行了Roe平均的向量),a_是声速
//出口:dt
double CFL_(double*** U_, double*** a_, double dx, double dy, double cfl) {
  int i, j;
  double u, v, rou, a, vel, maxvel;
  maxvel = 1e-100;
  for (i = 0; i <= Nx + 6; i++)
    for (j = 0; j <= Ny + 6; j++)
      if (U_[i][j][0] != 0) {
        rou = U_[i][j][0];
        u = U_[i][j][1];
        v = U_[i][j][2];
        a = a_[i][j][0];
        vel = a + sqrt(u * u + v * v);
        if (vel >= maxvel) maxvel = vel;
      }
  return cfl * MIN(dx, dy) / maxvel;
}

//边界函数
//入口:U
//出口:无
//作用:将左右边界以及虚拟节点赋值,并且处理特殊的凹凸角点以及虚拟节点,另外对边界上面的法相速度强制赋零
void bound(double*** U, double dx, double dy) {
  int i, j, k;
  double rou1 = 1, u1 = 0, v1 = 0, a1 = 1, p1 = 0.71429;
  double rou2 = 3.85714, p2 = 7.381, u2 = 2.22223, v2 = 0;
  for (i = 0; i <= 2; i++)  //左边界条件
    for (j = 0; j <= int(0.5 / dy) + 6; j++) {
      U[i][j][0] = rou2;
      U[i][j][1] = rou2 * u2;
      U[i][j][2] = rou2 * v2;
      U[i][j][3] = p2 / (GAMA - 1) + 0.5 * rou2 * (u2 * u2 + v2 * v2);
    }
  for (i = Nx + 4; i <= Nx + 6; i++)  //右边界条件
    for (j = 0; j <= int(0.5 / dy) + 6; j++) {
      U[i][j][0] = U[i - 1][j][0];
      U[i][j][1] = U[i - 1][j][1];
      U[i][j][2] = U[i - 1][j][2];
      U[i][j][3] = U[i - 1][j][3];
    }
  for (i = 3; i <= Nx + 3; i++)  //虚拟节点的处理,这里采用镜面反射,标量和切向速度取相同值,法相速度取相反值
    for (j = 0; j <= 2; j++) {
      U[i][j][0] = U[i][6 - j][0];
      U[i][j][1] = U[i][6 - j][1];
      U[i][j][2] = -U[i][6 - j][2];
      U[i][j][3] = U[i][6 - j][3];
    }
  for (i = 3; i <= int(1.0 / dx) + 2; i++)
    for (j = int(0.5 / dy) + 4; j <= int(0.5 / dy) + 6; j++) {
      U[i][j][0] = U[i][-j + 2 * int(0.5 / dy) + 6][0];
      U[i][j][1] = U[i][-j + 2 * int(0.5 / dy) + 6][1];
      U[i][j][2] = -U[i][-j + 2 * int(0.5 / dy) + 6][2];
      U[i][j][3] = U[i][-j + 2 * int(0.5 / dy) + 6][3];
    }
  for (i = int(1.0 / dx); i <= int(2.0 / dx) + 6; i++)
    for (j = Ny + 4; j <= Ny + 6; j++) {
      U[i][j][0] = U[i][-j + 2 * int(1.0 / dy) + 6][0];
      U[i][j][1] = U[i][-j + 2 * int(1.0 / dy) + 6][1];
      U[i][j][2] = -U[i][-j + 2 * int(1.0 / dy) + 6][2];
      U[i][j][3] = U[i][-j + 2 * int(1.0 / dy) + 6][3];
    }
  for (i = int(2.0 / dx) + 4; i <= Nx + 3; i++)
    for (j = int(0.5 / dy) + 4; j <= int(0.5 / dy) + 6; j++) {
      U[i][j][0] = U[i][-j + 2 * int(0.5 / dy) + 6][0];
      U[i][j][1] = U[i][-j + 2 * int(0.5 / dy) + 6][1];
      U[i][j][2] = -U[i][-j + 2 * int(0.5 / dy) + 6][2];
      U[i][j][3] = U[i][-j + 2 * int(0.5 / dy) + 6][3];
    }
  for (i = int(1.0 / dx); i <= int(1.0 / dx) + 2; i++)
    for (j = int(0.5 / dy) + 4; j <= Ny + 6; j++) {
      U[i][j][0] = U[2 * int(1.0 / dx) + 6 - i][j][0];
      U[i][j][1] = -U[2 * int(1.0 / dx) + 6 - i][j][1];
      U[i][j][2] = U[2 * int(1.0 / dx) + 6 - i][j][2];
      U[i][j][3] = U[2 * int(1.0 / dx) + 6 - i][j][3];
    }
  for (i = int(2.0 / dx) + 4; i <= int(2.0 / dx) + 6; i++)
    for (j = int(0.5 / dy) + 4; j <= Ny + 6; j++) {
      U[i][j][0] = U[2 * int(2.0 / dx) + 6 - i][j][0];
      U[i][j][1] = -U[2 * int(2.0 / dx) + 6 - i][j][1];
      U[i][j][2] = U[2 * int(2.0 / dx) + 6 - i][j][2];
      U[i][j][3] = U[2 * int(2.0 / dx) + 6 - i][j][3];
    }
  //角点的虚拟网格特殊处理
  int a, b;
  a = int(1.0 / dx) + 2;  //角凸点(a,b)以及其附近虚拟节点处理,处理方法参考了教材相关内容
  b = int(0.5 / dy) + 4;
  for (k = 0; k <= 3; k++) {
    U[a][b][k] = 0.5 * (U[a][b - 2][k] + U[a + 2][b][k]);
    U[a - 1][b][k] = U[a - 1][b - 2][k];
    U[a - 2][b][k] = U[a - 2][b - 2][k];
    U[a][b + 1][k] = U[a + 2][b + 1][k];
    U[a][b + 2][k] = U[a + 2][b + 2][k];
    U[a - 1][b + 1][k] = 0.5 * (U[a - 1][b - 3][k] + U[a + 3][b + 1][k]);
    U[a - 2][b + 1][k] = U[a - 2][b - 3][k];
    U[a - 1][b + 2][k] = U[a + 3][b + 2][k];
    U[a - 2][b + 2][k] = 0.5 * (U[a - 2][b - 4][k] + U[a + 4][b + 2][k]);
  }
  U[a - 1][b][2] = -U[a - 1][b - 2][2];
  U[a - 2][b][2] = -U[a - 2][b - 2][2];
  U[a - 2][b + 1][2] = -U[a - 2][b - 3][2];
  U[a][b + 1][1] = -U[a + 2][b + 1][1];
  U[a][b + 2][1] = -U[a + 2][b + 2][1];
  U[a - 1][b + 2][1] = -U[a + 3][b + 2][1];
  a = int(2.0 / dx) + 4;
  b = int(0.5 / dy) + 4;
  for (k = 0; k <= 3; k++) {
    U[a][b][k] = 0.5 * (U[a][b - 2][k] + U[a - 2][b][k]);
    U[a + 1][b][k] = U[a + 1][b - 2][k];
    U[a + 2][b][k] = U[a + 2][b - 2][k];
    U[a][b + 1][k] = U[a - 2][b + 1][k];
    U[a][b + 2][k] = U[a - 2][b + 2][k];
    U[a + 1][b + 1][k] = 0.5 * (U[a + 1][b - 3][k] + U[a - 3][b + 1][k]);
    U[a + 2][b + 1][k] = U[a + 2][b - 3][k];
    U[a + 1][b + 2][k] = U[a - 3][b + 2][k];
    U[a + 2][b + 2][k] = 0.5 * (U[a + 2][b - 4][k] + U[a - 4][b + 2][k]);
  }
  U[a + 1][b][2] = -U[a + 1][b - 2][2];
  U[a + 2][b][2] = -U[a + 2][b - 2][2];
  U[a + 2][b + 1][2] = -U[a + 2][b - 3][2];
  U[a][b + 1][1] = -U[a - 2][b + 1][1];
  U[a][b + 2][1] = -U[a - 2][b + 2][1];
  U[a + 1][b + 2][1] = -U[a - 3][b + 2][1];
  //边界上强行赋法向为零
  for (i = 3; i <= Nx + 3; i++) {
    U[i][3][2] = 0;
  }
  for (i = 3; i <= int(1.0 / dx) + 2; i++) {
    U[i][int(0.5 / dy) + 3][2] = 0;
  }
  for (i = int(2.0 / dx) + 4; i <= Nx + 3; i++) {
    U[i][int(0.5 / dy) + 3][2] = 0;
  }
  for (i = int(1.0 / dx) + 3; i <= int(2.0 / dx) + 3; i++) {
    U[i][Ny + 3][2] = 0;
  }
  for (j = int(0.5 / dy) + 4; j <= Ny + 3; j++) {
    U[int(1.0 / dx) + 3][j][1] = 0;
  }
  for (j = int(0.5 / dy) + 4; j <= Ny + 3; j++) {
    U[int(2.0 / dx) + 3][j][1] = 0;
  }
}
void bound_Tensor(Tensor& U, double dx, double dy) {
  int i, j, k;
  double rou1 = 1, u1 = 0, v1 = 0, a1 = 1, p1 = 0.71429;
  double rou2 = 3.85714, p2 = 7.381, u2 = 2.22223, v2 = 0;
  for (i = 0; i <= 2; i++)  //左边界条件
    for (j = 0; j <= int(0.5 / dy) + 6; j++) {
      U.indexing(i, j, 0) = rou2;
      U.indexing(i, j, 1) = rou2 * u2;
      U.indexing(i, j, 2) = rou2 * v2;
      U.indexing(i, j, 3) = p2 / (GAMA - 1) + 0.5 * rou2 * (u2 * u2 + v2 * v2);
    }
  for (i = Nx + 4; i <= Nx + 6; i++)  //右边界条件
    for (j = 0; j <= int(0.5 / dy) + 6; j++) {
      U.indexing(i, j, 0) = U.indexing(i - 1, j, 0);
      U.indexing(i, j, 1) = U.indexing(i - 1, j, 1);
      U.indexing(i, j, 2) = U.indexing(i - 1, j, 2);
      U.indexing(i, j, 3) = U.indexing(i - 1, j, 3);
    }
  for (i = 3; i <= Nx + 3; i++)  //虚拟节点的处理,这里采用镜面反射,标量和切向速度取相同值,法相速度取相反值
    for (j = 0; j <= 2; j++) {
      U.indexing(i, j, 0) = U.indexing(i, 6 - j, 0);
      U.indexing(i, j, 1) = U.indexing(i, 6 - j, 1);
      U.indexing(i, j, 2) = -U.indexing(i, 6 - j, 2);
      U.indexing(i, j, 3) = U.indexing(i, 6 - j, 3);
    }
  for (i = 3; i <= int(1.0 / dx) + 2; i++)
    for (j = int(0.5 / dy) + 4; j <= int(0.5 / dy) + 6; j++) {
      U.indexing(i, j, 0) = U.indexing(i, -j + 2 * int(0.5 / dy) + 6, 0);
      U.indexing(i, j, 1) = U.indexing(i, -j + 2 * int(0.5 / dy) + 6, 1);
      U.indexing(i, j, 2) = -U.indexing(i, -j + 2 * int(0.5 / dy) + 6, 2);
      U.indexing(i, j, 3) = U.indexing(i, -j + 2 * int(0.5 / dy) + 6, 3);
    }
  for (i = int(1.0 / dx); i <= int(2.0 / dx) + 6; i++)
    for (j = Ny + 4; j <= Ny + 6; j++) {
      U.indexing(i, j, 0) = U.indexing(i, -j + 2 * int(1.0 / dy) + 6, 0);
      U.indexing(i, j, 1) = U.indexing(i, -j + 2 * int(1.0 / dy) + 6, 1);
      U.indexing(i, j, 2) = -U.indexing(i, -j + 2 * int(1.0 / dy) + 6, 2);
      U.indexing(i, j, 3) = U.indexing(i, -j + 2 * int(1.0 / dy) + 6, 3);
    }
  for (i = int(2.0 / dx) + 4; i <= Nx + 3; i++)
    for (j = int(0.5 / dy) + 4; j <= int(0.5 / dy) + 6; j++) {
      U.indexing(i, j, 0) = U.indexing(i, -j + 2 * int(0.5 / dy) + 6, 0);
      U.indexing(i, j, 1) = U.indexing(i, -j + 2 * int(0.5 / dy) + 6, 1);
      U.indexing(i, j, 2) = -U.indexing(i, -j + 2 * int(0.5 / dy) + 6, 2);
      U.indexing(i, j, 3) = U.indexing(i, -j + 2 * int(0.5 / dy) + 6, 3);
    }
  for (i = int(1.0 / dx); i <= int(1.0 / dx) + 2; i++)
    for (j = int(0.5 / dy) + 4; j <= Ny + 6; j++) {
      U.indexing(i, j, 0) = U.indexing(2 * int(1.0 / dx) + 6 - i, j, 0);
      U.indexing(i, j, 1) = -U.indexing(2 * int(1.0 / dx) + 6 - i, j, 1);
      U.indexing(i, j, 2) = U.indexing(2 * int(1.0 / dx) + 6 - i, j, 2);
      U.indexing(i, j, 3) = U.indexing(2 * int(1.0 / dx) + 6 - i, j, 3);
    }
  for (i = int(2.0 / dx) + 4; i <= int(2.0 / dx) + 6; i++)
    for (j = int(0.5 / dy) + 4; j <= Ny + 6; j++) {
      U.indexing(i, j, 0) = U.indexing(2 * int(2.0 / dx) + 6 - i, j, 0);
      U.indexing(i, j, 1) = -U.indexing(2 * int(2.0 / dx) + 6 - i, j, 1);
      U.indexing(i, j, 2) = U.indexing(2 * int(2.0 / dx) + 6 - i, j, 2);
      U.indexing(i, j, 3) = U.indexing(2 * int(2.0 / dx) + 6 - i, j, 3);
    }
  //角点的虚拟网格特殊处理
  int a, b;
  a = int(1.0 / dx) + 2;  //角凸点(a,b)以及其附近虚拟节点处理,处理方法参考了教材相关内容
  b = int(0.5 / dy) + 4;
  for (k = 0; k <= 3; k++) {
    U.indexing(a, b, k) = 0.5 * (U.indexing(a, b - 2, k) + U.indexing(a + 2, b, k));
    U.indexing(a - 1, b, k) = U.indexing(a - 1, b - 2, k);
    U.indexing(a - 2, b, k) = U.indexing(a - 2, b - 2, k);
    U.indexing(a, b + 1, k) = U.indexing(a + 2, b + 1, k);
    U.indexing(a, b + 2, k) = U.indexing(a + 2, b + 2, k);
    U.indexing(a - 1, b + 1, k) = 0.5 * (U.indexing(a - 1, b - 3, k) + U.indexing(a + 3, b + 1, k));
    U.indexing(a - 2, b + 1, k) = U.indexing(a - 2, b - 3, k);
    U.indexing(a - 1, b + 2, k) = U.indexing(a + 3, b + 2, k);
    U.indexing(a - 2, b + 2, k) = 0.5 * (U.indexing(a - 2, b - 4, k) + U.indexing(a + 4, b + 2, k));
  }
  U.indexing(a - 1, b, 2) = -U.indexing(a - 1, b - 2, 2);
  U.indexing(a - 2, b, 2) = -U.indexing(a - 2, b - 2, 2);
  U.indexing(a - 2, b + 1, 2) = -U.indexing(a - 2, b - 3, 2);
  U.indexing(a, b + 1, 1) = -U.indexing(a + 2, b + 1, 1);
  U.indexing(a, b + 2, 1) = -U.indexing(a + 2, b + 2, 1);
  U.indexing(a - 1, b + 2, 1) = -U.indexing(a + 3, b + 2, 1);
  a = int(2.0 / dx) + 4;
  b = int(0.5 / dy) + 4;
  for (k = 0; k <= 3; k++) {
    U.indexing(a, b, k) = 0.5 * (U.indexing(a, b - 2, k) + U.indexing(a - 2, b, k));
    U.indexing(a + 1, b, k) = U.indexing(a + 1, b - 2, k);
    U.indexing(a + 2, b, k) = U.indexing(a + 2, b - 2, k);
    U.indexing(a, b + 1, k) = U.indexing(a - 2, b + 1, k);
    U.indexing(a, b + 2, k) = U.indexing(a - 2, b + 2, k);
    U.indexing(a + 1, b + 1, k) = 0.5 * (U.indexing(a + 1, b - 3, k) + U.indexing(a - 3, b + 1, k));
    U.indexing(a + 2, b + 1, k) = U.indexing(a + 2, b - 3, k);
    U.indexing(a + 1, b + 2, k) = U.indexing(a - 3, b + 2, k);
    U.indexing(a + 2, b + 2, k) = 0.5 * (U.indexing(a + 2, b - 4, k) + U.indexing(a - 4, b + 2, k));
  }
  U.indexing(a + 1, b, 2) = -U.indexing(a + 1, b - 2, 2);
  U.indexing(a + 2, b, 2) = -U.indexing(a + 2, b - 2, 2);
  U.indexing(a + 2, b + 1, 2) = -U.indexing(a + 2, b - 3, 2);
  U.indexing(a, b + 1, 1) = -U.indexing(a - 2, b + 1, 1);
  U.indexing(a, b + 2, 1) = -U.indexing(a - 2, b + 2, 1);
  U.indexing(a + 1, b + 2, 1) = -U.indexing(a - 3, b + 2, 1);
  //边界上强行赋法向为零
  for (i = 3; i <= Nx + 3; i++) {
    U.indexing(i, 3, 2) = 0;
  }
  for (i = 3; i <= int(1.0 / dx) + 2; i++) {
    U.indexing(i, int(0.5 / dy) + 3, 2) = 0;
  }
  for (i = int(2.0 / dx) + 4; i <= Nx + 3; i++) {
    U.indexing(i, int(0.5 / dy) + 3, 2) = 0;
  }
  for (i = int(1.0 / dx) + 3; i <= int(2.0 / dx) + 3; i++) {
    U.indexing(i, Ny + 3, 2) = 0;
  }
  for (j = int(0.5 / dy) + 4; j <= Ny + 3; j++) {
    U.indexing(int(1.0 / dx) + 3, j, 1) = 0;
  }
  for (j = int(0.5 / dy) + 4; j <= Ny + 3; j++) {
    U.indexing(int(2.0 / dx) + 3, j, 1) = 0;
  }
}

//求符号函数
//出口:大于零返回1,小于零返回-1,等于零返回0
double sign(double va) {
  if (va > 0) {
    return 1;
  }
  if (va < 0) {
    return -1;
  } else {
    return 0;
  }
}
//返回3个数中的较少者
double MIN3(double a1, double a2, double a3) {
  double re;
  if (a1 <= a2 && a1 <= a3)
    re = a1;
  else {
    if (a2 <= a3 && a2 <= a1)
      re = a2;
    else
      re = a3;
  }
  return re;
}
// minmod
double minmod(double w1, double w2) {
  double result;
  if (w1 * w2 > 0) {
    result = sign(w1) * MIN(fabs(w1), fabs(w2));
  } else
    result = 0;
  return result;
}
//对称TVD中用到的3个数的minmod
double minmod3(double w1, double w2, double w3) {
  double result;
  if (w1 * w2 * w3 > 0) {
    result = sign(w1) * MIN3(fabs(w1), fabs(w2), fabs(w3));
  } else
    result = 0;
  return result;
}
//通过U求F
void U2F(double U[4], double F[4]) {
  double u, v, p, rou;
  rou = U[0];
  u = U[1] / U[0];
  v = U[2] / U[0];
  p = (GAMA - 1) * (U[3] - 0.5 * rou * (u * u + v * v));
  F[0] = rou * u;
  F[1] = rou * u * u + p;
  F[2] = rou * u * v;
  F[3] = (U[3] + p) * u;
}
void U2G(double U[4], double G[4]) {
  double u, v, p, rou;
  rou = U[0];
  u = U[1] / U[0];
  v = U[2] / U[0];
  p = (GAMA - 1) * (U[3] - 0.5 * rou * (u * u + v * v));
  G[0] = rou * v;
  G[1] = rou * u * v;
  G[2] = rou * v * v + p;
  G[3] = (U[3] + p) * v;
}

void virtual_clear(double*** U, double dx, double dy) {
  int i, j, k;
  for (i = 3; i <= int(1.0 / dx) + 2; i++)  //为了显示正确,首先将虚拟节点赋0
    for (j = int(0.5 / dy) + 4; j <= int(0.5 / dy) + 6; j++)
      for (k = 0; k <= 3; k++) U[i][j][k] = 0;

  for (i = int(2.0 / dx) + 4; i <= Nx + 6; i++)
    for (j = int(0.5 / dy) + 4; j <= int(0.5 / dy) + 6; j++)
      for (k = 0; k <= 3; k++) U[i][j][k] = 0;

  for (i = int(1.0 / dx); i <= int(1.0 / dx) + 2; i++)
    for (j = int(0.5 / dy) + 4; j <= Ny + 6; j++)
      for (k = 0; k <= 3; k++) U[i][j][k] = 0;

  for (i = int(2.0 / dx) + 4; i <= int(2.0 / dx) + 6; i++)
    for (j = int(0.5 / dy) + 4; j <= Ny + 6; j++)
      for (k = 0; k <= 3; k++) U[i][j][k] = 0;
}

void virtual_clear_Tensor(Tensor& U, double dx, double dy) {
  int i, j, k;
  for (i = 3; i <= int(1.0 / dx) + 2; i++)  //为了显示正确,首先将虚拟节点赋0
    for (j = int(0.5 / dy) + 4; j <= int(0.5 / dy) + 6; j++)
      for (k = 0; k <= 3; k++) U.indexing(i, j, k) = 0;

  for (i = int(2.0 / dx) + 4; i <= Nx + 6; i++)
    for (j = int(0.5 / dy) + 4; j <= int(0.5 / dy) + 6; j++)
      for (k = 0; k <= 3; k++) U.indexing(i, j, k) = 0;

  for (i = int(1.0 / dx); i <= int(1.0 / dx) + 2; i++)
    for (j = int(0.5 / dy) + 4; j <= Ny + 6; j++)
      for (k = 0; k <= 3; k++) U.indexing(i, j, k) = 0;

  for (i = int(2.0 / dx) + 4; i <= int(2.0 / dx) + 6; i++)
    for (j = int(0.5 / dy) + 4; j <= Ny + 6; j++)
      for (k = 0; k <= 3; k++) U.indexing(i, j, k) = 0;
}

struct timer {
  std::unordered_map<const char*, double> elapsed_time_ms;
  std::unordered_map<const char*, int> count;
  std::unordered_map<const char*, timeval> time_point;

  void start(const char* func) {
    timeval s;
    gettimeofday(&s, NULL);
    time_point[func] = s;
  }
  void stop(const char* func) {
    timeval e;
    gettimeofday(&e, NULL);
    count[func]++;
    timeval s = time_point[func];
    elapsed_time_ms[func] += get_elapsed_time_ms(s, e);
  }
  void show_all() {
    for (auto it = elapsed_time_ms.begin(); it != elapsed_time_ms.end(); ++it) {
      auto func = it->first;
      double t = it->second;
      int c = count[func];
      printf("%s: %lf ms / %d count, avg %lf ms\n", func, t, c, t / c);
    }
  }
};

#endif