// MAC-Chorin2D.cpp : 定义控制台应用程序的入口点。
/*-----------------------------------------------------------------------------
利用MAC算法和Chorin压力迭代解法求解二维不可压缩黏性流动问题(C语言)   
-----------------------------------------------------------------------------*/
#include "stdio.h"
#include "cmath"
#include <string>

using namespace std;

#define Jx 400         //x方向网格点数
#define Jy 200          //y方向网格点数
#define Fx1 100        //第一方块左下角点外x位置
#define Fy1 80         //方块左下角点外y位置
#define Fy2 40         //第二方块（下）左下角点外y位置
#define Fy3 120         //第三方块（上）左下角点外y位置
#define Fx2 200        //第二、三方块左下角点外x位置
#define Kx 20         //方块X跨度
#define Ky 40         //方块y跨度
#define Re 100.0         //Reynolds数
#define dt 0.001      //时间步长
#define EPS 1e-4       //收敛限
//全局变量
double u[Jx+2][Jy+2],v[Jx+2][Jy+2],p[Jx+2][Jy+2],tu[Jx+2][Jy+2],tv[Jx+2][Jy+2];
/*---------------------------------------------------------
速度边界条件
入口: 给定值，自由输入条件；
出口: 自由输出条件；
上下边界及方柱为无滑移固壁，
-----------------------------------------------------------*/
void bounduv(double u[Jx+2][Jy+2],double v[Jx+2][Jy+2])
{
 int i,j;
 double dx,dy;
 dx=4.0/Jx;    //x方向网格间距
 dy=2.0/Jy;    //y方向网格间距
//初始速度分布
for(j=1;j<=Jy;j++)
  {
	if(j>=1&&j<=10)
		u[0][j]=10*(j*dy-0.5*dy);
		
	else if (j>=11&&j<=189)
		u[0][j]=1.0;

	else 
		u[0][j]=10*(2.0-(j*dy-0.5*dy));
  } 
//上下壁面边界
for(i=0;i<=Jx+1;i++)
 {
  u[i][0]=-u[i][1];
  u[i][Jy+1]=-u[i][Jy];
 }

for(i=0;i<=Jx+1;i++)
 {
  v[i][0]=0;
  v[i][Jy]=0;
 }
//进、出口条件 
for(j=0;j<=Jy+1;j++)
{
  u[1][j]=u[0][j];
  u[Jx+1][j]=u[Jx-1][j];
  u[Jx][j]=u[Jx-1][j];
}
for(j=0;j<=Jy+1;j++)
 {
  v[1][j]=0;
  v[0][j]=v[1][j];
  v[Jx+1][j]=v[Jx][j];
 }
 //方块上下边界
 for(i=Fx1+1;i<=Fx1+Kx-1;i++)
 {
	 u[i][Fy1+1]=-u[i][Fy1];
	 u[i][Fy1+Ky]=-u[i][Fy1+Ky+1];
 }
 for(i=Fx1+1;i<=Fx1+Kx;i++)
 {
	 v[i][Fy1]=0;
	 v[i][Fy1+Ky]=0;
 }

for(i=Fx2+1;i<=Fx2+Kx-1;i++)
 {
	 u[i][Fy2+1]=-u[i][Fy2];
	 u[i][Fy2+Ky]=-u[i][Fy2+Ky+1];
	 u[i][Fy3+1]=-u[i][Fy3];
	 u[i][Fy3+Ky]=-u[i][Fy3+Ky+1];
 }
 for(i=Fx2+1;i<=Fx2+Kx;i++)
 {
	 v[i][Fy2]=0;
	 v[i][Fy2+Ky]=0;
	 v[i][Fy3]=0;
	 v[i][Fy3+Ky]=0;
 }
//方块左右边界
 for(j=Fy1+1;j<=Fy1+Ky;j++)
 {
	 u[Fx1][j]=0;
	 u[Fx1+Kx][j]=0;
 }
 for(j=Fy1+1;j<=Fy1+Ky-1;j++)
 {
	 v[Fx1+1][j]=-v[Fx1][j];
	 v[Fx1+Kx][j]=-v[Fx1+Kx+1][j];
 }

for(j=Fy2+1;j<=Fy2+Ky;j++)
 {
	 u[Fx2][j]=0;
	 u[Fx2+Kx][j]=0;
 }
 for(j=Fy2+1;j<=Fy2+Ky-1;j++)
 {
	 v[Fx2+1][j]=-v[Fx2][j];
	 v[Fx2+Kx][j]=-v[Fx2+Kx+1][j];
 }

 for(j=Fy3+1;j<=Fy3+Ky;j++)
 {
	 u[Fx2][j]=0;
	 u[Fx2+Kx][j]=0;
 }
 for(j=Fy3+1;j<=Fy3+Ky-1;j++)
 {
	 v[Fx2+1][j]=-v[Fx2][j];
	 v[Fx2+Kx][j]=-v[Fx2+Kx+1][j];
 }
//对称边界条件
for(i=0;i<=Jx+1;i++)
	v[i][Jy/2]=0;
for(i=0;i<=Jx+1;i++)
for(j=0;j<=Jy/2;j++)
	u[i][j]=u[i][-j+1+Jy];
for(i=0;i<=Jx+1;i++)
for(j=0;j<=Jy/2-1;j++)
	v[i][j]=-v[i][-j+Jy];	
}
/*----------------------------------------------------------
压力边界条件
-----------------------------------------------------------*/ 
void boundp(double p[Jx+2][Jy+2])
{
 int i,j;
 //左右
 for(j=1;j<=Jy;j++)
 {
  p[1][j]=1.0;
  p[0][j]=p[1][j];
  p[Jx+1][j]=p[Jx][j];
  if(j>=Fy1+2&&j<=Fy1+Ky-1)
  {
	p[Fx1+1][j]=p[Fx1][j];
	p[Fx1+Kx][j]=p[Fx1+Kx+1][j];
  }
  if((j>=Fy2+2&&j<=Fy2+Ky-1)||(j>=Fy3+2&&j<=Fy3+Ky-1))
  {
	p[Fx2+1][j]=p[Fx2][j];
	p[Fx2+Kx][j]=p[Fx2+Kx+1][j];
  }
 }
 //上下
 for(i=0;i<=Jx+1;i++)
 {
    p[i][0]=p[i][1];
    p[i][Jy+1]=p[i][Jy];
  if(i>=Fx1+2&&i<=Fx1+Kx-1)
  {
    p[i][Fy1+1]=p[i][Fy1];
    p[i][Fy1+Ky]=p[i][Fy1+Ky+1];
  }
  if(i>=Fx2+2&&i<=Fx2+Kx-1)
  {
    p[i][Fy2+1]=p[i][Fy2];
    p[i][Fy2+Ky]=p[i][Fy2+Ky+1];
	p[i][Fy3+1]=p[i][Fy3];
    p[i][Fy3+Ky]=p[i][Fy3+Ky+1];
  }
 }
}
/*----------------------------------------------------------
方柱角点的压力处理
-----------------------------------------------------------*/ 
void boundjd(double p[Jx+2][Jy+2])	
{
	p[Fx1+1][Fy1+1]=(p[Fx1+1][Fy1]+p[Fx1][Fy1+1])/2;
	p[Fx1+1][Fy1+Ky]=(p[Fx1][Fy1+Ky]+p[Fx1+1][Fy1+Ky+1])/2;
	p[Fx1+Kx][Fy1+1]=(p[Fx1+Kx][Fy1]+p[Fx1+Kx+1][Fy1+1])/2;
	p[Fx1+Kx][Fy1+Ky]=(p[Fx1+Kx][Fy1+Ky+1]+p[Fx1+Kx+1][Fy1+Ky])/2;

	p[Fx2+1][Fy2+1]=(p[Fx2][Fy2+1]+p[Fx2+1][Fy2])/2;
	p[Fx2+1][Fy2+Ky]=(p[Fx2][Fy2+Ky]+p[Fx2+1][Fy2+Ky+1])/2;
	p[Fx2+Kx][Fy2+1]=(p[Fx2+Kx][Fy2]+p[Fx2+Kx+1][Fy2+1])/2;
	p[Fx2+Kx][Fy2+Ky]=(p[Fx2+Kx+1][Fy2+Ky]+p[Fx2+Kx][Fy2+Ky+1])/2;

	p[Fx2+1][Fy3+1]=(p[Fx2][Fy3+1]+p[Fx2+1][Fy3])/2;
	p[Fx2+1][Fy3+Ky]=(p[Fx2][Fy3+Ky]+p[Fx2+1][Fy3+Ky+1])/2;
	p[Fx2+Kx][Fy3+1]=(p[Fx2+Kx][Fy3]+p[Fx2+Kx+1][Fy3+1])/2;
	p[Fx2+Kx][Fy3+Ky]=(p[Fx2+Kx+1][Fy3+Ky]+p[Fx2+Kx][Fy3+Ky+1])/2;
}
/*-------------------------------------------------------
初始化：初始速度为0，初始压力为1.0。
---------------------------------------------------------*/
void init(double u[Jx+2][Jy+2],double v[Jx+2][Jy+2],double p[Jx+2][Jy+2],double dx,double dy)
{
 int i,j;
 for(i=0;i<=Jx+1;i++)
  for(j=0;j<=Jy+1;j++)
{
   u[i][j]=0;
   v[i][j]=0;
   p[i][j]=1.0;
 }
 bounduv(u,v);
 boundp(p);
 boundjd(p);
}
/*-------------------------------------------------------
根据动量方程求解u、v
---------------------------------------------------------*/
void solveuv(double u[Jx+2][Jy+2],double v[Jx+2][Jy+2],double p[Jx+2][Jy+2],double dx,double dy,double tu[Jx+2][Jy+2],double tv[Jx+2][Jy+2])
{
 double vav,uav,adv,vis,prs;
 int i,j;
 dx=4.0/Jx;    //x方向网格间距
 dy=2.0/Jy;    //y方向网格间距
 //第一部分
 for(i=1;i<=Fx1-1;i++)
  for(j=1;j<=Jy;j++)
  {
   vav=(v[i][j]+v[i+1][j]+v[i][j-1]+v[i+1][j-1])/4.0;
   adv=-u[i][j]*(u[i+1][j]-u[i-1][j])/dx/2.0-vav*(u[i][j+1]-u[i][j-1])/dy/2.0;
   vis=((u[i+1][j]-2*u[i][j]+u[i-1][j])/dx/dx+(u[i][j+1]-2*u[i][j]+u[i][j-1])/dy/dy)/Re;
   prs=-(p[i+1][j]-p[i][j])/dx;
   tu[i][j]=(adv+vis+prs)*dt+u[i][j];
  }
 for(i=1;i<=Fx1;i++)
  for(j=1;j<=Jy-1;j++)
  {
   uav=(u[i-1][j+1]+u[i][j+1]+u[i-1][j]+u[i][j])/4.0;
   adv=-uav*(v[i+1][j]-v[i-1][j])/dx/2.0-v[i][j]*(v[i][j+1]-v[i][j-1])/dy/2.0;
   vis=((v[i+1][j]-2*v[i][j]+v[i-1][j])/dx/dx+(v[i][j+1]-2*v[i][j]+v[i][j-1])/dy/dy)/Re;
   prs=-(p[i][j+1]-p[i][j])/dy;
   tv[i][j]=(adv+vis+prs)*dt+v[i][j];
  }

  for(i=1;i<=Fx1-1;i++)
  for(j=1;j<=Jy;j++)
   u[i][j]=tu[i][j];
  for(i=1;i<=Fx1;i++)
  for(j=1;j<=Jy-1;j++)
   v[i][j]=tv[i][j];
//第二部分下
  for(i=Fx1;i<=Fx1+Kx;i++)
   for(j=1;j<=Fy1;j++)
  {
   vav=(v[i][j]+v[i+1][j]+v[i][j-1]+v[i+1][j-1])/4.0;
   adv=-u[i][j]*(u[i+1][j]-u[i-1][j])/dx/2.0-vav*(u[i][j+1]-u[i][j-1])/dy/2.0;
   vis=((u[i+1][j]-2*u[i][j]+u[i-1][j])/dx/dx+(u[i][j+1]-2*u[i][j]+u[i][j-1])/dy/dy)/Re;
   prs=-(p[i+1][j]-p[i][j])/dx;
   tu[i][j]=(adv+vis+prs)*dt+u[i][j];
  }
  for(i=Fx1+1;i<=Fx1+Kx;i++)
  for(j=1;j<=Fy1-1;j++)
  {
   uav=(u[i-1][j+1]+u[i][j+1]+u[i-1][j]+u[i][j])/4.0;
  adv=-uav*(v[i+1][j]-v[i-1][j])/dx/2.0-v[i][j]*(v[i][j+1]-v[i][j-1])/dy/2.0;
   vis=((v[i+1][j]-2*v[i][j]+v[i-1][j])/dx/dx+(v[i][j+1]-2*v[i][j]+v[i][j-1])/dy/dy)/Re;
   prs=-(p[i][j+1]-p[i][j])/dy;
   tv[i][j]=(adv+vis+prs)*dt+v[i][j];
  } 

  for(i=Fx1;i<=Fx1+Kx;i++)
  for(j=1;j<=Fy1;j++)
   u[i][j]=tu[i][j];
  for(i=Fx1+1;i<=Fx1+Kx;i++)
  for(j=1;j<=Fy1-1;j++)
   v[i][j]=tv[i][j];
//第二部分上
  for(i=Fx1;i<=Fx1+Kx;i++)
   for(j=Fy1+Ky+1;j<=Jy;j++)
  {
   vav=(v[i][j]+v[i+1][j]+v[i][j-1]+v[i+1][j-1])/4.0;
  adv=-u[i][j]*(u[i+1][j]-u[i-1][j])/dx/2.0-vav*(u[i][j+1]-u[i][j-1])/dy/2.0;
   vis=((u[i+1][j]-2*u[i][j]+u[i-1][j])/dx/dx+(u[i][j+1]-2*u[i][j]+u[i][j-1])/dy/dy)/Re;
   prs=-(p[i+1][j]-p[i][j])/dx;
   tu[i][j]=(adv+vis+prs)*dt+u[i][j];
  }
  for(i=Fx1+1;i<=Fx1+Kx;i++)
  for(j=Fy1+Ky+1;j<=Jy-1;j++)
  {
   uav=(u[i-1][j+1]+u[i][j+1]+u[i-1][j]+u[i][j])/4.0;
   adv=-uav*(v[i+1][j]-v[i-1][j])/dx/2.0-v[i][j]*(v[i][j+1]-v[i][j-1])/dy/2.0;
   vis=((v[i+1][j]-2*v[i][j]+v[i-1][j])/dx/dx+(v[i][j+1]-2*v[i][j]+v[i][j-1])/dy/dy)/Re;
   prs=-(p[i][j+1]-p[i][j])/dy;
   tv[i][j]=(adv+vis+prs)*dt+v[i][j];
  } 

 for(i=Fx1;i<=Fx1+Kx;i++)
   for(j=Fy1+Ky+1;j<=Jy;j++)
   u[i][j]=tu[i][j];
  for(i=Fx1+1;i<=Fx1+Kx;i++)
  for(j=Fy1+Ky+1;j<=Jy-1;j++)
   v[i][j]=tv[i][j];
  //第三部分
  for(i=Fx1+Kx+1;i<=Fx2-1;i++)
  for(j=1;j<=Jy;j++)
  {
   vav=(v[i][j]+v[i+1][j]+v[i][j-1]+v[i+1][j-1])/4.0;
   adv=-u[i][j]*(u[i+1][j]-u[i-1][j])/dx/2.0-vav*(u[i][j+1]-u[i][j-1])/dy/2.0;
   vis=((u[i+1][j]-2*u[i][j]+u[i-1][j])/dx/dx+(u[i][j+1]-2*u[i][j]+u[i][j-1])/dy/dy)/Re;
   prs=-(p[i+1][j]-p[i][j])/dx;
   tu[i][j]=(adv+vis+prs)*dt+u[i][j];
  }
 for(i=Fx1+Kx+1;i<=Fx2;i++)
  for(j=1;j<=Jy-1;j++)
  {
   uav=(u[i-1][j+1]+u[i][j+1]+u[i-1][j]+u[i][j])/4.0;
  adv=-uav*(v[i+1][j]-v[i-1][j])/dx/2.0-v[i][j]*(v[i][j+1]-v[i][j-1])/dy/2.0;
   vis=((v[i+1][j]-2*v[i][j]+v[i-1][j])/dx/dx+(v[i][j+1]-2*v[i][j]+v[i][j-1])/dy/dy)/Re;
   prs=-(p[i][j+1]-p[i][j])/dy;
   tv[i][j]=(adv+vis+prs)*dt+v[i][j];
  } 

  for(i=Fx1+Kx+1;i<=Fx2-1;i++)
  for(j=1;j<=Jy;j++)
   u[i][j]=tu[i][j];
  for(i=Fx1+Kx+1;i<=Fx2;i++)
  for(j=1;j<=Jy-1;j++)
   v[i][j]=tv[i][j];
//第四部分下
for(i=Fx2;i<=Fx2+Kx;i++)
  for(j=1;j<=Fy2;j++)
  {
   vav=(v[i][j]+v[i+1][j]+v[i][j-1]+v[i+1][j-1])/4.0;
  adv=-u[i][j]*(u[i+1][j]-u[i-1][j])/dx/2.0-vav*(u[i][j+1]-u[i][j-1])/dy/2.0;
   vis=((u[i+1][j]-2*u[i][j]+u[i-1][j])/dx/dx+(u[i][j+1]-2*u[i][j]+u[i][j-1])/dy/dy)/Re;
   prs=-(p[i+1][j]-p[i][j])/dx;
   tu[i][j]=(adv+vis+prs)*dt+u[i][j];
  }
 for(i=Fx2+1;i<=Fx2+Kx;i++)
  for(j=1;j<=Fy2-1;j++)
  {
   uav=(u[i-1][j+1]+u[i][j+1]+u[i-1][j]+u[i][j])/4.0;
   adv=-uav*(v[i+1][j]-v[i-1][j])/dx/2.0-v[i][j]*(v[i][j+1]-v[i][j-1])/dy/2.0;
   vis=((v[i+1][j]-2*v[i][j]+v[i-1][j])/dx/dx+(v[i][j+1]-2*v[i][j]+v[i][j-1])/dy/dy)/Re;
   prs=-(p[i][j+1]-p[i][j])/dy;
   tv[i][j]=(adv+vis+prs)*dt+v[i][j];
  } 

 for(i=Fx2;i<=Fx2+Kx;i++)
  for(j=1;j<=Fy2;j++)
   u[i][j]=tu[i][j];
 for(i=Fx2+1;i<=Fx2+Kx;i++)
  for(j=1;j<=Fy2-1;j++)
   v[i][j]=tv[i][j];
  //第四部分中
  for(i=Fx2;i<=Fx2+Kx;i++)
  for(j=Fy2+Ky+1;j<=Fy3;j++)
  {
   vav=(v[i][j]+v[i+1][j]+v[i][j-1]+v[i+1][j-1])/4.0;
   adv=-u[i][j]*(u[i+1][j]-u[i-1][j])/dx/2.0-vav*(u[i][j+1]-u[i][j-1])/dy/2.0;
   vis=((u[i+1][j]-2*u[i][j]+u[i-1][j])/dx/dx+(u[i][j+1]-2*u[i][j]+u[i][j-1])/dy/dy)/Re;
   prs=-(p[i+1][j]-p[i][j])/dx;
   tu[i][j]=(adv+vis+prs)*dt+u[i][j];
  }
 for(i=Fx2+1;i<=Fx2+Kx;i++)
  for(j=Fy2+Ky+1;j<=Fy3-1;j++)
  {
   uav=(u[i-1][j+1]+u[i][j+1]+u[i-1][j]+u[i][j])/4.0;
   adv=-uav*(v[i+1][j]-v[i-1][j])/dx/2.0-v[i][j]*(v[i][j+1]-v[i][j-1])/dy/2.0;
   vis=((v[i+1][j]-2*v[i][j]+v[i-1][j])/dx/dx+(v[i][j+1]-2*v[i][j]+v[i][j-1])/dy/dy)/Re;
   prs=-(p[i][j+1]-p[i][j])/dy;
   tv[i][j]=(adv+vis+prs)*dt+v[i][j];
  } 

 for(i=Fx2;i<=Fx2+Kx;i++)
  for(j=Fy2+Ky+1;j<=Fy3;j++)
   u[i][j]=tu[i][j];
 for(i=Fx2+1;i<=Fx2+Kx;i++)
  for(j=Fy2+Ky+1;j<=Fy3-1;j++)
   v[i][j]=tv[i][j];
//第四部分上
  for(i=Fx2;i<=Fx2+Kx;i++)
  for(j=Fy3+Ky+1;j<=Jy;j++)
  {
   vav=(v[i][j]+v[i+1][j]+v[i][j-1]+v[i+1][j-1])/4.0;
   adv=-u[i][j]*(u[i+1][j]-u[i-1][j])/dx/2.0-vav*(u[i][j+1]-u[i][j-1])/dy/2.0;
   vis=((u[i+1][j]-2*u[i][j]+u[i-1][j])/dx/dx+(u[i][j+1]-2*u[i][j]+u[i][j-1])/dy/dy)/Re;
   prs=-(p[i+1][j]-p[i][j])/dx;
   tu[i][j]=(adv+vis+prs)*dt+u[i][j];
  }
 for(i=Fx2+1;i<=Fx2+Kx;i++)
  for(j=Fy3+Ky+1;j<=Jy-1;j++)
  {
   uav=(u[i-1][j+1]+u[i][j+1]+u[i-1][j]+u[i][j])/4.0;
   adv=-uav*(v[i+1][j]-v[i-1][j])/dx/2.0-v[i][j]*(v[i][j+1]-v[i][j-1])/dy/2.0;
   vis=((v[i+1][j]-2*v[i][j]+v[i-1][j])/dx/dx+(v[i][j+1]-2*v[i][j]+v[i][j-1])/dy/dy)/Re;
   prs=-(p[i][j+1]-p[i][j])/dy;
   tv[i][j]=(adv+vis+prs)*dt+v[i][j];
  } 

 for(i=Fx2;i<=Fx2+Kx;i++)
  for(j=Fy3+Ky+1;j<=Jy;j++)
   u[i][j]=tu[i][j];
 for(i=Fx2+1;i<=Fx2+Kx;i++)
  for(j=Fy3+Ky+1;j<=Jy-1;j++)
   v[i][j]=tv[i][j];
//第五部分
for(i=Fx2+Kx+1;i<=Jx-1;i++)
  for(j=1;j<=Jy;j++)
  {
   vav=(v[i][j]+v[i+1][j]+v[i][j-1]+v[i+1][j-1])/4.0;
   adv=-u[i][j]*(u[i+1][j]-u[i-1][j])/dx/2.0-vav*(u[i][j+1]-u[i][j-1])/dy/2.0;
   vis=((u[i+1][j]-2*u[i][j]+u[i-1][j])/dx/dx+(u[i][j+1]-2*u[i][j]+u[i][j-1])/dy/dy)/Re;
   prs=-(p[i+1][j]-p[i][j])/dx;
   tu[i][j]=(adv+vis+prs)*dt+u[i][j];
  }
 for(i=Fx2+Kx+1;i<=Jx;i++)
  for(j=1;j<=Jy-1;j++)
  {
   uav=(u[i-1][j+1]+u[i][j+1]+u[i-1][j]+u[i][j])/4.0;
   adv=-uav*(v[i+1][j]-v[i-1][j])/dx/2.0-v[i][j]*(v[i][j+1]-v[i][j-1])/dy/2.0;
   vis=((v[i+1][j]-2*v[i][j]+v[i-1][j])/dx/dx+(v[i][j+1]-2*v[i][j]+v[i][j-1])/dy/dy)/Re;
   prs=-(p[i][j+1]-p[i][j])/dy;
   tv[i][j]=(adv+vis+prs)*dt+v[i][j];
  } 

 for(i=Fx2+Kx+1;i<=Jx-1;i++)
  for(j=1;j<=Jy;j++)
   u[i][j]=tu[i][j];
 for(i=Fx2+Kx+1;i<=Jx;i++)
  for(j=1;j<=Jy-1;j++)
   v[i][j]=tv[i][j];
  
 bounduv(u,v);
}
/*-------------------------------------------------------
利用Chorin压力迭代解法求解压力p
--------------------------------------------------------*/
void solvep(double p[Jx+2][Jy+2],double u[Jx+2][Jy+2],double v[Jx+2][Jy+2],double dx,double dy,double lambda)
{
 int i,j;
 dx=4.0/Jx;    //x方向网格间距
 dy=2.0/Jy;    //y方向网格间距
 //第一部分
  for(i=1;i<=Fx1;i++)
  for(j=1;j<=Jy;j++)
   p[i][j]=p[i][j]-lambda*((u[i][j]-u[i-1][j])/dx+(v[i][j]-v[i][j-1])/dy);
//第二部分
  for(i=Fx1+1;i<=Fx1+Kx;i++)
  for(j=1;j<=Fy1;j++)
   p[i][j]=p[i][j]-lambda*((u[i][j]-u[i-1][j])/dx+(v[i][j]-v[i][j-1])/dy);

  for(i=Fx1+1;i<=Fx1+Kx;i++)
  for(j=Fy1+Ky+1;j<=Jy;j++)
   p[i][j]=p[i][j]-lambda*((u[i][j]-u[i-1][j])/dx+(v[i][j]-v[i][j-1])/dy);
//第三部分
  for(i=Fx1+Kx+1;i<=Fx2;i++)
  for(j=1;j<=Jy;j++)
   p[i][j]=p[i][j]-lambda*((u[i][j]-u[i-1][j])/dx+(v[i][j]-v[i][j-1])/dy);
//第四部分
  for(i=Fx2+1;i<=Fx2+Kx;i++)
  for(j=1;j<=Fy2;j++)
   p[i][j]=p[i][j]-lambda*((u[i][j]-u[i-1][j])/dx+(v[i][j]-v[i][j-1])/dy);

  for(i=Fx2+1;i<=Fx2+Kx;i++)
  for(j=Fy2+Ky+1;j<=Fy3;j++)
   p[i][j]=p[i][j]-lambda*((u[i][j]-u[i-1][j])/dx+(v[i][j]-v[i][j-1])/dy);

  for(i=Fx2+1;i<=Fx2+Kx;i++)
  for(j=Fy3+Ky+1;j<=Jy;j++)
   p[i][j]=p[i][j]-lambda*((u[i][j]-u[i-1][j])/dx+(v[i][j]-v[i][j-1])/dy);
//第五部分
  for(i=Fx2+Kx+1;i<=Jx;i++)
  for(j=1;j<=Jy;j++)
   p[i][j]=p[i][j]-lambda*((u[i][j]-u[i-1][j])/dx+(v[i][j]-v[i][j-1])/dy);

  boundp(p);
 boundjd(p);
}
/*-------------------------------------------------------
判断收敛:连续方程与0的误差小于误差限 
--------------------------------------------------------*/
int conv(double u[Jx+2][Jy+2],double v[Jx+2][Jy+2],double dx, double dy, double &dm)
{
 int i,j;
 double er;
 dx=4.0/Jx;    //x方向网格间距
 dy=2.0/Jy;    //y方向网格间距
 dm=0;
  //第一部分
  for(i=1;i<=Fx1;i++)
  for(j=1;j<=Jy;j++)
  { 
   er=fabs((u[i][j]-u[i-1][j])/dx+(v[i][j]-v[i][j-1])/dy);
   if(dm<er)dm=er;
  } 
  //第二部分
  for(i=Fx1+1;i<=Fx1+Kx;i++)
  for(j=1;j<=Fy1;j++)
  { 
   er=fabs((u[i][j]-u[i-1][j])/dx+(v[i][j]-v[i][j-1])/dy);
   if(dm<er)dm=er;
  } 

  for(i=Fx1+1;i<=Fx1+Kx;i++)
  for(j=Fy1+Ky+1;j<=Jy;j++)
  { 
   er=fabs((u[i][j]-u[i-1][j])/dx+(v[i][j]-v[i][j-1])/dy);
   if(dm<er)dm=er;
  } 
  //第三部分
  for(i=Fx1+Kx+1;i<=Fx2;i++)
  for(j=1;j<=Jy;j++)
  { 
   er=fabs((u[i][j]-u[i-1][j])/dx+(v[i][j]-v[i][j-1])/dy);
   if(dm<er)dm=er;
  } 
//第四部分
  for(i=Fx2+1;i<=Fx2+Kx;i++)
  for(j=1;j<=Fy2;j++)
  { 
   er=fabs((u[i][j]-u[i-1][j])/dx+(v[i][j]-v[i][j-1])/dy);
   if(dm<er)dm=er;
  } 

  for(i=Fx2+1;i<=Fx2+Kx;i++)
  for(j=Fy2+Ky+1;j<=Fy3;j++)
  { 
   er=fabs((u[i][j]-u[i-1][j])/dx+(v[i][j]-v[i][j-1])/dy);
   if(dm<er)dm=er;
  } 

  for(i=Fx2+1;i<=Fx2+Kx;i++)
  for(j=Fy3+Ky+1;j<=Jy;j++)
  { 
   er=fabs((u[i][j]-u[i-1][j])/dx+(v[i][j]-v[i][j-1])/dy);
   if(dm<er)dm=er;
  } 
//第五部分
  for(i=Fx2+Kx+1;i<=Jx;i++)
  for(j=1;j<=Jy;j++)
  { 
   er=fabs((u[i][j]-u[i-1][j])/dx+(v[i][j]-v[i][j-1])/dy);
   if(dm<er)dm=er;
  } 
  
 if(dm<EPS)return 1;
 else return 0;
}
/*-------------------------------------------------------
输出全场计算结果
-------------------------------------------------------*/
void output(const char *sN, double u[Jx+2][Jy+2],double v[Jx+2][Jy+2],double p[Jx+2][Jy+2],double dx, double dy)
{
 int i,j;
 dx=4.0/Jx;    //x方向网格间距
 dy=2.0/Jy;    //y方向网格间距
 FILE *fp;
 
 fp=fopen(sN,"w");
 fprintf(fp,"Title = \"computational Results\" \n");
 fprintf(fp,"VARIABLES = \"x\", \"y\", \"u\",\"v\",\"p\" \n");
 //1 
 fprintf(fp,"ZONE T= \"Zone 1\", I=%d, J=%d,F=POINT\n",Jy,Fx1); 
 for(i=1;i<=Fx1;i++)
  for(j=1;j<=Jy;j++)
   fprintf(fp,"%16f%16f%20e%20e%20e\n",i*dx-0.5*dx,j*dy-0.5*dy,(u[i][j]+u[i-1][j])/2,(v[i][j-1]+v[i][j])/2.0,p[i][j]);  
//2
 fprintf(fp,"ZONE T= \"Zone 2(1)\", I=%d, J=%d,F=POINT\n",Fy1,Kx+1); 
  for(i=Fx1;i<=Fx1+Kx;i++)
  for(j=1;j<=Fy1;j++)
   fprintf(fp,"%16f%16f%20e%20e%20e\n",i*dx-0.5*dx,j*dy-0.5*dy,(u[i][j]+u[i-1][j])/2,(v[i][j-1]+v[i][j])/2.0,p[i][j]); 

  fprintf(fp,"ZONE T= \"Zone 2(2)\", I=%d, J=%d,F=POINT\n",Jy-Fy1-Ky,Kx+1); 
  for(i=Fx1;i<=Fx1+Kx;i++)
  for(j=Fy1+Ky+1;j<=Jy;j++)
   fprintf(fp,"%16f%16f%20e%20e%20e\n",i*dx-0.5*dx,j*dy-0.5*dy,(u[i][j]+u[i-1][j])/2,(v[i][j-1]+v[i][j])/2.0,p[i][j]);  
 //3
fprintf(fp,"ZONE T= \"Zone 3\", I=%d, J=%d,F=POINT\n",Jy,Fx2-Fx1-Kx+1); 
  for(i=Fx1+Kx;i<=Fx2;i++)
  for(j=1;j<=Jy;j++)
   fprintf(fp,"%16f%16f%20e%20e%20e\n",i*dx-0.5*dx,j*dy-0.5*dy,(u[i][j]+u[i-1][j])/2,(v[i][j-1]+v[i][j])/2.0,p[i][j]);  
//4
fprintf(fp,"ZONE T= \"Zone 4(1)\", I=%d, J=%d,F=POINT\n",Fy2,Kx+1); 
  for(i=Fx2;i<=Fx2+Kx;i++)
  for(j=1;j<=Fy2;j++) 
    fprintf(fp,"%16f%16f%20e%20e%20e\n",i*dx-0.5*dx,j*dy-0.5*dy,(u[i][j]+u[i-1][j])/2,(v[i][j-1]+v[i][j])/2.0,p[i][j]); 

  fprintf(fp,"ZONE T= \"Zone 4(2)\", I=%d, J=%d,F=POINT\n",Fy3-Fy2-Ky,Kx+1); 
  for(i=Fx2;i<=Fx2+Kx;i++)
  for(j=Fy2+Ky+1;j<=Fy3;j++)
    fprintf(fp,"%16f%16f%20e%20e%20e\n",i*dx-0.5*dx,j*dy-0.5*dy,(u[i][j]+u[i-1][j])/2,(v[i][j-1]+v[i][j])/2.0,p[i][j]); 
  
  fprintf(fp,"ZONE T= \"Zone 4(3)\", I=%d, J=%d,F=POINT\n",Jy-Fy3-Ky,Kx+1); 
  for(i=Fx2;i<=Fx2+Kx;i++)
  for(j=Fy3+Ky+1;j<=Jy;j++)
    fprintf(fp,"%16f%16f%20e%20e%20e\n",i*dx-0.5*dx,j*dy-0.5*dy,(u[i][j]+u[i-1][j])/2,(v[i][j-1]+v[i][j])/2.0,p[i][j]); 
//5
fprintf(fp,"ZONE T= \"Zone 5\", I=%d, J=%d,F=POINT\n",Jy,Jx-Fx2-Kx+1); 
  for(i=Fx2+Kx;i<=Jx;i++)
  for(j=1;j<=Jy;j++)
    fprintf(fp,"%16f%16f%20e%20e%20e\n",i*dx-0.5*dx,j*dy-0.5*dy,(u[i][j]+u[i-1][j])/2,(v[i][j-1]+v[i][j])/2.0,p[i][j]); 

fclose(fp);
}
/*-------------------------------------------------------
求解流场的主函数
-------------------------------------------------------*/
int main()
{
 int n,isconv;
 double lambda;
 double dm;
 double dx;
 double dy;
init(u,v,p,dx,dy);
 /*-------------------------------------------------------------------------
 lambda常数与时间步长、网格间距和Reynolds数都有关。
 关系到收敛性，应仔细选择！
 ----------------------------------------------------------------------------*/
 lambda=0.03;
 n=0;
do 
 {
  solvep(p,u,v,dx,dy,lambda);
  solveuv(u,v,p,dx,dy,tu,tv);
  isconv=conv(u,v,dx,dy,dm);
  n++;
  if(n%10==0)
  {
   printf("%5d Step Max vel divergence=%16e\n",n,dm);
  }
 }while(!isconv);
  string output_str = string("result.plt");
 output(output_str.c_str(),u,v,p,dx,dy);
}
