#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath> 
#include <string.h>

using namespace std; 

double const GAMA=1.4;

#define MIN(x,y) (((x)<(y))?(x):(y))
#define MAX(x,y) (((x)>(y))?(x):(y))

#define Lx 4.0
#define Ly 1.0

#define l_x 1.0 //前台阶位置
#define l_y 0.2

#define Jx 200
#define Jy 60

#define TT 15 //总时间
#define Sf 0.1 //时间步长因子
#define Ebusila 0.05 //人工粘性项小参数


//全局变量
int N_x,N_y;

double U[Jx+5][Jy+5][4],F[Jx+5][Jy+5][4],G[Jx+5][Jy+5][4];

/*---------------------------------------------------
计算时间步长
入口：U，当前物理量
	  dx，dy，网格宽度
返回：时间步长
----------------------------------------------------*/
double CFL(double U[Jx+5][Jy+5][4],double dx,double dy)
{
	int i,j;
	double u1,u2,v1,v2,p1,p2,h1,h2,vel_x,vel_y;
	double Mx,My,u_temp,v_temp,h_temp,a_temp,Ratio_rou;
	Mx=1e-100;
	My=1e-100;
	for(i=2;i<=Jx+2;i++)
		for(j=2;j<=Jy+2;j++)
		{	
			if((i>N_x)&&(j<N_y))continue;
	
	//x方向
			u1=U[i][j][1]/U[i][j][0];
			v1=U[i][j][2]/U[i][j][0];
			p1=(GAMA-1)*(U[i][j][3]-0.5*U[i][j][0]*(u1*u1+v1*v1));
			h1=GAMA*p1/((GAMA-1)*U[i][j][0])+0.5*(u1*u1+v1*v1);

			u2=U[i+1][j][1]/U[i+1][j][0];
			v2=U[i+1][j][2]/U[i+1][j][0];
			p2=(GAMA-1)*(U[i+1][j][3]-0.5*U[i+1][j][0]*(u2*u2+v2*v2));
			h2=GAMA*p2/((GAMA-1)*U[i+1][j][0])+0.5*(u2*u2+v2*v2);

			Ratio_rou=sqrt(U[i+1][j][0]/U[i][j][0]);
			u_temp=(u1+u2*Ratio_rou)/(1.0+Ratio_rou);
			v_temp=(v1+v2*Ratio_rou)/(1.0+Ratio_rou);
			h_temp=(h1+h2*Ratio_rou)/(1.0+Ratio_rou);
			a_temp=sqrt(fabs((GAMA-1)*(h_temp-0.5*(u_temp*u_temp+v_temp*v_temp))));
			vel_x=a_temp+fabs(u_temp);
			if(vel_x>Mx)Mx=vel_x;			
	
	//y方向
			u1=U[i][j][1]/U[i][j][0];
			v1=U[i][j][2]/U[i][j][0];
			p1=(GAMA-1)*(U[i][j][3]-0.5*U[i][j][0]*(u1*u1+v1*v1));
			h1=GAMA*p1/((GAMA-1)*U[i][j][0])+0.5*(u1*u1+v1*v1);

			u2=U[i][j+1][1]/U[i][j+1][0];
			v2=U[i][j+1][2]/U[i][j+1][0];
			p2=(GAMA-1)*(U[i][j+1][3]-0.5*U[i][j+1][0]*(u2*u2+v2*v2));
			h2=GAMA*p2/((GAMA-1)*U[i][j+1][0])+0.5*(u2*u2+v2*v2);

			Ratio_rou=sqrt(U[i][j+1][0]/U[i][j][0]);
			u_temp=(u1+u2*Ratio_rou)/(1.0+Ratio_rou);
			v_temp=(v1+v2*Ratio_rou)/(1.0+Ratio_rou);
			h_temp=(h1+h2*Ratio_rou)/(1.0+Ratio_rou);
			a_temp=sqrt(fabs((GAMA-1)*(h_temp-0.5*(u_temp*u_temp+v_temp*v_temp))));
			vel_y=a_temp+fabs(v_temp);
			if(vel_y>My)My=vel_y;			
		}

	return Sf*MIN(dx,dy)/MAX(Mx,My);
}
/*--------------------------------------------------------
计算数值格式粘性函数Q_z
入口：lamda（特征值）
出口：返回Q_z值
--------------------------------------------------------*/
double Q_z(double x)
{
	double result;
	if(fabs(x)>=Ebusila)
		result=fabs(x);
	else
		result=(x*x/Ebusila+Ebusila)/2.0;
	return result;
}
/*-------------------------------------------------------
计算人工粘性修正项的minmod
入口：alfa[i],alfa[i-1],alfa[i+1]等
出口：返回minmod的值
--------------------------------------------------------*/
double minmod(double x,double y,double z)
{
	double result,temp;

	if((x<0&&y<0&&z<0)||(x>0.0&&y>0&&z>0))		
	{
		if(fabs(x)>fabs(y))
		{
			temp=x;
			x=y;
			y=temp;
		}
		if(fabs(x)>fabs(z))
		{
			temp=x;
			x=z;
			z=temp;
		}
		result=x;
	}
	else result=0.0;
	return result;	
}
/*-------------------------------------------------------
根据U计算F
入口：U，当前U矢量
出口：F计算得到F矢量
U,F定义参看Euler方程
---------------------------------------------------------*/
void U2F(double U[4],double F[4])
{
	double u,v,p;
	u=U[1]/U[0];
	v=U[2]/U[0];
	p=(GAMA-1)*(U[3]-0.5*U[0]*(u*u+v*v));
	F[0]=U[1];
	F[1]=U[0]*u*u+p;
	F[2]=U[0]*u*v;
	F[3]=(U[3]+p)*u;
}
/*-------------------------------------------------------
根据U计算G
入口：U，当前U矢量
出口：G计算得到G矢量
U,G定义参看Euler方程
---------------------------------------------------------*/
void U2G(double U[4],double G[4])
{
	double u,v,p;
	u=U[1]/U[0];
	v=U[2]/U[0];
	p=(GAMA-1)*(U[3]-0.5*U[0]*(u*u+v*v));
	G[0]=U[2];
	G[1]=U[0]*u*v;
	G[2]=U[0]*v*v+p;
	G[3]=(U[3]+p)*v;
}
/*-------------------------------------------------
初始化
入口：无
出口：U，已经给定初始值
	  dx，dy网格宽度
---------------------------------------------------*/
void Init(double U[Jx+5][Jy+5][4],double& dx,double& dy)
{
	int i,j;
	double rou1=1.0,u1=3.0,v1=0,p1=0.71429;

	dx=Lx/(Jx);
	dy=Ly/(Jy);
	N_x=(int)(l_x/dx)+2;
	N_y=(int)(l_y/dy)+2;

	for(i=0;i<Jx+5;i++)
		for(j=0;j<=Jy+4;j++)
		{
			U[i][j][0]=rou1;
			U[i][j][1]=rou1*u1;
			U[i][j][2]=rou1*v1;
			U[i][j][3]=p1/(GAMA-1)+rou1*(u1*u1+v1*v1)/2;
		}
}
/*------------------------------------------------------
边界条件
入口：dx，dy，网格宽度
出口：U，已经给定边界
--------------------------------------------------------*/

void bound(double U[Jx+5][Jy+5][4],double dx,double dy)
{
	int i,j,k;
	double rou1=1.0,u1=3.0,v1=0,p1=0.71429;

//左边界	
	for(j=0;j<=Jy+4;j++)
	{
		U[0][j][0]=rou1;
		U[0][j][1]=rou1*u1;
		U[0][j][2]=rou1*v1;
		U[0][j][3]=p1/(GAMA-1)+0.5*rou1*(u1*u1+v1*v1);

		U[1][j][0]=rou1;
		U[1][j][1]=rou1*u1;
		U[1][j][2]=rou1*v1;
		U[1][j][3]=p1/(GAMA-1)+0.5*rou1*(u1*u1+v1*v1);
	}

//上边界
	for(i=0;i<=Jx+4;i++)
	{
		U[i][Jy+3][0]=U[i][Jy+1][0];
		U[i][Jy+3][1]=U[i][Jy+1][1];
		U[i][Jy+3][2]=-U[i][Jy+1][2];
		U[i][Jy+3][3]=U[i][Jy+1][3];

		U[i][Jy+4][0]=U[i][Jy][0];
		U[i][Jy+4][1]=U[i][Jy][1];
		U[i][Jy+4][2]=-U[i][Jy][2];
		U[i][Jy+4][3]=U[i][Jy][3];
	}

//右边界
	for(j=N_y;j<=Jy+4;j++)
		for(k=0;k<4;k++)
		{
			U[Jx+3][j][k]=U[Jx+2][j][k];
			U[Jx+4][j][k]=U[Jx+3][j][k];                                                                                     
		}	

//左下边界
	for(i=0;i<=N_x;i++)
	{
		U[i][0][0]=U[i][4][0];
		U[i][0][1]=U[i][4][1];
		U[i][0][2]=-U[i][4][2];
		U[i][0][3]=U[i][4][3];

		U[i][1][0]=U[i][3][0];
		U[i][1][1]=U[i][3][1];
		U[i][1][2]=-U[i][3][2];
		U[i][1][3]=U[i][3][3];

	
	}

//右下边界
	for(i=N_x+3;i<=Jx+2;i++)
	{
		U[i][N_y-2][0]=U[i][N_y+2][0];
		U[i][N_y-2][1]=U[i][N_y+2][1];
		U[i][N_y-2][2]=-U[i][N_y+2][2];
		U[i][N_y-2][3]=U[i][N_y+2][3];

		U[i][N_y-1][0]=U[i][N_y+1][0];
		U[i][N_y-1][1]=U[i][N_y+1][1];
		U[i][N_y-1][2]=-U[i][N_y+1][2];
		U[i][N_y-1][3]=U[i][N_y+1][3];
	}

//台阶y方向的边界处理
	for(j=2;j<N_y-2;j++)
	{		
		U[N_x+2][j][0]=U[N_x-2][j][0];
		U[N_x+2][j][1]=-U[N_x-2][j][1];
		U[N_x+2][j][2]=U[N_x-2][j][2];
		U[N_x+2][j][3]=U[N_x-2][j][3];	

		U[N_x+1][j][0]=U[N_x-1][j][0];
		U[N_x+1][j][1]=-U[N_x-1][j][1];
		U[N_x+1][j][2]=U[N_x-1][j][2];
		U[N_x+1][j][3]=U[N_x-1][j][3];
	}

//下角点处理
	
	for(k=0;k<=3;k++)
	{
		U[N_x+1][0][k]=U[N_x-1][4][k];
		U[N_x+1][1][k]=U[N_x-1][3][k];
		U[N_x+2][1][k]=U[N_x-2][3][k];
		U[N_x+2][0][k]=U[N_x-2][4][k];	
	}

//上角点处理
		U[N_x+2][N_y-1][0]=U[N_x+2][N_y+1][0];
		U[N_x+2][N_y-1][1]=U[N_x+2][N_y+1][1];
		U[N_x+2][N_y-1][2]=-U[N_x+2][N_y+1][2];
		U[N_x+2][N_y-1][3]=U[N_x+2][N_y+1][3];

		U[N_x+1][N_y-2][0]=U[N_x-1][N_y-2][0];
		U[N_x+1][N_y-2][1]=-U[N_x-1][N_y-2][1];
		U[N_x+1][N_y-2][2]=U[N_x-1][N_y-2][2];
		U[N_x+1][N_y-2][3]=U[N_x-1][N_y-2][3];

		U[N_x+2][N_y-2][0]=U[N_x][N_y][0];
		U[N_x+2][N_y-2][1]=-U[N_x][N_y][1];
		U[N_x+2][N_y-2][2]=-U[N_x][N_y][2];
		U[N_x+2][N_y-2][3]=U[N_x][N_y][3];

		U[N_x+1][N_y-1][0]=0.5*(U[N_x-1][N_y-1][0]+U[N_x+1][N_y+1][0]);
		U[N_x+1][N_y-1][1]=0.5*(U[N_x-1][N_y-1][1]+U[N_x+1][N_y+1][1]);
		U[N_x+1][N_y-1][2]=0.5*(U[N_x-1][N_y-1][2]+U[N_x+1][N_y+1][2]);
		U[N_x+1][N_y-1][3]=0.5*(U[N_x-1][N_y-1][3]+U[N_x+1][N_y+1][3]);

		j=Jy+2;
		for(i=2;i<=Jx+2;i++)
	{
		U[i][j][2]=0;
	}//上边界壁面法向速度强制赋0

		j=2;
		for(i=2;i<=N_x;i++)
	{
		U[i][j][2]=0;
	}//下边界壁面法向速度强制赋0
		
		j=N_y;
		for(i=N_x;i<=Jx+2;i++)
	{
		U[i][j][2]=0;
	}//下边界壁面法向速度强制赋0

		i=N_x;
		for(j=2;j<=N_y;j++)
	{
		U[i][j][1]=0;
	}//台阶垂直壁面法向速度强制赋0

}
/*-------------------------------------------------------
计算x方向的特征值lamda_x
入口：U矢量
出口：lamda_x矢量
--------------------------------------------------------*/
void lamda_x(double lamda_x[4],double u,double v,double a)
{
	lamda_x[0]=u-a;
	lamda_x[1]=u;
	lamda_x[2]=u+a;
	lamda_x[3]=u;
}
/*-------------------------------------------------------
计算x方向的左右特征矢量矩阵L_x，R_x
入口：U矢量,L,R
出口：L_x、R_x矢量
--------------------------------------------------------*/
void LR_x(double R[4][4],double L[4][4],double u,double v,double a )
{
	double b1,b2;
	//计算x右特征矢量R
			R[0][0]=1.0;
			R[0][1]=1.0;
			R[0][2]=1.0;
			R[0][3]=0.0;

			R[1][0]=u-a;
			R[1][1]=u;
			R[1][2]=u+a;
			R[1][3]=0;

			R[2][0]=v;
			R[2][1]=v;
			R[2][2]=v;
			R[2][3]=1.0;

			R[3][0]=a*a/(GAMA-1.0)+0.5*(u*u+v*v)-u*a;
			R[3][1]=0.5*(u*u+v*v);
			R[3][2]=a*a/(GAMA-1.0)+0.5*(u*u+v*v)+u*a;
			R[3][3]=v;

	//计算x左特征矢量L
			b2=(GAMA-1)/(a*a);
			b1=0.5*b2*(u*u+v*v);

			L[0][0]=0.5*(b1+u/a);
			L[0][1]=0.5*(-b2*u-1.0/a);
			L[0][2]=0.5*(-b2*v);
			L[0][3]=0.5*b2;

			L[1][0]=1.0-b1;
			L[1][1]=b2*u;
			L[1][2]=b2*v;
			L[1][3]=-b2;

			L[2][0]=0.5*(b1-u/a);
			L[2][1]=0.5*(-b2*u+1.0/a);
			L[2][2]=0.5*(-b2*v);
			L[2][3]=0.5*b2;

			L[3][0]=-v;
			L[3][1]=0;
			L[3][2]=1;
			L[3][3]=0;
}
/*-------------------------------------------------------
计算y方向的特征值lamda_y
--------------------------------------------------------*/
void lamda_y(double lamda_y[4],double u,double v,double a)
{
	lamda_y[0]=v-a;
	lamda_y[1]=v;
	lamda_y[2]=v+a;
	lamda_y[3]=v;
}
/*-------------------------------------------------------
计算y方向的左右特征矢量矩阵L_y，R_y
入口：U矢量,L,R,(u[i][j]+u[i][j+1])/2,(v[i][j]+v[i][j+1])/2,(a[i][j]+a[i][j+1])/2
出口：L_y、R_y矢量
--------------------------------------------------------*/
void LR_y(double R[4][4],double L[4][4],double u,double v,double a)
{
	double b1,b2;
	//计算y右特征矢量R
			R[0][0]=1.0;
			R[0][1]=1.0;
			R[0][2]=1.0;
			R[0][3]=0.0;

			R[1][0]=u;
			R[1][1]=u;
			R[1][2]=u;
			R[1][3]=1;

			R[2][0]=v-a;
			R[2][1]=v;
			R[2][2]=v+a;
			R[2][3]=0.0;

			R[3][0]=a*a/(GAMA-1.0)+0.5*(u*u+v*v)-v*a;
			R[3][1]=0.5*(v*v+u*u);
			R[3][2]=a*a/(GAMA-1.0)+0.5*(u*u+v*v)+v*a;
			R[3][3]=u;

	//计算y左特征矢量L[i][j]
			b2=(GAMA-1)/(a*a);
			b1=0.5*b2*(u*u+v*v);

			L[0][0]=0.5*(b1+v/a);
			L[0][1]=0.5*(-b2*u);
			L[0][2]=0.5*(-b2*v-1/a);
			L[0][3]=0.5*b2;

			L[1][0]=1-b1;
			L[1][1]=b2*u;
			L[1][2]=b2*v;
			L[1][3]=-b2;

			L[2][0]=0.5*(b1-v/a);
			L[2][1]=0.5*(-b2*u);
			L[2][2]=0.5*(-b2*v+1/a);
			L[2][3]=0.5*b2;

			L[3][0]=-u;
			L[3][1]=1.0;
			L[3][2]=0.0;
			L[3][3]=0.0;
}

/*-------------------------------------------------------
对称TVD格式，计算x方向的差分算子
入口：U矢量,dx,dt
出口：U_x矢量
--------------------------------------------------------*/
void LLx(double U[Jx+5][Jy+5][4],double dx,double dt)
{
	int i,j,k,l;
	double u1,v1,p1,h1,u2,v2,p2,h2;
	double u_temp,v_temp,a_temp,h_temp,r;
	double Ratio_rou;
	double lamda[Jx+5][4],R[Jx+5][4][4],L[Jx+5][4][4];
	double Fai[Jx+5][4],alfa[Jx+5][4],g_b[Jx+5][4];

	r=dt/dx;

	for(i=0;i<=Jx+4;i++)
		for(j=0;j<=Jy+4;j++)
		{
				U2F(U[i][j],F[i][j]);
		}
			
	for(j=0;j<=Jy+4;j++)
	{
		for(i=0;i<=Jx+3;i++)
		{
			u1=U[i][j][1]/U[i][j][0];
			v1=U[i][j][2]/U[i][j][0];
			p1=(GAMA-1)*(U[i][j][3]-0.5*U[i][j][0]*(u1*u1+v1*v1));
			h1=GAMA*p1/((GAMA-1)*U[i][j][0])+0.5*(u1*u1+v1*v1);

			u2=U[i+1][j][1]/U[i+1][j][0];
			v2=U[i+1][j][2]/U[i+1][j][0];
			p2=(GAMA-1)*(U[i+1][j][3]-0.5*U[i+1][j][0]*(u2*u2+v2*v2));
			h2=GAMA*p2/((GAMA-1)*U[i+1][j][0])+0.5*(u2*u2+v2*v2);

			Ratio_rou=sqrt(U[i+1][j][0]/U[i][j][0]);
			u_temp=(u1+u2*Ratio_rou)/(1.0+Ratio_rou);
			v_temp=(v1+v2*Ratio_rou)/(1.0+Ratio_rou);
			h_temp=(h1+h2*Ratio_rou)/(1.0+Ratio_rou);			
			a_temp=sqrt(fabs((GAMA-1)*(h_temp-0.5*(u_temp*u_temp+v_temp*v_temp))));
			lamda_x(lamda[i],u_temp,v_temp,a_temp);
			LR_x(R[i],L[i],u_temp,v_temp,a_temp);			
		}
	
//初始化Fai[Jx+5][4],alfa[Jx+5][4],g_b[Jx+5][4];
		for(i=0;i<=Jx+4;i++)
			for(k=0;k<=3;k++)
			{
				alfa[i][k]=0.0;
				Fai[i][k]=0.0;
				g_b[i][k]=0.0;
			}

//计算alfa[i][k]
		for(i=0;i<=Jx+3;i++)			
			for(k=0;k<=3;k++)
				for(l=0;l<=3;l++)	
					alfa[i][k]=alfa[i][k]+L[i][k][l]*(U[i+1][j][l]-U[i][j][l]);
		
//计算g_b[Jx+5][Jy+5][4]
		for(i=1;i<=Jx+3;i++)		
			for(k=0;k<=3;k++)
				g_b[i][k]=minmod(alfa[i-1][k],alfa[i][k],alfa[i+1][k]);		

//计算Fai[Jx+5][4]
		for(i=0;i<=Jx+4;i++)		
			for(k=0;k<=3;k++)
				for(l=0;l<=3;l++)
					Fai[i][k]=Fai[i][k]-R[i][k][l]*(pow(r*lamda[i][l],2)*g_b[i][l]+Q_z(r*lamda[i][l])*(alfa[i][l]-g_b[i][l]))/r;
//计算U_x
		double fl,fr;
		for(i=2;i<=Jx+2;i++)
			for(k=0;k<=3;k++)
			{
				if(i>N_x+1&&j<N_y-1)
					continue;
				fl=0.5*(F[i][j][k]+F[i-1][j][k]+Fai[i-1][k]);
				fr=0.5*(F[i][j][k]+F[i+1][j][k]+Fai[i][k]);
				U[i][j][k]=U[i][j][k]-r*(fr-fl);
			}		
	}

}
/*-------------------------------------------------------
对称TVD格式，计算y方向的差分算子
入口：U矢量,dy,dt
出口：U_y矢量
--------------------------------------------------------*/
void LLy(double U[Jx+5][Jy+5][4],double dy,double dt)
{
	int i,j,k,l;
	double u1,v1,p1,h1,u2,v2,p2,h2;
	double Ratio_rou,u_temp,v_temp,a_temp,h_temp,r;
	double lamda[Jy+5][4],R[Jy+5][4][4],L[Jy+5][4][4];
	double Fai[Jy+5][4],alfa[Jy+5][4],g_b[Jy+5][4];

	r=dt/dy;
	for(i=0;i<=Jx+4;i++)
		for(j=0;j<=Jy+4;j++)		
			U2G(U[i][j],G[i][j]);
			
	for(i=0;i<=Jx+4;i++)
	{
		for(j=0;j<=Jy+3;j++)
		{						
			u1=U[i][j][1]/U[i][j][0];
			v1=U[i][j][2]/U[i][j][0];
			p1=(GAMA-1)*(U[i][j][3]-0.5*U[i][j][0]*(u1*u1+v1*v1));
			h1=GAMA*p1/((GAMA-1)*U[i][j][0])+0.5*(u1*u1+v1*v1);

			u2=U[i][j+1][1]/U[i][j+1][0];
			v2=U[i][j+1][2]/U[i][j+1][0];
			p2=(GAMA-1)*(U[i][j+1][3]-0.5*U[i][j+1][0]*(u2*u2+v2*v2));
			h2=GAMA*p2/((GAMA-1)*U[i][j+1][0])+0.5*(u2*u2+v2*v2);

			Ratio_rou=sqrt(U[i][j+1][0]/U[i][j][0]);
			u_temp=(u1+u2*Ratio_rou)/(1.0+Ratio_rou);
			v_temp=(v1+v2*Ratio_rou)/(1.0+Ratio_rou);
			h_temp=(h1+h2*Ratio_rou)/(1.0+Ratio_rou);
			a_temp=sqrt(fabs((GAMA-1)*(h_temp-0.5*(u_temp*u_temp+v_temp*v_temp))));
			lamda_y(lamda[j],u_temp,v_temp,a_temp);
			LR_y(R[j],L[j],u_temp,v_temp,a_temp);				
		}
	
//初始化Fai[Jy+5][4],alfa[Jy+5][4],g_b[Jy+5][4];
		for(j=0;j<=Jy+4;j++)
			for(k=0;k<=3;k++)
			{
				alfa[j][k]=0.0;
				Fai[j][k]=0.0;
				g_b[j][k]=0.0;
			}
//计算alfa[j][k]
		for(j=0;j<=Jy+3;j++)
			for(k=0;k<=3;k++)
				for(l=0;l<=3;l++)
					alfa[j][k]=alfa[j][k]+L[j][k][l]*(U[i][j+1][l]-U[i][j][l]);	
//计算g_b[Jy+5][4]
		for(j=1;j<=Jy+3;j++)			
			for(k=0;k<=3;k++)
				g_b[j][k]=minmod(alfa[j-1][k],alfa[j][k],alfa[j+1][k]);
		
//计算Fai[Jy+5][4]
		for(j=0;j<=Jy+4;j++)
			for(k=0;k<=3;k++)
				for(l=0;l<=3;l++)
					Fai[j][k]=Fai[j][k]-R[j][k][l]*(pow(r*lamda[j][l],2)*g_b[j][l]+Q_z(r*lamda[j][l])*(alfa[j][l]-g_b[j][l]))/r;

//计算U_y
		double fl,fr;
		for(j=2;j<=Jy+2;j++)
			for(k=0;k<=3;k++)
			{
				if(i>N_x+1&&j<N_y-1)
					continue;
				fl=0.5*(G[i][j][k]+G[i][j-1][k]+Fai[j-1][k]);
				fr=0.5*(G[i][j][k]+G[i][j+1][k]+Fai[j][k]);
				U[i][j][k]=U[i][j][k]-r*(fr-fl);			
			}
	}
}
/*-------------------------------------------------------
*对称型TVD格式
--------------------------------------------------------*/
void symmetry_TVD(double U[Jx+5][Jy+5][4],double dx,double dy,double dt)
{
	bound(U,dx,dy);	
	LLx(U,dx,dt/2.0);
	LLy(U,dy,dt/2.0);

	bound(U,dx,dy);
	LLy(U,dy,dt/2.0);
	LLx(U,dx,dt/2.0);

}

/*----------------------------------------------------------------
保存前一时间步的U矢量
入口：无
出口：无
----------------------------------------------------------------*/
void save_U_last(double U_old[Jx+5][Jy+5][4],double U[Jx+5][Jy+5][4])
{
	int i,j,k;
	for(i=0;i<=Jx+4;i++)
		for(j=0;j<=Jy+4;j++)
			for(k=0;k<=3;k++)
				U_old[i][j][k]=U[i][j][k];
}

/*-----------------------------------------------------
求解误差
-------------------------------------------------------*/
double error(double U1[Jx+5][Jy+5][4],double U2[Jx+5][Jy+5][4])
{
	int i,j,k;
	double dd_variable;
	double err=1e-10;
	for(i=2;i<=Jx+2;i++)
		for(j=2;j<=Jy+2;j++)
			for(k=0;k<=3;k++)
		{	
			if((i>N_x)&&(j<N_y))
				continue;
			if(U2[i][j][k]!=0)
				dd_variable=fabs((U2[i][j][k]-U1[i][j][k]));
			else
				dd_variable=0;
			if(dd_variable>err)err=dd_variable;
		}

	return err;
}
/*---------------------------------------------------------------
输出telplot文件格式数据
入口：U,当前时刻U矢量
	  dx，x方向上的网格宽度
	  dy，y方向上的网格宽度
出口：无
----------------------------------------------------------------*/
void Output(double U[Jx+5][Jy+5][4],double dx,double dy,double T)
{
	int i,j;
	double rou,u,v,p,a,Ma;

	FILE *fp1,*fp2;
	char head1[200]={"TITLE	=\"Dataset\"\nVARIABLES =\"X\"\n\"y\"\n\"rou\"\n\"u\"\n\"v\"\n\"p\"\n\"Ma\"\n"};
	char head2[60]={"ZONE T=\"Rectangular zone1\"\n"};
	char head3[200]={"ZONETYPE=Ordered\nDATAPACKING=POINT\nDT=(SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE )"};
	char filename1[100],filename_tp[256]={"result_T_"},filename2[10]={".plt"},filename_tx[256]={"result_T_"},filename3[10]={".txt"};
	
	gcvt(T,3,filename1); 
	strcat(filename_tp,filename1);
	strcat(filename_tp,filename2);
	strcat(filename_tx,filename1);
	strcat(filename_tx,filename3);

	fp1=fopen(filename_tp,"w");

	fprintf(fp1,"%s%s\nI=%d,J=%d,k=%d,%s\n",head1,head2,Jy+1,Jx+1,1,head3);

		
	for(i=2;i<=Jx+2;i++)
		for(j=2;j<=Jy+2;j++)
		{
			if(i>N_x&&j<N_y)
			{	rou=0;
				u=0;
				v=0;
				p=0;
				a=0;
				Ma=0;
			}
			else
			{
				rou=U[i][j][0];
				u=U[i][j][1]/rou;
				v=U[i][j][2]/rou;
				p=(GAMA-1)*(U[i][j][3]-0.5*rou*(u*u+v*v));
				a=sqrt(GAMA*p/rou);
				Ma=sqrt(u*u+v*v)/a;
			}
		
			fprintf(fp1,"%20f\t%20.10e\t%20.10e\t%20.10e\t%20.10e\t%20.10e\t%20.10e\n",(i-2)*dx,(j-2)*dy,rou,u,v,p,Ma);
		}	
	fclose(fp1);

	fp2=fopen(filename_tx,"w");
	j=Jy/2+2;
	for(i=2;i<=Jx+1;i++)
	{
		rou=U[i][j][0];
		u=U[i][j][1]/rou;
		v=U[i][j][2]/rou;
		p=(GAMA-1)*(U[i][j][3]-0.5*rou*(u*u+v*v));
		a=sqrt(GAMA*p/rou);
		Ma=sqrt(u*u+v*v)/a;

		fprintf(fp2,"%20f\t%20.10e\t%20.10e\t%20.10e\t%20.10e\t%20.10e\n",(i-2)*dx,rou,u,v,p,Ma);
	}
	fclose(fp2);
}
//主函数//
int main()
{
	double T,dx,dy,dt,err;
	double U_old[Jx+5][Jy+5][4];
	double n=0;
	FILE *fp;
	fp=fopen("err.txt","w");
	fprintf(fp,"n\t,err\n");
	err=1e-2;
	Init(U,dx,dy);
	bound(U,dx,dy);
	T=0;
	while(T<TT&&err>1e-5&&n<30000)
	{
		n+=1;
		dt=CFL(U,dx,dy);
		T+=dt;
		printf("T=%10g	dt=%10g	err=%f	rou=%10f\n",T,dt,err,U[N_x][N_y][0]);
		if(fmod(n,400)==0){printf("n=%f\n",n);	Output(U,dx,dy,T);}
		save_U_last(U_old,U);
		symmetry_TVD(U,dx,dy,dt);
		err=error(U_old,U);
		fprintf(fp,"%d\t%10.10e\n",n,err);
	}
	fclose(fp);
	Output(U,dx,dy,T);
}
