/*..............................
MUSCL格式
..............................*/
#include<cstdlib>
#include<iostream>
#include<fstream>
#include<math.h>
#include<iomanip>

using namespace std;

#define Nx 400    //x方向网格数
#define Ny 200   //y方向网格数
#define Lx 4.0    //计算区域x方向长度
#define Ly 2.0    //计算区域y方向长度

#define GAMA 1.4  //气体常数

#define TT 2   //计算总时间
#define Step 50
#define CFL 0.1

#define MIN(x,y)(((x)<(y))?(x):(y))

//全局变量 
double U[Nx+5][Ny+5][4],FL[Nx+5][Ny+5][4], FR[Nx+5][Ny+5][4], GL[Nx+5][Ny+5][4],GR[Nx+5][Ny+5][4]; 
double k=-1, b=1;   
double ULX[Nx+5][Ny+5][4],URX[Nx+5][Ny+5][4],ULY[Nx+5][Ny+5][4],URY[Nx+5][Ny+5][4];

/*定义初始条件和边界条件*/
void Initial()              
{
	int i,j;
	double rou1=1.0,u1=0,v1=0,a1=1,p1=0.71429,dx,dy;
	dx=Lx/Nx;
	dy=Ly/Ny;

	for(i=0;i<=Nx+4;i++)
	{
		for(j=0;j<=Ny+4;j++)  //全部赋值为rou1,u1,v1,p1
		{		
				U[i][j][0]=rou1;
                U[i][j][1]=rou1*u1;
                U[i][j][2]=rou1*v1;
	            U[i][j][3]=p1/(GAMA-1)+rou1*(u1*u1+v1*v1)/2;
		}
	}

	return;
}



void bound(double U[Nx+5][Ny+5][4])
{
	int i,j,a,b,s;
	double rou1=1,u1=0,v1=0,a1=1,p1=0.71429;
	double rou2=3.85714,p2=7.381,u2=2.22223,v2=0;
	double dx,dy;
	dx=Lx/Nx;
	dy=Ly/Ny;

	for(i=0;i<=2;i++)//左边界入口条件
	{
		for(j=0;j<=Ny+4;j++)
		{
			U[i][j][0]=rou2;
			U[i][j][1]=rou2*u2;
			U[i][j][2]=rou2*v2;
			U[i][j][3]=p2/(GAMA-1)+0.5*rou2*(u2*u2+v2*v2);
		}
	}
	for(i=Nx+3;i<=Nx+4;i++)//右边界出口条件
	{
		for(j=0;j<=Ny+4;j++)
		{
			U[i][j][0]=U[i-1][j][0];
			U[i][j][1]=U[i-1][j][1];
			U[i][j][2]=U[i-1][j][2];
			U[i][j][3]=U[i-1][j][3];
		}
	}
//虚拟节点的处理,这里采用镜面反射,标量和切向速度取相同值,法相速度取相反值
	
	for(i=3;i<=Nx+2;i++)//上下壁面
	{
	for(j=0;j<=1;j++)
		{
			U[i][j][0]=U[i][-j+4][0];
			U[i][j][1]=U[i][-j+4][1];
			U[i][j][2]=-U[i][-j+4][2];
			U[i][j][3]=U[i][-j+4][3];
		}
	for(j=Ny+3;j<=Ny+4;j++)
		{
			U[i][j][0]=U[i][-j+2*Ny+4][0];
			U[i][j][1]=U[i][-j+2*Ny+4][1];
			U[i][j][2]=-U[i][-j+2*Ny+4][2];
			U[i][j][3]=U[i][-j+2*Ny+4][3];
		}
	}

	for(i=2;i<=Nx+2;i++)//上下壁面边界法向强行赋值为零
		{
		U[i][2][2]=0;
		U[i][Ny+2][2]=0;
		}

	//左方柱

	for(j=int(0.8/dy)+5;j<=int(1.2/dy)-1;j++)//左方柱左右表面虚拟点
	{
		for(i=int(1.0/dx)+3;i<=int(1.0/dx)+4;i++)
		{
			U[i][j][0]=U[2*int(1.0/dx)+4-i][j][0];
			U[i][j][1]=-U[2*int(1.0/dx)+4-i][j][1];
			U[i][j][2]=U[2*int(1.0/dx)+4-i][j][2];
			U[i][j][3]=U[2*int(1.0/dx)+4-i][j][3];
		}
		for(i=int(1.2/dx);i<=int(1.2/dx)+1;i++)
			{		
			U[i][j][0]=U[2*int(1.2/dx)+4-i][j][0];
			U[i][j][1]=-U[2*int(1.2/dx)+4-i][j][1];
			U[i][j][2]=U[2*int(1.2/dx)+4-i][j][2];
			U[i][j][3]=U[2*int(1.2/dx)+4-i][j][3];	
			}
	}	
	for(i=int(1.0/dx)+5;i<=int(1.2/dx)-1;i++)//左方柱上下表面虚拟点
	{
		for(j=int(1.2/dy);j<=int(1.2/dy)+1;j++)
		{
			U[i][j][0]=U[i][-j+2*int(1.2/dy)+4][0];
			U[i][j][1]=U[i][-j+2*int(1.2/dy)+4][1];
			U[i][j][2]=-U[i][-j+2*int(1.2/dy)+4][2];
			U[i][j][3]=U[i][-j+2*int(1.2/dy)+4][3];
		}
		for(j=int(0.8/dy)+3;j<=int(0.8/dy)+4;j++)
		{
			U[i][j][0]=U[i][-j+2*int(0.8/dy)+4][0];
			U[i][j][1]=U[i][-j+2*int(0.8/dy)+4][1];
			U[i][j][2]=-U[i][-j+2*int(0.8/dy)+4][2];
			U[i][j][3]=U[i][-j+2*int(0.8/dy)+4][3];
		}
	}
	//左方柱 角点的虚拟网格处理
	{
	a=int(1.2/dx)+2;    //方柱右下角角凸点(a,b)
	b=int(0.8/dy)+2;

		U[a][b][0]=0.5*(U[a][b-1][0]+U[a+1][b][0]);
		U[a][b][1]=0;
		U[a][b][2]=0;
		U[a][b][3]=0.5*(U[a][b-1][3]+U[a+1][b][3]);

		U[a-1][b+2][0]=U[a+1][b+2][0];
		U[a-1][b+2][1]=-U[a+1][b+2][1];
		U[a-1][b+2][2]=U[a+1][b+2][2];
		U[a-1][b+2][3]=U[a+1][b+2][3];

		U[a-2][b+1][0]=U[a-2][b-1][0];
		U[a-2][b+1][1]=U[a-2][b-1][1];
		U[a-2][b+1][2]=-U[a-2][b-1][2];
		U[a-2][b+1][3]=U[a-2][b-1][3];

		U[a-1][b+1][0]=0.5*(U[a-1][b-1][0]+U[a+1][b+1][0]);
		U[a-1][b+1][1]=0.5*(U[a-1][b-1][1]+U[a+1][b+1][1]);
		U[a-1][b+1][2]=0.5*(U[a-1][b-1][2]+U[a+1][b+1][2]);
		U[a-1][b+1][3]=0.5*(U[a-1][b-1][3]+U[a+1][b+1][3]);

		U[a-2][b+2][0]=U[a][b][0];
		U[a-2][b+2][1]=-U[a][b][1];
		U[a-2][b+2][2]=-U[a][b][2];
		U[a-2][b+2][3]=U[a][b][3];
	
	a=int(1.0/dx)+2;//方柱左下角
	b=int(0.8/dy)+2;

		U[a][b][0]=0.5*(U[a-1][b][0]+U[a][b-1][0]);
		U[a][b][1]=0;
		U[a][b][2]=0;
		U[a][b][3]=0.5*(U[a-1][b][3]+U[a][b-1][3]);
	
		U[a+1][b+2][0]=U[a-1][b+2][0];
		U[a+1][b+2][1]=-U[a-1][b+2][1];
		U[a+1][b+2][2]=U[a-1][b+2][2];
		U[a+1][b+2][3]=U[a-1][b+2][3];

		U[a+2][b+1][0]=U[a+2][b-1][0];
		U[a+2][b+1][1]=U[a+2][b-1][1];
		U[a+2][b+1][2]=-U[a+2][b-1][2];
		U[a+2][b+1][3]=U[a+2][b-1][3];

		U[a+1][b+1][0]=0.5*(U[a-1][b+1][0]+U[a+1][b-1][0]);
		U[a+1][b+1][1]=0.5*(U[a-1][b+1][1]+U[a+1][b-1][1]);
		U[a+1][b+1][2]=0.5*(U[a-1][b+1][2]+U[a+1][b-1][2]);
		U[a+1][b+1][3]=0.5*(U[a-1][b+1][3]+U[a+1][b-1][3]);

		U[a+2][b+2][0]=U[a][b][0];
		U[a+2][b+2][1]=-U[a][b][1];
		U[a+2][b+2][2]=-U[a][b][2];
		U[a+2][b+2][3]=U[a][b][3];
	
	a=int(1.0/dx)+2;//方柱左上角
	b=int(1.2/dy)+2;

		U[a][b][0]=0.5*(U[a][b+1][0]+U[a-1][b][0]);
		U[a][b][1]=0;
		U[a][b][2]=0;
		U[a][b][3]=0.5*(U[a][b+1][3]+U[a-1][b][3]);
		
		U[a+1][b-2][0]=U[a-1][b-2][0];
		U[a+1][b-2][1]=-U[a-1][b-2][1];
		U[a+1][b-2][2]=U[a-1][b-2][2];
		U[a+1][b-2][3]=U[a-1][b-2][3];

		U[a+2][b-1][0]=U[a+2][b+1][0];
		U[a+2][b-1][1]=U[a+2][b+1][1];
		U[a+2][b-1][2]=-U[a+2][b+1][2];
		U[a+2][b-1][3]=U[a+2][b+1][3];

		U[a+1][b-1][0]=0.5*(U[a-1][b-1][0]+U[a+1][b+1][0]);
		U[a+1][b-1][1]=0.5*(U[a-1][b-1][1]+U[a+1][b+1][1]);
		U[a+1][b-1][2]=0.5*(U[a-1][b-1][2]+U[a+1][b+1][2]);
		U[a+1][b-1][3]=0.5*(U[a-1][b-1][3]+U[a+1][b+1][3]);

		U[a+2][b-2][0]=U[a][b][0];
		U[a+2][b-2][1]=-U[a][b][1];
		U[a+2][b-2][2]=-U[a][b][2];
		U[a+2][b-2][3]=U[a][b][3];

	a=int(1.2/dx)+2;//方柱右上角
	b=int(1.2/dy)+2;

		U[a][b][0]=0.5*(U[a][b+1][0]+U[a+1][b][0]);
		U[a][b][1]=0;
		U[a][b][2]=0;
		U[a][b][3]=0.5*(U[a][b+1][3]+U[a+1][b][3]);
	
		U[a-1][b-2][0]=U[a+1][b-2][0];
		U[a-1][b-2][1]=-U[a+1][b-2][1];
		U[a-1][b-2][2]=U[a+1][b-2][2];
		U[a-1][b-2][3]=U[a+1][b-2][3];

		U[a-2][b-1][0]=U[a-2][b+1][0];
		U[a-2][b-1][1]=U[a-2][b+1][1];
		U[a-2][b-1][2]=-U[a-2][b+1][2];
		U[a-2][b-1][3]=U[a-2][b+1][3];

		U[a-1][b-1][0]=0.5*(U[a-1][b+1][0]+U[a+1][b-1][0]);
		U[a-1][b-1][1]=0.5*(U[a-1][b+1][1]+U[a+1][b-1][1]);
		U[a-1][b-1][2]=0.5*(U[a-1][b+1][2]+U[a+1][b-1][2]);
		U[a-1][b-1][3]=0.5*(U[a-1][b+1][3]+U[a+1][b-1][3]);

		U[a-2][b-2][0]=U[a][b][0];
		U[a-2][b-2][1]=-U[a][b][1];
		U[a-2][b-2][2]=-U[a][b][2];
		U[a-2][b-2][3]=U[a][b][3];
	}

	//边界上强行赋法向为零

	
	for(i=int(1.0/dx)+2;i<=int(1.2/dx)+2;i++)//左方柱上下壁面边界
		{
		U[i][int(1.2/dy)+2][2]=0;
		U[i][int(0.8/dy)+2][2]=0;
		}
	for(j=int(0.8/dy)+2;j<=int(1.2/dy)+2;j++)//左方柱左右边界
		{
		U[int(1.0/dx)+2][j][1]=0;
		U[int(1.2/dx)+2][j][1]=0;
		}
	//左方柱内部赋值为零

	for(i=int(1.0/dx)+5;i<=int(1.2/dy)-1;i++)
			for(j=int(0.8/dy)+5;j<=int(1.2/dy)-1;j++)
						for(s=0;s<=3;s++)
						U[i][j][s]=0;


//右方柱
	for(j=int(0.8/dy)+5;j<=int(1.2/dy)-1;j++) //右方柱左右表面虚拟点
	{
		for(i=int(2.0/dx)+3;i<=int(2.0/dx)+4;i++)
		{
			U[i][j][0]=U[2*int(2.0/dx)+4-i][j][0];
			U[i][j][1]=-U[2*int(2.0/dx)+4-i][j][1];
			U[i][j][2]=U[2*int(2.0/dx)+4-i][j][2];
			U[i][j][3]=U[2*int(2.0/dx)+4-i][j][3];
		}
		for(i=int(2.2/dx);i<=int(2.2/dx)+1;i++)
			{		
			U[i][j][0]=U[2*int(2.2/dx)+4-i][j][0];
			U[i][j][1]=-U[2*int(2.2/dx)+4-i][j][1];
			U[i][j][2]=U[2*int(2.2/dx)+4-i][j][2];
			U[i][j][3]=U[2*int(2.2/dx)+4-i][j][3];	
			}
	}	
	for(i=int(2.0/dx)+5;i<=int(2.2/dx)-1;i++) //方柱上下表面虚拟点
	{
		for(j=int(1.2/dy);j<=int(1.2/dy)+1;j++)
		{
			U[i][j][0]=U[i][-j+2*int(1.2/dy)+4][0];
			U[i][j][1]=U[i][-j+2*int(1.2/dy)+4][1];
			U[i][j][2]=-U[i][-j+2*int(1.2/dy)+4][2];
			U[i][j][3]=U[i][-j+2*int(1.2/dy)+4][3];
		}
		for(j=int(0.8/dy)+3;j<=int(0.8/dy)+4;j++)
		{
			U[i][j][0]=U[i][-j+2*int(0.8/dy)+4][0];
			U[i][j][1]=U[i][-j+2*int(0.8/dy)+4][1];
			U[i][j][2]=-U[i][-j+2*int(0.8/dy)+4][2];
			U[i][j][3]=U[i][-j+2*int(0.8/dy)+4][3];
		}
	}
	//右角点的虚拟网格处理
	{
	a=int(2.2/dx)+2;    //方柱右下角角凸点(a,b)以及其附近虚拟节点处理
	b=int(0.8/dy)+2;

		U[a][b][0]=0.5*(U[a][b-1][0]+U[a+1][b][0]);
		U[a][b][1]=0;
		U[a][b][2]=0;
		U[a][b][3]=0.5*(U[a][b-1][3]+U[a+1][b][3]);

		U[a-1][b+2][0]=U[a+1][b+2][0];
		U[a-1][b+2][1]=-U[a+1][b+2][1];
		U[a-1][b+2][2]=U[a+1][b+2][2];
		U[a-1][b+2][3]=U[a+1][b+2][3];

		U[a-2][b+1][0]=U[a-2][b-1][0];
		U[a-2][b+1][1]=U[a-2][b-1][1];
		U[a-2][b+1][2]=-U[a-2][b-1][2];
		U[a-2][b+1][3]=U[a-2][b-1][3];

		U[a-1][b+1][0]=0.5*(U[a-1][b-1][0]+U[a+1][b+1][0]);
		U[a-1][b+1][1]=0.5*(U[a-1][b-1][1]+U[a+1][b+1][1]);
		U[a-1][b+1][2]=0.5*(U[a-1][b-1][2]+U[a+1][b+1][2]);
		U[a-1][b+1][3]=0.5*(U[a-1][b-1][3]+U[a+1][b+1][3]);

		U[a-2][b+2][0]=U[a][b][0];
		U[a-2][b+2][1]=-U[a][b][1];
		U[a-2][b+2][2]=-U[a][b][2];
		U[a-2][b+2][3]=U[a][b][3];
	
	a=int(2.0/dx)+2; //方柱左下角
	b=int(0.8/dy)+2;

		U[a][b][0]=0.5*(U[a-1][b][0]+U[a][b-1][0]);
		U[a][b][1]=0;
		U[a][b][2]=0;
		U[a][b][3]=0.5*(U[a-1][b][3]+U[a][b-1][3]);
	
		U[a+1][b+2][0]=U[a-1][b+2][0];
		U[a+1][b+2][1]=-U[a-1][b+2][1];
		U[a+1][b+2][2]=U[a-1][b+2][2];
		U[a+1][b+2][3]=U[a-1][b+2][3];

		U[a+2][b+1][0]=U[a+2][b-1][0];
		U[a+2][b+1][1]=U[a+2][b-1][1];
		U[a+2][b+1][2]=-U[a+2][b-1][2];
		U[a+2][b+1][3]=U[a+2][b-1][3];

		U[a+1][b+1][0]=0.5*(U[a-1][b+1][0]+U[a+1][b-1][0]);
		U[a+1][b+1][1]=0.5*(U[a-1][b+1][1]+U[a+1][b-1][1]);
		U[a+1][b+1][2]=0.5*(U[a-1][b+1][2]+U[a+1][b-1][2]);
		U[a+1][b+1][3]=0.5*(U[a-1][b+1][3]+U[a+1][b-1][3]);

		U[a+2][b+2][0]=U[a][b][0];
		U[a+2][b+2][1]=-U[a][b][1];
		U[a+2][b+2][2]=-U[a][b][2];
		U[a+2][b+2][3]=U[a][b][3];
	
	a=int(2.0/dx)+2;//方柱左上角
	b=int(1.2/dy)+2;

		U[a][b][0]=0.5*(U[a][b+1][0]+U[a-1][b][0]);
		U[a][b][1]=0;
		U[a][b][2]=0;
		U[a][b][3]=0.5*(U[a][b+1][3]+U[a-1][b][3]);
		
		U[a+1][b-2][0]=U[a-1][b-2][0];
		U[a+1][b-2][1]=-U[a-1][b-2][1];
		U[a+1][b-2][2]=U[a-1][b-2][2];
		U[a+1][b-2][3]=U[a-1][b-2][3];

		U[a+2][b-1][0]=U[a+2][b+1][0];
		U[a+2][b-1][1]=U[a+2][b+1][1];
		U[a+2][b-1][2]=-U[a+2][b+1][2];
		U[a+2][b-1][3]=U[a+2][b+1][3];

		U[a+1][b-1][0]=0.5*(U[a-1][b-1][0]+U[a+1][b+1][0]);
		U[a+1][b-1][1]=0.5*(U[a-1][b-1][1]+U[a+1][b+1][1]);
		U[a+1][b-1][2]=0.5*(U[a-1][b-1][2]+U[a+1][b+1][2]);
		U[a+1][b-1][3]=0.5*(U[a-1][b-1][3]+U[a+1][b+1][3]);

		U[a+2][b-2][0]=U[a][b][0];
		U[a+2][b-2][1]=-U[a][b][1];
		U[a+2][b-2][2]=-U[a][b][2];
		U[a+2][b-2][3]=U[a][b][3];

	a=int(2.2/dx)+2;//方柱右上角
	b=int(1.2/dy)+2;

		U[a][b][0]=0.5*(U[a][b+1][0]+U[a+1][b][0]);
		U[a][b][1]=0;
		U[a][b][2]=0;
		U[a][b][3]=0.5*(U[a][b+1][3]+U[a+1][b][3]);
	
		U[a-1][b-2][0]=U[a+1][b-2][0];
		U[a-1][b-2][1]=-U[a+1][b-2][1];
		U[a-1][b-2][2]=U[a+1][b-2][2];
		U[a-1][b-2][3]=U[a+1][b-2][3];

		U[a-2][b-1][0]=U[a-2][b+1][0];
		U[a-2][b-1][1]=U[a-2][b+1][1];
		U[a-2][b-1][2]=-U[a-2][b+1][2];
		U[a-2][b-1][3]=U[a-2][b+1][3];

		U[a-1][b-1][0]=0.5*(U[a-1][b+1][0]+U[a+1][b-1][0]);
		U[a-1][b-1][1]=0.5*(U[a-1][b+1][1]+U[a+1][b-1][1]);
		U[a-1][b-1][2]=0.5*(U[a-1][b+1][2]+U[a+1][b-1][2]);
		U[a-1][b-1][3]=0.5*(U[a-1][b+1][3]+U[a+1][b-1][3]);

		U[a-2][b-2][0]=U[a][b][0];
		U[a-2][b-2][1]=-U[a][b][1];
		U[a-2][b-2][2]=-U[a][b][2];
		U[a-2][b-2][3]=U[a][b][3];
	}
	//边界上强行赋法向为零
	
	for(i=int(2.0/dx)+2;i<=int(2.2/dx)+2;i++)//方柱上下壁面边界
		{
		U[i][int(1.2/dy)+2][2]=0;
		U[i][int(0.8/dy)+2][2]=0;
		}
	for(j=int(0.8/dy)+2;j<=int(1.2/dy)+2;j++)//方柱左右边界
		{
		U[int(2.0/dx)+2][j][1]=0;
		U[int(2.2/dx)+2][j][1]=0;
		}

	//方柱内部赋值为零
	for(i=int(2.0/dx)+5;i<=int(2.2/dy)-1;i++)
			for(j=int(0.8/dy)+5;j<=int(1.2/dy)-1;j++)
						for(s=0;s<=3;s++)
						U[i][j][s]=0;			
}
//时间步长
double TIME(double U[Nx+5][Ny+5][4],double dx, double dy)
{
	int i,j;
	double u,v,rou,a,vel,p,maxvel;
	maxvel=maxvel=1e-100;
	for(i=0;i<=Nx+4;i++)
	{
		for(j=2;j<=Ny+2;j++)
		{
			if(U[i][j][0]!=0)
			{
				rou=U[i][j][0];
				u=U[i][j][1]/U[i][j][0];
				v=U[i][j][2]/U[i][j][0];
				p=(GAMA-1)*(U[i][j][3]-0.5*rou*(u*u+v*v));
				a=pow(GAMA*p/rou,0.5);
				vel=a+pow(u*u+v*v,0.5);
				if(vel>=maxvel)
					maxvel=vel;
			}
		}
	}
	return CFL*MIN(dx,dy)/maxvel;
}
//Sign函数
double sign(double va)
{
	if(va>0)
	{return 1;}
	if(va<0)
	{return -1;}
	else
	{return 0;}
}
//mimmod 函数
double minmod(double w1,double w2)
{
	double result;
	if(w1*w2>0)
	{
		result=sign(w1)*MIN(fabs(w1),fabs(w2));
	}
	else result=0;
	return result;
}


/*计算f+*/
void U2FL(double U[Nx+5][Ny+5][4],double FL[Nx+5][Ny+5][4],int i,int j)
{
	double rou,u,v,p,a,h,f1,f2,f3;
	int s=0;
	for(s=0;s<=3;s++)
	{
        ULX[i][j][s]=U[i][j][s]+0.25*(1-k)*(minmod((U[i][j][s]-U[i-1][j][s]),(b*(U[i+1][j][s]-U[i][j][s]))))
			+0.25*(1+k)*(minmod((U[i+1][j][s]-U[i][j][s]),(b*(U[i][j][s]-U[i-1][j][s]))));
	}
    rou=ULX[i][j][0];
	u=ULX[i][j][1]/ULX[i][j][0];
	v=ULX[i][j][2]/ULX[i][j][0];
	p=(GAMA-1)*(ULX[i][j][3]-0.5*(u*u+v*v)*rou);   	 
	a=sqrt(GAMA*p/rou);
	h=p/rou*GAMA/(GAMA-1)+(u*u+v*v)/2;

    f1=0.5*(u+sqrt(u*u+0.0001*0.0001));
	f2=0.5*(u-a+sqrt((u-a)*(u-a)+0.0001*0.0001));
	f3=0.5*(u+a+sqrt((u+a)*(u+a)+0.0001*0.0001));
	
	FL[i][j][0]=rou/2/GAMA*(2*(GAMA-1)*f1+f2+f3);	
	FL[i][j][1]=rou/2/GAMA*(2*(GAMA-1)*u*f1+(u-a)*f2+(u+a)*f3);
	FL[i][j][2]=rou/2/GAMA*(2*(GAMA-1)*f1+f2+f3)*v;
	FL[i][j][3]=rou/2/GAMA*((GAMA-1)*(u*u+v*v)*f1+(h-a*u)*f2+(h+a*u)*f3);
}

/*计算f-*/
void U2FR(double U[Nx+5][Ny+5][4],double FR[Nx+5][Ny+5][4],int i,int j)
{
	double rou,u,v,p,a,h,f1,f2,f3;
	int s=0;
	for(s=0;s<=3;s++)
	{
        URX[i][j][s]=U[i+1][j][s]-0.25*(1-k)*(minmod((U[i+2][j][s]-U[i+1][j][s]),(b*(U[i+1][j][s]-U[i][j][s]))))
			-0.25*(1+k)*(minmod((U[i+1][j][s]-U[i][j][s]),(b*(U[i+2][j][s]-U[i+1][j][s]))));
	} 
	rou=URX[i][j][0];
	u=URX[i][j][1]/URX[i][j][0];
	v=URX[i][j][2]/URX[i][j][0];
	p=(GAMA-1)*(URX[i][j][3]-0.5*rou*(u*u+v*v));
	a=sqrt(GAMA*p/rou);
	h=p/rou*GAMA/(GAMA-1)+(u*u+v*v)/2;

    f1=0.5*(u-sqrt(u*u+0.0001*0.0001));
	f2=0.5*(u-a-sqrt((u-a)*(u-a)+0.0001*0.0001));
	f3=0.5*(u+a-sqrt((u+a)*(u+a)+0.0001*0.0001));
				
	FR[i][j][0]=rou/2/GAMA*(2*(GAMA-1)*f1+f2+f3);	
	FR[i][j][1]=rou/2/GAMA*(2*(GAMA-1)*u*f1+(u-a)*f2+(u+a)*f3);	
	FR[i][j][2]=rou/2/GAMA*(2*(GAMA-1)*f1+f2+f3)*v;	
	FR[i][j][3]=rou/2/GAMA*((GAMA-1)*(u*u+v*v)*f1+(h-a*u)*f2+(h+a*u)*f2);	
}

/*计算g+*/
void U2GL(double U[Nx+5][Ny+5][4],double GL[Nx+5][Ny+5][4],int i,int j)
 {
	double rou,u,v,p,a,h,g1,g2,g3;
	int s=0;
	for(s=0;s<=3;s++)
	{
        ULY[i][j][s]=U[i][j][s]+0.25*(1-k)*(minmod((U[i][j][s]-U[i][j-1][s]),(b*(U[i][j+1][s]-U[i][j][s]))))
			 +0.25*(1+k)*(minmod((U[i][j+1][s]-U[i][j][s]),(b*(U[i][j][s]-U[i][j-1][s]))));
	}
	rou=ULY[i][j][0];
	u=ULY[i][j][1]/ULY[i][j][0];
	v=ULY[i][j][2]/ULY[i][j][0];
	p=(GAMA-1)*(ULY[i][j][3]-0.5*(u*u+v*v)*rou);
	a=sqrt(GAMA*p/rou);
	h=p/rou*GAMA/(GAMA-1)+(u*u+v*v)/2;
	
	g1=0.5*(v+sqrt(v*v+0.0001*0.0001));
	g2=0.5*(v-a+sqrt((v-a)*(v-a)+0.0001*0.0001));
	g3=0.5*(v+a+sqrt((v+a)*(v+a)+0.0001*0.0001));
	
    GL[i][j][0]=rou/2/GAMA*(2*(GAMA-1)*g1+g2+g3);
    GL[i][j][1]=rou/2/GAMA*(2*(GAMA-1)*g1+g2+g3)*u;
    GL[i][j][2]=rou/2/GAMA*(2*(GAMA-1)*g1*v+(v-a)*g2+(v+a)*g3);
	GL[i][j][3]=rou/2/GAMA*((GAMA-1)*(u*u+v*v)*g1+(h-a*v)*g2+(h+a*v)*g3);	
}

/*计算g-*/
void U2GR(double U[Nx+5][Ny+5][4],double GR[Nx+5][Ny+5][4],int i,int j)
{
	double rou,u,v,p,a,h,g1,g2,g3;
	int s=0;
	for(s=0;s<=3;s++)
	{
		URY[i][j][s]=U[i][j+1][s]-0.25*(1-k)*(minmod((U[i][j+2][s]-U[i][j+1][s]),(b*(U[i][j+1][s]-U[i][j][s]))))
			-0.25*(1+k)*(minmod((U[i][j+1][s]-U[i][j][s]),(b*(U[i][j+2][s]-U[i][j+1][s]))));
	}
	rou=URY[i][j][0];
	u=URY[i][j][1]/URY[i][j][0];
	v=URY[i][j][2]/URY[i][j][0];
	p=(GAMA-1)*(URY[i][j][3]-0.5*(u*u+v*v)*rou);	
	a=sqrt(GAMA*p/rou);
	h=p/rou*GAMA/(GAMA-1)+(u*u+v*v)/2;
	
    g1=0.5*(v-sqrt(v*v+0.0001));
	g2=0.5*(v-a-sqrt((v-a)*(v-a)+0.0001*0.0001));
	g3=0.5*(v+a-sqrt((v+a)*(v+a)+0.0001*0.0001));
	
    GR[i][j][0]=rou/2/GAMA*(2*(GAMA-1)*g1+g2+g3);
    GR[i][j][1]=rou/2/GAMA*(2*(GAMA-1)*g1+g2+g3)*u;
    GR[i][j][2]=rou/2/GAMA*(2*(GAMA-1)*g1*v+(v-a)*g2+(v+a)*g3);
	GR[i][j][3]=rou/2/GAMA*((GAMA-1)*(u*u+v*v)*g1+(h-a*v)*g2+(h+a*v)*g3);
}
//MUSCL
void MUSCL(double U[Nx+5][Ny+5][4],double dx,double dy,double dt)
{
	int i,j,s;
	for(i=0;i<=Nx+4;i++)
		{
			for(j=0;(j<=Ny+4);j++)
			{
				if(U[i][j][0]!=0)
				{
					U2FL(U,FL,i,j); 
					U2FR(U,FR,i,j); 
					U2GL(U,GL,i,j);
					U2GR(U,GR,i,j);
				}
            }
		}	
		for(i=2;i<=Nx+2;i++)
		{
			for(j=2;(j<=Ny+2);j++)
			{
				if(i>=int(1.0/dx)+5&&i<=int(1.2/dx)-1&&j>=int(0.8/dy)+5&&j<=int(1.2/dy)-1)
				{
						for(s=0;s<=3;s++)
						{
						U[i][j][s]=0;
						}
				}
				else if(i>=int(2.0/dx)+5&&i<=int(2.2/dx)-1&&j>=int(0.8/dy)+5&&j<=int(1.2/dy)-1)
				{
						for(s=0;s<=3;s++)
						{
						U[i][j][s]=0;
						}
				}
				else 
					{
						for(s=0;s<=3;s++)
						{   
					   U[i][j][s]=U[i][j][s]-dt/dx*((FL[i][j][s]+FR[i][j][s])-(FL[i-1][j][s]+FR[i-1][j][s]))-dt/dy*((GL[i][j][s]+GR[i][j][s])-(GL[i][j-1][s]+GR[i][j-1][s]));
						}
					}
			 }
		 }
}
/*输出函数*/ 
int Output(double U[Nx+5][Ny+5][4],int m,int n)
{
	int i,j;
	double rou,u,v,p,dx,dy;
	char Name[]={"CFD000.plt"};
	FILE *fp;
	dx=Lx/Nx;
	dy=Ly/Ny;
	if(n%Step==0)
	{
		Name[5]=char(m%10+48);
		Name[4]=char((m/10)%10+48);
		Name[3]=char(m/100+48);
		fp=fopen(Name,"w");
		fprintf(fp,"TITLE     = \"Dataset\"\nVARIABLES = \"x\" \"y\" \"rou\" \"u\" \"v\" \"p\" ");
		fprintf(fp,"ZONE T=\"Zone 1\"\nI=%d J=%d K=%d ZONETYPE=Ordered\n",Ny+1,Nx+1,1);
		fprintf(fp,"DATAPACKING=POINT\n");
	for(i=2;i<=Nx+2;i++)
		for(j=2;j<=Ny+2;j++)
		{			
			if(i>int(1.0/dx)+2&&i<int(1.2/dx)+2&&j>int(0.8/dy)+2&&j<int(1.2/dy)+2)
				{
					rou=0;
					u=0;
					v=0;
					p=0;
				}
			else if(i>int(2.0/dx)+2&&i<int(2.2/dx)+2&&j>int(0.8/dy)+2&&j<int(1.2/dy)+2)
				{
					rou=0;
					u=0;
					v=0;
					p=0;
				}
			else 
				{
					rou=U[i][j][0];
					u=U[i][j][1]/U[i][j][0];
					v=U[i][j][2]/U[i][j][0];
					p=(GAMA-1)*(U[i][j][3]-0.5*rou*(u*u+v*v));
				}				
		fprintf(fp,"%f  %f  %f  %f   %f    %f\n",(i-2)*dx,(j-2)*dy,rou,u,v,p);
		}
	fclose(fp);
		return 1;
	}
	else
		return 0;
}

int main()
{
	
	double T=0,dt=0;
	int m1=0,n1=0;
	double dx,dy;
	dx=Lx/Nx;
	dy=Ly/Ny;
	Initial();
	Output(U,m1,n1);

	while(T<TT)  	
	{
	   bound(U);
       dt=TIME(U,dx,dy);
	   T=dt+T;
	MUSCL( U, dx,dy,dt);

		bound(U);
		printf("T=%10g    dt=%10g\n",T,dt);
        n1++;
        m1=m1+Output(U,m1,n1);
    } 
	return 0;
}

