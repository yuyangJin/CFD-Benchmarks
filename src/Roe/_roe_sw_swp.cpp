//-----------------Roe,MUSCL,迎风TVD,对称TVD,NND,ENO,WENO以及紧致格式8种算法计算激波通过突然扩腰的管道流动Cpp源程序-----------

//-------------各种算法的结果保存在D:\CFD\的相应文件夹下。这里各种格式对应的结果保存文件夹如下：
//Roe格式---Roe_      MUSCL格式---MSCL    迎风TVD格式-----uTVD   对称TVD格式----sTVD  
//NND格式-----NND_    ENO格式-----ENO_    WENO格式---WENO      紧致格式---CMPT




#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip> 
#include <string>

using namespace std;

#define GAMA 1.4
#define PI 3.1415926
#define MIN(x,y)(((x)<(y))?(x):(y))
#define MAX(x,y)(((x)>(y))?(x):(y))
#define Lx 4.0                  //x方向的无量纲长度
#define Ly 1.0                
#define TT 10                  //总时间
#define Nx 400		          //x方向的分的控制单元数	
#define Ny 100            //x方向的分的控制单元数
//各个方法的CFL数
#define ROECFL 0.4
#define MUSCLCFL 0.4
#define UPWINDTVDCFL 0.4
#define SYMTVDCFL 0.3
#define NNDCFL 0.4
#define ENOCFL 0.1
#define WENOCFL 0.2
#define COMPACTCFL 0.1
//各个方法中用到的常数

//MUSCL
#define K -1        //MUSCL算法中用到的两个常数 
#define BETA 1
//迎风TVD
#define DELTA 0.5    //Q修正系数
#define OMEGA 0.3       //g修正系数
//WENO
#define P 2             //模板不光滑性的放大系数
#define OPDT 4	

//全局变量
//首先是一些通用的全局变量 其中后缀p表示plus 后缀d表示decrease 下标_一般表示i+1/2点的值或者Roe平均以后的i+1/2点处的值
double U[Nx+7][Ny+7][4],U_[Nx+7][Ny+7][4];                        //其中U_表示Roe平均以后的U,U_的三个分量于U不同
double Ut[Nx+7][Ny+7][4],U1[Nx+7][Ny+7][4],U2[Nx+7][Ny+7][4];     //Ut是临时变量,U1和U2是Runge-Kutta方法中的临时变量u1和u2
double L_[Nx+7][Ny+7][4][4],R_[Nx+7][Ny+7][4][4];                 //L_和R_是Roe平均以后的左右特征向量
double a_[Nx+7][Ny+7][1],LAMDA_[Nx+7][Ny+7][4][4];                 //a_是Roe平均以后的声速,LAMDA_是Roe平均以后的特征值
double G[Nx+7][Ny+7][4],G_[Nx+7][Ny+7][4],F[Nx+7][Ny+7][4],F_[Nx+7][Ny+7][4];   //F,G是方程中的FG,F_表示离散后的F通量
//Roe
double alpha_[Nx+7][Ny+7][4],theta[Nx+7][Ny+7][4];       //alpha_表示Roe方法中的alpha=L*(U[i+1]-U[i]).theta表示特征向量叠加的和
//MUSCL
double U_L[Nx+7][Ny+7][4],U_R[Nx+7][Ny+7][4];            //U_L和U_R分别表示MUSCL方法中的ULR
double Fp[Nx+7][Ny+7][4],Fd[Nx+7][Ny+7][4];            //Fp和Fd表示对F作Sterger-Warming分裂以后的正负F 其中Fp=Fplus Fd=Fdecrease
double Gp[Nx+7][Ny+7][4],Gd[Nx+7][Ny+7][4];           //G表示y方向的量
//迎风TVD
double g_[Nx+7][Ny+7][4];  //这里alpha_,g_,g,gama_,Q_,theta分别表示迎风TVD算法中用到的中间变量,其意义见教材
double g[Nx+7][Ny+7][4],gama_[Nx+7][Ny+7][4],Q_[Nx+7][Ny+7][4];  //其中theta定义于Roe方法中的意义不同
//对称TVD---其用到的变量都在之前被包含了

//NND---------其用到的变量都在之前被包含了

//ENO
double q3p[Nx+7][Ny+7][4][3],q3d[Nx+7][Ny+7][4][3];     //这里q3p和q3d 表示模板的正负差值f
double F_p[Nx+7][Ny+7][4],F_d[Nx+7][Ny+7][4];        //这里F_p表示正的流通量F[i+1/2]
double G_p[Nx+7][Ny+7][4],G_d[Nx+7][Ny+7][4];
//WENO
double ISd[Nx+7][Ny+7][4][3],ISp[Nx+7][Ny+7][4][3];  //这里ISd和omegap以及alphap的定义见教材
double omegap[Nx+7][Ny+7][4][3],omegad[Nx+7][Ny+7][4][3]; 
double alphap[Nx+7][Ny+7][4][3],alphad[Nx+7][Ny+7][4][3];
//紧致格式
double fp[Nx+7][Ny+7][4],fd[Nx+7][Ny+7][4];    //紧致格式中用小写的f表示F
double gp[Nx+7][Ny+7][4],gd[Nx+7][Ny+7][4],p[Nx+7][Ny+7];  //p是压力

std::string output_dir = string("_result/");

//初始化函数
//作用:将全场赋初值
void Initial(double U[Nx+7][Ny+7][4],double &dx,double &dy)
{
	int i,j;
	dx=Lx/Nx;
	dy=Ly/Ny;
	double rou1=1,u1=0,v1=0,a1=1,p1=0.71429;
	double rou2=3.85714,p2=7.381,u2=2.22223,v2=0;
	for(i=0;i<=Nx+6;i++)    //初始条件  这里先把所有区域都赋值u1 v1 p1 rou1
		for(j=0;j<=Ny+6;j++)
		{
			U[i][j][0]=rou1;
			U[i][j][1]=rou1*u1;
			U[i][j][2]=rou1*v1;
			U[i][j][3]=p1/(GAMA-1)+0.5*rou1*(u1*u1+v1*v1);			
		}
    //边界以外赋零    
	for(i=0;i<=int(1.0/dx)-1;i++)          
		for(j=int(0.5/dy)+7;j<=Ny+6;j++)
		{
			U[i][j][0]=0;
			U[i][j][1]=0;
			U[i][j][2]=0;
			U[i][j][3]=0;	
		}
	for(i=int(2.0/dx)+7;i<=Nx+6;i++)          
		for(j=int(0.5/dy)+7;j<=Ny+6;j++)
		{
			U[i][j][0]=0;
			U[i][j][1]=0;
			U[i][j][2]=0;
			U[i][j][3]=0;			
		}
}
//CFL稳定性条件
//入口:U(没有进行过Roe平均的U向量)
//出口:dt
double CFL(double U[Nx+7][Ny+7][4],double dx,double dy,double cfl)
{
	int i,j;
	double u,v,rou,a,vel,maxvel,p;
	maxvel=pow(10,-100);
	for(i=0;i<=Nx+6;i++)    
		for(j=0;j<=Ny+6;j++)
			if(U[i][j][0]!=0)
			{
				rou=U[i][j][0];
				u=U[i][j][1]/U[i][j][0];
				v=U[i][j][2]/U[i][j][0];
				p=(GAMA-1)*(U[i][j][3]-0.5*rou*(u*u+v*v));
				a=pow(p*GAMA/rou,0.5);
				vel=a+fabs(u);
				if(vel>=maxvel)maxvel=vel;
				vel=a+fabs(v);
				if(vel>=maxvel)maxvel=vel;
			}
	return cfl*MIN(dx,dy)/maxvel;
}
//CFL_稳定性条件
//入口:U_是进行了Roe平均的向量),a_是声速
//出口:dt
double CFL_(double U_[Nx+7][Ny+7][4],double a_[Nx+7][Ny+7][1],double dx,double dy,double cfl)
{
	int i,j;
	double u,v,rou,a,vel,maxvel;
	maxvel=pow(10,-100);
	for(i=0;i<=Nx+6;i++)    
		for(j=0;j<=Ny+6;j++)
			if(U[i][j][0]!=0)
			{
				rou=U_[i][j][0];
				u=U_[i][j][1];
				v=U_[i][j][2];
				a=a_[i][j][0];
				vel=a+pow((u*u+v*v),0.5);
				if(vel>=maxvel)maxvel=vel;
			}
	return cfl*MIN(dx,dy)/maxvel;
}
//边界函数
//入口:U
//出口:无
//作用:将左右边界以及虚拟节点赋值,并且处理特殊的凹凸角点以及虚拟节点,另外对边界上面的法相速度强制赋零
void bound(double U[Nx+7][Ny+7][4],double dx,double dy)
{
	int i,j,k;
	double rou1=1,u1=0,v1=0,a1=1,p1=0.71429;
	double rou2=3.85714,p2=7.381,u2=2.22223,v2=0;	
	for(i=0;i<=2;i++)         //左边界条件
		for(j=0;j<=int(0.5/dy)+6;j++)
		{
			U[i][j][0]=rou2;
			U[i][j][1]=rou2*u2;
			U[i][j][2]=rou2*v2;
			U[i][j][3]=p2/(GAMA-1)+0.5*rou2*(u2*u2+v2*v2);
		}
	for(i=Nx+4;i<=Nx+6;i++)         //右边界条件
		for(j=0;j<=int(0.5/dy)+6;j++)
		{
			U[i][j][0]=U[i-1][j][0];
			U[i][j][1]=U[i-1][j][1];
			U[i][j][2]=U[i-1][j][2];
			U[i][j][3]=U[i-1][j][3];
		}
	for(i=3;i<=Nx+3;i++)       //虚拟节点的处理,这里采用镜面反射,标量和切向速度取相同值,法相速度取相反值
		for(j=0;j<=2;j++)
		{
			U[i][j][0]=U[i][6-j][0];
			U[i][j][1]=U[i][6-j][1];
			U[i][j][2]=-U[i][6-j][2];
			U[i][j][3]=U[i][6-j][3];
		}
	for(i=3;i<=int(1.0/dx)+2;i++)
		for(j=int(0.5/dy)+4;j<=int(0.5/dy)+6;j++)
		{
			U[i][j][0]=U[i][-j+2*int(0.5/dy)+6][0];
			U[i][j][1]=U[i][-j+2*int(0.5/dy)+6][1];
			U[i][j][2]=-U[i][-j+2*int(0.5/dy)+6][2];
			U[i][j][3]=U[i][-j+2*int(0.5/dy)+6][3];
		}
	for(i=int(1.0/dx);i<=int(2.0/dx)+6;i++)
		for(j=Ny+4;j<=Ny+6;j++)
		{
			U[i][j][0]=U[i][-j+2*int(1.0/dy)+6][0];
			U[i][j][1]=U[i][-j+2*int(1.0/dy)+6][1];
			U[i][j][2]=-U[i][-j+2*int(1.0/dy)+6][2];
			U[i][j][3]=U[i][-j+2*int(1.0/dy)+6][3];
		}
	for(i=int(2.0/dx)+4;i<=Nx+3;i++)
		for(j=int(0.5/dy)+4;j<=int(0.5/dy)+6;j++)
		{
			U[i][j][0]=U[i][-j+2*int(0.5/dy)+6][0];
			U[i][j][1]=U[i][-j+2*int(0.5/dy)+6][1];
			U[i][j][2]=-U[i][-j+2*int(0.5/dy)+6][2];
			U[i][j][3]=U[i][-j+2*int(0.5/dy)+6][3];
		}
	for(i=int(1.0/dx);i<=int(1.0/dx)+2;i++)
		for(j=int(0.5/dy)+4;j<=Ny+6;j++)
		{
			U[i][j][0]=U[2*int(1.0/dx)+6-i][j][0];
			U[i][j][1]=-U[2*int(1.0/dx)+6-i][j][1];
			U[i][j][2]=U[2*int(1.0/dx)+6-i][j][2];
			U[i][j][3]=U[2*int(1.0/dx)+6-i][j][3];
		}
	for(i=int(2.0/dx)+4;i<=int(2.0/dx)+6;i++)
		for(j=int(0.5/dy)+4;j<=Ny+6;j++)
		{
			U[i][j][0]=U[2*int(2.0/dx)+6-i][j][0];
			U[i][j][1]=-U[2*int(2.0/dx)+6-i][j][1];
			U[i][j][2]=U[2*int(2.0/dx)+6-i][j][2];
			U[i][j][3]=U[2*int(2.0/dx)+6-i][j][3];
		} 
	//角点的虚拟网格特殊处理
	int a,b;
	a=int(1.0/dx)+2;       //角凸点(a,b)以及其附近虚拟节点处理,处理方法参考了教材相关内容
	b=int(0.5/dy)+4;
    for(k=0;k<=3;k++)
	{
		U[a][b][k]=0.5*(U[a][b-2][k]+U[a+2][b][k]);
		U[a-1][b][k]=U[a-1][b-2][k];
		U[a-2][b][k]=U[a-2][b-2][k];
		U[a][b+1][k]=U[a+2][b+1][k];
		U[a][b+2][k]=U[a+2][b+2][k];
		U[a-1][b+1][k]=0.5*(U[a-1][b-3][k]+U[a+3][b+1][k]);
		U[a-2][b+1][k]=U[a-2][b-3][k];
		U[a-1][b+2][k]=U[a+3][b+2][k];
		U[a-2][b+2][k]=0.5*(U[a-2][b-4][k]+U[a+4][b+2][k]);
	}
	U[a-1][b][2]=-U[a-1][b-2][2];
	U[a-2][b][2]=-U[a-2][b-2][2];
	U[a-2][b+1][2]=-U[a-2][b-3][2];
	U[a][b+1][1]=-U[a+2][b+1][1];
	U[a][b+2][1]=-U[a+2][b+2][1];
	U[a-1][b+2][1]=-U[a+3][b+2][1];
	a=int(2.0/dx)+4;
	b=int(0.5/dy)+4;
	for(k=0;k<=3;k++)
	{
		U[a][b][k]=0.5*(U[a][b-2][k]+U[a-2][b][k]);
		U[a+1][b][k]=U[a+1][b-2][k];
		U[a+2][b][k]=U[a+2][b-2][k];
		U[a][b+1][k]=U[a-2][b+1][k];
		U[a][b+2][k]=U[a-2][b+2][k];
		U[a+1][b+1][k]=0.5*(U[a+1][b-3][k]+U[a-3][b+1][k]);
		U[a+2][b+1][k]=U[a+2][b-3][k];
		U[a+1][b+2][k]=U[a-3][b+2][k];
		U[a+2][b+2][k]=0.5*(U[a+2][b-4][k]+U[a-4][b+2][k]);
	}
	U[a+1][b][2]=-U[a+1][b-2][2];
	U[a+2][b][2]=-U[a+2][b-2][2];
	U[a+2][b+1][2]=-U[a+2][b-3][2];
	U[a][b+1][1]=-U[a-2][b+1][1];
	U[a][b+2][1]=-U[a-2][b+2][1];
	U[a+1][b+2][1]=-U[a-3][b+2][1];
	//边界上强行赋法向为零
	for(i=3;i<=Nx+3;i++)
	{
		U[i][3][2]=0;
	}
	for(i=3;i<=int(1.0/dx)+2;i++)
	{
		U[i][int(0.5/dy)+3][2]=0;
	}
	for(i=int(2.0/dx)+4;i<=Nx+3;i++)
	{
		U[i][int(0.5/dy)+3][2]=0;
	}
	for(i=int(1.0/dx)+3;i<=int(2.0/dx)+3;i++)
	{
		U[i][Ny+3][2]=0;
	}
	for(j=int(0.5/dy)+4;j<=Ny+3;j++)
	{
		U[int(1.0/dx)+3][j][1]=0;
	}
	for(j=int(0.5/dy)+4;j<=Ny+3;j++)
	{
		U[int(2.0/dx)+3][j][1]=0;
	}
}
//求符号函数
//出口:大于零返回1,小于零返回-1,等于零返回0
double sign(double va)
{
	if(va>0)
	{return 1;}
	if(va<0)
	{return -1;}
	else
	{return 0;}
}
//返回3个数中的较少者
double MIN3(double a1,double a2,double a3)
{
	double re;
	if(a1<=a2&&a1<=a3)re=a1;
	else
	{
		if(a2<=a3&&a2<=a1)re=a2;
		else re=a3;
	}
	return re;
}
//minmod
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
//对称TVD中用到的3个数的minmod
double minmod3(double w1,double w2,double w3)
{
	double result;
	if(w1*w2*w3>0)
	{
		result=sign(w1)*MIN3(fabs(w1),fabs(w2),fabs(w3));
	}
	else result=0;
	return result;
}
//通过U求F
void U2F(double U[4],double F[4])
{
	double u,v,p,rou;
	rou=U[0];
	u=U[1]/U[0];
	v=U[2]/U[0];
	p=(GAMA-1)*(U[3]-0.5*rou*(u*u+v*v));
	F[0]=rou*u;
	F[1]=rou*u*u+p;
	F[2]=rou*u*v;
	F[3]=(U[3]+p)*u;
}
void U2G(double U[4],double G[4])
{
	double u,v,p,rou;
	rou=U[0];
	u=U[1]/U[0];
	v=U[2]/U[0];
	p=(GAMA-1)*(U[3]-0.5*rou*(u*u+v*v));
	G[0]=rou*v;
	G[1]=rou*u*v;
	G[2]=rou*v*v+p;
	G[3]=(U[3]+p)*v;
}














//---------------------------------Roe方法--------------------------
//其中用到的各个量的解释:
//U是待求向量,U_是Roe平均后的U向量,U_的四个分量分别是Roe平均以后的rou_,u_,v_,H=a*a/(GAMA-1)+0.5*(u*u+v*v)
//LAMDA_,L_和R_,a_分别是Roe平均以后的特征值以及左右特征向量和声速.
//alpha=L*(U[i+1]-U[i]),theta是各个特征值局部特征空间的叠加量.F和G分别是原方程中的FG
void A1_Roe_QuickLocate()
{
}

//x方向的Roe平均
//作用:将U做Roe平均以后返回U_和声速a_
void RoeAVG_x(double U[Nx+7][Ny+7][4],double U_[Nx+7][Ny+7][4],double a_[Nx+7][Ny+7][1])
{
	int i,j;
	double uL,uR,vL,vR,rouL,rouR,pL,pR,D,HL,HR;        //这里uL和uR分别表示u(i+1/2)点左右的速度
	for(i=0;i<=Nx+5;i++)    
		for(j=0;j<=Ny+5;j++)
			if(U[i][j][0]!=0&&U[i+1][j][0]!=0)
			{
				rouL=U[i][j][0];
				rouR=U[i+1][j][0];
				D=pow(rouR/rouL,0.5);
				uL=U[i][j][1]/U[i][j][0];
				uR=U[i+1][j][1]/U[i+1][j][0];
				vL=U[i][j][2]/U[i][j][0];
				vR=U[i+1][j][2]/U[i+1][j][0];
				pL=(GAMA-1)*(U[i][j][3]-0.5*rouL*(uL*uL+vL*vL));
				pR=(GAMA-1)*(U[i+1][j][3]-0.5*rouR*(uR*uR+vR*vR));
				HL=pL*GAMA/(GAMA-1)/rouL+(uL*uL+vL*vL)/2;
				HR=pR*GAMA/(GAMA-1)/rouR+(uR*uR+vR*vR)/2;
				U_[i][j][0]=rouL*(1+D)*(1+D)/4;        //U_表示roe平均以后的向量其4个分量分别是rou,u,v,H				                                       
				U_[i][j][1]=(uL+D*uR)/(1+D);             //注意它和U分量不同,其中U_[i]表示U[i+0.5]处roe平均的结果
				U_[i][j][2]=(vL+D*vR)/(1+D);
				U_[i][j][3]=(D*HR+HL)/(1+D);
				a_[i][j][0]=pow((GAMA-1)*(U_[i][j][3]-0.5*(U_[i][j][1]*U_[i][j][1]+U_[i][j][2]*U_[i][j][2])),0.5);
			}
}
void RoeAVG_y(double U[Nx+7][Ny+7][4],double U_[Nx+7][Ny+7][4],double a_[Nx+7][Ny+7][1])  // 这里U_[i][j][k]为在i+0.5处的Roe平均值
{
	int i,j;
	double uL,uR,vL,vR,rouL,rouR,pL,pR,D,HL,HR;
	for(i=0;i<=Nx+5;i++)    
		for(j=0;j<=Ny+5;j++)
			if(U[i][j][0]!=0&&U[i][j+1][0]!=0)
			{
				rouL=U[i][j][0];
				rouR=U[i][j+1][0];
				D=pow(rouR/rouL,0.5);
				uL=U[i][j][1]/U[i][j][0];
				uR=U[i][j+1][1]/U[i][j+1][0];
				vL=U[i][j][2]/U[i][j][0];
				vR=U[i][j+1][2]/U[i][j+1][0];
				pL=(GAMA-1)*(U[i][j][3]-0.5*rouL*(uL*uL+vL*vL));
				pR=(GAMA-1)*(U[i][j+1][3]-0.5*rouR*(uR*uR+vR*vR));
				HL=pL*GAMA/(GAMA-1)/rouL+(uL*uL+vL*vL)/2;
				HR=pR*GAMA/(GAMA-1)/rouR+(uR*uR+vR*vR)/2;
				U_[i][j][0]=rouL*(1+D)*(1+D)/4;          //U_表示roe平均以后的向量其3个分量分别是rou,u,v,H注意它和U分量不同
				U_[i][j][1]=(uL+D*uR)/(1+D);                 //其中U[i]表示U[i+0.5]处roe平均的结果
				U_[i][j][2]=(vL+D*vR)/(1+D);
				U_[i][j][3]=(D*HR+HL)/(1+D);
				a_[i][j][0]=pow((GAMA-1)*(U_[i][j][3]-0.5*(U_[i][j][1]*U_[i][j][1]+U_[i][j][2]*U_[i][j][2])),0.5);
			}
}
//作用:用Roe平均以后的U值求特征值与特征向量
void Roe_Eig_x(double U[Nx+7][Ny+7][4],double U_[Nx+7][Ny+7][4],double LAMDA_[Nx+7][Ny+7][4][4],double a_[Nx+7][Ny+7][1],
		   double R_[Nx+7][Ny+7][4][4],double L_[Nx+7][Ny+7][4][4])       //利用Roe平均以后的值计算x方向的特征值与特征向量
{
	int i,j,k,m,n;
	for(i=0;i<=Nx+5;i++)          //先将LAMDA_全赋值0
		for(j=0;j<=Ny+5;j++)
			if(U[i][j][0]!=0)
				for(m=0;m<=3;m++)
					for(n=0;n<=3;n++)
						LAMDA_[i][j][m][n]=0;

	for(i=0;i<=Nx+5;i++)        //对LAMDA_赋值
		for(j=0;j<=Ny+5;j++)
			if(U[i][j][0]!=0)
			{
				LAMDA_[i][j][0][0]=U_[i][j][1]-a_[i][j][0];
				LAMDA_[i][j][1][1]=U_[i][j][1];
				LAMDA_[i][j][2][2]=U_[i][j][1]+a_[i][j][0];
				LAMDA_[i][j][3][3]=U_[i][j][1];
				for(k=0;k<=3;k++)          //这里对lambda做一个修正,防止声速0点出现的问题
				{
					if(LAMDA_[i][j][k][k]>=0)
					{
						LAMDA_[i][j][k][k]=0.5*(LAMDA_[i][j][k][k]+pow(LAMDA_[i][j][k][k]*LAMDA_[i][j][k][k]+pow(10,-8),0.5));
					}
					else LAMDA_[i][j][k][k]=0.5*(LAMDA_[i][j][k][k]-pow(LAMDA_[i][j][k][k]*LAMDA_[i][j][k][k]+pow(10,-8),0.5));
				}
			}
//求右特征向量	
	double u,v,a,rou,H;
	for(i=0;i<=Nx+5;i++)
		for(j=0;j<=Ny+5;j++)
			if(U[i][j][0]!=0)
			{
				rou=U_[i][j][0];
				u=U_[i][j][1];
				v=U_[i][j][2];
				a=a_[i][j][0];
				H=U_[i][j][3];
				R_[i][j][0][0]=1;
				R_[i][j][0][1]=1;
				R_[i][j][0][2]=1;
				R_[i][j][0][3]=0;
				R_[i][j][1][0]=-a+u;
				R_[i][j][1][1]=u;
				R_[i][j][1][2]=u+a;
				R_[i][j][1][3]=0;
				R_[i][j][2][0]=v;
				R_[i][j][2][1]=v;
				R_[i][j][2][2]=v;
				R_[i][j][2][3]=1;
				R_[i][j][3][0]=H-a*u;
				R_[i][j][3][1]=0.5*(u*u+v*v);
				R_[i][j][3][2]=H+u*a;
				R_[i][j][3][3]=v;
			}
//求左特征向量
	double b1,b2;
	for(i=0;i<=Nx+5;i++)
		for(j=0;j<=Ny+5;j++)
			if(U[i][j][0]!=0)
			{
				rou=U_[i][j][0];
				u=U_[i][j][1];
				v=U_[i][j][2];
				a=a_[i][j][0];
				b2=(GAMA-1)/a/a;
				b1=b2*(u*u+v*v)/2;
				L_[i][j][0][0]=(b1+u/a)/2;
				L_[i][j][0][1]=-(b2*u+1/a)/2;
				L_[i][j][0][2]=-b2*v/2;
				L_[i][j][0][3]=b2/2;
				L_[i][j][1][0]=1-b1;
				L_[i][j][1][1]=b2*u;
				L_[i][j][1][2]=b2*v;
				L_[i][j][1][3]=-b2;
				L_[i][j][2][0]=0.5*(b1-u/a);
				L_[i][j][2][1]=0.5*(1/a-b2*u);
				L_[i][j][2][2]=-0.5*b2*v;
				L_[i][j][2][3]=b2/2;
				L_[i][j][3][0]=-v;
				L_[i][j][3][1]=0;
				L_[i][j][3][2]=1;
				L_[i][j][3][3]=0;
			}
}
//用Roe平均以后的U值求特征值与特征向量
void Roe_Eig_y(double U[Nx+7][Ny+7][4],double U_[Nx+7][Ny+7][4],double LAMDA_[Nx+7][Ny+7][4][4],double a_[Nx+7][Ny+7][1],
		   double R_[Nx+7][Ny+7][4][4],double L_[Nx+7][Ny+7][4][4])
{
	int i,j,k,m,n;
	for(i=0;i<=Nx+5;i++)          //先将LAMDA_全赋值0
		for(j=0;j<=Ny+5;j++)
			if(U[i][j][0]!=0)
				for(m=0;m<=3;m++)
					for(n=0;n<=3;n++)
						LAMDA_[i][j][m][n]=0;

	for(i=0;i<=Nx+5;i++)        //对LAMDA_赋值
		for(j=0;j<=Ny+5;j++)
			if(U[i][j][0]!=0)
			{
				LAMDA_[i][j][0][0]=U_[i][j][2]-a_[i][j][0];
				LAMDA_[i][j][1][1]=U_[i][j][2];
				LAMDA_[i][j][2][2]=U_[i][j][2]+a_[i][j][0];
				LAMDA_[i][j][3][3]=U_[i][j][2];
				for(k=0;k<=3;k++)
				{
					if(LAMDA_[i][j][k][k]>=0)
					{
						LAMDA_[i][j][k][k]=0.5*(LAMDA_[i][j][k][k]+pow(LAMDA_[i][j][k][k]*LAMDA_[i][j][k][k]+pow(10,-8),0.5));
					}
					else LAMDA_[i][j][k][k]=0.5*(LAMDA_[i][j][k][k]-pow(LAMDA_[i][j][k][k]*LAMDA_[i][j][k][k]+pow(10,-8),0.5));
				}
			}
//求右特征向量	
	double u,v,a,rou,H;
	for(i=0;i<=Nx+5;i++)
		for(j=0;j<=Ny+5;j++)
			if(U[i][j][0]!=0)
			{
				rou=U_[i][j][0];
				u=U_[i][j][1];
				v=U_[i][j][2];
				a=a_[i][j][0];
				H=U_[i][j][3];
				R_[i][j][0][0]=1;
				R_[i][j][0][1]=1;
				R_[i][j][0][2]=1;
				R_[i][j][0][3]=0;
				R_[i][j][1][0]=u;
				R_[i][j][1][1]=u;
				R_[i][j][1][2]=u;
				R_[i][j][1][3]=1;
				R_[i][j][2][0]=v-a;
				R_[i][j][2][1]=v;
				R_[i][j][2][2]=v+a;
				R_[i][j][2][3]=0;
				R_[i][j][3][0]=H-a*v;
				R_[i][j][3][1]=0.5*(u*u+v*v);
				R_[i][j][3][2]=H+v*a;
				R_[i][j][3][3]=u;
			}
//求左特征向量
	double b1,b2;
	for(i=0;i<=Nx+5;i++)
		for(j=0;j<=Ny+5;j++)
			if(U[i][j][0]!=0)
			{
				rou=U_[i][j][0];
				u=U_[i][j][1];
				v=U_[i][j][2];
				a=a_[i][j][0];
				b2=(GAMA-1)/a/a;
				b1=b2*(u*u+v*v)/2;
				L_[i][j][0][0]=(b1+v/a)/2;
				L_[i][j][0][1]=b2*u/2;
				L_[i][j][0][2]=-(b2*v+1/a)/2;
				L_[i][j][0][3]=b2/2;
				L_[i][j][1][0]=1-b1;
				L_[i][j][1][1]=b2*u;
				L_[i][j][1][2]=b2*v;
				L_[i][j][1][3]=-b2;
				L_[i][j][2][0]=0.5*(b1-v/a);
				L_[i][j][2][1]=-0.5*(b2*u);
				L_[i][j][2][2]=0.5*(1/a-b2*v);
				L_[i][j][2][3]=b2/2;
				L_[i][j][3][0]=-u;
				L_[i][j][3][1]=1;
				L_[i][j][3][2]=0;
				L_[i][j][3][3]=0;
			}
}
//---------------------Roe_x方向算子
void Roe_x(double U[Nx+7][Ny+7][4],double L_[Nx+7][Ny+7][4][4],double R_[Nx+7][Ny+7][4][4],double F_[Nx+7][Ny+7][4],
		   double F[Nx+7][Ny+7][4],double LAMDA_[Nx+7][Ny+7][4][4],double alpha_[Nx+7][Ny+7][4],
		   double theta[Nx+7][Ny+7][4],double dx,double dy,double &dt)
{
	int i,j,k,l;
	double r;
	//先进行Roe平均
	RoeAVG_x(U,U_,a_);
	//求dt
	dt=CFL_(U_,a_,dx,dy,ROECFL);  
	r=dt/dx;
	//求特征向量和特征值
	Roe_Eig_x(U,U_,LAMDA_,a_,R_,L_);
	//计算alpha_
	for(i=2;i<=Nx+3;i++)    
		for(j=2;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
				{
					alpha_[i][j][k]=0;
					theta[i][j][k]=0;
				}
	for(i=2;i<=Nx+3;i++)    
		for(j=2;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					for(l=0;l<=3;l++)  
						alpha_[i][j][k]+=L_[i][j][k][l]*(U[i+1][j][l]-U[i][j][l]);						

	//计算theta项
	for(i=2;i<=Nx+3;i++)    
		for(j=2;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					for(l=0;l<=3;l++)
						theta[i][j][k]+=fabs(LAMDA_[i][j][l][l])*alpha_[i][j][l]*R_[i][j][k][l];
	//U2F
	for(i=2;i<=Nx+4;i++)    
		for(j=2;j<=Ny+4;j++)
			if(U[i][j][0]!=0)
				U2F(U[i][j],F[i][j]);
	//计算F_
	for(i=2;i<=Nx+3;i++)    
		for(j=2;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					F_[i][j][k]=0.5*(F[i][j][k]+F[i+1][j][k]-theta[i][j][k]);			
    //分区计算U
	for(i=3;i<=Nx+3;i++)    
		for(j=3;j<=int(0.5/dy)+3;j++)
			for(k=0;k<=3;k++)
				U[i][j][k]=U[i][j][k]-r*(F_[i][j][k]-F_[i-1][j][k]);

	for(i=int(1.0/dx)+3;i<=int(2.0/dx)+3;i++)    
		for(j=int(0.5/dy)+4;j<=Ny+3;j++)
			for(k=0;k<=3;k++)
				U[i][j][k]=U[i][j][k]-r*(F_[i][j][k]-F_[i-1][j][k]);
}
//-------------------Roe_y方向算子
void Roe_y(double U[Nx+7][Ny+7][4],double L_[Nx+7][Ny+7][4][4],double R_[Nx+7][Ny+7][4][4],double G_[Nx+7][Ny+7][4],
		   double G[Nx+7][Ny+7][4],double LAMDA_[Nx+7][Ny+7][4][4],double alpha_[Nx+7][Ny+7][4],
		   double theta[Nx+7][Ny+7][4],double dx,double dy,double &dt)
{
	int i,j,k,l;
	double r;
	//先进行Roe平均
	RoeAVG_y(U,U_,a_);
	dt=CFL_(U_,a_,dx,dy,ROECFL);  
	r=dt/dy;
	//求特征向量和特征值
	Roe_Eig_y(U,U_,LAMDA_,a_,R_,L_);
	//计算alpha_
	for(i=2;i<=Nx+3;i++)    
		for(j=2;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
				{
					alpha_[i][j][k]=0;
					theta[i][j][k]=0;
				}
	for(i=2;i<=Nx+3;i++)    
		for(j=2;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					for(l=0;l<=3;l++)  
						alpha_[i][j][k]+=L_[i][j][k][l]*(U[i][j+1][l]-U[i][j][l]);						
	//计算theta项
	for(i=2;i<=Nx+3;i++)    
		for(j=2;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					for(l=0;l<=3;l++)
						theta[i][j][k]+=fabs(LAMDA_[i][j][l][l])*alpha_[i][j][l]*R_[i][j][k][l];
	//U2G
	for(i=2;i<=Nx+4;i++)    
		for(j=2;j<=Ny+4;j++)
			if(U[i][j][0]!=0)
				U2G(U[i][j],G[i][j]);

	//计算G_
	for(i=2;i<=Nx+3;i++)    
		for(j=2;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					G_[i][j][k]=0.5*(G[i][j][k]+G[i][j+1][k]-theta[i][j][k]);
					
    //分区计算U
	for(i=3;i<=Nx+3;i++)    
		for(j=3;j<=int(0.5/dy)+3;j++)
			for(k=0;k<=3;k++)
				U[i][j][k]=U[i][j][k]-r*(G_[i][j][k]-G_[i][j-1][k]);

	for(i=int(1.0/dx)+3;i<=int(2.0/dx)+3;i++)    
		for(j=int(0.5/dy)+4;j<=Ny+3;j++)
			for(k=0;k<=3;k++)
				U[i][j][k]=U[i][j][k]-r*(G_[i][j][k]-G_[i][j-1][k]);
}
void Roe_Solver(double U[Nx+7][Ny+7][4],double U_[Nx+7][Ny+7][4],double LAMDA_[Nx+7][Ny+7][4][4],
				double L_[Nx+7][Ny+7][4][4],double R_[Nx+7][Ny+7][4][4],double alpha_[Nx+7][Ny+7][4],
				double theta[Nx+7][Ny+7][4],double F[Nx+7][Ny+7][4],double F_[Nx+7][Ny+7][4],
				double G[Nx+7][Ny+7][4],double G_[Nx+7][Ny+7][4],double dx,double dy,
			    double &dt,double a_[Nx+7][Ny+7][1])   //Roe求解器
{
	bound(U,dx,dy);                            
	Roe_x(U,L_,R_,F_,F,LAMDA_,alpha_,theta,dx,dy,dt);
	bound(U,dx,dy); 
	Roe_y(U,L_,R_,G_,G,LAMDA_,alpha_,theta,dx,dy,dt);
	bound(U,dx,dy);
	Roe_y(U,L_,R_,G_,G,LAMDA_,alpha_,theta,dx,dy,dt);
	bound(U,dx,dy);                                       
	Roe_x(U,L_,R_,F_,F,LAMDA_,alpha_,theta,dx,dy,dt);
	bound(U,dx,dy);
}














//-------------------------MUSCL格式------------------
//变量解释
//U_L和U_R的定义见教材3.39,Fp和Fd分别是正负通量.F_是F[i+1/2]通量
void A1_MUSCL_QuickLocation()
{
}
//作用:求Fp和Fd
void U2FpFd_AVG(double U_L[4],double U_R[4],double Fp[4],double Fd[4])
{
	double H,rou,u,v,p,a,lambda1,lambda2,lambda3,lambda4,lambda1p,lambda2p,lambda3p,lambda4p,lambda1d,lambda2d,lambda3d,lambda4d;
	//先算Fp即F+
	rou=U_L[0];
	u=U_L[1]/U_L[0];
	v=U_L[2]/U_L[0];
	p=(GAMA-1)*(U_L[3]-0.5*rou*(u*u+v*v));
	a=pow(GAMA*p/rou,0.5);
	H=a*a/(GAMA-1)+0.5*(u*u+v*v);
	lambda1=u;
	lambda2=u;
	lambda3=u-a;
	lambda4=u+a;
	//计算正特征值
	lambda1p=0.5*(fabs(lambda1)+lambda1);
	lambda2p=0.5*(fabs(lambda2)+lambda2);
	lambda3p=0.5*(fabs(lambda3)+lambda3);
	lambda4p=0.5*(fabs(lambda4)+lambda4);
	//特征值修正
	lambda1p=0.5*(lambda1p+pow(lambda1p*lambda1p+pow(10,-16),0.5));
	lambda2p=0.5*(lambda2p+pow(lambda2p*lambda2p+pow(10,-16),0.5));
	lambda3p=0.5*(lambda3p+pow(lambda3p*lambda3p+pow(10,-16),0.5));
	lambda4p=0.5*(lambda4p+pow(lambda4p*lambda4p+pow(10,-16),0.5));
	Fp[0]=rou/2.0/GAMA*(2*(GAMA-1)*lambda1p+lambda3p+lambda4p);
	Fp[1]=rou/2.0/GAMA*(2*u*(GAMA-1)*lambda1p+(u-a)*lambda3p+(u+a)*lambda4p);
	Fp[2]=rou/2.0/GAMA*(v*2*(GAMA-1)*lambda1p+v*lambda3p+v*lambda4p);
	Fp[3]=rou/2.0/GAMA*(2*((GAMA-1)*H-a*a)*lambda1p+(H-a*u)*lambda3p+(H+a*u)*lambda4p);
	//现在就算Fd即F-
	rou=U_R[0];
	u=U_R[1]/U_R[0];
	v=U_R[2]/U_R[0];
	p=(GAMA-1)*(U_R[3]-0.5*rou*(u*u+v*v));
	a=pow(GAMA*p/rou,0.5);
	H=a*a/(GAMA-1)+0.5*(u*u+v*v);
	lambda1=u;
	lambda2=u;
	lambda3=u-a;
	lambda4=u+a;
	lambda1d=-0.5*(fabs(lambda1)-lambda1);
	lambda2d=-0.5*(fabs(lambda2)-lambda2);
	lambda3d=-0.5*(fabs(lambda3)-lambda3);
	lambda4d=-0.5*(fabs(lambda4)-lambda4);
	lambda1d=0.5*(lambda1d-pow(lambda1d*lambda1d+pow(10,-16),0.5));
	lambda2d=0.5*(lambda2d-pow(lambda2d*lambda2d+pow(10,-16),0.5));
	lambda3d=0.5*(lambda3d-pow(lambda3d*lambda3d+pow(10,-16),0.5));
	lambda4d=0.5*(lambda4d-pow(lambda4d*lambda4d+pow(10,-16),0.5));
	Fd[0]=rou/2.0/GAMA*(2*(GAMA-1)*lambda1d+lambda3d+lambda4d);
	Fd[1]=rou/2.0/GAMA*(2*u*(GAMA-1)*lambda1d+(u-a)*lambda3d+(u+a)*lambda4d);
	Fd[2]=rou/2.0/GAMA*(v*2*(GAMA-1)*lambda1d+v*lambda3d+v*lambda4d);
	Fd[3]=rou/2.0/GAMA*(2*((GAMA-1)*H-a*a)*lambda1d+(H-a*u)*lambda3d+(H+a*u)*lambda4d);
}
void MUSCL_x(double U[Nx+7][Ny+7][4],double U_L[Nx+7][Ny+7][4],double U_R[Nx+7][Ny+7][4],
			 double F_[Nx+7][Ny+7][4],double Fp[Nx+7][Ny+7][4],double Fd[Nx+7][Ny+7][4],double dx,double dy,double dt)
{
	int i,j,k;
	double r=dt/dx;
	//计算U_L和U_R
	for(i=2;i<=Nx+3;i++)
		for(j=2;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
				{
					U_L[i][j][k]=U[i][j][k]+0.25*(1-K)*minmod(U[i][j][k]-U[i-1][j][k],BETA*(U[i+1][j][k]-U[i][j][k]))
						+0.25*(1+K)*minmod(U[i+1][j][k]-U[i][j][k],BETA*(U[i][j][k]-U[i-1][j][k]));
					U_R[i][j][k]=U[i+1][j][k]-0.25*(1-K)*minmod(U[i+2][j][k]-U[i+1][j][k],BETA*(U[i+1][j][k]-U[i][j][k]))
						-0.25*(1+K)*minmod(U[i+1][j][k]-U[i][j][k],BETA*(U[i+2][j][k]-U[i+1][j][k]));
				}
	//计算Fp和Fd
	for(i=2;i<=Nx+3;i++)
		for(j=2;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				U2FpFd_AVG(U_L[i][j],U_R[i][j],Fp[i][j],Fd[i][j]);
	//计算F_
	for(i=2;i<=Nx+3;i++)
		for(j=2;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					F_[i][j][k]=Fp[i][j][k]+Fd[i][j][k];
	//计算U
	for(i=3;i<=Nx+3;i++)
		for(j=3;j<=int(0.5/dy)+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					U[i][j][k]=U[i][j][k]-r*(F_[i][j][k]-F_[i-1][j][k]);

	for(i=int(1.0/dx)+3;i<=int(2.0/dx)+3;i++)       
		for(j=int(0.5/dy)+4;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					U[i][j][k]=U[i][j][k]-r*(F_[i][j][k]-F_[i-1][j][k]);
}
//算Gp和Gd
void U2GpGd_AVG(double U_L[4],double U_R[4],double Gp[4],double Gd[4])
{
	double H,rou,u,v,p,a,lambda1,lambda2,lambda3,lambda4,lambda1p,lambda2p,lambda3p,lambda4p;
	double lambda1d,lambda2d,lambda3d,lambda4d;
	//计算Gp即G+的函数
	rou=U_L[0];
	u=U_L[1]/U_L[0];
	v=U_L[2]/U_L[0];
	p=(GAMA-1)*(U_L[3]-0.5*rou*(u*u+v*v));
	a=pow(GAMA*p/rou,0.5);
	H=a*a/(GAMA-1)+0.5*(u*u+v*v);
	lambda1=v;
	lambda2=v;
	lambda3=v-a;
	lambda4=v+a;
	lambda1p=0.5*(fabs(lambda1)+lambda1);
	lambda2p=0.5*(fabs(lambda2)+lambda2);
	lambda3p=0.5*(fabs(lambda3)+lambda3);
	lambda4p=0.5*(fabs(lambda4)+lambda4);
	lambda1p=0.5*(lambda1p+pow(lambda1p*lambda1p+pow(10,-16),0.5));
	lambda2p=0.5*(lambda2p+pow(lambda2p*lambda2p+pow(10,-16),0.5));
	lambda3p=0.5*(lambda3p+pow(lambda3p*lambda3p+pow(10,-16),0.5));
	lambda4p=0.5*(lambda4p+pow(lambda4p*lambda4p+pow(10,-16),0.5));
	Gp[0]=rou/2.0/GAMA*(2*(GAMA-1)*lambda1p+lambda3p+lambda4p);
	Gp[1]=rou/2.0/GAMA*(u*2*(GAMA-1)*lambda1p+u*lambda3p+u*lambda4p);
	Gp[2]=rou/2.0/GAMA*(2*v*(GAMA-1)*lambda1p+(v-a)*lambda3p+(v+a)*lambda4p);
	Gp[3]=rou/2.0/GAMA*(2*((GAMA-1)*H-a*a)*lambda1p+(H-a*v)*lambda3p+(H+a*v)*lambda4p);
	//计算Gd即G-的函数
	rou=U_R[0];
	u=U_R[1]/U_R[0];
	v=U_R[2]/U_R[0];
	p=(GAMA-1)*(U_R[3]-0.5*rou*(u*u+v*v));
	a=pow(GAMA*p/rou,0.5);
	H=a*a/(GAMA-1)+0.5*(u*u+v*v);
	lambda1=v;
	lambda2=v;
	lambda3=v-a;
	lambda4=v+a;
	lambda1d=-0.5*(fabs(lambda1)-lambda1);
	lambda2d=-0.5*(fabs(lambda2)-lambda2);
	lambda3d=-0.5*(fabs(lambda3)-lambda3);
	lambda4d=-0.5*(fabs(lambda4)-lambda4);
	lambda1d=0.5*(lambda1d-pow(lambda1d*lambda1d+pow(10,-16),0.5));
	lambda2d=0.5*(lambda2d-pow(lambda2d*lambda2d+pow(10,-16),0.5));
	lambda3d=0.5*(lambda3d-pow(lambda3d*lambda3d+pow(10,-16),0.5));
	lambda4d=0.5*(lambda4d-pow(lambda4d*lambda4d+pow(10,-16),0.5));
	Gd[0]=rou/2.0/GAMA*(2*(GAMA-1)*lambda1d+lambda3d+lambda4d);
	Gd[1]=rou/2.0/GAMA*(u*2*(GAMA-1)*lambda1d+u*lambda3d+u*lambda4d);
	Gd[2]=rou/2.0/GAMA*(2*v*(GAMA-1)*lambda1d+(v-a)*lambda3d+(v+a)*lambda4d);
	Gd[3]=rou/2.0/GAMA*(2*((GAMA-1)*H-a*a)*lambda1d+(H-a*v)*lambda3d+(H+a*v)*lambda4d);
}
void MUSCL_y(double U[Nx+7][Ny+7][4],double U_L[Nx+7][Ny+7][4],double U_R[Nx+7][Ny+7][4],
			 double G_[Nx+7][Ny+7][4],double Gp[Nx+7][Ny+7][4],double Gd[Nx+7][Ny+7][4],double dx,double dy,double dt)
{
	int i,j,k;
	double r=dt/dy;
	//计算U_L和U_R
	for(i=2;i<=Nx+3;i++)
		for(j=2;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
				{
					U_L[i][j][k]=U[i][j][k]+0.25*(1-K)*minmod(U[i][j][k]-U[i][j-1][k],BETA*(U[i][j+1][k]-U[i][j][k]))
						+0.25*(1+K)*minmod(U[i][j+1][k]-U[i][j][k],BETA*(U[i][j][k]-U[i][j-1][k]));
					U_R[i][j][k]=U[i][j+1][k]-0.25*(1-K)*minmod(U[i][j+2][k]-U[i][j+1][k],BETA*(U[i][j+1][k]-U[i][j][k]))
						-0.25*(1+K)*minmod(U[i][j+1][k]-U[i][j][k],BETA*(U[i][j+2][k]-U[i][j+1][k]));
				}
	//计算Gp和Gd
	for(i=2;i<=Nx+3;i++)
		for(j=2;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				U2GpGd_AVG(U_L[i][j],U_R[i][j],Gp[i][j],Gd[i][j]);
	//计算G_
	for(i=2;i<=Nx+3;i++)
		for(j=2;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					G_[i][j][k]=Gp[i][j][k]+Gd[i][j][k];
	//计算U
	for(i=3;i<=Nx+3;i++)
		for(j=3;j<=int(0.5/dy)+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					U[i][j][k]=U[i][j][k]-r*(G_[i][j][k]-G_[i][j-1][k]);

	for(i=int(1.0/dx)+3;i<=int(2.0/dx)+3;i++)       
		for(j=int(0.5/dy)+4;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					U[i][j][k]=U[i][j][k]-r*(G_[i][j][k]-G_[i][j-1][k]);
}

void MUSCL_Solver(double U[Nx+7][Ny+7][4],double U_L[Nx+7][Ny+7][4],double U_R[Nx+7][Ny+7][4],double Fp[Nx+7][Ny+7][4],
							  double Fd[Nx+7][Ny+7][4],double F_[Nx+7][Ny+7][4],double Gp[Nx+7][Ny+7][4],double Gd[Nx+7][Ny+7][4],
							  double G_[Nx+7][Ny+7][4],double dx,double dy,double &dt)
{
	bound(U,dx,dy);
	dt=CFL(U,dx,dy,MUSCLCFL);
	MUSCL_x(U,U_L,U_R,F_,Fp,Fd,dx,dy,dt/2.0);
	bound(U,dx,dy);
	dt=CFL(U,dx,dy,MUSCLCFL);
	MUSCL_y(U,U_L,U_R,G_,Gp,Gd,dx,dy,dt/2.0);
	bound(U,dx,dy);
	dt=CFL(U,dx,dy,MUSCLCFL);
	MUSCL_y(U,U_L,U_R,G_,Gp,Gd,dx,dy,dt/2.0);
	bound(U,dx,dy);
	dt=CFL(U,dx,dy,MUSCLCFL);
	MUSCL_x(U,U_L,U_R,F_,Fp,Fd,dx,dy,dt/2.0);
	bound(U,dx,dy);
}











//---------------------------迎风TVD-----------------------------
//变量解释:
//g_,g,gama_,Q_,的定义见教材3.66  theta是翻扩散黏性项在特征空间的叠加之和
void A1_UpWind_TVD_QuickLocation()
{
}
void UpWindTVD_x(double U[Nx+7][Ny+7][4],double U_[Nx+7][Ny+7][4],double LAMDA_[Nx+7][Ny+7][4][4],double L_[Nx+7][Ny+7][4][4],double R_[Nx+7][Ny+7][4][4],
			double alpha_[Nx+7][Ny+7][4],double g_[Nx+7][Ny+7][4],double g[Nx+7][Ny+7][4],double gama_[Nx+7][Ny+7][4],
			double Q_[Nx+7][Ny+7][4],double theta[Nx+7][Ny+7][4],double F[Nx+7][Ny+7][4],double F_[Nx+7][Ny+7][4],
			double dx,double dy,double &dt)       //此函数将之前计算的特征之余特征向量通过TVD算法计算U
{
	int i,j,k,l;
	double r,z;
	RoeAVG_x(U,U_,a_);          //求x方向的Roe平均
	dt=CFL_(U_,a_,dx,dy,UPWINDTVDCFL);             //用Roe平均以后的流场计算dt
	Roe_Eig_x(U,U_,LAMDA_,a_,R_,L_);
	r=dt/dx;
	//计算alpha_
	for(i=0;i<=Nx+5;i++)    
		for(j=0;j<=Ny+5;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					alpha_[i][j][k]=0;

	for(i=0;i<=Nx+5;i++)    
		for(j=0;j<=Ny+5;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					for(l=0;l<=3;l++)  
						alpha_[i][j][k]+=L_[i][j][k][l]*(U[i+1][j][l]-U[i][j][l]);						
	//计算g_
	for(i=0;i<=Nx+5;i++)    
		for(j=0;j<=Ny+5;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
				{
					z=LAMDA_[i][j][k][k];
					if(fabs(z)<DELTA)
					{											
						g_[i][j][k]=alpha_[i][j][k]*0.5*(0.5*(z*z/DELTA+DELTA)+r*z*z);
					}
					else g_[i][j][k]=alpha_[i][j][k]*0.5*(fabs(z)+r*z*z);
				}
	//用g_通过minmod来计算g
	for(i=2;i<=Nx+4;i++)    
		for(j=2;j<=Ny+4;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					g[i][j][k]=minmod(g_[i][j][k],g_[i-1][j][k]);
	//对g进行修正
	for(i=2;i<=Nx+4;i++)    
		for(j=2;j<=Ny+4;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
				{
					if(alpha_[i][j][k]!=0||alpha_[i-1][j][k]!=0)
					{
						g[i][j][k]=(1+OMEGA*fabs(alpha_[i][j][k]-alpha_[i-1][j][k])/(fabs(alpha_[i][j][k])+fabs(alpha_[i-1][j][k])))*g[i][j][k];
					}
					else g[i][j][k]=0;
				}
	//计算gama_
	for(i=2;i<=Nx+3;i++)    
		for(j=2;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
				{
					if(alpha_[i][j][k]!=0)
					{
						gama_[i][j][k]=(g[i+1][j][k]-g[i][j][k])/alpha_[i][j][k];
					}
					else gama_[i][j][k]=0;
				}
	//用gama_和lambda_来计算Q_
	for(i=2;i<=Nx+3;i++)    
		for(j=2;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
				{
					z=LAMDA_[i][j][k][k]+gama_[i][j][k];
					if(fabs(z)<DELTA)
					{											
						Q_[i][j][k]=0.5*(z*z/DELTA+DELTA);
					}
					else Q_[i][j][k]=fabs(z);
				}
	//计算反扩散粘性项theta=(g[i]+g[i+1]-Q(a+gama_)*alpha_)*R
	for(i=2;i<=Nx+3;i++)    
		for(j=2;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					theta[i][j][k]=0;

	for(i=2;i<=Nx+3;i++)    
		for(j=2;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
				{
					for(l=0;l<=3;l++)
					{
						theta[i][j][k]+=R_[i][j][k][l]*(g[i][j][l]+g[i+1][j][l]-Q_[i][j][l]*alpha_[i][j][l]);
					}
				}
	//计算F
	for(i=0;i<=Nx+6;i++)    
		for(j=0;j<=Ny+6;j++)
			if(U[i][j][0]!=0)
				U2F(U[i][j],F[i][j]);
	//计算F_
	for(i=2;i<=Nx+3;i++)    
		for(j=2;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					F_[i][j][k]=0.5*(F[i][j][k]+F[i+1][j][k])+0.5*theta[i][j][k];
    //分区计算U
	for(i=3;i<=Nx+3;i++)    
		for(j=3;j<=int(0.5/dy)+3;j++)	
			for(k=0;k<=3;k++)
				U[i][j][k]=U[i][j][k]-r*(F_[i][j][k]-F_[i-1][j][k]);

	for(i=int(1.0/dx)+3;i<=int(2.0/dx)+3;i++)    
		for(j=int(0.5/dy)+4;j<=Ny+3;j++)
			for(k=0;k<=3;k++)
				U[i][j][k]=U[i][j][k]-r*(F_[i][j][k]-F_[i-1][j][k]);
}
void UpWindTVD_y(double U[Nx+7][Ny+7][4],double U_[Nx+7][Ny+7][4],double LAMDA_[Nx+7][Ny+7][4][4],
				 double L_[Nx+7][Ny+7][4][4],double R_[Nx+7][Ny+7][4][4],
				 double alpha_[Nx+7][Ny+7][4],double g_[Nx+7][Ny+7][4],double g[Nx+7][Ny+7][4],double gama_[Nx+7][Ny+7][4],
				 double Q_[Nx+7][Ny+7][4],double theta[Nx+7][Ny+7][4],double G[Nx+7][Ny+7][4],double G_[Nx+7][Ny+7][4],
				 double dx,double dy,double &dt)
{
	int i,j,k,l;
	double r,z;
	RoeAVG_y(U,U_,a_);          
	dt=CFL_(U_,a_,dx,dy,UPWINDTVDCFL);            
	Roe_Eig_y(U,U_,LAMDA_,a_,R_,L_);
	r=dt/dy;
	//计算alpha_
	for(i=0;i<=Nx+5;i++)    
		for(j=0;j<=Ny+5;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					alpha_[i][j][k]=0;

	for(i=0;i<=Nx+5;i++)    
		for(j=0;j<=Ny+5;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					for(l=0;l<=3;l++)  
						alpha_[i][j][k]+=L_[i][j][k][l]*(U[i][j+1][l]-U[i][j][l]);						
	//计算g_
	for(i=0;i<=Nx+5;i++)    
		for(j=0;j<=Ny+5;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
				{
					z=LAMDA_[i][j][k][k];
					if(fabs(z)<DELTA)
					{											
						g_[i][j][k]=alpha_[i][j][k]*0.5*(0.5*(z*z/DELTA+DELTA)+r*z*z);
					}
					else g_[i][j][k]=alpha_[i][j][k]*0.5*(fabs(z)+r*z*z);
				}
	//计算g
	for(i=2;i<=Nx+4;i++)    
		for(j=2;j<=Ny+4;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					g[i][j][k]=minmod(g_[i][j][k],g_[i][j-1][k]);
	//对g做修正
	for(i=2;i<=Nx+4;i++)    
		for(j=2;j<=Ny+4;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
				{
					if(alpha_[i][j][k]!=0||alpha_[i][j-1][k]!=0)
					{
						g[i][j][k]=(1+OMEGA*fabs(alpha_[i][j][k]-alpha_[i][j-1][k])/(fabs(alpha_[i][j][k])+fabs(alpha_[i][j-1][k])))*g[i][j][k];
					}
					else g[i][j][k]=0;
				}
	//计算gama_
	for(i=2;i<=Nx+3;i++)    
		for(j=2;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
				{
					if(alpha_[i][j][k]!=0)
					{
						gama_[i][j][k]=(g[i][j+1][k]-g[i][j][k])/alpha_[i][j][k];
					}
					else gama_[i][j][k]=0;
				}
	//计算Q_
	for(i=2;i<=Nx+3;i++)    
		for(j=2;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
				{
					z=LAMDA_[i][j][k][k]+gama_[i][j][k];
					if(fabs(z)<DELTA)
					{											
						Q_[i][j][k]=0.5*(z*z/DELTA+DELTA);
					}
					else Q_[i][j][k]=fabs(z);
				}
	//计算粘性项theta
	for(i=2;i<=Nx+3;i++)    
		for(j=2;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					theta[i][j][k]=0;

	for(i=2;i<=Nx+3;i++)    
		for(j=2;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					for(l=0;l<=3;l++)
						theta[i][j][k]+=R_[i][j][k][l]*(g[i][j][l]+g[i][j+1][l]-Q_[i][j][l]*alpha_[i][j][l]);
	//计算G
	for(i=0;i<=Nx+6;i++)    
		for(j=0;j<=Ny+6;j++)
			if(U[i][j][0]!=0)
				U2G(U[i][j],G[i][j]);
	//计算G_
	for(i=2;i<=Nx+3;i++)    
		for(j=2;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					G_[i][j][k]=0.5*(G[i][j][k]+G[i][j+1][k])+0.5*theta[i][j][k];
    //分区计算U
	for(i=3;i<=Nx+3;i++)    
		for(j=3;j<=int(0.5/dy)+3;j++)
			for(k=0;k<=3;k++)
				U[i][j][k]=U[i][j][k]-r*(G_[i][j][k]-G_[i][j-1][k]);

	for(i=int(1.0/dx)+3;i<=int(2.0/dx)+3;i++)    
		for(j=int(0.5/dy)+4;j<=Ny+3;j++)
			for(k=0;k<=3;k++)
				U[i][j][k]=U[i][j][k]-r*(G_[i][j][k]-G_[i][j-1][k]);

}
void UpWindTVD_Solver(double U[Nx+7][Ny+7][4],double U_[Nx+7][Ny+7][4],double LAMDA_[Nx+7][Ny+7][4][4],
					double L_[Nx+7][Ny+7][4][4],double R_[Nx+7][Ny+7][4][4],double alpha_[Nx+7][Ny+7][4],
					 double g_[Nx+7][Ny+7][4],
				  double g[Nx+7][Ny+7][4],double gama_[Nx+7][Ny+7][4],double Q_[Nx+7][Ny+7][4],
				  double theta[Nx+7][Ny+7][4],double F[Nx+7][Ny+7][4],double F_[Nx+7][Ny+7][4],
				  double G[Nx+7][Ny+7][4],double G_[Nx+7][Ny+7][4],double dx,double dy,
				  double &dt,double a_[Nx+7][Ny+7][1])   //TVD求解器
{
	bound(U,dx,dy);                 //先对边界以及虚拟节点赋值
	UpWindTVD_x(U,U_,LAMDA_,L_,R_,alpha_,g_,g,gama_,Q_,theta,F,F_,dx,dy,dt);  //TVD算法求解U
	bound(U,dx,dy);
	UpWindTVD_y(U,U_,LAMDA_,L_,R_,alpha_,g_,g,gama_,Q_,theta,G,G_,dx,dy,dt);
	bound(U,dx,dy);
	UpWindTVD_y(U,U_,LAMDA_,L_,R_,alpha_,g_,g,gama_,Q_,theta,G,G_,dx,dy,dt);
	bound(U,dx,dy);
	UpWindTVD_x(U,U_,LAMDA_,L_,R_,alpha_,g_,g,gama_,Q_,theta,F,F_,dx,dy,dt);
	bound(U,dx,dy);
}










//---------------------------对称TVD----------------------------
//变量解释:
//其中各个变量的定义同迎风TVD
void A1_SymTVD_QuickLocation()
{
}
void SymTVD_x(double U[Nx+7][Ny+7][4],double U_[Nx+7][Ny+7][4],double LAMDA_[Nx+7][Ny+7][4][4],
			  double L_[Nx+7][Ny+7][4][4],double R_[Nx+7][Ny+7][4][4],
			double alpha_[Nx+7][Ny+7][4],double g_[Nx+7][Ny+7][4],double gama_[Nx+7][Ny+7][4],
			double Q_[Nx+7][Ny+7][4],double theta[Nx+7][Ny+7][4],double F[Nx+7][Ny+7][4],double F_[Nx+7][Ny+7][4],
			double dx,double dy,double &dt)       //此函数将之前计算的特征之余特征向量通过TVD算法计算U
{
	int i,j,k,l;
	double r,z;
	RoeAVG_x(U,U_,a_);          //求x方向的Roe平均
	dt=CFL_(U_,a_,dx,dy,SYMTVDCFL);             //用Roe平均以后的流场计算dt
	Roe_Eig_x(U,U_,LAMDA_,a_,R_,L_);
	r=dt/dx;
	//计算alpha_
	for(i=0;i<=Nx+5;i++)    
		for(j=0;j<=Ny+5;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					alpha_[i][j][k]=0;

	for(i=0;i<=Nx+5;i++)    
		for(j=0;j<=Ny+5;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					for(l=0;l<=3;l++)  
						alpha_[i][j][k]+=L_[i][j][k][l]*(U[i+1][j][l]-U[i][j][l]);						
	//计算g_
	for(i=2;i<=Nx+3;i++)    
		for(j=2;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					g_[i][j][k]=minmod3(alpha_[i-1][j][k],alpha_[i][j][k],alpha_[i+1][j][k]);
	//计算Q_
	for(i=2;i<=Nx+3;i++)    
		for(j=2;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
				{
					z=LAMDA_[i][j][k][k]*r;
					if(fabs(z)<DELTA)
					{											
						Q_[i][j][k]=0.5*(z*z/DELTA+DELTA);
					}
					else Q_[i][j][k]=fabs(z);
				}
	//计算反扩散粘性项theta=(g[i]+g[i+1]-Q(a+gama_)*alpha_)*R
	for(i=2;i<=Nx+3;i++)    
		for(j=2;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					theta[i][j][k]=0;

	for(i=2;i<=Nx+3;i++)    
		for(j=2;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					for(l=0;l<=3;l++)
						theta[i][j][k]+=-R_[i][j][k][l]/r*(g_[i][j][l]*r*r*LAMDA_[i][j][l][l]*LAMDA_[i][j][l][l]
							+Q_[i][j][l]*(alpha_[i][j][l]-g_[i][j][l]));

	//计算F
	for(i=2;i<=Nx+4;i++)    
		for(j=2;j<=Ny+4;j++)
			if(U[i][j][0]!=0)
				U2F(U[i][j],F[i][j]);
	//计算F_
	for(i=2;i<=Nx+3;i++)    
		for(j=2;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					F_[i][j][k]=0.5*(F[i][j][k]+F[i+1][j][k])+0.5*theta[i][j][k];
    //分区计算U
	for(i=3;i<=Nx+3;i++)    
		for(j=3;j<=int(0.5/dy)+3;j++)	
			for(k=0;k<=3;k++)
				U[i][j][k]=U[i][j][k]-r*(F_[i][j][k]-F_[i-1][j][k]);

	for(i=int(1.0/dx)+3;i<=int(2.0/dx)+3;i++)    
		for(j=int(0.5/dy)+4;j<=Ny+3;j++)
			for(k=0;k<=3;k++)
				U[i][j][k]=U[i][j][k]-r*(F_[i][j][k]-F_[i-1][j][k]);
}
void SymTVD_y(double U[Nx+7][Ny+7][4],double U_[Nx+7][Ny+7][4],
			  double LAMDA_[Nx+7][Ny+7][4][4],double L_[Nx+7][Ny+7][4][4],double R_[Nx+7][Ny+7][4][4],
			double alpha_[Nx+7][Ny+7][4],double g_[Nx+7][Ny+7][4],double gama_[Nx+7][Ny+7][4],
			double Q_[Nx+7][Ny+7][4],double theta[Nx+7][Ny+7][4],double G[Nx+7][Ny+7][4],double G_[Nx+7][Ny+7][4],
			double dx,double dy,double &dt)       //此函数将之前计算的特征之余特征向量通过TVD算法计算U
{
	int i,j,k,l;
	double r,z;
	RoeAVG_y(U,U_,a_);          //求x方向的Roe平均
	dt=CFL_(U_,a_,dx,dy,SYMTVDCFL);             //用Roe平均以后的流场计算dt
	Roe_Eig_y(U,U_,LAMDA_,a_,R_,L_);
	r=dt/dy;
	//计算alpha_
	for(i=0;i<=Nx+5;i++)    
		for(j=0;j<=Ny+5;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					alpha_[i][j][k]=0;

	for(i=0;i<=Nx+5;i++)    
		for(j=0;j<=Ny+5;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					for(l=0;l<=3;l++)  
						alpha_[i][j][k]+=L_[i][j][k][l]*(U[i][j+1][l]-U[i][j][l]);						
	//计算g_
	for(i=2;i<=Nx+3;i++)    
		for(j=2;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					g_[i][j][k]=minmod3(alpha_[i][j-1][k],alpha_[i][j][k],alpha_[i][j+1][k]);
	//计算Q_
	for(i=2;i<=Nx+3;i++)    
		for(j=2;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
				{
					z=LAMDA_[i][j][k][k]*r;
					if(fabs(z)<DELTA)
					{											
						Q_[i][j][k]=0.5*(z*z/DELTA+DELTA);
					}
					else Q_[i][j][k]=fabs(z);
				}
	//计算反扩散粘性项theta=(g[i]+g[i+1]-Q(a+gama_)*alpha_)*R
	for(i=2;i<=Nx+3;i++)    
		for(j=2;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					theta[i][j][k]=0;

	for(i=2;i<=Nx+3;i++)    
		for(j=2;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					for(l=0;l<=3;l++)
						theta[i][j][k]+=-R_[i][j][k][l]/r*(g_[i][j][l]*r*r*LAMDA_[i][j][l][l]*LAMDA_[i][j][l][l]
							+Q_[i][j][l]*(alpha_[i][j][l]-g_[i][j][l]));
	//计算G
	for(i=0;i<=Nx+6;i++)    
		for(j=0;j<=Ny+6;j++)
			if(U[i][j][0]!=0)
				U2G(U[i][j],G[i][j]);
	//计算G_
	for(i=2;i<=Nx+3;i++)    
		for(j=2;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					G_[i][j][k]=0.5*(G[i][j][k]+G[i][j+1][k])+0.5*theta[i][j][k];
    //分区计算U
	for(i=3;i<=Nx+3;i++)    
		for(j=3;j<=int(0.5/dy)+3;j++)
			for(k=0;k<=3;k++)
				U[i][j][k]=U[i][j][k]-r*(G_[i][j][k]-G_[i][j-1][k]);

	for(i=int(1.0/dx)+3;i<=int(2.0/dx)+3;i++)    
		for(j=int(0.5/dy)+4;j<=Ny+3;j++)
			for(k=0;k<=3;k++)
				U[i][j][k]=U[i][j][k]-r*(G_[i][j][k]-G_[i][j-1][k]);
}
void SymTVD_Solver(double U[Nx+7][Ny+7][4],double U_[Nx+7][Ny+7][4],double LAMDA_[Nx+7][Ny+7][4][4],
				  double L_[Nx+7][Ny+7][4][4],double R_[Nx+7][Ny+7][4][4],double alpha_[Nx+7][Ny+7][4],
				  double g_[Nx+7][Ny+7][4],
				  double Q_[Nx+7][Ny+7][4],
				  double theta[Nx+7][Ny+7][4],double F[Nx+7][Ny+7][4],double F_[Nx+7][Ny+7][4],
				  double G[Nx+7][Ny+7][4],double G_[Nx+7][Ny+7][4],double dx,double dy,
				  double &dt,double a_[Nx+7][Ny+7][1])   //TVD求解器
{
	bound(U,dx,dy);             
	SymTVD_x(U,U_,LAMDA_,L_,R_,alpha_,g_,gama_,Q_,theta,F,F_,dx,dy,dt);
	bound(U,dx,dy);
	SymTVD_y(U,U_,LAMDA_,L_,R_,alpha_,g_,gama_,Q_,theta,G,G_,dx,dy,dt);
	bound(U,dx,dy);
	SymTVD_y(U,U_,LAMDA_,L_,R_,alpha_,g_,gama_,Q_,theta,G,G_,dx,dy,dt);
	bound(U,dx,dy);               
	SymTVD_x(U,U_,LAMDA_,L_,R_,alpha_,g_,gama_,Q_,theta,F,F_,dx,dy,dt);
}













//----------------------------NND-----------------------
//变量定义:
//Fp=Fplus=F+ Fd=F-
void A1_NND_QuickLocation()
{
}

void U2FpFd(double U[4],double Fp[4],double Fd[4])
{
	double H,rou,u,v,p,a,lambda1,lambda2,lambda3,lambda4,lambda1p,lambda2p,lambda3p,lambda4p,lambda1d,lambda2d,lambda3d,lambda4d;
	//先算Fp即F+
	rou=U[0];
	u=U[1]/U[0];
	v=U[2]/U[0];
	p=(GAMA-1)*(U[3]-0.5*rou*(u*u+v*v));
	a=pow(GAMA*p/rou,0.5);
	H=a*a/(GAMA-1)+0.5*(u*u+v*v);
	lambda1=u;
	lambda2=u;
	lambda3=u-a;
	lambda4=u+a;
	//计算正特征值
	lambda1p=0.5*(fabs(lambda1)+lambda1);
	lambda2p=0.5*(fabs(lambda2)+lambda2);
	lambda3p=0.5*(fabs(lambda3)+lambda3);
	lambda4p=0.5*(fabs(lambda4)+lambda4);
	//特征值修正
	lambda1p=0.5*(lambda1p+pow(lambda1p*lambda1p+pow(10,-16),0.5));
	lambda2p=0.5*(lambda2p+pow(lambda2p*lambda2p+pow(10,-16),0.5));
	lambda3p=0.5*(lambda3p+pow(lambda3p*lambda3p+pow(10,-16),0.5));
	lambda4p=0.5*(lambda4p+pow(lambda4p*lambda4p+pow(10,-16),0.5));
	Fp[0]=rou/2/GAMA*(2*(GAMA-1)*lambda1p+lambda3p+lambda4p);
	Fp[1]=rou/2/GAMA*(2*u*(GAMA-1)*lambda1p+(u-a)*lambda3p+(u+a)*lambda4p);
	Fp[2]=rou/2/GAMA*(v*2*(GAMA-1)*lambda1p+v*lambda3p+v*lambda4p);
	Fp[3]=rou/2/GAMA*(2*((GAMA-1)*H-a*a)*lambda1p+(H-a*u)*lambda3p+(H+a*u)*lambda4p);
	//现在就算Fd即F-
	lambda1d=-0.5*(fabs(lambda1)-lambda1);
	lambda2d=-0.5*(fabs(lambda2)-lambda2);
	lambda3d=-0.5*(fabs(lambda3)-lambda3);
	lambda4d=-0.5*(fabs(lambda4)-lambda4);
	lambda1d=0.5*(lambda1d-pow(lambda1d*lambda1d+pow(10,-16),0.5));
	lambda2d=0.5*(lambda2d-pow(lambda2d*lambda2d+pow(10,-16),0.5));
	lambda3d=0.5*(lambda3d-pow(lambda3d*lambda3d+pow(10,-16),0.5));
	lambda4d=0.5*(lambda4d-pow(lambda4d*lambda4d+pow(10,-16),0.5));
	Fd[0]=rou/2/GAMA*(2*(GAMA-1)*lambda1d+lambda3d+lambda4d);
	Fd[1]=rou/2/GAMA*(2*u*(GAMA-1)*lambda1d+(u-a)*lambda3d+(u+a)*lambda4d);
	Fd[2]=rou/2/GAMA*(v*2*(GAMA-1)*lambda1d+v*lambda3d+v*lambda4d);
	Fd[3]=rou/2/GAMA*(2*((GAMA-1)*H-a*a)*lambda1d+(H-a*u)*lambda3d+(H+a*u)*lambda4d);
}
void U2GpGd(double U[4],double Gp[4],double Gd[4])
{
	double H,rou,u,v,p,a,lambda1,lambda2,lambda3,lambda4,lambda1p,lambda2p,lambda3p,lambda4p,lambda1d,lambda2d,lambda3d,lambda4d;
	//计算Gp即G+的函数
	rou=U[0];
	u=U[1]/U[0];
	v=U[2]/U[0];
	p=(GAMA-1)*(U[3]-0.5*rou*(u*u+v*v));
	a=pow(GAMA*p/rou,0.5);
	H=a*a/(GAMA-1)+0.5*(u*u+v*v);
	lambda1=v;
	lambda2=v;
	lambda3=v-a;
	lambda4=v+a;
	lambda1p=0.5*(fabs(lambda1)+lambda1);
	lambda2p=0.5*(fabs(lambda2)+lambda2);
	lambda3p=0.5*(fabs(lambda3)+lambda3);
	lambda4p=0.5*(fabs(lambda4)+lambda4);
	lambda1p=0.5*(lambda1p+pow(lambda1p*lambda1p+pow(10,-16),0.5));
	lambda2p=0.5*(lambda2p+pow(lambda2p*lambda2p+pow(10,-16),0.5));
	lambda3p=0.5*(lambda3p+pow(lambda3p*lambda3p+pow(10,-16),0.5));
	lambda4p=0.5*(lambda4p+pow(lambda4p*lambda4p+pow(10,-16),0.5));
	Gp[0]=rou/2/GAMA*(2*(GAMA-1)*lambda1p+lambda3p+lambda4p);
	Gp[1]=rou/2/GAMA*(u*2*(GAMA-1)*lambda1p+u*lambda3p+u*lambda4p);
	Gp[2]=rou/2/GAMA*(2*v*(GAMA-1)*lambda1p+(v-a)*lambda3p+(v+a)*lambda4p);
	Gp[3]=rou/2/GAMA*(2*((GAMA-1)*H-a*a)*lambda1p+(H-a*v)*lambda3p+(H+a*v)*lambda4p);
	//计算Gd即G-的函数
	lambda1d=-0.5*(fabs(lambda1)-lambda1);
	lambda2d=-0.5*(fabs(lambda2)-lambda2);
	lambda3d=-0.5*(fabs(lambda3)-lambda3);
	lambda4d=-0.5*(fabs(lambda4)-lambda4);
	lambda1d=0.5*(lambda1d-pow(lambda1d*lambda1d+pow(10,-16),0.5));
	lambda2d=0.5*(lambda2d-pow(lambda2d*lambda2d+pow(10,-16),0.5));
	lambda3d=0.5*(lambda3d-pow(lambda3d*lambda3d+pow(10,-16),0.5));
	lambda4d=0.5*(lambda4d-pow(lambda4d*lambda4d+pow(10,-16),0.5));
	Gd[0]=rou/2/GAMA*(2*(GAMA-1)*lambda1d+lambda3d+lambda4d);
	Gd[1]=rou/2/GAMA*(u*2*(GAMA-1)*lambda1d+u*lambda3d+u*lambda4d);
	Gd[2]=rou/2/GAMA*(2*v*(GAMA-1)*lambda1d+(v-a)*lambda3d+(v+a)*lambda4d);
	Gd[3]=rou/2/GAMA*(2*((GAMA-1)*H-a*a)*lambda1d+(H-a*v)*lambda3d+(H+a*v)*lambda4d);
}

void NND_x(double U[Nx+7][Ny+7][4],double Fp[Nx+7][Ny+7][4],double Fd[Nx+7][Ny+7][4],double dx,double dy,double &dt)
{
	int i,j,k;
	double r;
	dt=CFL(U,dx,dy,NNDCFL);
	r=dt/dx;
	for(i=0;i<=Nx+6;i++)
		for(j=0;j<=Ny+6;j++)
			if(U[i][j][0]!=0)
				U2FpFd(U[i][j],Fp[i][j],Fd[i][j]);
	//分区计算U
	for(i=3;i<=Nx+3;i++)
		for(j=3;j<=int(0.5/dy)+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
				{
					U[i][j][k]=U[i][j][k]-r*(Fp[i][j][k]+0.5*minmod(Fp[i][j][k]-Fp[i-1][j][k],Fp[i+1][j][k]-Fp[i][j][k])
						+Fd[i+1][j][k]-0.5*minmod(Fd[i+1][j][k]-Fd[i][j][k],Fd[i+2][j][k]-Fd[i+1][j][k])
						-Fp[i-1][j][k]-0.5*minmod(Fp[i-1][j][k]-Fp[i-2][j][k],Fp[i][j][k]-Fp[i-1][j][k])
						-Fd[i][j][k]+0.5*minmod(Fd[i][j][k]-Fd[i-1][j][k],Fd[i+1][j][k]-Fd[i][j][k]));
				}

	for(i=int(1.0/dx)+3;i<=int(2.0/dx)+3;i++)
		for(j=int(0.5/dy)+4;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
				{
					U[i][j][k]=U[i][j][k]-r*(Fp[i][j][k]+0.5*minmod(Fp[i][j][k]-Fp[i-1][j][k],Fp[i+1][j][k]-Fp[i][j][k])
						+Fd[i+1][j][k]-0.5*minmod(Fd[i+1][j][k]-Fd[i][j][k],Fd[i+2][j][k]-Fd[i+1][j][k])
						-Fp[i-1][j][k]-0.5*minmod(Fp[i-1][j][k]-Fp[i-2][j][k],Fp[i][j][k]-Fp[i-1][j][k])
						-Fd[i][j][k]+0.5*minmod(Fd[i][j][k]-Fd[i-1][j][k],Fd[i+1][j][k]-Fd[i][j][k]));
				}
}
void NND_y(double U[Nx+7][Ny+7][4],double Gp[Nx+7][Ny+7][4],double Gd[Nx+7][Ny+7][4],double dx,double dy,double &dt)
{
	int i,j,k;
	double r;
	dt=CFL(U,dx,dy,NNDCFL);
	r=dt/dy;
	for(i=0;i<=Nx+6;i++)
		for(j=0;j<=Ny+6;j++)
			if(U[i][j][0]!=0)
				U2GpGd(U[i][j],Gp[i][j],Gd[i][j]);
	//计算U
	for(i=3;i<=Nx+3;i++)
		for(j=3;j<=int(0.5/dy)+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
				{
					U[i][j][k]=U[i][j][k]-r*(Gp[i][j][k]+0.5*minmod(Gp[i][j][k]-Gp[i][j-1][k],Gp[i][j+1][k]-Gp[i][j][k])
						+Gd[i][j+1][k]-0.5*minmod(Gd[i][j+1][k]-Gd[i][j][k],Gd[i][j+2][k]-Gd[i][j+1][k])
						-Gp[i][j-1][k]-0.5*minmod(Gp[i][j-1][k]-Gp[i][j-2][k],Gp[i][j][k]-Gp[i][j-1][k])
						-Gd[i][j][k]+0.5*minmod(Gd[i][j][k]-Gd[i][j-1][k],Gd[i][j+1][k]-Gd[i][j][k]));
				}

	for(i=int(1.0/dx)+3;i<=int(2.0/dx)+3;i++)
		for(j=int(0.5/dy)+4;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
				{
					U[i][j][k]=U[i][j][k]-r*(Gp[i][j][k]+0.5*minmod(Gp[i][j][k]-Gp[i][j-1][k],Gp[i][j+1][k]-Gp[i][j][k])
						+Gd[i][j+1][k]-0.5*minmod(Gd[i][j+1][k]-Gd[i][j][k],Gd[i][j+2][k]-Gd[i][j+1][k])
						-Gp[i][j-1][k]-0.5*minmod(Gp[i][j-1][k]-Gp[i][j-2][k],Gp[i][j][k]-Gp[i][j-1][k])
						-Gd[i][j][k]+0.5*minmod(Gd[i][j][k]-Gd[i][j-1][k],Gd[i][j+1][k]-Gd[i][j][k]));
				}
}

void NND_Solver(double U[Nx+7][Ny+7][4],double Fp[Nx+7][Ny+7][4],double Fd[Nx+7][Ny+7][4],
				double Gp[Nx+7][Ny+7][4],double Gd[Nx+7][Ny+7][4],double dx,double dy,double &dt)
{
	bound(U,dx,dy);
	NND_x(U,Fp,Fd,dx,dy,dt);
	bound(U,dx,dy);
	NND_y(U,Gp,Gd,dx,dy,dt);
	bound(U,dx,dy);
	NND_y(U,Gp,Gd,dx,dy,dt);
	bound(U,dx,dy);
	NND_x(U,Fp,Fd,dx,dy,dt);
}












//-----------------------------ENO---------------------
//这里使用的是Lax-Friedrichs数值通量分裂和三阶Runge-Kutta时间离散
//变量解释:
//F_p表示F+[i+1/2],q3p表示F+[i+1/2]的三点模板差值,同理q3d表示F-的三点模板差值
//Gp=G+  U1和U2是Runge-Kutta中用到的中间变量
void A1_ENO_QuickLocation()
{
}
void LF_x(double U[Nx+7][Ny+7][4],double LAMDA_[Nx+7][Ny+7][4][4],double F[Nx+7][Ny+7][4],double Fp[Nx+7][Ny+7][4],
		   double Fd[Nx+7][Ny+7][4])         //利用Roe平均以后的值计算x方向的特征值与特征向量
{
	int i,j,k;
	double rou,u,v,p,a,maxlamda=pow(10,-100);
	for(i=0;i<=Nx+6;i++)        //对LAMDA_赋值
		for(j=0;j<=Ny+6;j++)
			if(U[i][j][0]!=0)
			{
				rou=U[i][j][0];
				u=U[i][j][1]/U[i][j][0];
				v=U[i][j][2]/U[i][j][0];
				p=(GAMA-1)*(U[i][j][3]-0.5*rou*(u*u+v*v));
				a=pow(p*GAMA/rou,0.5);
				LAMDA_[i][j][0][0]=u-a;
				LAMDA_[i][j][1][1]=u;
				LAMDA_[i][j][2][2]=u+a;
				LAMDA_[i][j][3][3]=u;
				for(k=0;k<=3;k++)
				{					
					if(fabs(LAMDA_[i][j][k][k])>=maxlamda)maxlamda=fabs(LAMDA_[i][j][k][k]);
				}
			}
	//U2F
	for(i=0;i<=Nx+6;i++)      
		for(j=0;j<=Ny+6;j++)
			if(U[i][j][0]!=0)
				U2F(U[i][j],F[i][j]);

	//计算Fpd
	for(i=0;i<=Nx+6;i++)      
		for(j=0;j<=Ny+6;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
				{
					Fp[i][j][k]=0.5*(F[i][j][k]+maxlamda*U[i][j][k]);
					Fd[i][j][k]=0.5*(F[i][j][k]-maxlamda*U[i][j][k]);
				}

}
void ENO_x(double U[Nx+7][Ny+7][4],double F[Nx+7][Ny+7][4],double Fp[Nx+7][Ny+7][4],double Fd[Nx+7][Ny+7][4],
		   double F_p[Nx+7][Ny+7][4],double F_d[Nx+7][Ny+7][4],double F_[Nx+7][Ny+7][4],
		   double q3p[Nx+7][Ny+7][4][3],double q3d[Nx+7][Ny+7][4][3],double dx,double dy,double dt)       //此函数将之前计算的特征之余特征向量通过TVD算法计算U
{
	int i,j,k;
	double r;
	dt=CFL(U,dx,dy,ENOCFL);
	r=dt/dx;
	//计算q3pd q5pd
	for(i=2;i<=Nx+3;i++)    
		for(j=2;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
				{
					q3p[i][j][k][0]=1.0/3.0*Fp[i-2][j][k]-7.0/6.0*Fp[i-1][j][k]+11.0/6.0*Fp[i][j][k];
					q3p[i][j][k][1]=-1.0/6.0*Fp[i-1][j][k]+5.0/6.0*Fp[i][j][k]+1.0/3.0*Fp[i+1][j][k];
					q3p[i][j][k][2]=1.0/3.0*Fp[i][j][k]+5.0/6.0*Fp[i+1][j][k]-1.0/6.0*Fp[i+2][j][k];

					q3d[i][j][k][0]=-1.0/6.0*Fd[i-1][j][k]+5.0/6.0*Fd[i][j][k]+1.0/3.0*Fd[i+1][j][k];
					q3d[i][j][k][1]=1.0/3.0*Fd[i][j][k]+5.0/6.0*Fd[i+1][j][k]-1.0/6.0*Fd[i+2][j][k];
					q3d[i][j][k][2]=11.0/6.0*Fd[i+1][j][k]-7.0/6.0*Fd[i+2][j][k]+1.0/3.0*Fd[i+3][j][k];

				}

		//判断差商大小并且赋值
	for(i=2;i<=Nx+3;i++)    
		for(j=2;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
				{
					if(fabs(U[i+1][j][k]-U[i][j][k])>fabs(U[i][j][k]-U[i-1][j][k])&&fabs(U[i+1][j][k]-U[i][j][k])>fabs(U[i-1][j][k]-U[i-2][j][k]))
						F_p[i][j][k]=q3p[i][j][k][0];
					else 
					{
						if(fabs(U[i][j][k]-U[i-1][j][k])>fabs(U[i+1][j][k]-U[i][j][k])&&fabs(U[i][j][k]-U[i-1][j][k])>fabs(U[i+2][j][k]-U[i+1][j][k]))
							F_p[i][j][k]=q3p[i][j][k][2];
						else
							F_p[i][j][k]=q3p[i][j][k][1];
					}

				}

	for(i=2;i<=Nx+3;i++)    
		for(j=2;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
				{
					if(fabs(U[i+2][j][k]-U[i+1][j][k])>fabs(U[i+1][j][k]-U[i][j][k])&&fabs(U[i+2][j][k]-U[i+1][j][k])>fabs(U[i][j][k]-U[i-1][j][k]))
						F_d[i][j][k]=q3d[i][j][k][0];
					else
					{
						if(fabs(U[i+1][j][k]-U[i][j][k])>fabs(U[i+2][j][k]-U[i+1][j][k])&&fabs(U[i+1][j][k]-U[i][j][k])>fabs(U[i+3][j][k]-U[i+2][j][k]))
							F_d[i][j][k]=q3d[i][j][k][2];
						else
							F_d[i][j][k]=q3d[i][j][k][1];
					}
				}

	//计算F_
	for(i=2;i<=Nx+3;i++)    
		for(j=2;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					F_[i][j][k]=F_p[i][j][k]+F_d[i][j][k];
	//分区计算U
	for(i=3;i<=Nx+3;i++)    
		for(j=3;j<=int(0.5/dy)+3;j++)	
			for(k=0;k<=3;k++)
				U[i][j][k]=U[i][j][k]-r*(F_[i][j][k]-F_[i-1][j][k]);

	for(i=int(1.0/dx)+3;i<=int(2.0/dx)+3;i++)    
		for(j=int(0.5/dy)+4;j<=Ny+3;j++)
			for(k=0;k<=3;k++)
				U[i][j][k]=U[i][j][k]-r*(F_[i][j][k]-F_[i-1][j][k]);
	
}

///////////-------------------------------y
void LF_y(double U[Nx+7][Ny+7][4],double LAMDA_[Nx+7][Ny+7][4][4],double G[Nx+7][Ny+7][4],double Gp[Nx+7][Ny+7][4],
		   double Gd[Nx+7][Ny+7][4])        //利用Roe平均以后的值计算x方向的特征值与特征向量
{
	int i,j,k;
	double rou,u,v,a,p,maxlamda=pow(10,-100);
	for(i=0;i<=Nx+5;i++)        //对LAMDA_赋值
		for(j=0;j<=Ny+5;j++)
			if(U[i][j][0]!=0)
			{
				rou=U[i][j][0];
				u=U[i][j][1]/U[i][j][0];
				v=U[i][j][2]/U[i][j][0];
				p=(GAMA-1)*(U[i][j][3]-0.5*rou*(u*u+v*v));
				a=pow(p*GAMA/rou,0.5);
				LAMDA_[i][j][0][0]=v-a;
				LAMDA_[i][j][1][1]=v;
				LAMDA_[i][j][2][2]=v+a;
				LAMDA_[i][j][3][3]=v;
				for(k=0;k<=3;k++)
				{					
					if(fabs(LAMDA_[i][j][k][k])>=maxlamda)maxlamda=fabs(LAMDA_[i][j][k][k]);
				}
			}

		//U2G
	for(i=0;i<=Nx+6;i++)      
		for(j=0;j<=Ny+6;j++)
			if(U[i][j][0]!=0)
				U2G(U[i][j],G[i][j]);

	//计算Gpd
	for(i=0;i<=Nx+6;i++)      
		for(j=0;j<=Ny+6;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
				{
					Gp[i][j][k]=0.5*(G[i][j][k]+maxlamda*U[i][j][k]);
					Gd[i][j][k]=0.5*(G[i][j][k]-maxlamda*U[i][j][k]);
				}

}
void ENO_y(double U[Nx+7][Ny+7][4],double G[Nx+7][Ny+7][4],double Gp[Nx+7][Ny+7][4],double Gd[Nx+7][Ny+7][4],
		   double G_p[Nx+7][Ny+7][4],double G_d[Nx+7][Ny+7][4],double G_[Nx+7][Ny+7][4],
		   double q3p[Nx+7][Ny+7][4][3],double q3d[Nx+7][Ny+7][4][3],double dx,double dy,double dt)    //此函数将之前计算的特征之余特征向量通过TVD算法计算U
{
	int i,j,k;
	double r;
	dt=CFL(U,dx,dy,ENOCFL);
	r=dt/dy;
	//计算q3pd q5pd
	for(i=2;i<=Nx+3;i++)    
		for(j=2;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
			{
				for(k=0;k<=3;k++)
				{
					q3p[i][j][k][0]=1.0/3.0*Gp[i][j-2][k]-7.0/6.0*Gp[i][j-1][k]+11.0/6.0*Gp[i][j][k];
					q3p[i][j][k][1]=-1.0/6.0*Gp[i][j-1][k]+5.0/6.0*Gp[i][j][k]+1.0/3.0*Gp[i][j+1][k];
					q3p[i][j][k][2]=1.0/3.0*Gp[i][j][k]+5.0/6.0*Gp[i][j+1][k]-1.0/6.0*Gp[i][j+2][k];

					q3d[i][j][k][0]=-1.0/6.0*Gd[i][j-1][k]+5.0/6.0*Gd[i][j][k]+1.0/3.0*Gd[i][j+1][k];
					q3d[i][j][k][1]=1.0/3.0*Gd[i][j][k]+5.0/6.0*Gd[i][j+1][k]-1.0/6.0*Gd[i][j+2][k];
					q3d[i][j][k][2]=11.0/6.0*Gd[i][j+1][k]-7.0/6.0*Gd[i][j+2][k]+1.0/3.0*Gd[i][j+3][k];
			
				}
			}

	//判断差商大小并且赋值
	for(i=2;i<=Nx+3;i++)    
		for(j=2;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
				{
					if(fabs(U[i][j+1][k]-U[i][j][k])>fabs(U[i][j][k]-U[i][j-1][k])&&fabs(U[i][j+1][k]-U[i][j][k])>fabs(U[i][j-1][k]-U[i][j-2][k]))
						G_p[i][j][k]=q3p[i][j][k][0];
					else 
					{
						if(fabs(U[i][j][k]-U[i][j-1][k])>fabs(U[i][j+1][k]-U[i][j][k])&&fabs(U[i][j][k]-U[i][j-1][k])>fabs(U[i][j+2][k]-U[i][j+1][k]))
							G_p[i][j][k]=q3p[i][j][k][2];
						else
							G_p[i][j][k]=q3p[i][j][k][1];
					}
				}

	for(i=2;i<=Nx+3;i++)    
		for(j=2;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
				{
					if(fabs(U[i][j+2][k]-U[i][j+1][k])>fabs(U[i][j+1][k]-U[i][j][k])&&fabs(U[i][j+2][k]-U[i][j+1][k])>fabs(U[i][j][k]-U[i][j-1][k]))
						G_d[i][j][k]=q3d[i][j][k][0];
					else
					{
						if(fabs(U[i][j+1][k]-U[i][j][k])>fabs(U[i][j+2][k]-U[i][j+1][k])&&fabs(U[i][j+1][k]-U[i][j][k])>fabs(U[i][j+3][k]-U[i][j+2][k]))
							G_d[i][j][k]=q3d[i][j][k][2];
						else
							G_d[i][j][k]=q3d[i][j][k][1];
					}

				}


	//计算G_
	for(i=2;i<=Nx+3;i++)    
		for(j=2;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					G_[i][j][k]=G_p[i][j][k]+G_d[i][j][k];

	//分区计算U
	for(i=3;i<=Nx+3;i++)    
		for(j=3;j<=int(0.5/dy)+3;j++)	
			for(k=0;k<=3;k++)
				U[i][j][k]=U[i][j][k]-r*(G_[i][j][k]-G_[i][j-1][k]);

	for(i=int(1.0/dx)+3;i<=int(2.0/dx)+3;i++)    
		for(j=int(0.5/dy)+4;j<=Ny+3;j++)
			for(k=0;k<=3;k++)
				U[i][j][k]=U[i][j][k]-r*(G_[i][j][k]-G_[i][j-1][k]);
	
}

void ENO_Solver(double U[Nx+7][Ny+7][4],double U1[Nx+7][Ny+7][4],double U2[Nx+7][Ny+7][4],
				double F[Nx+7][Ny+7][4],double Fp[Nx+7][Ny+7][4],double Fd[Nx+7][Ny+7][4],
		   double F_p[Nx+7][Ny+7][4],double F_d[Nx+7][Ny+7][4],double F_[Nx+7][Ny+7][4],
		   double G[Nx+7][Ny+7][4],double Gp[Nx+7][Ny+7][4],double Gd[Nx+7][Ny+7][4],
		   double G_p[Nx+7][Ny+7][4],double G_d[Nx+7][Ny+7][4],double G_[Nx+7][Ny+7][4],
		   double LAMDA_[Nx+7][Ny+7][4][4],
		   double q3p[Nx+7][Ny+7][4][3],double q3d[Nx+7][Ny+7][4][3],double dx,double dy,double &dt)
{
	bound(U,dx,dy);
	LF_x(U,LAMDA_,F,Fp,Fd);
	ENO_x(U,F,Fp,Fd,F_p,F_d,F_,q3p,q3d,dx,dy,dt);
	bound(U,dx,dy);
	LF_y(U,LAMDA_,G,Gp,Gd);
	ENO_y(U,G,Gp,Gd,G_p,G_d,G_,q3p,q3d,dx,dy,dt);
	bound(U,dx,dy);
	LF_y(U,LAMDA_,G,Gp,Gd);
	ENO_y(U,G,Gp,Gd,G_p,G_d,G_,q3p,q3d,dx,dy,dt);
	bound(U,dx,dy);
	LF_x(U,LAMDA_,F,Fp,Fd);
	ENO_x(U,F,Fp,Fd,F_p,F_d,F_,q3p,q3d,dx,dy,dt);
	bound(U,dx,dy);
}












//------------------------------------WENO----------------------
//同ENO这里WENO用的也是Lax-Friedrichs数值通量分裂
//变量解释:
//ISp=IS+ 其中IS,alpha,omega的定义见教材5.13
void A1_WENO_QuickLocation()
{
}

void WENO_x(double U[Nx+7][Ny+7][4],double ISp[Nx+7][Ny+7][4][3],double ISd[Nx+7][Ny+7][4][3],double omegap[Nx+7][Ny+7][4][3],
			double omegad[Nx+7][Ny+7][4][3],double alphap[Nx+7][Ny+7][4][3],double alphad[Nx+7][Ny+7][4][3],
			double q3p[Nx+7][Ny+7][4][3],double q3d[Nx+7][Ny+7][4][3],
			double Fp[Nx+7][Ny+7][4],double Fd[Nx+7][Ny+7][4],
			double F_p[Nx+7][Ny+7][4],double F_d[Nx+7][Ny+7][4],double F_[Nx+7][Ny+7][4],
			double dx,double dy,double dt)       //此函数将之前计算的特征之余特征向量通过TVD算法计算U
{
	int i,j,k,l;
	dt=CFL(U,dx,dy,WENOCFL);
	double r=dt/dx;
	//计算ISpd
	for(i=2;i<=Nx+3;i++)    
		for(j=2;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
				{
					ISp[i][j][k][0]=13.0/12.0*pow(Fp[i-2][j][k]-2*Fp[i-1][j][k]+Fp[i][j][k],2)
						+0.25*pow(Fp[i-2][j][k]-4*Fp[i-1][j][k]+3*Fp[i][j][k],2);
					ISp[i][j][k][1]=13.0/12.0*pow(Fp[i-1][j][k]-2*Fp[i][j][k]+Fp[i+1][j][k],2)
						+0.25*pow(Fp[i-1][j][k]-Fp[i+1][j][k],2);
					ISp[i][j][k][2]=13.0/12.0*pow(Fp[i][j][k]-2*Fp[i+1][j][k]+Fp[i+2][j][k],2)
						+0.25*pow(3*Fp[i][j][k]-4*Fp[i+1][j][k]+Fp[i+2][j][k],2);

					ISd[i][j][k][0]=13.0/12.0*pow(Fd[i+1][j][k]-2*Fd[i][j][k]+Fd[i-1][j][k],2)
						+0.25*pow(3*Fd[i+1][j][k]-4*Fd[i][j][k]+Fd[i-1][j][k],2);
					ISd[i][j][k][1]=13.0/12.0*pow(Fd[i+2][j][k]-2*Fd[i+1][j][k]+Fd[i][j][k],2)
						+0.25*pow(Fd[i+2][j][k]-Fd[i][j][k],2);
					ISd[i][j][k][2]=13.0/12.0*pow(Fd[i+3][j][k]-2*Fd[i+2][j][k]+Fd[i+1][j][k],2)
						+0.25*pow(Fd[i+3][j][k]-4*Fd[i+2][j][k]+3*Fd[i+1][j][k],2);
				}
	//计算alphapd以及omegapd
	for(i=2;i<=Nx+3;i++)    
		for(j=2;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
				{
					alphap[i][j][k][0]=0.1/pow(pow(10,-16)+ISp[i][j][k][0],P);
					alphap[i][j][k][1]=0.6/pow(pow(10,-16)+ISp[i][j][k][1],P);
					alphap[i][j][k][2]=0.3/pow(pow(10,-16)+ISp[i][j][k][2],P);

					alphad[i][j][k][0]=0.3/pow(pow(10,-16)+ISd[i][j][k][0],P);
					alphad[i][j][k][1]=0.6/pow(pow(10,-16)+ISd[i][j][k][1],P);
					alphad[i][j][k][2]=0.1/pow(pow(10,-16)+ISd[i][j][k][2],P);
					for(l=0;l<=2;l++)
					{
						omegap[i][j][k][l]=alphap[i][j][k][l]/(alphap[i][j][k][0]+alphap[i][j][k][1]+alphap[i][j][k][2]);
						omegad[i][j][k][l]=alphad[i][j][k][l]/(alphad[i][j][k][0]+alphad[i][j][k][1]+alphad[i][j][k][2]);						
					}
				}
	//计算q3pd q5pd
	for(i=2;i<=Nx+3;i++)    
		for(j=2;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
				{
					q3p[i][j][k][0]=1.0/3.0*Fp[i-2][j][k]-7.0/6.0*Fp[i-1][j][k]+11.0/6.0*Fp[i][j][k];
					q3p[i][j][k][1]=-1.0/6.0*Fp[i-1][j][k]+5.0/6.0*Fp[i][j][k]+1.0/3.0*Fp[i+1][j][k];
					q3p[i][j][k][2]=1.0/3.0*Fp[i][j][k]+5.0/6.0*Fp[i+1][j][k]-1.0/6.0*Fp[i+2][j][k];

					q3d[i][j][k][0]=-1.0/6.0*Fd[i-1][j][k]+5.0/6.0*Fd[i][j][k]+1.0/3.0*Fd[i+1][j][k];
					q3d[i][j][k][1]=1.0/3.0*Fd[i][j][k]+5.0/6.0*Fd[i+1][j][k]-1.0/6.0*Fd[i+2][j][k];
					q3d[i][j][k][2]=11.0/6.0*Fd[i+1][j][k]-7.0/6.0*Fd[i+2][j][k]+1.0/3.0*Fd[i+3][j][k];
				}

	//计算LF_pd
	for(i=2;i<=Nx+3;i++)    
		for(j=2;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
				{
					F_p[i][j][k]=0;
					F_d[i][j][k]=0;
				}

	for(i=2;i<=Nx+3;i++)    
		for(j=2;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					for(l=0;l<=2;l++)
					{
						F_p[i][j][k]+=omegap[i][j][k][l]*q3p[i][j][k][l];
						F_d[i][j][k]+=omegad[i][j][k][l]*q3d[i][j][k][l];
					}
	//计算F_
	for(i=2;i<=Nx+3;i++)    
		for(j=2;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					F_[i][j][k]=F_p[i][j][k]+F_d[i][j][k];

	
	//分区计算U
	for(i=3;i<=Nx+3;i++)    
		for(j=3;j<=int(0.5/dy)+3;j++)	
			for(k=0;k<=3;k++)
				U[i][j][k]=U[i][j][k]-r*(F_[i][j][k]-F_[i-1][j][k]);

	for(i=int(1.0/dx)+3;i<=int(2.0/dx)+3;i++)    
		for(j=int(0.5/dy)+4;j<=Ny+3;j++)
			for(k=0;k<=3;k++)
				U[i][j][k]=U[i][j][k]-r*(F_[i][j][k]-F_[i-1][j][k]);


}


void WENO_y(double U[Nx+7][Ny+7][4],double ISp[Nx+7][Ny+7][4][3],
			double ISd[Nx+7][Ny+7][4][3],double omegap[Nx+7][Ny+7][4][3],
			double omegad[Nx+7][Ny+7][4][3],double alphap[Nx+7][Ny+7][4][3],double alphad[Nx+7][Ny+7][4][3],
			double q3p[Nx+7][Ny+7][4][3],double q3d[Nx+7][Ny+7][4][3],
			double Gp[Nx+7][Ny+7][4],double Gd[Nx+7][Ny+7][4],
			double G_p[Nx+7][Ny+7][4],double G_d[Nx+7][Ny+7][4],double G_[Nx+7][Ny+7][4],
			double dx,double dy,double dt)       //此函数将之前计算的特征之余特征向量通过TVD算法计算U
{
	int i,j,k,l;
	dt=CFL(U,dx,dy,WENOCFL);
	double r=dt/dy;
	//计算ISpd
	for(i=2;i<=Nx+3;i++)    
		for(j=2;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
				{
					ISp[i][j][k][0]=13.0/12.0*pow(Gp[i][j-2][k]-2*Gp[i][j-1][k]+Gp[i][j][k],2)
						+0.25*pow(Gp[i][j-2][k]-4*Gp[i][j-1][k]+3*Gp[i][j][k],2);
					ISp[i][j][k][1]=13.0/12.0*pow(Gp[i][j-1][k]-2*Gp[i][j][k]+Gp[i][j+1][k],2)
						+0.25*pow(Gp[i][j-1][k]-Gp[i][j+1][k],2);
					ISp[i][j][k][2]=13.0/12.0*pow(Gp[i][j][k]-2*Gp[i][j+1][k]+Gp[i][j+2][k],2)
						+0.25*pow(3*Gp[i][j][k]-4*Gp[i][j+1][k]+Gp[i][j+2][k],2);

					ISd[i][j][k][0]=13.0/12.0*pow(Gd[i][j+1][k]-2*Gd[i][j][k]+Gd[i][j-1][k],2)
						+0.25*pow(3*Gd[i][j+1][k]-4*Gd[i][j][k]+Gd[i][j-1][k],2);
					ISd[i][j][k][1]=13.0/12.0*pow(Gd[i][j+2][k]-2*Gd[i][j+1][k]+Gd[i][j][k],2)
						+0.25*pow(Gd[i][j+2][k]-Gd[i][j][k],2);
					ISd[i][j][k][2]=13.0/12.0*pow(Gd[i][j+3][k]-2*Gd[i][j+2][k]+Gd[i][j+1][k],2)
						+0.25*pow(Gd[i][j+3][k]-4*Gd[i][j+2][k]+3*Gd[i][j+1][k],2);
				}

	//计算alphapd以及omegapd
	for(i=2;i<=Nx+3;i++)    
		for(j=2;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
				{
					alphap[i][j][k][0]=0.1/pow(pow(10,-16)+ISp[i][j][k][0],P);
					alphap[i][j][k][1]=0.6/pow(pow(10,-16)+ISp[i][j][k][1],P);
					alphap[i][j][k][2]=0.3/pow(pow(10,-16)+ISp[i][j][k][2],P);

					alphad[i][j][k][0]=0.3/pow(pow(10,-16)+ISd[i][j][k][0],P);
					alphad[i][j][k][1]=0.6/pow(pow(10,-16)+ISd[i][j][k][1],P);
					alphad[i][j][k][2]=0.1/pow(pow(10,-16)+ISd[i][j][k][2],P);
					for(l=0;l<=2;l++)
					{
						omegap[i][j][k][l]=alphap[i][j][k][l]/(alphap[i][j][k][0]+alphap[i][j][k][1]+alphap[i][j][k][2]);
						omegad[i][j][k][l]=alphad[i][j][k][l]/(alphad[i][j][k][0]+alphad[i][j][k][1]+alphad[i][j][k][2]);						
					}
				}

	//计算q3pd q5pd
	for(i=2;i<=Nx+3;i++)    
		for(j=2;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
				{
					q3p[i][j][k][0]=1.0/3.0*Gp[i][j-2][k]-7.0/6.0*Gp[i][j-1][k]+11.0/6.0*Gp[i][j][k];
					q3p[i][j][k][1]=-1.0/6.0*Gp[i][j-1][k]+5.0/6.0*Gp[i][j][k]+1.0/3.0*Gp[i][j+1][k];
					q3p[i][j][k][2]=1.0/3.0*Gp[i][j][k]+5.0/6.0*Gp[i][j+1][k]-1.0/6.0*Gp[i][j+2][k];

					q3d[i][j][k][0]=-1.0/6.0*Gd[i][j-1][k]+5.0/6.0*Gd[i][j][k]+1.0/3.0*Gd[i][j+1][k];
					q3d[i][j][k][1]=1.0/3.0*Gd[i][j][k]+5.0/6.0*Gd[i][j+1][k]-1.0/6.0*Gd[i][j+2][k];
					q3d[i][j][k][2]=11.0/6.0*Gd[i][j+1][k]-7.0/6.0*Gd[i][j+2][k]+1.0/3.0*Gd[i][j+3][k];
				}
	
	//计算LG_pd
	for(i=2;i<=Nx+3;i++)    
		for(j=2;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
				{
					G_p[i][j][k]=0;
					G_d[i][j][k]=0;
				}

	for(i=2;i<=Nx+3;i++)    
		for(j=2;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					for(l=0;l<=2;l++)
					{
						G_p[i][j][k]+=omegap[i][j][k][l]*q3p[i][j][k][l];
						G_d[i][j][k]+=omegad[i][j][k][l]*q3d[i][j][k][l];
					}

	//计算G_
	for(i=2;i<=Nx+3;i++)    
		for(j=2;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					G_[i][j][k]=G_p[i][j][k]+G_d[i][j][k];

	//分区计算U
	for(i=3;i<=Nx+3;i++)    
		for(j=3;j<=int(0.5/dy)+3;j++)
			for(k=0;k<=3;k++)
				U[i][j][k]=U[i][j][k]-r*(G_[i][j][k]-G_[i][j-1][k]);

	for(i=int(1.0/dx)+3;i<=int(2.0/dx)+3;i++)    
		for(j=int(0.5/dy)+4;j<=Ny+3;j++)
			for(k=0;k<=3;k++)
				U[i][j][k]=U[i][j][k]-r*(G_[i][j][k]-G_[i][j-1][k]);

}


void WENO_Solver(double U[Nx+7][Ny+7][4],double U1[Nx+7][Ny+7][4],double U2[Nx+7][Ny+7][4],
				 double ISp[Nx+7][Ny+7][4][3],double ISd[Nx+7][Ny+7][4][3],double omegap[Nx+7][Ny+7][4][3],
			double omegad[Nx+7][Ny+7][4][3],double alphap[Nx+7][Ny+7][4][3],double alphad[Nx+7][Ny+7][4][3],
			double q3p[Nx+7][Ny+7][4][3],double q3d[Nx+7][Ny+7][4][3],double LAMDA_[Nx+7][Ny+7][4][4],
			double Fp[Nx+7][Ny+7][4],double Fd[Nx+7][Ny+7][4],
			double F_p[Nx+7][Ny+7][4],double F_d[Nx+7][Ny+7][4],double F[Nx+7][Ny+7][4],double F_[Nx+7][Ny+7][4],
			double G[Nx+7][Ny+7][4],double G_[Nx+7][Ny+7][4],double Gp[Nx+7][Ny+7][4],double Gd[Nx+7][Ny+7][4],
			double G_p[Nx+7][Ny+7][4],double G_d[Nx+7][Ny+7][4],
			double dx,double dy,double &dt)
{
	bound(U,dx,dy);
	LF_x(U,LAMDA_,F,Fp,Fd);
	WENO_x(U,ISp,ISd,omegap,omegad,alphap,alphad,q3p,q3d,Fp,Fd,F_p,F_d,F_,dx,dy,dt);
	bound(U,dx,dy);
	LF_y(U,LAMDA_,G,Gp,Gd);
	WENO_y(U,ISp,ISd,omegap,omegad,alphap,alphad,q3p,q3d,Gp,Gd,G_p,G_d,G_,dx,dy,dt);
	bound(U,dx,dy);
	LF_y(U,LAMDA_,G,Gp,Gd);
	WENO_y(U,ISp,ISd,omegap,omegad,alphap,alphad,q3p,q3d,Gp,Gd,G_p,G_d,G_,dx,dy,dt);
	bound(U,dx,dy);
	LF_x(U,LAMDA_,F,Fp,Fd);
	WENO_x(U,ISp,ISd,omegap,omegad,alphap,alphad,q3p,q3d,Fp,Fd,F_p,F_d,F_,dx,dy,dt);
	bound(U,dx,dy);
}














//--------------------------紧致格式------------------------------
//这里使用的是人工滤波法,选择的是Harten提出的开关函数
//变量解释:
//注意这里大写的F表示f的导数的逼近,而f表示的是原方程的变量.
void A1_Compact_QuickLocation()
{
}
void CPT_fpfd(double U[4],double fp[4],double fd[4])
{
	double H,rou,u,v,p,a,lambda1,lambda2,lambda3,lambda4,lambda1p,lambda2p,lambda3p,lambda4p,lambda1d,lambda2d,lambda3d,lambda4d;
	//先算fp即f+
	rou=U[0];
	u=U[1]/U[0];
	v=U[2]/U[0];
	p=(GAMA-1)*(U[3]-0.5*rou*(u*u+v*v));
	a=pow(GAMA*p/rou,0.5);
	H=a*a/(GAMA-1)+0.5*(u*u+v*v);
	lambda1=u;
	lambda2=u;
	lambda3=u-a;
	lambda4=u+a;
	//计算正特征值
	lambda1p=0.5*(fabs(lambda1)+lambda1);
	lambda2p=0.5*(fabs(lambda2)+lambda2);
	lambda3p=0.5*(fabs(lambda3)+lambda3);
	lambda4p=0.5*(fabs(lambda4)+lambda4);
	//特征值修正
	lambda1p=0.5*(lambda1p+pow(lambda1p*lambda1p+pow(10,-16),0.5));
	lambda2p=0.5*(lambda2p+pow(lambda2p*lambda2p+pow(10,-16),0.5));
	lambda3p=0.5*(lambda3p+pow(lambda3p*lambda3p+pow(10,-16),0.5));
	lambda4p=0.5*(lambda4p+pow(lambda4p*lambda4p+pow(10,-16),0.5));
	fp[0]=rou/2.0/GAMA*(2*(GAMA-1)*lambda1p+lambda3p+lambda4p);
	fp[1]=rou/2.0/GAMA*(2*u*(GAMA-1)*lambda1p+(u-a)*lambda3p+(u+a)*lambda4p);
	fp[2]=rou/2.0/GAMA*(v*2*(GAMA-1)*lambda1p+v*lambda3p+v*lambda4p);
	fp[3]=rou/2.0/GAMA*(2*((GAMA-1)*H-a*a)*lambda1p+(H-a*u)*lambda3p+(H+a*u)*lambda4p);
	//现在就算fd即f-
	lambda1d=-0.5*(fabs(lambda1)-lambda1);
	lambda2d=-0.5*(fabs(lambda2)-lambda2);
	lambda3d=-0.5*(fabs(lambda3)-lambda3);
	lambda4d=-0.5*(fabs(lambda4)-lambda4);
	lambda1d=0.5*(lambda1d-pow(lambda1d*lambda1d+pow(10,-16),0.5));
	lambda2d=0.5*(lambda2d-pow(lambda2d*lambda2d+pow(10,-16),0.5));
	lambda3d=0.5*(lambda3d-pow(lambda3d*lambda3d+pow(10,-16),0.5));
	lambda4d=0.5*(lambda4d-pow(lambda4d*lambda4d+pow(10,-16),0.5));
	fd[0]=rou/2.0/GAMA*(2*(GAMA-1)*lambda1d+lambda3d+lambda4d);
	fd[1]=rou/2.0/GAMA*(2*u*(GAMA-1)*lambda1d+(u-a)*lambda3d+(u+a)*lambda4d);
	fd[2]=rou/2.0/GAMA*(v*2*(GAMA-1)*lambda1d+v*lambda3d+v*lambda4d);
	fd[3]=rou/2.0/GAMA*(2*((GAMA-1)*H-a*a)*lambda1d+(H-a*u)*lambda3d+(H+a*u)*lambda4d);
}


void Compact_x(double U[Nx+7][Ny+7][4],double Ut[Nx+7][Ny+7][4],double Fp[Nx+7][Ny+7][4],double Fd[Nx+7][Ny+7][4],
			   double fp[Nx+7][Ny+7][4],double fd[Nx+7][Ny+7][4],
			   double dx,double dy,double dt,double p[Nx+7][Ny+7])
{
	int i,j,k;
	double r=dt/dx;
	double q,rou,u,v;
	double nu;
	for(i=0;i<=Nx+6;i++)
		for(j=0;j<=Ny+6;j++)
			if(U[i][j][0]!=0)
			{
				rou=U[i][j][0];
				u=U[i][j][1]/U[i][j][0];
				v=U[i][j][2]/U[i][j][0];
				p[i][j]=(GAMA-1)*(U[i][j][3]-0.5*rou*(u*u+v*v));
			}

	nu=5*r*(1-5*r);
	//人工黏性
	for(i=3;i<=Nx+3;i++)    
		for(j=3;j<=int(0.5/dy)+3;j++)
		{
			q=fabs(fabs(U[i+1][j][0]-U[i][j][0])-fabs(U[i][j][0]-U[i-1][j][0]))/(fabs(U[i+1][j][0]-U[i][j][0])+fabs(U[i][j][0]-U[i-1][j][0])+pow(10,-100));			
			for(k=0;k<=3;k++)
			{
				Ut[i][j][k]=U[i][j][k]+0.5*nu*q*(U[i+1][j][k]-2*U[i][j][k]+U[i-1][j][k]);
			}
		}

	for(i=int(1.0/dx)+3;i<=int(2.0/dx)+3;i++)    
		for(j=int(0.5/dy)+4;j<=Ny+3;j++)
		{
			q=fabs(fabs(U[i+1][j][0]-U[i][j][0])-fabs(U[i][j][0]-U[i-1][j][0]))/(fabs(U[i+1][j][0]-U[i][j][0])+fabs(U[i][j][0]-U[i-1][j][0])+pow(10,-100));	
			for(k=0;k<=3;k++)
			{
				Ut[i][j][k]=U[i][j][k]+0.5*nu*q*(U[i+1][j][k]-2*U[i][j][k]+U[i-1][j][k]);
			}
		}

	for(i=3;i<=Nx+3;i++)    
		for(j=3;j<=int(0.5/dy)+3;j++)
		{
			for(k=0;k<=3;k++)
			{
				U[i][j][k]=Ut[i][j][k];
			}
		}

	for(i=int(1.0/dx)+3;i<=int(2.0/dx)+3;i++)    
		for(j=int(0.5/dy)+4;j<=Ny+3;j++)
			for(k=0;k<=3;k++)
				U[i][j][k]=Ut[i][j][k];

	//计算fp和fd
	for(i=0;i<=Nx+6;i++)
		for(j=0;j<=Ny+6;j++)
			if(U[i][j][0]!=0)
				CPT_fpfd(U[i][j],fp[i][j],fd[i][j]);

	//分区计算F+的左值
	for(j=3;j<=int(0.5/dy)+3;j++)
		for(k=0;k<=3;k++)
			Fp[2][j][k]=0.5*(3*(fp[3][j][k]-fp[2][j][k])-(fp[4][j][k]-fp[3][j][k]));

	for(j=int(0.5/dy)+4;j<=Ny+3;j++)
		for(k=0;k<=3;k++)
			Fp[int(1.0/dx)+2][j][k]=0.5*(3*(fp[int(1.0/dx)+3][j][k]-fp[int(1.0/dx)+2][j][k])-(fp[int(1.0/dx)+4][j][k]-fp[int(1.0/dx)+3][j][k]));

	//分区计算全场的F+
	for(i=3;i<=Nx+3;i++)
		for(j=3;j<=int(0.5/dy)+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					Fp[i][j][k]=-0.5*Fp[i-1][j][k]+5.0/4.0*(fp[i][j][k]-fp[i-1][j][k])+0.25*(fp[i+1][j][k]-fp[i][j][k]);

	for(i=int(1.0/dx)+3;i<=int(2.0/dx)+3;i++)
		for(j=int(0.5/dy)+4;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					Fp[i][j][k]=-0.5*Fp[i-1][j][k]+5.0/4.0*(fp[i][j][k]-fp[i-1][j][k])+0.25*(fp[i+1][j][k]-fp[i][j][k]);

	//分区计算F-的右值
	for(j=3;j<=int(0.5/dy)+3;j++)
		for(k=0;k<=3;k++)
			Fd[Nx+4][j][k]=0.5*(3*(fd[Nx+4][j][k]-fd[Nx+3][j][k])-(fd[Nx+3][j][k]-fd[Nx+2][j][k]));

	for(j=int(0.5/dy)+4;j<=Ny+3;j++)
		for(k=0;k<=3;k++)
			Fd[int(2.0/dx)+4][j][k]=0.5*(3*(fd[int(2.0/dx)+4][j][k]-fd[int(2.0/dx)+3][j][k])-(fd[int(2.0/dx)+3][j][k]-fd[int(2.0/dx)+2][j][k]));

	//分区计算全场的F-值
	for(i=Nx+3;i>=3;i--)
		for(j=3;j<=int(0.5/dy)+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					Fd[i][j][k]=-0.5*Fd[i+1][j][k]+5.0/4.0*(fd[i+1][j][k]-fd[i][j][k])+0.25*(fd[i][j][k]-fd[i-1][j][k]);

	for(i=int(2.0/dx)+3;i>=int(1.0/dx)+3;i--)
		for(j=int(0.5/dy)+4;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					Fd[i][j][k]=-0.5*Fd[i+1][j][k]+5.0/4.0*(fd[i+1][j][k]-fd[i][j][k])+0.25*(fd[i][j][k]-fd[i-1][j][k]);

}

void Compact_RK_x(double U[Nx+7][Ny+7][4],double Ut[Nx+7][Ny+7][4],double U1[Nx+7][Ny+7][4],double U2[Nx+7][Ny+7][4],
		  double Fp[Nx+7][Ny+7][4],double Fd[Nx+7][Ny+7][4],
			   double fp[Nx+7][Ny+7][4],double fd[Nx+7][Ny+7][4],
			   double dx,double dy,double dt,double p[Nx+7][Ny+7],double LAMDA_[Nx+7][Ny+7][4][4])
{
	int i,j,k;
	dt=CFL(U,dx,dy,COMPACTCFL);
	double r=dt/dx;
	//分区计算U1
	for(i=3;i<=Nx+3;i++)
		for(j=3;j<=int(0.5/dy)+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					U1[i][j][k]=U[i][j][k]-r*(Fp[i][j][k]+Fd[i][j][k]);

	for(i=int(1.0/dx)+3;i<=int(2.0/dx)+3;i++)       
		for(j=int(0.5/dy)+4;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					U1[i][j][k]=U[i][j][k]-r*(Fp[i][j][k]+Fd[i][j][k]);

	//U1求边界条件
	bound(U1,dx,dy);
	//计算dt
	dt=CFL(U1,dx,dy,COMPACTCFL);
	r=dt/dx;
	for(i=0;i<=Nx+6;i++)
		for(j=0;j<=Ny+6;j++)
			if(U[i][j][0]!=0)
				CPT_fpfd(U1[i][j],fp[i][j],fd[i][j]);
	//用U1算Fpd
	Compact_x(U1,Ut,Fp,Fd,fp,fd,dx,dy,dt,p);
	//RK算U2
	for(i=3;i<=Nx+3;i++)
		for(j=3;j<=int(0.5/dy)+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					U2[i][j][k]=0.75*U[i][j][k]+0.25*U1[i][j][k]-r*(Fp[i][j][k]+Fd[i][j][k]);

	for(i=int(1.0/dx)+3;i<=int(2.0/dx)+3;i++)       
		for(j=int(0.5/dy)+4;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					U2[i][j][k]=0.75*U[i][j][k]+0.25*U1[i][j][k]-r*(Fp[i][j][k]+Fd[i][j][k]);

	bound(U2,dx,dy);
	dt=CFL(U2,dx,dy,COMPACTCFL);
	r=dt/dx;
	for(i=0;i<=Nx+6;i++)
		for(j=0;j<=Ny+6;j++)
			if(U[i][j][0]!=0)
				CPT_fpfd(U2[i][j],fp[i][j],fd[i][j]);
	Compact_x(U2,Ut,Fp,Fd,fp,fd,dx,dy,dt,p);
	//RK算U3-即U
	for(i=3;i<=Nx+3;i++)
		for(j=3;j<=int(0.5/dy)+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					U[i][j][k]=1.0/3.0*U[i][j][k]+2.0/3.0*U2[i][j][k]-2.0/3.0*r*(Fp[i][j][k]+Fd[i][j][k]);

	for(i=int(1.0/dx)+3;i<=int(2.0/dx)+3;i++)       
		for(j=int(0.5/dy)+4;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					U[i][j][k]=1.0/3.0*U[i][j][k]+2.0/3.0*U2[i][j][k]-2.0/3.0*r*(Fp[i][j][k]+Fd[i][j][k]);


}
///////////-------------------------------y
void CPT_gpgd(double U[4],double gp[4],double gd[4])
{
	double H,rou,u,v,p,a,lambda1,lambda2,lambda3,lambda4,lambda1p,lambda2p,lambda3p,lambda4p,lambda1d,lambda2d,lambda3d,lambda4d;
	//计算Gp即G+的函数
	rou=U[0];
	u=U[1]/U[0];
	v=U[2]/U[0];
	p=(GAMA-1)*(U[3]-0.5*rou*(u*u+v*v));
	a=pow(GAMA*p/rou,0.5);
	H=a*a/(GAMA-1)+0.5*(u*u+v*v);
	lambda1=v;
	lambda2=v;
	lambda3=v-a;
	lambda4=v+a;
	lambda1p=0.5*(fabs(lambda1)+lambda1);
	lambda2p=0.5*(fabs(lambda2)+lambda2);
	lambda3p=0.5*(fabs(lambda3)+lambda3);
	lambda4p=0.5*(fabs(lambda4)+lambda4);
	lambda1p=0.5*(lambda1p+pow(lambda1p*lambda1p+pow(10,-16),0.5));
	lambda2p=0.5*(lambda2p+pow(lambda2p*lambda2p+pow(10,-16),0.5));
	lambda3p=0.5*(lambda3p+pow(lambda3p*lambda3p+pow(10,-16),0.5));
	lambda4p=0.5*(lambda4p+pow(lambda4p*lambda4p+pow(10,-16),0.5));
	gp[0]=rou/2.0/GAMA*(2*(GAMA-1)*lambda1p+lambda3p+lambda4p);
	gp[1]=rou/2.0/GAMA*(u*2*(GAMA-1)*lambda1p+u*lambda3p+u*lambda4p);
	gp[2]=rou/2.0/GAMA*(2*v*(GAMA-1)*lambda1p+(v-a)*lambda3p+(v+a)*lambda4p);
	gp[3]=rou/2.0/GAMA*(2*((GAMA-1)*H-a*a)*lambda1p+(H-a*v)*lambda3p+(H+a*v)*lambda4p);
	//计算Gd即G-的函数
	lambda1d=-0.5*(fabs(lambda1)-lambda1);
	lambda2d=-0.5*(fabs(lambda2)-lambda2);
	lambda3d=-0.5*(fabs(lambda3)-lambda3);
	lambda4d=-0.5*(fabs(lambda4)-lambda4);
	lambda1d=0.5*(lambda1d-pow(lambda1d*lambda1d+pow(10,-16),0.5));
	lambda2d=0.5*(lambda2d-pow(lambda2d*lambda2d+pow(10,-16),0.5));
	lambda3d=0.5*(lambda3d-pow(lambda3d*lambda3d+pow(10,-16),0.5));
	lambda4d=0.5*(lambda4d-pow(lambda4d*lambda4d+pow(10,-16),0.5));
	gd[0]=rou/2.0/GAMA*(2*(GAMA-1)*lambda1d+lambda3d+lambda4d);
	gd[1]=rou/2.0/GAMA*(u*2*(GAMA-1)*lambda1d+u*lambda3d+u*lambda4d);
	gd[2]=rou/2.0/GAMA*(2*v*(GAMA-1)*lambda1d+(v-a)*lambda3d+(v+a)*lambda4d);
	gd[3]=rou/2.0/GAMA*(2*((GAMA-1)*H-a*a)*lambda1d+(H-a*v)*lambda3d+(H+a*v)*lambda4d);
}

void Compact_y(double U[Nx+7][Ny+7][4],double Ut[Nx+7][Ny+7][4],double Gp[Nx+7][Ny+7][4],double Gd[Nx+7][Ny+7][4],
			   double gp[Nx+7][Ny+7][4],double gd[Nx+7][Ny+7][4],double dx,double dy,double dt,double p[Nx+7][Ny+7])
{
	int i,j,k;
	double r=dt/dy;
	double q,rou,u,v;
	double nu;
	for(i=0;i<=Nx+6;i++)
		for(j=0;j<=Ny+6;j++)
			if(U[i][j][0]!=0)
			{
				rou=U[i][j][0];
			    u=U[i][j][1]/U[i][j][0];
			    v=U[i][j][2]/U[i][j][0];
			    p[i][j]=(GAMA-1)*(U[i][j][3]-0.5*rou*(u*u+v*v));
			}

	nu=5*r*(1-5*r);
	//人工黏性
	for(i=3;i<=Nx+3;i++)    
		for(j=3;j<=int(0.5/dy)+3;j++)
		{
			q=fabs(fabs(U[i][j+1][0]-U[i][j][0])-fabs(U[i][j][0]-U[i][j-1][0]))/(fabs(U[i][j+1][0]-U[i][j][0])+fabs(U[i][j][0]-U[i][j-1][0])+pow(10,-100));		
			for(k=0;k<=3;k++)
			{
				Ut[i][j][k]=U[i][j][k]+0.5*nu*q*(U[i][j+1][k]-2*U[i][j][k]+U[i][j-1][k]);
			}
		}
	for(i=int(1.0/dx)+3;i<=int(2.0/dx)+3;i++)    
		for(j=int(0.5/dy)+4;j<=Ny+3;j++)
		{
			q=fabs(fabs(U[i][j+1][0]-U[i][j][0])-fabs(U[i][j][0]-U[i][j-1][0]))/(fabs(U[i][j+1][0]-U[i][j][0])+fabs(U[i][j][0]-U[i][j-1][0])+pow(10,-100));		
			for(k=0;k<=3;k++)
			{
				Ut[i][j][k]=U[i][j][k]+0.5*nu*q*(U[i][j+1][k]-2*U[i][j][k]+U[i][j-1][k]);
			}
		}

	for(i=3;i<=Nx+3;i++)    
		for(j=3;j<=int(0.5/dy)+3;j++)
			for(k=0;k<=3;k++)
				U[i][j][k]=Ut[i][j][k];

	for(i=int(1.0/dx)+3;i<=int(2.0/dx)+3;i++)    
		for(j=int(0.5/dy)+4;j<=Ny+3;j++)
			for(k=0;k<=3;k++)
				U[i][j][k]=Ut[i][j][k];

	//计算gp和gd
	for(i=0;i<=Nx+6;i++)
		for(j=0;j<=Ny+6;j++)
			if(U[i][j][0]!=0)
				CPT_gpgd(U[i][j],gp[i][j],gd[i][j]);

	//分区计算G+的下值
	for(i=3;i<=Nx+3;i++)
		for(k=0;k<=3;k++)
			Gp[i][2][k]=0.5*(3*(gp[i][3][k]-gp[i][2][k])-(gp[i][4][k]-gp[i][3][k]));

	//分区计算全场的G+
	for(j=3;j<=Ny+3;j++)
		for(i=3;i<=Nx+3;i++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					Gp[i][j][k]=-0.5*Gp[i][j-1][k]+5.0/4.0*(gp[i][j][k]-gp[i][j-1][k])+0.25*(gp[i][j+1][k]-gp[i][j][k]);

	//分区计算G-的上值
	for(i=3;i<=int(1.0/dx)+2;i++)
		for(k=0;k<=3;k++)
			Gd[i][int(0.5/dy)+4][k]=0.5*(3*(gd[i][int(0.5/dy)+4][k]-gd[i][int(0.5/dy)+3][k])-(gd[i][int(0.5/dy)+3][k]-gd[i][int(0.5/dy)+2][k]));

	for(i=int(1.0/dx)+3;i<=int(2.0/dx)+3;i++)
		for(k=0;k<=3;k++)
			Gd[i][int(1.0/dy)+4][k]=0.5*(3*(gd[i][int(1.0/dy)+4][k]-gd[i][int(1.0/dy)+3][k])-(gd[i][int(1.0/dy)+3][k]-gd[i][int(1.0/dy)+2][k]));

	for(i=int(2.0/dx)+4;i<=Nx+3;i++)
		for(k=0;k<=3;k++)
			Gd[i][int(0.5/dy)+4][k]=0.5*(3*(gd[i][int(0.5/dy)+4][k]-gd[i][int(0.5/dy)+3][k])-(gd[i][int(0.5/dy)+3][k]-gd[i][int(0.5/dy)+2][k]));

	//分区计算全场的F-值
	for(j=int(0.5/dy)+3;j>=3;j--)
		for(i=3;i<=int(1.0/dx)+2;i++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					Gd[i][j][k]=-0.5*Gd[i][j+1][k]+5.0/4.0*(gd[i][j+1][k]-gd[i][j][k])+0.25*(gd[i][j][k]-gd[i][j-1][k]);

	for(j=int(0.5/dy)+3;j>=3;j--)
		for(i=int(2.0/dx)+4;i<=Nx+3;i++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					Gd[i][j][k]=-0.5*Gd[i][j+1][k]+5.0/4.0*(gd[i][j+1][k]-gd[i][j][k])+0.25*(gd[i][j][k]-gd[i][j-1][k]);

	for(j=int(1.0/dy)+3;j>=3;j--)
		for(i=int(1.0/dx)+3;i<=int(2.0/dx)+3;i++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					Gd[i][j][k]=-0.5*Gd[i][j+1][k]+5.0/4.0*(gd[i][j+1][k]-gd[i][j][k])+0.25*(gd[i][j][k]-gd[i][j-1][k]);

}
void Compact_RK_y(double U[Nx+7][Ny+7][4],double Ut[Nx+7][Ny+7][4],double U1[Nx+7][Ny+7][4],double U2[Nx+7][Ny+7][4],
		  double Gp[Nx+7][Ny+7][4],double Gd[Nx+7][Ny+7][4],
			   double gp[Nx+7][Ny+7][4],double gd[Nx+7][Ny+7][4],double dx,double dy,double dt,double p[Nx+7][Ny+7],double LAMDA_[Nx+7][Ny+7][4][4])
{
	int i,j,k;
	dt=CFL(U,dx,dy,COMPACTCFL);
	double r=dt/dy;
	//分区计算U
	for(i=3;i<=Nx+3;i++)
		for(j=3;j<=int(0.5/dy)+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					U1[i][j][k]=U[i][j][k]-r*(Gp[i][j][k]+Gd[i][j][k]);

	for(i=int(1.0/dx)+3;i<=int(2.0/dx)+3;i++)       
		for(j=int(0.5/dy)+4;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					U1[i][j][k]=U[i][j][k]-r*(Gp[i][j][k]+Gd[i][j][k]);

	bound(U1,dx,dy);
	dt=CFL(U1,dx,dy,COMPACTCFL);
	r=dt/dy;
	for(i=0;i<=Nx+6;i++)
		for(j=0;j<=Ny+6;j++)
			if(U[i][j][0]!=0)
				CPT_gpgd(U1[i][j],gp[i][j],gd[i][j]);
	Compact_y(U1,Ut,Gp,Gd,gp,gd,dx,dy,dt,p);
	for(i=3;i<=Nx+3;i++)
		for(j=3;j<=int(0.5/dy)+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					U2[i][j][k]=0.75*U[i][j][k]+0.25*U1[i][j][k]-r*(Gp[i][j][k]+Gd[i][j][k]);

	for(i=int(1.0/dx)+3;i<=int(2.0/dx)+3;i++)       
		for(j=int(0.5/dy)+4;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					U2[i][j][k]=0.75*U[i][j][k]+0.25*U1[i][j][k]-r*(Gp[i][j][k]+Gd[i][j][k]);

	bound(U2,dx,dy);
	dt=CFL(U2,dx,dy,COMPACTCFL);
	r=dt/dy;
	for(i=0;i<=Nx+6;i++)
		for(j=0;j<=Ny+6;j++)
			if(U[i][j][0]!=0)
				CPT_gpgd(U2[i][j],gp[i][j],gd[i][j]);
	Compact_y(U2,Ut,Gp,Gd,gp,gd,dx,dy,dt,p);
	for(i=3;i<=Nx+3;i++)
		for(j=3;j<=int(0.5/dy)+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					U[i][j][k]=1.0/3.0*U[i][j][k]+2.0/3.0*U2[i][j][k]-2.0/3.0*r*(Gp[i][j][k]+Gd[i][j][k]);

	for(i=int(1.0/dx)+3;i<=int(2.0/dx)+3;i++)       
		for(j=int(0.5/dy)+4;j<=Ny+3;j++)
			if(U[i][j][0]!=0)
				for(k=0;k<=3;k++)
					U[i][j][k]=1.0/3.0*U[i][j][k]+2.0/3.0*U2[i][j][k]-2.0/3.0*r*(Gp[i][j][k]+Gd[i][j][k]);

}

void Compact_Solver(double U[Nx+7][Ny+7][4],double Ut[Nx+7][Ny+7][4],double U1[Nx+7][Ny+7][4],double U2[Nx+7][Ny+7][4],
					double Fp[Nx+7][Ny+7][4],double Fd[Nx+7][Ny+7][4],
					double fp[Nx+7][Ny+7][4],double fd[Nx+7][Ny+7][4],
					double Gp[Nx+7][Ny+7][4],double Gd[Nx+7][Ny+7][4],
					double gp[Nx+7][Ny+7][4],double gd[Nx+7][Ny+7][4],
					double dx,double dy,double dt,double p[Nx+7][Ny+7])
{
	bound(U,dx,dy);
	Compact_x(U,Ut,Fp,Fd,fp,fd,dx,dy,dt,p);
	Compact_RK_x(U,Ut,U1,U2,Fp,Fd,fp,fd,dx,dy,dt,p,LAMDA_);
	bound(U,dx,dy);
	Compact_y(U,Ut,Gp,Gd,gp,gd,dx,dy,dt,p);
	Compact_RK_y(U,Ut,U1,U2,Gp,Gd,gp,gd,dx,dy,dt,p,LAMDA_);
	bound(U,dx,dy);
	Compact_y(U,Ut,Gp,Gd,gp,gd,dx,dy,dt,p);
	Compact_RK_y(U,Ut,U1,U2,Gp,Gd,gp,gd,dx,dy,dt,p,LAMDA_);
	bound(U,dx,dy);
	Compact_x(U,Ut,Fp,Fd,fp,fd,dx,dy,dt,p);
	Compact_RK_x(U,Ut,U1,U2,Fp,Fd,fp,fd,dx,dy,dt,p,LAMDA_);
}

//动画plt文件输出器函数
//入口:U,dx,dy,n,m   n表示调用求解器计算的次数,m表示输出的plt文件个数,它们在main函数中初值都为零.
//出口:如果成功输出plt文件则返回1,否则返回0
//作用:通过判断(n/4==m)来决定是否输出plt文件,这里表示每4次调用求解器输出一个plt文件,如果输出了plt文件则返回1否则返回0
//在main函数中如果返回1,则	m+=Animation(U,dx,dy,n,m);	则m的值加1表示plt文件的个数加了1
int Animation(double U[Nx+7][Ny+7][4],double dx,double dy,int n,int m,int method) 
{
	int i,j,k;
	double u,v,p,rou;
	char *name;
	//对文件存放的路径进行判断
	char path1[]="Roe000.plt\0";
	char path2[]="MUSCL000.plt\0";
	char path3[]="uTVD000.plt\0";
	char path4[]="sTVD000.plt\0";
	char path5[]="NND000.plt\0";
	char path6[]="ENO000.plt\0";
	char path7[]="WENO000.plt\0";
	char path8[]="CMPT000.plt\0";
	if(method==1)name=path1;
	if(method==2)name=path2;
	if(method==3)name=path3;
	if(method==4)name=path4;
	if(method==5)name=path5;
	if(method==6)name=path6;
	if(method==7)name=path7;
	if(method==8)name=path8;
	for(i=3;i<=int(1.0/dx)+2;i++)     //为了显示正确,首先将虚拟节点赋0
		for(j=int(0.5/dy)+4;j<=int(0.5/dy)+6;j++)
			for(k=0;k<=3;k++)
				U[i][j][k]=0;

	for(i=int(2.0/dx)+4;i<=Nx+6;i++)     
		for(j=int(0.5/dy)+4;j<=int(0.5/dy)+6;j++)
			for(k=0;k<=3;k++)
				U[i][j][k]=0;

	for(i=int(1.0/dx);i<=int(1.0/dx)+2;i++)     
		for(j=int(0.5/dy)+4;j<=Ny+6;j++)
			for(k=0;k<=3;k++)
				U[i][j][k]=0;


	for(i=int(2.0/dx)+4;i<=int(2.0/dx)+6;i++)     
		for(j=int(0.5/dy)+4;j<=Ny+6;j++)
			for(k=0;k<=3;k++)
				U[i][j][k]=0;

	if(n/OPDT==m)         //其中n表示的是调用FinalHomework_TVD_Solver求解器的次数,m表示已经输出的plt文件的个数 OPDT=Output Per Dt
	{		                  //这里表示每调用OPDT次求解器的时间间隔输出一个plt文件,并返回值1,否则返回0
		if(m<100)
		{
			*(name+17)=char(m%10+48);            //对plt文件名进行操作  0对应的ASCII码是48
			*(name+16)=char(m/10+48);
		}
		if(m>=100)
		{
			*(name+17)=char(m%10+48);
			*(name+16)=char((m/10)%10+48);
			*(name+15)=char(m/100+48);
		}
		ofstream tec(name);
	    tec<<"title=\"Finalwork\"\n";
	    tec<<"variables=\"x\",\"y\",\"rou\",\"u\",\"v\",\"p\"\n";
	    tec<<"zone i="<<Nx+1<<", j="<<2*Ny+1<<"\n";
	    for(j=-Ny;j<=Ny;j++)                          //由于计算时只取了一半,这里将所有量做轴对称显示即为所有流场数据
		{
		    for(i=3;i<=Nx+3;i++)
			{		
			    if(j>=0)
				{
				    if(U[i][j+3][0]!=0)
					{
					    rou=U[i][j+3][0];
					    u=U[i][j+3][1]/U[i][j+3][0];			
				        v=U[i][j+3][2]/U[i][j+3][0];			
				        p=(GAMA-1)*(U[i][j+3][3]-rou*(u*u+v*v)/2);			
				        tec<<(i-3)*dx<<" "<<j*dy<<" "<<rou<<" "<<u<<" "<<v<<" "<<p<<"\n";
					}
			    	else tec<<(i-3)*dx<<" "<<j*dy<<" 0 0 0 0\n";
				}
			
			    if(j<0)
				{
			     	if(U[i][-j+3][0]!=0)
					{
			     		rou=U[i][-j+3][0];
				    	u=U[i][-j+3][1]/U[i][-j+3][0];			
				        v=-U[i][-j+3][2]/U[i][-j+3][0];			
				        p=(GAMA-1)*(U[i][-j+3][3]-rou*(u*u+v*v)/2);			
				        tec<<(i-3)*dx<<" "<<j*dy<<" "<<rou<<" "<<u<<" "<<v<<" "<<p<<"\n";
					}
				    else tec<<(i-3)*dx<<" "<<j*dy<<" 0 0 0 0\n";
				}
			}
		}
		return 1;       //返回值1表示成功输出了一个plt文件
	}
	else return 0;
}

int main()
{
	int m=0,n=0,method;      //n和m分别表示调用求解器的次数和输出plt文件的个数
	double dx,dy,dt=0,T;
	// system("md D:\\CFD\\Roe_");
	// system("md D:\\CFD\\MSCL");
	// system("md D:\\CFD\\uTVD");
	// system("md D:\\CFD\\sTVD");
	// system("md D:\\CFD\\NND_");
	// system("md D:\\CFD\\ENO_");
	// system("md D:\\CFD\\WENO");
	// system("md D:\\CFD\\CMPT");
	// system("cls");
	Initial(U,dx,dy);
	T=0;
	cout<<"请注意，生成的计算结果保存在D:\\CFD相应的目录下，具体见本程序最开始的说明\n\n";
	cout<<"请输入相应的数字选择方法\n\n";
	cout<<"1-Roe    2-MUSCL   3-迎风TVD   4-对称TVD   5-NND   6-ENO   7-WENO   8-紧致格式\n";
	cout<<"\n";
	cin>>method;
	if(method==1)
	{
		while(T<=TT)
		{
			dt=CFL(U,dx,dy,ROECFL);
			Roe_Solver(U,U_,LAMDA_,L_,R_,alpha_,theta,F,F_,G,G_,dx,dy,dt,a_);
			T+=dt;
			n++;              
			m+=Animation(U,dx,dy,n,m,method);	
		}
	}
	if(method==2)
	{
		while(T<=TT)
		{
			dt=CFL(U,dx,dy,MUSCLCFL);
			MUSCL_Solver(U,U_L,U_R,Fp,Fd,F_,Gp,Gd,G_,dx,dy,dt);
			T+=dt;
			n++;              
			m+=Animation(U,dx,dy,n,m,method);	
		}
	}
	if(method==3)
	{
		while(T<=TT)
		{
			dt=CFL(U,dx,dy,UPWINDTVDCFL);
			UpWindTVD_Solver(U,U_,LAMDA_,L_,R_,alpha_,g_,g,gama_,Q_,theta,F,F_,G,G_,dx,dy,dt,a_);
			T+=dt;
			n++;              
			m+=Animation(U,dx,dy,n,m,method);	
		}
	}
	if(method==4)
	{
		while(T<=TT)
		{
			dt=CFL(U,dx,dy,SYMTVDCFL);
			SymTVD_Solver(U,U_,LAMDA_,L_,R_,alpha_,g_,Q_,theta,F,F_,G,G_,dx,dy,dt,a_);
			T+=dt;
			n++;              
			m+=Animation(U,dx,dy,n,m,method);	
		}
	}
	if(method==5)
	{
		while(T<=TT)
		{
			dt=CFL(U,dx,dy,NNDCFL);
			NND_Solver(U,Fp,Fd,Gp,Gd,dx,dy,dt);
			T+=dt;
			n++;              
			m+=Animation(U,dx,dy,n,m,method);	
		}
	}
	if(method==6)
	{
		while(T<=TT)
		{
			ENO_Solver(U,U1,U2,F,Fp,Fd,F_p,F_d,F_,G,Gp,Gd,F_p,F_d,F_,LAMDA_,q3p,q3d,dx,dy,dt);
			T+=dt;
			n++;              
			m+=Animation(U,dx,dy,n,m,method);	
		}
	}
	if(method==7)
	{
		while(T<=TT)
		{
			WENO_Solver(U,U1,U2,ISp,ISd,omegap,omegad,alphap,alphad,q3p,q3d,LAMDA_,Fp,Fd,F_p,F_d,F,F_,G,G_,Gp,Gd,G_p,G_d,dx,dy,dt);
			T+=dt;
			n++;              
			m+=Animation(U,dx,dy,n,m,method);	
		}
	}
	if(method==8)
	{
		while(T<=TT)
		{
			dt=CFL(U,dx,dy,COMPACTCFL);
			Compact_Solver(U,Ut,U1,U2,Fp,Fd,fp,fd,Gp,Gd,gp,gd,dx,dy,dt,p);
			T+=dt;
			n++;              
			m+=Animation(U,dx,dy,n,m,method);	
		}
	}
	else cout<<"您的输入非法\n";
}