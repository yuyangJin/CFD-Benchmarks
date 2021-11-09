#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

using namespace std;

#define GAMA 1.4 //气体常数
#define PI 3.1415926

#define MIN(x,y) (((x)<(y))?(x):(y))
#define MAX(x,y) (((x)>(y))?(x):(y))


#define Lx 2.0 //计算区域
#define Ly 2.0

#define TT 0.6 //总时间
#define Sf 0.2 //时间步长因子
#define epsilon 0.12 //人工粘性项小参数

#define Jx 100 
#define Jy 100 //网格数

//全局变量
double U[Jx+2][Jy+2][4],U_old[Jx+2][Jy+2][4], F[Jx+2][Jy+2][4],G[Jx+2][Jy+2][4];
double UF[Jx+2][Jy+2][4],UG[Jx+2][Jy+2][4];
double rou[Jx+2][Jy+2],u[Jx+2][Jy+2],v[Jx+2][Jy+2],h[Jx+2][Jy+2],p[Jx+2][Jy+2];


/*-------------------------------------------------
初始化
---------------------------------------------------*/
void Init(double U[Jx+2][Jy+2][4],double& dx,double& dy)
{
	int i,j;
	double rou1=1.5,p1=1.0,u1=-3.0,v1=1.0;
	double rou2=0.5,p2=1.0,u2=-3.0,v2=-1.0;
	double rou3=0.5,p3=1.0,u3=3.0,v3=1.0;
	double rou4=1.5,p4=1.0,u4=3.0,v4=-1.0;
	dx=Lx/(Jx);
	dy=Ly/(Jy);
	for(i=0;i<=Jx/2;i++)
	{
		for(j=0;j<=Jy/2;j++)
		{
			U[i][j][0]=rou3;
			U[i][j][1]=rou3*u3;
			U[i][j][2]=rou3*v3;
			U[i][j][3]=p3/(GAMA-1)+0.5*rou3*(u3*u3+v3*v3);
		}
	}
	for(i=Jx/2+1;i<=Jx+1;i++)
	{
		for(j=0;j<=Jy/2;j++)
		{
			U[i][j][0]=rou4;
			U[i][j][1]=rou4*u4;
			U[i][j][2]=rou4*v4;
			U[i][j][3]=p4/(GAMA-1)+0.5*rou4*(u4*u4+v4*v4);
		}
	}
	for(i=0;i<=Jx/2;i++)
	{
		for(j=Jy/2+1;j<=Jy+1;j++)
		{
			U[i][j][0]=rou1;
			U[i][j][1]=rou1*u1;
			U[i][j][2]=rou1*v1;
			U[i][j][3]=p1/(GAMA-1)+0.5*rou1*(u1*u1+v1*v1);
		}
	}
	for(i=Jx/2+1;i<=Jx+1;i++)
	{
		for(j=Jy/2+1;j<=Jy+1;j++)
		{
			U[i][j][0]=rou2;
			U[i][j][1]=rou2*u2;
			U[i][j][2]=rou2*v2;
			U[i][j][3]=p2/(GAMA-1)+0.5*rou2*(u2*u2+v2*v2);
		}
	}
}

/*------------------------------------------------------
边界条件
入口：dx，dy，网格宽度
出口：U，已经给定边界
--------------------------------------------------------*/
void bound(double U[Jx+2][Jy+2][4])
{
	int i,j;
	for(j=0;j<=Jx;j++)
		for(i=0;i<=3;i++)
			U[0][j][i]=U[1][j][i];//下边界
	for(j=0;j<=Jx;j++)
		for(i=0;i<=3;i++)
			U[Jy][j][i]=U[Jy-1][j][i];//上边界
	for(j=0;j<=Jx;j++)
		for(i=0;i<=3;i++)
			U[j][0][i]=U[j][1][i];//左边界
	for(j=0;j<=Jx;j++)
		for(i=0;i<=3;i++)
			U[j][Jx][i]=U[j][Jx-1][i];//右边界
}

/*---------------------------------------------------
计算时间步长
----------------------------------------------------*/
double CFL(double U[Jx+2][Jy+2][4],double dx,double dy)
{
	int i,j;
	double p,u,v,vel;
	double Mx,My,maxvel;
	Mx=1e-100;
	My=1e-100;
	maxvel=1e-100;
	for(i=1;i<=Jx;i++)
		for(j=1;j<=Jy;j++)
		{
			u=U[i][j][1]/U[i][j][0];
			v=U[i][j][2]/U[i][j][0];
			p=(GAMA-1.0)*(U[i][j][3]-0.5*U[i][j][0]*(u*u+v*v));
			vel=sqrt(GAMA*p/U[i][j][0])+sqrt(u*u+v*v);		
			if(vel>maxvel)
				maxvel=vel;			
		}

	return Sf*MIN(dx,dy)/maxvel;
}
/*-------------------------------------------------------
计算数值粘性函数Q_z
入口：lamda
出口：返回Q_z值
--------------------------------------------------------*/
double Q_z(double x)
{
	double result;
	if(fabs(x)>=epsilon)
		result=fabs(x);
	else
		result=(x*x/epsilon+epsilon)/2.0;
	return result;
}
/*-------------------------------------------------------
计算delta
入口：lamda,r=dx/dt
出口：返回delta
--------------------------------------------------------*/
double delta(double x,double r)
{
	double result;
	result=0.5*(Q_z(x)-r*x*x);
	return result;
}
/*-------------------------------------------------------
计算人工粘性修正项的minmod
--------------------------------------------------------*/
double minmod(double x,double y,double z)
{
	double result,temp;

	if((x<0&&y<0&&z<0)||(x>0&&y>0&&z>0))		
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
		if(fabs(y)>fabs(z))
		{
			temp=y;
			y=z;
			z=temp;
		}
		result=x;
	}
	else result=0.0;
	return result;	
}
/*-------------------------------------------------------
根据U计算F
---------------------------------------------------------*/
void U2F(double UF[4],double F[4])
{
	double u,v,p;
	u=UF[1]/UF[0];
	v=UF[2]/UF[0];
	p=(GAMA-1)*(UF[3]-0.5*UF[0]*(u*u+v*v));
	F[0]=UF[1];
	F[1]=UF[0]*u*u+p;
	F[2]=UF[0]*u*v;
	F[3]=(UF[3]+p)*u;
}
/*-------------------------------------------------------
根据U计算G
---------------------------------------------------------*/
void U2G(double UG[4],double G[4])
{
	double u,v,p;
	u=UG[1]/UG[0];
	v=UG[2]/UG[0];
	p=(GAMA-1)*(UG[3]-0.5*UG[0]*(u*u+v*v));
	G[0]=UG[2];
	G[1]=UG[0]*u*v;
	G[2]=UG[0]*v*v+p;
	G[3]=(UG[3]+p)*v;
}
/*------------------------------------------------------
计算流场参数
-------------------------------------------------------*/
void Calculate_field_p(double U[Jx+2][Jy+2][4],double rou[Jx+2][Jy+2],double u[Jx+2][Jy+2],double v[Jx+2][Jy+2],double h[Jx+2][Jy+2],double p[Jx+2][Jy+2])
{
	int i,j;
	for(i=0;i<Jx+1;i++)
		for(j=0;j<Jy+1;j++)
		{	
			rou[i][j]=U[i][j][0];
			u[i][j]=U[i][j][1]/rou[i][j];
			v[i][j]=U[i][j][2]/rou[i][j];			
			p[i][j]=(GAMA-1.0)*(U[i][j][3]-0.5*rou[i][j]*(u[i][j]*u[i][j]+v[i][j]*v[i][j]));
			h[i][j]=GAMA*p[i][j]/((GAMA-1)*rou[i][j])+0.5*(u[i][j]*u[i][j]+v[i][j]*v[i][j]);
		}
}
/*-------------------------------------------------------
计算x方向的特征值lamda_x
入口：U矢量,(u[i][j]+u[i+1][j])/2,(a[i][j]+a[i+1][j])/2
出口：lamda_x矢量
--------------------------------------------------------*/
void lamda_x(double lamda_x[4],double u,double v,double a)
{
	lamda_x[0]=u-a;
	lamda_x[1]=u;
	lamda_x[2]=(u+a);
	lamda_x[3]=u;
}
/*-------------------------------------------------------
计算x方向的左右特征矢量矩阵L_x，R_x
入口：U矢量,L,R,(u[i][j]+u[i+1][j])/2,(v[i][j]+v[i+1][j])/2,(a[i][j]+a[i+1][j])/2
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

	//计算左x特征矢量L
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
入口：U矢量,(u[i][j]+u[i][j+1])/2,(a[i][j]+a[i][j+1])/2
出口：lamda_y矢量
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

	//计算左y特征矢量L[i][j]
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
/*----------------------------------------------------------------
保存前一时间步的U矢量
入口：无
出口：无
----------------------------------------------------------------*/
void save_U_last()
{
	int i,j,k;
	for(i=0;i<=Jx+1;i++)
		for(j=0;j<=Jy+1;j++)
			for(k=0;k<=3;k++)
				U_old[i][j][k]=U[i][j][k];
}
void symmetry_TVD(double U[Jx+2][Jy+2][4],double dx,double dy,double dt)
{
	bound(U);
	int i,j,k,l;
	double u_temp,v_temp,a_temp,h_temp,rx,ry,rou_temp,p_temp,E_temp;
	double Ratio_rou;
	double lamda[Jx+2][4],R[Jx+2][4][4],L[Jx+2][4][4];
	double Fai[Jx+2][4],alfa[Jx+2][4],g_b[Jx+2][4];

	rx=dt/dx;
	ry=dt/dy;

	for(j=0;j<=Jy;j++)
	{
		for(i=0;i<=Jx;i++)
			{
				Ratio_rou=sqrt(rou[i+1][j]/rou[i][j]);
				u_temp=(u[i][j]+u[i+1][j]*Ratio_rou)/(1.0+Ratio_rou);
				v_temp=(v[i][j]+v[i+1][j]*Ratio_rou)/(1.0+Ratio_rou);
				h_temp=(h[i][j]+h[i+1][j]*Ratio_rou)/(1.0+Ratio_rou);
				a_temp=sqrt((GAMA-1)*(h_temp-0.5*(u_temp*u_temp+v_temp*v_temp)));
				rou_temp=rou[i][j]*(0.5*(1+Ratio_rou))*(0.5*(1+Ratio_rou));
				p_temp=rou_temp*a_temp*a_temp/GAMA;
				E_temp=rou_temp*h_temp-p_temp;
				lamda_x(lamda[i],u_temp,v_temp,a_temp);
				LR_x(R[i],L[i],u_temp,v_temp,a_temp);
				UF[i][j][0]=rou_temp;
				UF[i][j][1]=rou_temp*u_temp;
				UF[i][j][2]=rou_temp*v_temp;
				UF[i][j][3]=E_temp;
				U2F(UF[i][j],F[i][j]);
			}
		//初始化Fai[Jx+2][4],alfa[Jx+2][4],g_b[Jx+2][4];
		for(i=0;i<=Jx+1;i++)
			for(k=0;k<=3;k++)
			{
				alfa[i][k]=0.0;
				Fai[i][k]=0.0;
				g_b[i][k]=0.0;
			}
		//计算alfa[i][k]
		for(i=0;i<Jx;i++)			
			for(k=0;k<=3;k++)
				for(l=0;l<=3;l++)	
					alfa[i][k]=alfa[i][k]+L[i][k][l]*(U[i+1][j][l]-U[i][j][l]);
		//计算g_b
		for(i=1;i<Jx;i++)		
			for(k=0;k<=3;k++)
				g_b[i][k]=minmod(alfa[i-1][k],alfa[i][k],alfa[i+1][k]);
		//计算Fai
		for(i=0;i<Jx;i++)		
			for(k=0;k<=3;k++)
				for(l=0;l<=3;l++)
					Fai[i][k]=Fai[i][k]-R[i][k][l]*(pow(rx*lamda[i][l],2)*g_b[i][l]+Q_z(rx*lamda[i][l])*(alfa[i][l]-g_b[i][l]))/rx;

		//计算修正后的F
		for(i=0;i<=Jx;i++)
			for(k=0;k<=3;k++)
			{
				F[i][j][k]=F[i][j][k]+0.5*Fai[i][k];
			}
	}


	for(i=0;i<=Jx;i++)
	{
		for(j=0;j<=Jy;j++)
			{
				Ratio_rou=sqrt(rou[i][j+1]/rou[i][j]);
				u_temp=(u[i][j+1]+u[i][j]*Ratio_rou)/(1.0+Ratio_rou);
				v_temp=(v[i][j+1]+v[i][j]*Ratio_rou)/(1.0+Ratio_rou);
				h_temp=(h[i][j+1]+h[i][j]*Ratio_rou)/(1.0+Ratio_rou);
				a_temp=sqrt((GAMA-1)*(h_temp-0.5*(u_temp*u_temp+v_temp*v_temp)));
				rou_temp=rou[i][j]*(0.5*(1+Ratio_rou))*(0.5*(1+Ratio_rou));
				p_temp=rou_temp*a_temp*a_temp/GAMA;
				E_temp=rou_temp*h_temp-p_temp;
				lamda_y(lamda[j],u_temp,v_temp,a_temp);
				LR_y(R[j],L[j],u_temp,v_temp,a_temp);
				UG[i][j][0]=rou_temp;
				UG[i][j][1]=rou_temp*u_temp;
				UG[i][j][2]=rou_temp*v_temp;
				UG[i][j][3]=E_temp;
				U2G(UG[i][j],G[i][j]);
			}
		for(j=0;j<=Jy+1;j++)
			for(k=0;k<=3;k++)
			{
				alfa[j][k]=0.0;
				Fai[j][k]=0.0;
				g_b[j][k]=0.0;
			}
		//计算alfa[j][k]
		for(j=0;j<=Jy;j++)
			for(k=0;k<=3;k++)
				for(l=0;l<=3;l++)
					alfa[j][k]=alfa[j][k]+L[j][k][l]*(U[i][j+1][l]-U[i][j][l]);	
		//计算g_b[Jy+1][4]
		for(j=1;j<=Jy;j++)			
			for(k=0;k<=3;k++)
				g_b[j][k]=minmod(alfa[j-1][k],alfa[j][k],alfa[j+1][k]);
		//计算Fai[Jy+1][4]
		for(j=0;j<=Jy;j++)
			for(k=0;k<=3;k++)
				for(l=0;l<=3;l++)
					Fai[j][k]=Fai[j][k]-R[j][k][l]*(pow(ry*lamda[j][l],2)*g_b[j][l]
								+Q_z(ry*lamda[j][l])*(alfa[j][l]-g_b[j][l]))/ry;
		//计算修正后的G
		for(j=0;j<=Jy;j++)
			for(k=0;k<=3;k++)
			{
				G[i][j][k]=G[i][j][k]+0.5*Fai[j][k];				
			}
	}
	for(i=0;i<=Jx;i++)
		for(j=0;j<=Jy;j++)
			for(k=0;k<=3;k++)
				U[i][j][k]=U[i][j][k]-rx*(F[i][j][k]-F[i-1][j][k])-ry*(G[i][j][k]-G[i][j-1][k]);//计算U
	bound(U);

}
/*----------------------------------------------------------------
计算一时间步的当前U矢量的相对误差
入口：无
出口：无
----------------------------------------------------------------*/
double error(double U1[Jx+2][Jy+2][4],double U2[Jx+2][Jy+2][4])
{
	int i,j,k;
	double dd_variable;
	double err=1e-10;
	for(i=1;i<=Jx;i++)
		for(j=1;j<=Jy;j++)
			for(k=0;k<=3;k++)
		{	
			if(U2[i][j][k]!=0)
				dd_variable=fabs((U2[i][j][k]-U1[i][j][k])/U2[i][j][k]);
			else
				dd_variable=0;
			if(dd_variable>err)err=dd_variable;
		}
	
	return err;
}
void Output(double U[Jx+2][Jy+2][4],double dx,double dy,double T)
{
	int i,j;
	double rou,u,v,p,a,Ma;

	FILE *fp1,*fp2;
	char head1[200]={"TITLE	=\"Dataset\"\nVARIABLES =\"X\"\n\"y\"\n\"rou\"\n\"u\"\n\"v\"\n\"p\"\n\"Ma\"\n"};
	char head2[200]={"ZONE T=\"Rectangular zone1\"\n"};
	char head3[200]={"ZONETYPE=Ordered\nDATAPACKING=POINT\nDT=(SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE )"};
	char filename1[100],filename_tp[256]={"result_T_"},filename2[10]={".plt"},filename_tx[256]={"result_T_"},filename3[10]={".txt"};
	
	gcvt(T,3,filename1); 
	strcat(filename_tp,filename1);
	strcat(filename_tp,filename2);
	strcat(filename_tx,filename1);
	strcat(filename_tx,filename3);

	fp1=fopen(filename_tp,"w");

	fprintf(fp1,"%s%s\nI=%d,J=%d,k=%d,%s\n",head1,head2,Jy+1,Jx+1,1,head3);

		
	for(i=0;i<=Jx;i++)
		for(j=0;j<=Jy;j++)
		{
			rou=U[i][j][0];
			u=U[i][j][1]/rou;
			v=U[i][j][2]/rou;
			p=(GAMA-1)*(U[i][j][3]-0.5*rou*(u*u+v*v));
			a=sqrt(GAMA*p/rou);
			Ma=sqrt(u*u+v*v)/a;
			fprintf(fp1,"%20f\t%20.10e\t%20.10e\t%20.10e\t%20.10e\t%20.10e\t%20.10e\n",(i-2)*dx,(j-2)*dy,rou,u,v,p,Ma);
		}	
	fclose(fp1);

	fp2=fopen("result.txt","w");
	j=Jy/2;
	i=Jx/2;
	for(i=1;i<=Jx;i++)
	//for(j=0;j<=Jy+1;j++)
	{
		rou=U[i][j][0];
		u=U[i][j][1]/rou;
		v=U[i][j][2]/rou;
		p=(GAMA-1)*(U[i][j][3]-0.5*rou*(u*u+v*v));
		a=sqrt(GAMA*p/rou);
		fprintf(fp2,"%20f\t%20.10e\t%20.10e\t%20.10e\t%20.10e\t%20.10e\n",i*dx,rou,u,v,p,a);

	}
	fclose(fp2);
}
int main()
{
	double T,dx,dy,dt,err;
	double n=0;
	err=1e-7;
	Init(U,dx,dy);
	T=0;
	while(T<TT&&err>1e-8)
	{
		n+=1;
		Calculate_field_p(U,rou,u,v,h,p);
		dt=CFL(U,dx,dy);
		T+=dt;
		printf("T=%10g	dt=%10g	err=%f\n",T,dt,err);
		if(fmod(n,10)==0)
		{
			printf("n=%10f\n",n);
			Output(U,dx,dy,T);
		}
		save_U_last();
		symmetry_TVD(U,dx,dy,dt);
		err=error(U_old,U);
	}
	Output(U,dx,dy,T);
}
