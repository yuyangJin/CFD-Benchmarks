// Finite Volume Symmetric TVD Based on Operator Splitting Method
#include "Riemann2d.h"

/*-------------------------------------------------------
对称TVD格式，计算x方向的差分算子
入口：U矢量,dx,dt
出口：U_x矢量
--------------------------------------------------------*/
void LLx(double U[Jx+2][Jy+2][4],double dx,double dt)
{
	int i,j,k,l;
	double u1,v1,p1,h1,u2,v2,p2,h2;
	double u_temp,v_temp,a_temp,h_temp,r,rou_temp,p_temp,E_temp;
	double Ratio_rou;
	double lamda[Jx+2][4],R[Jx+2][4][4],L[Jx+2][4][4];
	double Fai[Jx+2][4],alfa[Jx+2][4],g_b[Jx+2][4],U_temp[Jx+2][Jy+2][4];
	r=dt/dx;
	for(j=0;j<=Jy+1;j++)
	{
		for(i=0;i<=Jx;i++)
		{
			u1=U[i][j][1]/U[i][j][0];
			v1=U[i][j][2]/U[i][j][0];
			p1=(GAMA-1)*(U[i][j][3]-0.5*U[i][j][0]*(u1*u1+v1*v1));
			h1=GAMA*p1/((GAMA-1)*U[i][j][0])+0.5*(u1*u1+v1*v1);

			u2=U[i+1][j][1]/U[i+1][j][0];
			v2=U[i+1][j][2]/U[i+1][j][0];
			p2=(GAMA-1)*(U[i+1][j][3]-0.5*U[i+1][j][0]*(u2*u2+v2*v2));
			h2=GAMA*p2/((GAMA-1)*U[i+1][j][0])+0.5*(u2*u2+v2*v2);
			//这里计算节点上状态参量的值有重复，影响计算时间，但是不知道为什么统一计算状态参量的时候总是出问题。
			Ratio_rou=sqrt(U[i+1][j][0]/U[i][j][0]);
			u_temp=(u1+u2*Ratio_rou)/(1.0+Ratio_rou);
			v_temp=(v1+v2*Ratio_rou)/(1.0+Ratio_rou);
			rou_temp=U[i][j][0]*(1+Ratio_rou)*(1+Ratio_rou)/4;
			h_temp=(h1+h2*Ratio_rou)/(1.0+Ratio_rou);			
			a_temp=sqrt(fabs((GAMA-1)*(h_temp-0.5*(u_temp*u_temp+v_temp*v_temp))));
			p_temp=a_temp*a_temp*rou_temp/GAMA;
			E_temp=rou_temp*h_temp-p_temp;
			
			U_temp[i][j][0]=rou_temp;
			U_temp[i][j][1]=rou_temp*u_temp;
			U_temp[i][j][2]=rou_temp*v_temp;
			U_temp[i][j][3]=E_temp;
			lamda_x(lamda[i],u_temp,v_temp,a_temp);
			LR_x(R[i],L[i],u_temp,v_temp,a_temp);
			U2F(U_temp[i][j],F[i][j]);
		}
		//初始化Fai,alfa,g_b;
		for(i=0;i<=Jx+1;i++)
			for(k=0;k<=3;k++)
			{
				alfa[i][k]=0.0;
				Fai[i][k]=0.0;
				g_b[i][k]=0.0;
			}
		//计算alfa
		for(i=0;i<=Jx;i++)			
			for(k=0;k<=3;k++)
				for(l=0;l<=3;l++)	
					alfa[i][k]=alfa[i][k]+L[i][k][l]*(U[i+1][j][l]-U[i][j][l]);
		//计算g_b
		for(i=1;i<=Jx;i++)		
			for(k=0;k<=3;k++)
				g_b[i][k]=minmod(alfa[i-1][k],alfa[i][k],alfa[i+1][k]);		
		//计算Fai
		for(i=0;i<=Jx;i++)		
			for(k=0;k<=3;k++)
				for(l=0;l<=3;l++)
Fai[i][k]=Fai[i][k]-R[i][k][l]*(pow(r*lamda[i][l],2)*g_b[i][l]+Q_z(r*lamda[i][l])*(alfa[i][l]-g_b[i][l]))/r;
		//计算U_x
		double fl,fr;
		for(i=2;i<=Jx;i++)
			for(k=0;k<=3;k++)
			{
			    fl=F[i-1][j][k]+0.5*Fai[i-1][k];
				fr=F[i][j][k]+0.5*Fai[i][k];
				U[i][j][k]=U[i][j][k]-r*(fr-fl);
			}		
	}

}
/*-------------------------------------------------------
对称TVD格式，计算y方向的差分算子
入口：U矢量,dy,dt
出口：U_y矢量
--------------------------------------------------------*/
void LLy(double U[Jx+2][Jy+2][4],double dy,double dt)
{
	int i,j,k,l;
	double u1,v1,p1,h1,u2,v2,p2,h2;
	double Ratio_rou,u_temp,v_temp,a_temp,h_temp,r,rou_temp,p_temp,E_temp;
	double lamda[Jy+2][4],R[Jy+2][4][4],L[Jy+2][4][4];
	double Fai[Jy+2][4],alfa[Jy+2][4],g_b[Jy+2][4],U_temp[Jx+2][Jy+2][4];
	r=dt/dy;
	for(i=0;i<=Jx+1;i++)
	{
		for(j=0;j<=Jy;j++)
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
			rou_temp=U[i][j][0]*(1+Ratio_rou)*(1+Ratio_rou)/4;
			h_temp=(h1+h2*Ratio_rou)/(1.0+Ratio_rou);
			a_temp=sqrt(fabs((GAMA-1)*(h_temp-0.5*(u_temp*u_temp+v_temp*v_temp))));
			p_temp=a_temp*a_temp*rou_temp/GAMA;
			E_temp=rou_temp*h_temp-p_temp;
			U_temp[i][j][0]=rou_temp;
			U_temp[i][j][1]=rou_temp*u_temp;
			U_temp[i][j][2]=rou_temp*v_temp;
			U_temp[i][j][3]=E_temp;
			lamda_y(lamda[j],u_temp,v_temp,a_temp);
			LR_y(R[j],L[j],u_temp,v_temp,a_temp);
			U2G(U_temp[i][j],G[i][j]);				
		}	
		//初始化Fai[Jy+5][4],alfa[Jy+5][4],g_b[Jy+5][4];
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
		//计算g_b[Jy+5][4]
		for(j=1;j<=Jy;j++)			
			for(k=0;k<=3;k++)
				g_b[j][k]=minmod(alfa[j-1][k],alfa[j][k],alfa[j+1][k]);
		//计算Fai[Jy+5][4]，这里Fai的含义是R*Fai
		for(j=0;j<=Jy;j++)
			for(k=0;k<=3;k++)
				for(l=0;l<=3;l++)
					Fai[j][k]=Fai[j][k]-R[j][k][l]*(pow(r*lamda[j][l],2)*g_b[j][l]
								+Q_z(r*lamda[j][l])*(alfa[j][l]-g_b[j][l]))/r;
		//计算U_y
		double fl,fr;
		for(j=2;j<=Jy;j++)
			for(k=0;k<=3;k++)
			{
				fl=G[i][j-1][k]+0.5*Fai[j-1][k];
				fr=G[i][j][k]+0.5*Fai[j][k];
				U[i][j][k]=U[i][j][k]-r*(fr-fl);				
			}
	}
}
/*-------------------------------------------------------
*对称型TVD格式
--------------------------------------------------------*/
void symmetry_TVD(double U[Jx+2][Jy+2][4],double dx,double dy,double dt)
{
	bound(U);	
	LLx(U,dx,dt/2.0);
	bound(U);
	LLy(U,dy,dt/2.0);
	bound(U);
	LLy(U,dy,dt/2.0);
	bound(U);
	LLx(U,dx,dt/2.0);
}
