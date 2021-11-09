#include "Riemann2d.h"

/*-------------------------------------------------------
对称TVD格式，计算x方向的差分算子
入口：U矢量,dx,dt
出口：U_x矢量
--------------------------------------------------------*/
void LLx(double U[Jx+2][Jy+2][4],double dx,double dt)
{
	int i,j,k,l;
	double u_temp,v_temp,a_temp,h_temp,r;
	double Ratio_rou;
	double lamda[Jx+2][4],R[Jx+2][4][4],L[Jx+2][4][4];
	double Fai[Jx+2][4],alfa[Jx+2][4],g_b[Jx+2][4];
	r=dt/dx;
	Calculate_field_p();
	for(i=0;i<=Jx+1;i++)
		for(j=0;j<=Jy+1;j++)
		U2F(U[i][j],F[i][j]);
			
	for(j=0;j<=Jy+1;j++)
		{
			for(i=0;i<=Jx;i++)
			{
				Ratio_rou=sqrt(rou[i+1][j]/rou[i][j]);
				u_temp=(u[i][j]+u[i+1][j]*Ratio_rou)/(1.0+Ratio_rou);
				v_temp=(v[i][j]+v[i+1][j]*Ratio_rou)/(1.0+Ratio_rou);
				h_temp=(h[i][j]+h[i+1][j]*Ratio_rou)/(1.0+Ratio_rou);
				a_temp=sqrt((GAMA-1)*(h_temp-0.5*(u_temp*u_temp+v_temp*v_temp)));
				lamda_x(lamda[i],u_temp,v_temp,a_temp);
				LR_x(R[i],L[i],u_temp,v_temp,a_temp);
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
			for(i=0;i<=Jx;i++)
				for(k=0;k<=3;k++)
					for(l=0;l<=3;l++)	
						alfa[i][k]=alfa[i][k]+L[i][k][l]*(U[i+1][j][l]-U[i][j][l]);
		//计算g_b[Jx+2][Jy+2][4]
			for(i=1;i<=Jx;i++)
				for(k=0;k<=3;k++)
					g_b[i][k]=minmod(alfa[i-1][k],alfa[i][k],alfa[i+1][k]);
		//计算Fai[Jx+2][4]
			for(i=0;i<=Jx;i++)
				for(k=0;k<=3;k++)
					for(l=0;l<=3;l++)
					Fai[i][k]=Fai[i][k]-R[i][k][l]*(r*pow(lamda[i][l],2)*g_b[i][l]
								+Q_z(lamda[i][l])*(alfa[i][l]-g_b[i][l]));
		//计算U_x
			double fl,fr;
			for(i=1;i<=Jx;i++)
				for(k=0;k<=3;k++)
				{
					fl=0.5*(F[i][j][k]+F[i-1][j][k]+Fai[i-1][k]);
					fr=0.5*(F[i][j][k]+F[i+1][j][k]+Fai[i][k]);
					U[i][j][k]=U[i][j][k]-r*(fr-fl);
					//printf("fl=%f  fr=%f \n",fl,fr);
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
	double Ratio_rou,u_temp,v_temp,a_temp,h_temp,r;
	double lamda[Jy+2][4],R[Jy+2][4][4],L[Jy+2][4][4];
	double Fai[Jy+2][4],alfa[Jy+2][4],g_b[Jy+2][4];
	r=dt/dy;
	Calculate_field_p();
	for(i=0;i<=Jx+1;i++)
		for(j=0;j<=Jy+1;j++)
		U2G(U[i][j],G[i][j]);
			
	for(i=0;i<=Jx+1;i++)
		{
			for(j=0;j<=Jy;j++)
			{
				Ratio_rou=sqrt(rou[i][j+1]/rou[i][j]);
				u_temp=(u[i][j+1]+u[i][j]*Ratio_rou)/(1.0+Ratio_rou);
				v_temp=(v[i][j+1]+v[i][j]*Ratio_rou)/(1.0+Ratio_rou);
				h_temp=(h[i][j+1]+h[i][j]*Ratio_rou)/(1.0+Ratio_rou);
				a_temp=sqrt((GAMA-1)*(h_temp-0.5*(u_temp*u_temp+v_temp*v_temp)));
				lamda_y(lamda[j],u_temp,v_temp,a_temp);
				LR_y(R[j],L[j],u_temp,v_temp,a_temp);
			}
		//初始化Fai[Jy+2][4],alfa[Jy+2][4],g_b[Jy+2][4];
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
		//计算g_b[Jy+2][4]
			for(j=1;j<=Jy;j++)
				for(k=0;k<=3;k++)
					g_b[j][k]=minmod(alfa[j-1][k],alfa[j][k],alfa[j+1][k]);
		//计算Fai[Jx+2][4]
			for(j=0;j<=Jy+1;j++)
				for(k=0;k<=3;k++)
					for(l=0;l<=3;l++)
						Fai[j][k]=Fai[j][k]-R[j][k][l]*(r*pow(lamda[j][l],2)*g_b[j][l]
									+Q_z(lamda[j][l])*(alfa[j][l]-g_b[j][l]));
		//计算U_y
			double fl,fr;
			for(j=1;j<=Jy;j++)
				for(k=0;k<=3;k++)
				{
					fl=0.5*(G[i][j][k]+G[i][j-1][k]+Fai[j-1][k]);
					fr=0.5*(G[i][j][k]+G[i][j+1][k]+Fai[j][k]);
					U[i][j][k]=U[i][j][k]-r*(fr-fl);
					//printf("U[%d][%d][%d]=%f\n",i,j,k,U[i][j][k]);
				}
		}

}
/*-------------------------------------------------------
*对称型TVD格式
--------------------------------------------------------*/
void symmetry_TVD(double U[Jx+2][Jy+2][4],double dx,double dy,double dt)
{
	bound(U);
	LLy(U,dy,dt/2.0);
	bound(U);
	LLx(U,dx,dt/2.0);
	bound(U);
	LLy(U,dy,dt/2.0);	
	LLx(U,dx,dt/2.0);
	bound(U);	
}
