#include "Riemann2d.h"

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