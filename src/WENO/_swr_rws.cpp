/****************************************************
可改动量：1.网格数Jx、Jy，在全局变量中改动
          2.马赫数Ms，在#define中改动
		  3.输出plt文件时间间隔T0，在#define中改动
****************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <cmath>

using namespace std;

#define GAMA 1.4 //气体常数
#define MIN(x,y) (((x)<(y))?(x):(y))
#define Lx 1.0 //计算区域
#define Ly 1.0
#define TT 0.25 //总时间
#define Sf 0.8 //时间步长因子
#define T0 0.01 //输出plt文件时间间隔
#define Jx 400 //网格数
#define Jy 400                                                   
#define P 3  //模板不光滑性的放大系数
#define Ms 4.0
//全局变量
double U[Jx+7][Jy+7][4],U1[Jx+7][Jy+7][4],U2[Jx+7][Jy+7][4],U3[Jx+7][Jy+7][4];
double lamda[Jx+7][Jy+7][4],lamdab[Jx+7][Jy+7][4];
double miu[Jx+7][Jy+7][4],miub[Jx+7][Jy+7][4];
double GP[Jx+7][Jy+7][4],GN[Jx+7][Jy+7][4];
double FP[Jx+7][Jy+7][4],FN[Jx+7][Jy+7][4];
double ISP[Jx+7][Jy+7][4][3],ISN[Jx+7][Jy+7][4][3];
double alphaP[Jx+7][Jy+7][4][3],alphaN[Jx+7][Jy+7][4][3];
double omegaP[Jx+7][Jy+7][4][3],omegaN[Jx+7][Jy+7][4][3];
double q3P[Jx+7][Jy+7][4][3],q3N[Jx+7][Jy+7][4][3];
double F_P[Jx+7][Jy+7][4],F_N[Jx+7][Jy+7][4];
double G_P[Jx+7][Jy+7][4],G_N[Jx+7][Jy+7][4];
double F_[Jx+7][Jy+7][4],G_[Jx+7][Jy+7][4];
//double Ms;

/*-------------------------------------------------------
计算时间步长
入口: U， 当前物理量，dx、dy， 网格宽度；
返回: 时间步长。
---------------------------------------------------------*/
double CFL(double Ut[Jx+7][Jy+7][4],double dx,double dy)
{
 int i,j,jx;
 double maxvel,p,u,v,vel;
 maxvel=1.0e-100;
 jx=(Jx*15)/100;
	for(i=3;i<=Jx+3;i++)                            
		for(j=3;j<=Jy+3;j++)
		{
			if(i<=j+jx)         				
              {
                  u=Ut[i][j][1]/Ut[i][j][0];
                  v=Ut[i][j][2]/Ut[i][j][0];
                  p=(GAMA-1.0)*(Ut[i][j][3]-0.5*Ut[i][j][0]*(u*u+v*v));  
                  vel=sqrt(GAMA*p/Ut[i][j][0])+sqrt(u*u+v*v);
                  if(vel>maxvel)maxvel=vel;
               }
		}
 return Sf*MIN(dx,dy)/maxvel;
}

/*-----------------------------------------------------------------
初始化
入口: 无；
出口: U， 已经给定的初始值，dx、dy， 网格宽度。
--------------------------------------------------------------------*/
void Init(double dx)
{
 int i,j;
 double rou1=1.0, u1=0, v1=0, p1=0.71429; 
 double rou2, u2, v2, p2;

rou2=rou1*(GAMA+1)*Ms*Ms/(2.0+(GAMA-1)*Ms*Ms);       
u2=2.0/(GAMA+1)*(Ms-1.0/Ms);                           
v2=0;
p2=p1*(2.0*GAMA/(GAMA+1)*Ms*Ms-(GAMA-1.0)/(GAMA+1.0));  
 
 for(j=0;j<=Jy+6;j++)
 for(i=0;i<=Jx+6;i++)
  {
    if(i*dx<=0.03)   
   {
    U[i][j][0]=rou2;
    U[i][j][1]=rou2*u2;
    U[i][j][2]=rou2*v2;
    U[i][j][3]=p2/(GAMA-1.0)+rou2*(u2*u2+v2*v2)/2.0;
   }
   else
   {
    U[i][j][0]=rou1;
    U[i][j][1]=rou1*u1;
    U[i][j][2]=rou1*v1;
	U[i][j][3]=p1/(GAMA-1.0)+rou1*(u1*u1+v1*v1)/2.0;  
   }
  }
}

/*-------------------------------------------------------
边界条件
入口: dx、dy， 网格宽度；
出口: U， 已经给定边界。
---------------------------------------------------------*/
void bound(double Ut[Jx+7][Jy+7][4])
{
 int i,j,k,jx;
 double rou1=1.0, u1=0, v1=0 ,p1=0.71429; 
 double rou2, u2, v2, p2;
 jx=(Jx*15)/100;

rou2=rou1*(GAMA+1)*Ms*Ms/(2.0+(GAMA-1)*Ms*Ms);     
u2=2.0/(GAMA+1)*(Ms-1.0/Ms);                           
v2=0;
p2=p1*(2.0*GAMA/(GAMA+1)*Ms*Ms-(GAMA-1.0)/(GAMA+1.0)); 
 
//左边界波后参数
for(i=0;i<=3;i++)
    for(j=3;j<=Jy+3;j++)
    {
      Ut[i][j][0]=rou2;
      Ut[i][j][1]=rou2*u2;
      Ut[i][j][2]=rou2*v2;
      Ut[i][j][3]=p2/(GAMA-1.0)+rou2*(u2*u2+v2*v2)/2.0;
    }

//下边界刚性壁面
 for(i=3;i<=3+jx;i++)
 {
	  for(j=0;j<=2;j++)		 
         {
          Ut[i][j][0]=Ut[i][6-j][0];      
          Ut[i][j][1]=Ut[i][6-j][1];      
          Ut[i][j][2]=-Ut[i][6-j][2];     
          Ut[i][j][3]=Ut[i][6-j][3];     
		 }
 }
		
for(i=3;i<=3+jx;i++)
     {	        
          Ut[i][3][0]=Ut[i][4][0];      
          Ut[i][3][1]=Ut[i][4][1];     
          Ut[i][3][2]=0;               
          Ut[i][3][3]=Ut[i][4][3];     	 
     }
	 

//斜面刚性壁面
for(i=4+jx;i<=Jx+4;i++)
    for(j=3;j<=Jy+3-jx;j++)   
		if(i==(j+jx+1))
		  {			 
			Ut[i][j][1]=Ut[i-1][j+1][2];			  
			Ut[i][j][2]=Ut[i-1][j+1][1]; 			  
			Ut[i][j][0]=Ut[i-1][j+1][0];			  
			Ut[i][j][3]=Ut[i-1][j+1][3]; 
		  } 
	 
	 for(i=4+jx;i<=Jx+5;i++)
		  for(j=2;j<=Jy+3-jx;j++)  
		  if(i==(j+jx+2))
		  {
			  Ut[i][j][1]=Ut[i-2][j+2][2];
			  Ut[i][j][2]=Ut[i-2][j+2][1]; 
			  Ut[i][j][0]=Ut[i-2][j+2][0];
			  Ut[i][j][3]=Ut[i-2][j+2][3];  
		  } 
	
	 for(i=4+jx;i<=Jx+6;i++)
		  for(j=1;j<=Jy+3-jx;j++)   
		  if(i==(j+jx+3))
		  {
			  Ut[i][j][1]=Ut[i-3][j+3][2];
			  Ut[i][j][2]=Ut[i-3][j+3][1]; 
			  Ut[i][j][0]=Ut[i-3][j+3][0];
			  Ut[i][j][3]=Ut[i-3][j+3][3]; 
		  } 

	 for(i=4+jx;i<=Jx+3;i++)
	     for(j=4;j<=Jy+3-jx;j++)       
		    if(i==(j+jx))
			{          
				Ut[i][j][0]=(Ut[i][j+1][0]+Ut[i-1][j][0])/2.0;     
				Ut[i][j][2]=(Ut[i][j+1][2]+Ut[i-1][j][1])/2.0;               
				Ut[i][j][1]=(Ut[i][j+1][2]+Ut[i-1][j][1])/2.0;               
				Ut[i][j][3]=(Ut[i][j+1][3]+Ut[i-1][j][3])/2.0;     
			}
		
//右边界自由输出
 	 for(i=Jx+4;i<=Jx+6;i++)
	     for(j=Jy+4-jx;j<=Jy+3;j++)  
             for(k=0;k<=3;k++)
                  { 
                     Ut[i][j][k]=Ut[Jx+3][j][k];     
                  }
 
//上边界自由输出
 for(i=3;i<=Jx+3;i++)
	 for(j=Jy+4;j<=Jy+6;j++)
		   for(k=0;k<=3;k++)
		   {
			   Ut[i][j][k]=Ut[i][Jy+3][k];
		   }
}

/*-------------------------------------------------------
由U求F的分裂后的量，FP和FN
入口: U，上一时刻U矢量
     dx，x方向的网格宽度，dt， 时间步长；
出口: FP和FN，计算得到的F分裂后的量。
---------------------------------------------------------*/
void FPFN(double Ut[Jx+7][Jy+7][4])
{
 int i,j,k;
 double delta,u,v,p,rou,E,a,h;
 delta=1.0e-10;

//正通量FP
for(i=0;i<=Jx+6;i++)
    for(j=0;j<=Jy+6;j++)
		{
			u=Ut[i][j][1]/Ut[i][j][0];
			v=Ut[i][j][2]/Ut[i][j][0];
			E=Ut[i][j][3];
			rou=Ut[i][j][0];
			p=(GAMA-1.0)*(E-0.5*rou*(u*u+v*v));
			a=sqrt(GAMA*p/rou);
			h=(p+E)/rou;          

			lamda[i][j][0]=u;                 
            lamda[i][j][1]=u;                                        
			lamda[i][j][2]=u-a;
            lamda[i][j][3]=u+a;

			for(k=0;k<=3;k++)
				lamdab[i][j][k]=0.5*(lamda[i][j][k]+sqrt(lamda[i][j][k]*lamda[i][j][k]+delta*delta));    

            FP[i][j][0]=rou/(2.0*GAMA)*(2.0*(GAMA-1.0)*lamdab[i][j][0]+lamdab[i][j][2]+lamdab[i][j][3]);
            FP[i][j][1]=rou/(2.0*GAMA)*(2.0*(GAMA-1.0)*u*lamdab[i][j][0]+(u-a)*lamdab[i][j][2]+(u+a)*lamdab[i][j][3]);
			FP[i][j][2]=rou/(2.0*GAMA)*(2.0*(GAMA-1.0)*v*lamdab[i][j][0]+v*lamdab[i][j][2]+v*lamdab[i][j][3]);
			FP[i][j][3]=rou/(2.0*GAMA)*((GAMA-1.0)*(u*u+v*v)*lamdab[i][j][0]+(h-a*u)*lamdab[i][j][2]+(h+a*u)*lamdab[i][j][3]);  
//负通量FN
			for(k=0;k<=3;k++)
				lamdab[i][j][k]=0.5*(lamda[i][j][k]-sqrt(lamda[i][j][k]*lamda[i][j][k]+delta*delta)); 

			FN[i][j][0]=rou/(2.0*GAMA)*(2.0*(GAMA-1.0)*lamdab[i][j][0]+lamdab[i][j][2]+lamdab[i][j][3]);
            FN[i][j][1]=rou/(2.0*GAMA)*(2.0*(GAMA-1.0)*u*lamdab[i][j][0]+(u-a)*lamdab[i][j][2]+(u+a)*lamdab[i][j][3]);
			FN[i][j][2]=rou/(2.0*GAMA)*(2.0*(GAMA-1.0)*v*lamdab[i][j][0]+v*lamdab[i][j][2]+v*lamdab[i][j][3]);
			FN[i][j][3]=rou/(2.0*GAMA)*((GAMA-1.0)*(u*u+v*v)*lamdab[i][j][0]+(h-a*u)*lamdab[i][j][2]+(h+a*u)*lamdab[i][j][3]);        
		} 
 }

/*-------------------------------------------------------
由U求G的分裂后的量，GP和GN
入口: U，上一时刻U矢量
     dy，y方向的网格宽度，dt，时间步长；
出口: GP和GN，计算得到的G分裂后的量。
---------------------------------------------------------*/
void GPGN(double Ut[Jx+7][Jy+7][4])
{
 int i,j,k;
 double delta,u,v,p,rou,E,a,h;
 delta=1.0e-10;

//正通量GP
for(i=0;i<=Jx+6;i++)
    for(j=0;j<=Jy+6;j++)   
		{
			u=Ut[i][j][1]/Ut[i][j][0];
			v=Ut[i][j][2]/Ut[i][j][0];
			E=Ut[i][j][3];
			rou=Ut[i][j][0];
			p=(GAMA-1.0)*(E-0.5*rou*(u*u+v*v));
			a=sqrt(GAMA*p/rou);
			h=(p+E)/rou;                  

			miu[i][j][0]=v;                 
            miu[i][j][1]=v;                                        
			miu[i][j][2]=v-a;
            miu[i][j][3]=v+a;

			for(k=0;k<=3;k++)
				miub[i][j][k]=0.5*(miu[i][j][k]+sqrt(miu[i][j][k]*miu[i][j][k]+delta*delta)); 

            GP[i][j][0]=rou/(2.0*GAMA)*(2.0*(GAMA-1.0)*miub[i][j][0]+miub[i][j][2]+miub[i][j][3]);
            GP[i][j][1]=rou/(2.0*GAMA)*(2.0*(GAMA-1.0)*u*miub[i][j][0]+u*miub[i][j][2]+u*miub[i][j][3]);
			GP[i][j][2]=rou/(2.0*GAMA)*(2.0*(GAMA-1.0)*v*miub[i][j][0]+(v-a)*miub[i][j][2]+(v+a)*miub[i][j][3]);
			GP[i][j][3]=rou/(2.0*GAMA)*((GAMA-1.0)*(u*u+v*v)*miub[i][j][0]+(h-a*u)*miub[i][j][2]+(h+a*u)*miub[i][j][3]);  
//负通量GN
			for(k=0;k<=3;k++)
				miub[i][j][k]=0.5*(miu[i][j][k]-sqrt(miu[i][j][k]*miu[i][j][k]+delta*delta)); 

            GN[i][j][0]=rou/(2.0*GAMA)*(2.0*(GAMA-1.0)*miub[i][j][0]+miub[i][j][2]+miub[i][j][3]);
            GN[i][j][1]=rou/(2.0*GAMA)*(2.0*(GAMA-1.0)*u*miub[i][j][0]+u*miub[i][j][2]+u*miub[i][j][3]);
			GN[i][j][2]=rou/(2.0*GAMA)*(2.0*(GAMA-1.0)*v*miub[i][j][0]+(v-a)*miub[i][j][2]+(v+a)*miub[i][j][3]);
			GN[i][j][3]=rou/(2.0*GAMA)*((GAMA-1.0)*(u*u+v*v)*miub[i][j][0]+(h-a*u)*miub[i][j][2]+(h+a*u)*miub[i][j][3]);  
		}  
 }

void WENO_x()     
{   int i,j,k,s,jx;
	double delta=1.0e-10;
	jx=(Jx*15)/100;
	
//计算ISPN
	for(i=2;i<=Jx+3;i++)                            
		for(j=2;j<=Jy+3;j++)
		{
			if(i<=j+jx)   
				for(k=0;k<=3;k++)
				{
					ISP[i][j][k][0]=13.0/12.0*pow(FP[i-2][j][k]-2.0*FP[i-1][j][k]+FP[i][j][k],2)+0.25*pow(FP[i-2][j][k]-4.0*FP[i-1][j][k]+3.0*FP[i][j][k],2);
					ISP[i][j][k][1]=13.0/12.0*pow(FP[i-1][j][k]-2.0*FP[i][j][k]+FP[i+1][j][k],2)+0.25*pow(FP[i-1][j][k]-FP[i+1][j][k],2);
					ISP[i][j][k][2]=13.0/12.0*pow(FP[i][j][k]-2.0*FP[i+1][j][k]+FP[i+2][j][k],2)+0.25*pow(3.0*FP[i][j][k]-4.0*FP[i+1][j][k]+FP[i+2][j][k],2);

					ISN[i][j][k][0]=13.0/12.0*pow(FN[i+1][j][k]-2.0*FN[i][j][k]+FN[i-1][j][k],2)+0.25*pow(3.0*FN[i+1][j][k]-4.0*FN[i][j][k]+FN[i-1][j][k],2);
					ISN[i][j][k][1]=13.0/12.0*pow(FN[i+2][j][k]-2.0*FN[i+1][j][k]+FN[i][j][k],2)+0.25*pow(FN[i+2][j][k]-FN[i][j][k],2);
					ISN[i][j][k][2]=13.0/12.0*pow(FN[i+3][j][k]-2.0*FN[i+2][j][k]+FN[i+1][j][k],2)+0.25*pow(FN[i+3][j][k]-4.0*FN[i+2][j][k]+3.0*FN[i+1][j][k],2);
				}
		}

//计算alphaPN以及omegaPN
	for(i=2;i<=Jx+3;i++)                            
		for(j=2;j<=Jy+3;j++)
		{
			if(i<=j+jx)        
				for(k=0;k<=3;k++)
				{
					alphaP[i][j][k][0]=0.1/pow(delta+ISP[i][j][k][0],P);
					alphaP[i][j][k][1]=0.6/pow(delta+ISP[i][j][k][1],P);
					alphaP[i][j][k][2]=0.3/pow(delta+ISP[i][j][k][2],P);

					alphaN[i][j][k][0]=0.3/pow(delta+ISN[i][j][k][0],P);
					alphaN[i][j][k][1]=0.6/pow(delta+ISN[i][j][k][1],P);
					alphaN[i][j][k][2]=0.1/pow(delta+ISN[i][j][k][2],P);
					for(s=0;s<=2;s++)
					{
						omegaP[i][j][k][s]=alphaP[i][j][k][s]/(alphaP[i][j][k][0]+alphaP[i][j][k][1]+alphaP[i][j][k][2]);
						omegaN[i][j][k][s]=alphaN[i][j][k][s]/(alphaN[i][j][k][0]+alphaN[i][j][k][1]+alphaN[i][j][k][2]);						
					}
				}
		}

	//计算q3PN
	for(i=2;i<=Jx+3;i++)                            
		for(j=2;j<=Jy+3;j++)
		{
			if(i<=j+jx)        
				for(k=0;k<=3;k++)
				{
					q3P[i][j][k][0]=1.0/3.0*FP[i-2][j][k]-7.0/6.0*FP[i-1][j][k]+11.0/6.0*FP[i][j][k];
					q3P[i][j][k][1]=-1.0/6.0*FP[i-1][j][k]+5.0/6.0*FP[i][j][k]+1.0/3.0*FP[i+1][j][k];
					q3P[i][j][k][2]=1.0/3.0*FP[i][j][k]+5.0/6.0*FP[i+1][j][k]-1.0/6.0*FP[i+2][j][k];

					q3N[i][j][k][0]=-1.0/6.0*FN[i-1][j][k]+5.0/6.0*FN[i][j][k]+1.0/3.0*FN[i+1][j][k];
					q3N[i][j][k][1]=1.0/3.0*FN[i][j][k]+5.0/6.0*FN[i+1][j][k]-1.0/6.0*FN[i+2][j][k];
					q3N[i][j][k][2]=11.0/6.0*FN[i+1][j][k]-7.0/6.0*FN[i+2][j][k]+1.0/3.0*FN[i+3][j][k];
				}
		}
	
	//计算F_PN
	for(i=2;i<=Jx+3;i++)                            
		for(j=2;j<=Jy+3;j++)
		{
			if(i<=j+jx)       
				for(k=0;k<=3;k++)
				{
					F_P[i][j][k]=0;
					F_N[i][j][k]=0;
				}
			
		}
	
	for(i=2;i<=Jx+3;i++)                            
		for(j=2;j<=Jy+3;j++)
		{
			if(i<=j+jx)      
				for(k=0;k<=3;k++)
				{
					for(s=0;s<=2;s++)
					{
						F_P[i][j][k]+=omegaP[i][j][k][s]*q3P[i][j][k][s];
						F_N[i][j][k]+=omegaN[i][j][k][s]*q3N[i][j][k][s];
					}
				}
			
		}

	//计算F_
	for(i=2;i<=Jx+3;i++)                            
		for(j=2;j<=Jy+3;j++)
		{
			if(i<=j+jx)        
				for(k=0;k<=3;k++)
				{
					F_[i][j][k]=F_P[i][j][k]+F_N[i][j][k];
				}
		}
}

void WENO_y()       
{
	int i,j,k,s,jx;
	double delta=1.0e-10;     
	jx=(Jx*15)/100;
	
//计算ISPN
	for(i=2;i<=Jx+3;i++)                            
		for(j=2;j<=Jy+3;j++)
		{
			if(i<=j+jx)       
				for(k=0;k<=3;k++)
				{
					ISP[i][j][k][0]=13.0/12.0*pow(GP[i][j-2][k]-2.0*GP[i][j-1][k]+GP[i][j][k],2)+0.25*pow(GP[i][j-2][k]-4.0*GP[i][j-1][k]+3.0*GP[i][j][k],2);
					ISP[i][j][k][1]=13.0/12.0*pow(GP[i][j-1][k]-2.0*GP[i][j][k]+GP[i][j+1][k],2)+0.25*pow(GP[i][j-1][k]-GP[i][j+1][k],2);
					ISP[i][j][k][2]=13.0/12.0*pow(GP[i][j][k]-2.0*GP[i][j+1][k]+GP[i][j+2][k],2)+0.25*pow(3.0*GP[i][j][k]-4.0*GP[i][j+1][k]+GP[i][j+2][k],2);

					ISN[i][j][k][0]=13.0/12.0*pow(GN[i][j+1][k]-2.0*GN[i][j][k]+GN[i][j-1][k],2)+0.25*pow(3.0*GN[i][j+1][k]-4.0*GN[i][j][k]+GN[i][j-1][k],2);
					ISN[i][j][k][1]=13.0/12.0*pow(GN[i][j+2][k]-2.0*GN[i][j+1][k]+GN[i][j][k],2)+0.25*pow(GN[i][j+2][k]-GN[i][j][k],2);
					ISN[i][j][k][2]=13.0/12.0*pow(GN[i][j+3][k]-2.0*GN[i][j+2][k]+GN[i][j+1][k],2)+0.25*pow(GN[i][j+3][k]-4.0*GN[i][j+2][k]+3.0*GN[i][j+1][k],2);
				}
		}
			
	//计算alphaPN以及omegaPN
	for(i=2;i<=Jx+3;i++)                            
		for(j=2;j<=Jy+3;j++)
		{
			if(i<=j+jx)        
				for(k=0;k<=3;k++)
				{
					alphaP[i][j][k][0]=0.1/pow(delta+ISP[i][j][k][0],P);
					alphaP[i][j][k][1]=0.6/pow(delta+ISP[i][j][k][1],P);
					alphaP[i][j][k][2]=0.3/pow(delta+ISP[i][j][k][2],P);

					alphaN[i][j][k][0]=0.3/pow(delta+ISN[i][j][k][0],P);
					alphaN[i][j][k][1]=0.6/pow(delta+ISN[i][j][k][1],P);
					alphaN[i][j][k][2]=0.1/pow(delta+ISN[i][j][k][2],P);
					for(s=0;s<=2;s++)
					{
						omegaP[i][j][k][s]=alphaP[i][j][k][s]/(alphaP[i][j][k][0]+alphaP[i][j][k][1]+alphaP[i][j][k][2]);
						omegaN[i][j][k][s]=alphaN[i][j][k][s]/(alphaN[i][j][k][0]+alphaN[i][j][k][1]+alphaN[i][j][k][2]);						
					}
				}
		}
			
	//计算q3PN
	for(i=2;i<=Jx+3;i++)                            
		for(j=2;j<=Jy+3;j++)
		{
			if(i<=j+jx)      
				for(k=0;k<=3;k++)
				{
					q3P[i][j][k][0]=1.0/3.0*GP[i][j-2][k]-7.0/6.0*GP[i][j-1][k]+11.0/6.0*GP[i][j][k];
					q3P[i][j][k][1]=-1.0/6.0*GP[i][j-1][k]+5.0/6.0*GP[i][j][k]+1.0/3.0*GP[i][j+1][k];
					q3P[i][j][k][2]=1.0/3.0*GP[i][j][k]+5.0/6.0*GP[i][j+1][k]-1.0/6.0*GP[i][j+2][k];

					q3N[i][j][k][0]=-1.0/6.0*GN[i][j-1][k]+5.0/6.0*GN[i][j][k]+1.0/3.0*GN[i][j+1][k];
					q3N[i][j][k][1]=1.0/3.0*GN[i][j][k]+5.0/6.0*GN[i][j+1][k]-1.0/6.0*GN[i][j+2][k];
					q3N[i][j][k][2]=11.0/6.0*GN[i][j+1][k]-7.0/6.0*GN[i][j+2][k]+1.0/3.0*GN[i][j+3][k];
				}
		}
	
	//计算G_PN
	for(i=2;i<=Jx+3;i++)                            
		for(j=2;j<=Jy+3;j++)
		{
			if(i<=j+jx)    
				for(k=0;k<=3;k++)
				{
					G_P[i][j][k]=0;
					G_N[i][j][k]=0;
				}
			
		}
	
	for(i=2;i<=Jx+3;i++)                            
		for(j=2;j<=Jy+3;j++)
		{
			if(i<=j+jx)         
				for(k=0;k<=3;k++)
				{
					for(s=0;s<=2;s++)
					{
						G_P[i][j][k]+=omegaP[i][j][k][s]*q3P[i][j][k][s];
						G_N[i][j][k]+=omegaN[i][j][k][s]*q3N[i][j][k][s];
					}
				}
			
		}

	//计算G_
	for(i=2;i<=Jx+3;i++)                            
		for(j=2;j<=Jy+3;j++)
		{
			if(i<=j+jx)      
				for(k=0;k<=3;k++)
				{
					G_[i][j][k]=G_P[i][j][k]+G_N[i][j][k];
				}
			
		}
}

void WENO_solver(double dx,double dy,double dt)
{
    int i,j,k,jx;
	jx=(Jx*15)/100;

    FPFN(U);
	WENO_x();
	GPGN(U);
    WENO_y();
	for(i=3;i<=Jx+3;i++)                             
		for(j=3;j<=Jy+3;j++)
		{
			if(i<=j+jx)     
			for(k=0;k<=3;k++)
			U1[i][j][k]=U[i][j][k]-dt/dx*(F_[i][j][k]-F_[i-1][j][k])-dt/dy*(G_[i][j][k]-G_[i][j-1][k]);
		}	
    bound(U1);

    FPFN(U1);
    WENO_x();
	GPGN(U1);
    WENO_y();  
	for(i=3;i<=Jx+3;i++)                             
		for(j=3;j<=Jy+3;j++)
		{
			if(i<=j+jx)    
				for(k=0;k<=3;k++)
            U2[i][j][k]=0.75*U[i][j][k]+0.25*U1[i][j][k]-0.25*dt/dx*(F_[i][j][k]-F_[i-1][j][k])-0.25*dt/dy*(G_[i][j][k]-G_[i][j-1][k]);
		}
    bound(U2);

    FPFN(U2);
	WENO_x();
	GPGN(U2);
    WENO_y();  
	for(i=3;i<=Jx+3;i++)                             
		for(j=3;j<=Jy+3;j++)
		{
			if(i<=j+jx)       
			for(k=0;k<=3;k++)
            U3[i][j][k]=1.0/3.0*U[i][j][k]+2.0/3.0*U2[i][j][k]-2.0/3.0*dt/dx*(F_[i][j][k]-F_[i-1][j][k])-2.0/3.0*dt/dy*(G_[i][j][k]-G_[i][j-1][k]);
		}

	for(i=3;i<=Jx+3;i++)                             
		for(j=3;j<=Jy+3;j++)
		{
			if(i<=j+jx)       
			for(k=0;k<=3;k++)
            U[i][j][k]=U3[i][j][k];                
			
		}
    bound(U);
}

/*-------------------------------------------------------
输出 文件格式数据画图
入口: U，当前时刻U矢量，
     dx， x方向的网格宽度，dy， y方向的网格宽度；
出口: 无。
---------------------------------------------------------*/
void Output(double Ut[Jx+7][Jy+7][4],double dx,double dy,double T)
{
 int i,j;
 FILE *fp;
 double rou,u,v,p;
 char filename[30];

 sprintf(filename,"%1.4f.plt",T); 
 fp=fopen(filename,"w");
 fprintf(fp,"TITLE     = \"Dataset\"\nVARIABLES = \"x\" \"y\" \"rou\" \"u\" \"v\" \"p\" \"E\"");
 fprintf(fp,"ZONE T=\"Zone 1\"\nI=%d J=%d K=%d ZONETYPE=Ordered\n",Jx+1,Jy+1,1);
 fprintf(fp,"DATAPACKING=POINT\n");
  
 for(i=3;i<=Jx+3;i++)
  for(j=3;j<=Jy+3;j++) 
  {
      if(j<i-int(0.15*Jx))
      {
         rou=-5;
         u=-5;
         v=-5;
         p=-5;
      }
      else
      {                      
         rou=Ut[i][j][0];
         u=Ut[i][j][1]/rou;
         v=Ut[i][j][2]/rou;
         p=(GAMA-1)*(Ut[i][j][3]-0.5*Ut[i][j][0]*(u*u+v*v));
      }           
     

fprintf(fp,"%20f%20f%20.10e%20.10e%20.10e%20.10e%20.10e\n",(i-3)*dx,(j-3)*dy,rou,u,v,p,Ut[i][j][3]); 

  }

 fclose(fp);

}

/*-------------------------------------------------------
主函数
入口: 无；
出口: 无。
---------------------------------------------------------*/
int main()
{
  double T,dx,dy,dt;
  int m;

  dx=Lx/Jx;
  dy=Ly/Jy;
  m=1;
  T=0;

  Init(dx);
  bound(U);
  Output(U,dx,dy,T);

  while(T<TT)
  {
  dt=CFL(U,dx,dy);
  T+=dt; 
  printf("T=%10g    dt=%10g\n",T,dt);
  WENO_solver(dx,dy,dt);
  if(m*T0<=T&&m*T0>T-dt)
	{
     Output(U,dx,dy,T);
     m++;
	}
  }
  Output(U,dx,dy,T);
}