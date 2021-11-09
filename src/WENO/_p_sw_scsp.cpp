#include <stdio.h>
#include <math.h>

#define GAMA 1.4  //气体常数
#define Ms 3.0

#define MIN(x,y) (((x)<(y))?(x):(y))
#define MAX(x,y) (((x)>(y))?(x):(y))

#define Ly 2.0  //计算区域
#define Lx 6.0

#define Sf 0.5  //时间步长因子
#define nu 0.1  //人工粘性系数
#define TT 4.0  //总运行时间
#define T0 0.01 //输出间隔时间

#define Ny  100   //网格数
#define Nx  300 

double U[Nx+5][Ny+5][4];//全局变量
double Ff[Nx+5][Ny+5][4],Gf[Nx+5][Ny+5][4];//粘性项TEMP 
///时间步长函数
double CFL(double U[Nx+5][Ny+5][4],double dx,double dy)
{
  int i,j;
  double maxvel,p,u,v,vel;
  maxvel=1e-100;
  for(i=2;i<=Nx+2;i++)
    for(j=2;j<=Ny+2;j++)
      {
        u=U[i][j][1]/U[i][j][0];
        v=U[i][j][2]/U[i][j][0];
        p=(GAMA-1)*(U[i][j][3]-0.5*U[i][j][0]*(u*u+v*v));  
        vel=sqrt(GAMA*p/U[i][j][0])+sqrt(u*u+v*v);
        if(vel>maxvel)maxvel=vel;
      }
  return Sf*MIN(dx,dy)/maxvel;
}

/////初值函数
 void Init(double U[Nx+5][Ny+5][4])
{
     int i,j;
     double rou=1.0,u=0,v=0,p=0.71429;

     for(i=0;i<=Nx+4;i++)
     for(j=0;j<=Ny+4;j++)
       {
         U[i][j][0]=rou;
         U[i][j][1]=rou*u;
         U[i][j][2]=rou*v;
         U[i][j][3]=p/(GAMA-1)+0.5*rou*(u*u+v*v);
       }
}
/////边界条件函数（含角点）
 void bound(double U[Nx+5][Ny+5][4])
{
     int i,j,k;
     double rou1=1.0,u1=0,v1=0,p1=0.71429;
     double rou2,u2,v2,p2;

     //波后参数
     rou2=rou1*Ms*Ms/(1+(GAMA-1)/(GAMA+1)*(Ms*Ms-1));
     u2=2/(GAMA+1)*(Ms-1/Ms)*sqrt(GAMA*p1/rou1);
     v2=0;
     p2=p1*(1+2*GAMA/(GAMA+1)*(Ms*Ms-1));

     //左边界
     for(i=0;i<=2;i++)
     for(j=2;j<=Ny+2;j++)
     {
      U[i][j][0]=rou2;
      U[i][j][1]=rou2*u2;
      U[i][j][2]=0;
      U[i][j][3]=p2/(GAMA-1)+0.5*rou2*(u2*u2+v2*v2);
     }
     //右边界
     for(i=Nx+3;i<=Nx+4;i++)for(j=2;j<=Ny+2;j++)for(k=0;k<=3;k++)U[i][j][k]=U[Nx+2][j][k];





      //固壁边界法向速度为零
      //X方向线段
      for(i=0;i<=Nx+4;i++)
      {
       U[i][Ny+2][2]=0;
       U[i][2][2]=0; 
      }
      for(i=Nx/6+3;i<=Nx/3+1;i++)U[i][Ny/4+2][2]=0;
      for(i=Nx/2+3;i<=Nx*2/3+1;i++)U[i][Ny*3/4+2][2]=0;

      //y方向线段
      for(j=Ny*3/4+3;j<=Ny+2;j++)
      {
       U[Nx/2+2][j][1]=0;
       U[Nx*2/3+2][j][1]=0;
      }
      for(j=2;j<=Ny/4+1;j++)
      {
       U[Nx/6+2][j][1]=0;
       U[Nx/3+2][j][1]=0;
      }      

      //虚拟边界的镜像处理
      //上边界
       for(j=Ny*3/4+5;j<=Ny+2;j++)//第二段 //
        for(i=Nx/2+3;i<=Nx/2+4;i++)
        {
         U[i][j][0]=U[Nx+4-i][j][0];
         U[i][j][1]=-U[Nx+4-i][j][1];
         U[i][j][2]=U[Nx+4-i][j][2];
         U[i][j][3]=U[Nx+4-i][j][3];
        }

      for(i=Nx/2+5;i<=Nx*2/3-1;i++)//第三段
       for(j=Ny*3/4+3;j<=Ny*3/4+4;j++)
       {
        U[i][j][0]=U[i][Ny*3/2+4-j][0];
        U[i][j][1]=U[i][Ny*3/2+4-j][1];
        U[i][j][2]=-U[i][Ny*3/2+4-j][2];
        U[i][j][3]=U[i][Ny*3/2+4-j][3];
       }        



       for(j=Ny*3/4+5;j<=Ny+2;j++)//第四段 //
        for(i=Nx*2/3;i<=Nx*2/3+1;i++)
        {
         U[i][j][0]=U[Nx*4/3+4-i][j][0];
         U[i][j][1]=-U[Nx*4/3+4-i][j][1];
         U[i][j][2]=U[Nx*4/3+4-i][j][2];
         U[i][j][3]=U[Nx*4/3+4-i][j][3];
        }          

      for(i=0;i<=Nx+4;i++)//第一段//
       for(j=Ny+3;j<=Ny+4;j++)
       {
        U[i][j][0]=U[i][2*Ny+4-j][0];
        U[i][j][1]=U[i][2*Ny+4-j][1];
        U[i][j][2]=-U[i][2*Ny+4-j][2];
        U[i][j][3]=U[i][2*Ny+4-j][3];
       }

        //下边界
       for(j=2;j<=Ny/4-1;j++)//第二段 /
        for(i=Nx/6+3;i<=Nx/6+4;i++)
        {
         U[i][j][0]=U[Nx/3+4-i][j][0];
         U[i][j][1]=-U[Nx/3+4-i][j][1];
         U[i][j][2]=U[Nx/3+4-i][j][2];
         U[i][j][3]=U[Nx/3+4-i][j][3];
        }
for(i=Nx/6+5;i<=Nx/3-1;i++)//第三段
       for(j=Ny/4;j<=Ny/4+1;j++)
       {
        U[i][j][0]=U[i][Ny/2+4-j][0];
        U[i][j][1]=U[i][Ny/2+4-j][1];
        U[i][j][2]=-U[i][Ny/2+4-j][2];
        U[i][j][3]=U[i][Ny/2+4-j][3];
       }        
for(j=2;j<=Ny/4-1;j++)//第四段
        for(i=Nx/3;i<=Nx/3+1;i++)
        {
         U[i][j][0]=U[Nx*2/3+4-i][j][0];
         U[i][j][1]=-U[Nx*2/3+4-i][j][1];
         U[i][j][2]=U[Nx*2/3+4-i][j][2];
         U[i][j][3]=U[Nx*2/3+4-i][j][3];
        }    

      for(i=0;i<=Nx+4;i++)//第一段
       for(j=0;j<=1;j++)
       {
        U[i][j][0]=U[i][4-j][0];
        U[i][j][1]=U[i][4-j][1];
        U[i][j][2]=-U[i][4-j][2];
        U[i][j][3]=U[i][4-j][3];
       }

       //非计算区域置零
      for(i=Nx/6+5;i<=Nx/3-1;i++)
       for(j=0;j<=Ny/4-1;j++)
       {
         U[i][j][0]=1;
         U[i][j][1]=0;
         U[i][j][2]=0;
         U[i][j][3]=0.71429/(GAMA-1);
       }

      for(i=Nx/2+5;i<=Nx*2/3-1;i++)
       for(j=Ny*3/4+5;j<=Ny+4;j++)
       {
         U[i][j][0]=1;
         U[i][j][1]=0;
         U[i][j][2]=0;
         U[i][j][3]=0.71429/(GAMA-1);
       }

      //凸角点处理共四个
      //上边界第二三段交界角点
      U[Nx/2+3][Ny*3/4+4][0]=U[Nx/2+1][Ny*3/4+4][0];
      U[Nx/2+3][Ny*3/4+4][1]=-U[Nx/2+1][Ny*3/4+4][1];
      U[Nx/2+3][Ny*3/4+4][2]=-U[Nx/2+1][Ny*3/4+4][2];
      U[Nx/2+3][Ny*3/4+4][3]=U[Nx/2+1][Ny*3/4+4][3];

      U[Nx/2+4][Ny*3/4+3][0]=U[Nx/2+4][Ny*3/4+1][0];
      U[Nx/2+4][Ny*3/4+3][1]=-U[Nx/2+4][Ny*3/4+1][1];
      U[Nx/2+4][Ny*3/4+3][2]=-U[Nx/2+4][Ny*3/4+1][2];
      U[Nx/2+4][Ny*3/4+3][3]=U[Nx/2+4][Ny*3/4+1][3];

      U[Nx/2+3][Ny*3/4+3][0]=0.5*(U[Nx/2+3][Ny*3/4+1][0]+U[Nx/2+1][Ny*3/4+3][0]);
      U[Nx/2+3][Ny*3/4+3][1]=0.5*(U[Nx/2+3][Ny*3/4+1][1]+U[Nx/2+1][Ny*3/4+3][1]);
      U[Nx/2+3][Ny*3/4+3][2]=0.5*(U[Nx/2+3][Ny*3/4+1][2]+U[Nx/2+1][Ny*3/4+3][2]);
      U[Nx/2+3][Ny*3/4+3][3]=0.5*(U[Nx/2+3][Ny*3/4+1][3]+U[Nx/2+1][Ny*3/4+3][3]); 

      U[Nx/2+2][Ny*3/4+2][0]=0.5*(U[Nx/2+1][Ny*3/4+2][0]+U[Nx/2+2][Ny*3/4+1][0]);
      U[Nx/2+2][Ny*3/4+2][1]=U[Nx/2+2][Ny*3/4+1][1];
      U[Nx/2+2][Ny*3/4+2][2]=U[Nx/2+1][Ny*3/4+2][2];
      U[Nx/2+2][Ny*3/4+2][3]=0.5*(U[Nx/2+1][Ny*3/4+2][3]+U[Nx/2+2][Ny*3/4+1][3]);

      U[Nx/2+4][Ny*3/4+4][0]=U[Nx/2+2][Ny*3/4+2][0];
      U[Nx/2+4][Ny*3/4+4][1]=-U[Nx/2+2][Ny*3/4+2][1];
      U[Nx/2+4][Ny*3/4+4][2]=-U[Nx/2+2][Ny*3/4+2][2];
      U[Nx/2+4][Ny*3/4+4][3]=U[Nx/2+2][Ny*3/4+2][3];


      //上边界第三四段交界角点
      U[Nx*2/3][Ny*3/4+3][0]=U[Nx*2/3][Ny*3/4+1][0];
      U[Nx*2/3][Ny*3/4+3][1]=-U[Nx*2/3][Ny*3/4+1][1];
      U[Nx*2/3][Ny*3/4+3][2]=-U[Nx*2/3][Ny*3/4+1][2];
      U[Nx*2/3][Ny*3/4+3][3]=U[Nx*2/3][Ny*3/4+1][3]; 

      U[Nx*2/3+1][Ny*3/4+4][0]=U[Nx*2/3+3][Ny*3/4+4][0];
      U[Nx*2/3+1][Ny*3/4+4][1]=-U[Nx*2/3+3][Ny*3/4+4][1];
      U[Nx*2/3+1][Ny*3/4+4][2]=-U[Nx*2/3+3][Ny*3/4+4][2];
      U[Nx*2/3+1][Ny*3/4+4][3]=U[Nx*2/3+3][Ny*3/4+4][3]; 

      U[Nx*2/3+1][Ny*3/4+3][0]=0.5*(U[Nx*2/3+1][Ny*3/4+1][0]+U[Nx*2/3+3][Ny*3/4+3][0]);
      U[Nx*2/3+1][Ny*3/4+3][1]=0.5*(U[Nx*2/3+1][Ny*3/4+1][1]+U[Nx*2/3+3][Ny*3/4+3][1]);
      U[Nx*2/3+1][Ny*3/4+3][2]=0.5*(U[Nx*2/3+1][Ny*3/4+1][2]+U[Nx*2/3+3][Ny*3/4+3][2]);
      U[Nx*2/3+1][Ny*3/4+3][3]=0.5*(U[Nx*2/3+1][Ny*3/4+1][3]+U[Nx*2/3+3][Ny*3/4+3][3]); 

      U[Nx*2/3+2][Ny*3/4+2][0]=0.5*(U[Nx*2/3+2][Ny*3/4+1][0]+U[Nx*2/3+3][Ny*3/4+2][0]);
      U[Nx*2/3+2][Ny*3/4+2][1]=U[Nx*2/3+2][Ny*3/4+1][1];
      U[Nx*2/3+2][Ny*3/4+2][2]=U[Nx*2/3+3][Ny*3/4+2][2];
      U[Nx*2/3+2][Ny*3/4+2][3]=0.5*(U[Nx*2/3+2][Ny*3/4+1][3]+U[Nx*2/3+3][Ny*3/4+2][3]);

      U[Nx*2/3][Ny*3/4+4][0]=U[Nx*2/3+2][Ny*3/4+2][0];
      U[Nx*2/3][Ny*3/4+4][1]=-U[Nx*2/3+2][Ny*3/4+2][1];
      U[Nx*2/3][Ny*3/4+4][2]=-U[Nx*2/3+2][Ny*3/4+2][2];
      U[Nx*2/3][Ny*3/4+4][3]=U[Nx*2/3+2][Ny*3/4+2][3]; 

//下边界第二三段交界角点
      U[Nx/6+3][Ny/4][0]=U[Nx/6+1][Ny/4][0];
      U[Nx/6+3][Ny/4][1]=-U[Nx/6+1][Ny/4][1];
      U[Nx/6+3][Ny/4][2]=-U[Nx/6+1][Ny/4][2];
      U[Nx/6+3][Ny/4][3]=U[Nx/6+1][Ny/4][3];

      U[Nx/6+4][Ny/4+1][0]=U[Nx/6+4][Ny/4+3][0];
      U[Nx/6+4][Ny/4+1][1]=-U[Nx/6+4][Ny/4+3][1];
      U[Nx/6+4][Ny/4+1][2]=-U[Nx/6+4][Ny/4+3][2];
      U[Nx/6+4][Ny/4+1][3]=U[Nx/6+4][Ny/4+3][3];

      U[Nx/6+3][Ny/4+1][0]=0.5*(U[Nx/6+1][Ny/4+1][0]+U[Nx/6+3][Ny/4+3][0]);
      U[Nx/6+3][Ny/4+1][1]=0.5*(U[Nx/6+1][Ny/4+1][1]+U[Nx/6+3][Ny/4+3][1]);
      U[Nx/6+3][Ny/4+1][2]=0.5*(U[Nx/6+1][Ny/4+1][2]+U[Nx/6+3][Ny/4+3][2]);
      U[Nx/6+3][Ny/4+1][3]=0.5*(U[Nx/6+1][Ny/4+1][3]+U[Nx/6+3][Ny/4+3][3]);

      U[Nx/6+2][Ny/4+2][0]=0.5*(U[Nx/6+1][Ny/4+2][0]+U[Nx/6+2][Ny/4+3][0]);
      U[Nx/6+2][Ny/4+2][1]=U[Nx/6+2][Ny/4+3][1];
      U[Nx/6+2][Ny/4+2][2]=U[Nx/6+1][Ny/4+2][2];
      U[Nx/6+2][Ny/4+2][3]=0.5*(U[Nx/6+1][Ny/4+2][3]+U[Nx/6+2][Ny/4+3][3]);

      U[Nx/6+4][Ny/4][0]=U[Nx/6+2][Ny/4+2][0];
      U[Nx/6+4][Ny/4][1]=-U[Nx/6+2][Ny/4+2][1];
      U[Nx/6+4][Ny/4][2]=-U[Nx/6+2][Ny/4+2][2];
      U[Nx/6+4][Ny/4][3]=U[Nx/6+2][Ny/4+2][3]; 
      //下边界第三四段交界角点
      U[Nx/3][Ny/4+1][0]=U[Nx/3][Ny/4+3][0];
      U[Nx/3][Ny/4+1][1]=-U[Nx/3][Ny/4+3][1];
      U[Nx/3][Ny/4+1][2]=-U[Nx/3][Ny/4+3][2];
      U[Nx/3][Ny/4+1][3]=U[Nx/3][Ny/4+3][3];

      U[Nx/3+1][Ny/4][0]=U[Nx/3+3][Ny/4][0];
      U[Nx/3+1][Ny/4][1]=-U[Nx/3+3][Ny/4][1];
      U[Nx/3+1][Ny/4][2]=-U[Nx/3+3][Ny/4][2];
      U[Nx/3+1][Ny/4][3]=U[Nx/3+3][Ny/4][3];

      U[Nx/3+1][Ny/4+1][0]=0.5*(U[Nx/3+3][Ny/4+1][0]+U[Nx/3+1][Ny/4+3][0]);
      U[Nx/3+1][Ny/4+1][1]=0.5*(U[Nx/3+3][Ny/4+1][1]+U[Nx/3+1][Ny/4+3][1]);
      U[Nx/3+1][Ny/4+1][2]=0.5*(U[Nx/3+3][Ny/4+1][2]+U[Nx/3+1][Ny/4+3][2]);
      U[Nx/3+1][Ny/4+1][3]=0.5*(U[Nx/3+3][Ny/4+1][3]+U[Nx/3+1][Ny/4+3][3]); 

      U[Nx/3+2][Ny/4+2][0]=0.5*(U[Nx/3+3][Ny/4+2][0]+U[Nx/3+2][Ny/4+3][0]);
      U[Nx/3+2][Ny/4+2][1]=U[Nx/3+2][Ny/4+3][1];
      U[Nx/3+2][Ny/4+2][2]=U[Nx/3+3][Ny/4+2][2];
      U[Nx/3+2][Ny/4+2][3]=0.5*(U[Nx/3+3][Ny/4+2][3]+U[Nx/3+2][Ny/4+3][3]);

      U[Nx/3][Ny/4][0]=U[Nx/3+2][Ny/4+2][0];
      U[Nx/3][Ny/4][1]=-U[Nx/3+2][Ny/4+2][1];
      U[Nx/3][Ny/4][2]=-U[Nx/3+2][Ny/4+2][2];
      U[Nx/3][Ny/4][3]=U[Nx/3+2][Ny/4+2][3]; 
}     
//ENO差分求解器

/*X方向流通矢量f(u)分裂
入口：U
出口： FP FN*/
void S_Wx(double U[Nx+5][Ny+5][4],double FP[Nx+5][Ny+5][4],double FN[Nx+5][Ny+5][4])
{
     int i,j,k;
     double lambda[4],lambda_p[4],lambda_n[4],x;
     x=1e-4;
     double u,v,p,a,h,q;
     //粘性处理
     for(i=2;i<=Nx+2;i++)
     for(j=2;j<=Ny+2;j++) 
      {
         q=fabs(fabs(U[i+1][j][0]-U[i][j][0])-fabs(U[i][j][0]-U[i-1][j][0]))/(fabs(U[i+1][j][0]-U[i][j][0])+fabs(U[i][j][0]-U[i-1][j][0])+1e-100);
         for(k=0;k<4;k++)
         Ff[i][j][k]=U[i][j][k]+0.5*nu*q*(U[i+1][j][k]-2*U[i][j][k]+U[i-1][j][k]); 
      }
     for(i=2;i<=Nx+2;i++)
     for(j=2;j<=Ny+2;j++)
     for(k=0;k<4;k++)U[i][j][k]=Ff[i][j][k];
     bound(U); 
     //粘性处理结束
     for(i=0;i<=Nx+4;i++)
      for(j=0;j<=Ny+4;j++)
      {
                         u=U[i][j][1]/U[i][j][0];
                         v=U[i][j][2]/U[i][j][0];
                         p=(GAMA-1)*(U[i][j][3]-0.5*U[i][j][0]*(u*u+v*v));
                         a=sqrt(GAMA*p/U[i][j][0]);
                         h=(U[i][j][3]+p)/U[i][j][0];

                         lambda[0]=u;
                         lambda[1]=u;
                         lambda[2]=u-a;
                         lambda[3]=u+a;

                         for(k=0;k<=3;k++)
                         {
                           lambda_p[k]=(lambda[k]+sqrt(lambda[k]*lambda[k]+x*x))/2;
                           lambda_n[k]=(lambda[k]-sqrt(lambda[k]*lambda[k]+x*x))/2;
                         }

                         FP[i][j][0]=U[i][j][0]/(2*GAMA)*(2*(GAMA-1)*lambda_p[0]+lambda_p[2]+lambda_p[3]);
                         FP[i][j][1]=U[i][j][0]/(2*GAMA)*(2*(GAMA-1)*u*lambda_p[0]+(u-a)*lambda_p[2]+(u+a)*lambda_p[3]);
                         FP[i][j][2]=U[i][j][0]/(2*GAMA)*(2*(GAMA-1)*v*lambda_p[1]+v*lambda_p[2]+v*lambda_p[3]);
                         FP[i][j][3]=U[i][j][0]/(2*GAMA)*((GAMA-1)*(u*u-v*v)*lambda_p[0]+2*(GAMA-1)*v*v*lambda_p[1]+(h-a*u)*lambda_p[2]+(h+a*u)*lambda_p[3]);

                         FN[i][j][0]=U[i][j][0]/(2*GAMA)*(2*(GAMA-1)*lambda_n[0]+lambda_n[2]+lambda_n[3]);
                         FN[i][j][1]=U[i][j][0]/(2*GAMA)*(2*(GAMA-1)*u*lambda_n[0]+(u-a)*lambda_n[2]+(u+a)*lambda_n[3]);
                         FN[i][j][2]=U[i][j][0]/(2*GAMA)*(2*(GAMA-1)*v*lambda_n[1]+v*lambda_n[2]+v*lambda_n[3]);
                         FN[i][j][3]=U[i][j][0]/(2*GAMA)*((GAMA-1)*(u*u-v*v)*lambda_n[0]+2*(GAMA-1)*v*v*lambda_n[1]+(h-a*u)*lambda_n[2]+(h+a*u)*lambda_n[3]);
       }
}

/*y方向流通矢量g(u)分裂
入口：U
出口：GP GN*/
void S_Wy(double U[Nx+5][Ny+5][4],double GP[Nx+5][Ny+5][4],double GN[Nx+5][Ny+5][4])
{
     int i,j,k;
     double mu[4],mu_p[4],mu_n[4],x;
     x=1e-4;
     double u,v,p,a,h,q;
     //粘性处理
     for(i=2;i<=Nx+2;i++)
     for(j=2;j<=Ny+2;j++) 
      {
      q=fabs(fabs(U[i][j+1][0]-U[i][j][0])-fabs(U[i][j][0]-U[i][j-1][0]))/(fabs(U[i][j+1][0]-U[i][j][0])+fabs(U[i][j][0]-U[i][j-1][0])+1e-100);
      for(k=0;k<4;k++)
      Gf[i][j][k]=U[i][j][k]+0.5*nu*q*(U[i][j+1][k]-2*U[i][j][k]+U[i][j-1][k]); 
      }
     for(i=2;i<=Nx+2;i++)
     for(j=2;j<=Ny+2;j++) 
     for(k=0;k<4;k++)U[i][j][k]=Gf[i][j][k]; 
     bound(U); 
     //粘性处理结束
     for(i=0;i<=Nx+4;i++)
     for(j=0;j<=Ny+4;j++)
     {
                u=U[i][j][1]/U[i][j][0];
                v=U[i][j][2]/U[i][j][0];
                p=(GAMA-1)*(U[i][j][3]-0.5*U[i][j][0]*(u*u+v*v));
                a=sqrt(GAMA*p/U[i][j][0]);
                h=(U[i][j][3]+p)/U[i][j][0];

                mu[0]=v;
                mu[1]=v;
                mu[2]=v-a;
                mu[3]=v+a;

                for(k=0;k<=3;k++)
                {
                 mu_p[k]=(mu[k]+sqrt(mu[k]*mu[k]+x*x))/2;
                 mu_n[k]=(mu[k]-sqrt(mu[k]*mu[k]+x*x))/2;
                }

                GP[i][j][0]=U[i][j][0]/(2*GAMA)*(2*(GAMA-1)*mu_p[1]+mu_p[2]+mu_p[3]);
                GP[i][j][1]=U[i][j][0]/(2*GAMA)*(2*(GAMA-1)*u*mu_p[0]+u*mu_p[2]+u*mu_p[3]);
                GP[i][j][2]=U[i][j][0]/(2*GAMA)*(2*(GAMA-1)*v*mu_p[1]+(v-a)*mu_p[2]+(v+a)*mu_p[3]);
                GP[i][j][3]=U[i][j][0]/(2*GAMA)*(2*(GAMA-1)*u*u*mu_p[0]+(GAMA-1)*(v*v-u*u)*mu_p[1]+(h-a*v)*mu_p[2]+(h+a*v)*mu_p[3]);

                GN[i][j][0]=U[i][j][0]/(2*GAMA)*(2*(GAMA-1)*mu_n[1]+mu_n[2]+mu_n[3]);
                GN[i][j][1]=U[i][j][0]/(2*GAMA)*(2*(GAMA-1)*u*mu_n[0]+u*mu_n[2]+u*mu_n[3]);
                GN[i][j][2]=U[i][j][0]/(2*GAMA)*(2*(GAMA-1)*v*mu_n[1]+(v-a)*mu_n[2]+(v+a)*mu_n[3]);
                GN[i][j][3]=U[i][j][0]/(2*GAMA)*(2*(GAMA-1)*u*u*mu_n[0]+(GAMA-1)*(v*v-u*u)*mu_n[1]+(h-a*v)*mu_n[2]+(h+a*v)*mu_n[3]);
     }
}

/*求解数值通量*/

//光滑模板判别式
double minmod(double a,double b,double c,double x,double y,double z)
{
       double l;
       if(x<=y&&x<=z) l=a;
       else if(y<=x&&y<=z) l=b;
       else if(z<=x&&z<=y) l=c;
       return l;
}

//半点上数值通量的计算
/*入口：FP,FN 
出口：F1*/
void ENOX(double F1[Nx+5][Ny+5][4],double FP[Nx+5][Ny+5][4],double FN[Nx+5][Ny+5][4])
{
     int i,j,k;
     double qp[3],ISp[3],qn[3],ISn[3];

     for(i=2;i<=Nx+2;i++)
     for(j=2;j<=Ny+2;j++)
     for(k=0;k<=3;k++)
     {
      qp[0]=1.0/3*FP[i-2][j][k]-7.0/6*FP[i-1][j][k]+11.0/6*FP[i][j][k];
      qp[1]=-1.0/6*FP[i-1][j][k]+5.0/6*FP[i][j][k]+1.0/3*FP[i+1][j][k];
      qp[2]=1.0/3*FP[i][j][k]+5.0/6*FP[i+1][j][k]-1.0/6*FP[i+2][j][k];

      ISp[0]=13.0/12*pow((FP[i-2][j][k]-2*FP[i-1][j][k]+FP[i][j][k]),2.0)+1.0/4*pow((FP[i-2][j][k]-4*FP[i-1][j][k]+3*FP[i][j][k]),2.0);
      ISp[1]=13.0/12*pow((FP[i-1][j][k]-2*FP[i][j][k]+FP[i+1][j][k]),2.0)+1.0/4*pow((FP[i-1][j][k]-FP[i+1][j][k]),2.0);
      ISp[2]=13.0/12*pow((FP[i][j][k]-2*FP[i+1][j][k]+FP[i+2][j][k]),2.0)+1.0/4*pow((3*FP[i][j][k]-4*FP[i+1][j][k]+FP[i+2][j][k]),2.0);

      qn[0]=-1.0/6*FN[i-1][j][k]+5.0/6*FN[i][j][k]+1.0/3*FN[i+1][j][k];
      qn[1]=1.0/3*FN[i][j][k]+5.0/6*FN[i+1][j][k]-1.0/6*FN[i+2][j][k];
      qn[2]=11.0/6*FN[i+1][j][k]-7.0/6*FN[i+2][j][k]+1.0/3*FN[i+3][j][k];

      ISn[0]=13.0/12*pow((FN[i-1][j][k]-2*FN[i][j][k]+FN[i+1][j][k]),2.0)+1.0/4*pow((FN[i-1][j][k]-4*FN[i][j][k]+3*FN[i+1][j][k]),2.0);
      ISn[1]=13.0/12*pow((FN[i][j][k]-2*FN[i+1][j][k]+FN[i+2][j][k]),2.0)+1.0/4*pow((FN[i+2][j][k]-FN[i][j][k]),2.0);
      ISn[2]=13.0/12*pow((FN[i+1][j][k]-2*FN[i+2][j][k]+FN[i+3][j][k]),2.0)+1.0/4*pow((3*FN[i+1][j][k]-4*FN[i+2][j][k]+FN[i+3][j][k]),2.0);

      F1[i][j][k]=minmod(qp[0],qp[1],qp[2],ISp[0],ISp[1],ISp[2])+minmod(qn[0],qn[1],qn[2],ISn[0],ISn[1],ISn[2]);//i点相当于半点(i+1/2)
      }
}

/*入口：GP,GN
出口：G1*/
void ENOY(double G1[Nx+5][Ny+5][4],double GP[Nx+5][Ny+5][4],double GN[Nx+5][Ny+5][4])
{
     int i,j,k;
     double qp[3],ISp[3],qn[3],ISn[3];

     for(j=2;j<=Ny+2;j++)
     for(i=2;i<=Nx+2;i++)
     for(k=0;k<=3;k++)
     {
      qp[0]=1.0/3*GP[i][j-2][k]-7.0/6*GP[i][j-1][k]+11.0/6*GP[i][j][k];
      qp[1]=-1.0/6*GP[i][j-1][k]+5.0/6*GP[i][j][k]+1.0/3*GP[i][j+1][k];
      qp[2]=1.0/3*GP[i][j][k]+5.0/6*GP[i][j+1][k]-1.0/6*GP[i][j+2][k];

      ISp[0]=13.0/12*pow((GP[i][j-2][k]-2*GP[i][j-1][k]+GP[i][j][k]),2.0)+1.0/4*pow((GP[i][j-2][k]-4*GP[i][j-1][k]+3*GP[i][j][k]),2.0);
      ISp[1]=13.0/12*pow((GP[i][j-1][k]-2*GP[i][j][k]+GP[i][j+1][k]),2.0)+1.0/4*pow((GP[i][j-1][k]-GP[i][j+1][k]),2.0);
      ISp[2]=13.0/12*pow((GP[i][j][k]-2*GP[i][j+1][k]+GP[i][j+2][k]),2.0)+1.0/4*pow((3*GP[i][j][k]-4*GP[i][j+1][k]+GP[i][j+2][k]),2.0);

      qn[0]=-1.0/6*GN[i][j-1][k]+5.0/6*GN[i][j][k]+1.0/3*GN[i][j+1][k];
      qn[1]=1.0/3*GN[i][j][k]+5.0/6*GN[i][j+1][k]-1.0/6*GN[i][j+2][k];
      qn[2]=11.0/6*GN[i][j+1][k]-7.0/6*GN[i][j+2][k]+1.0/3*GN[i][j+3][k];

      ISn[0]=13.0/12*pow((GN[i][j-1][k]-2*GN[i][j][k]+GN[i][j+1][k]),2.0)+1.0/4*pow((GN[i][j-1][k]-4*GN[i][j][k]+3*GN[i][j+1][k]),2.0);
      ISn[1]=13.0/12*pow((GN[i][j][k]-2*GN[i][j+1][k]+GN[i][j+2][k]),2.0)+1.0/4*pow((GN[i][j+2][k]-GN[i][j][k]),2.0);
      ISn[2]=13.0/12*pow((GN[i][j+1][k]-2*GN[i][j+2][k]+GN[i][j+3][k]),2.0)+1.0/4*pow((3*GN[i][j+1][k]-4*GN[i][j+2][k]+GN[i][j+3][k]),2.0);

      G1[i][j][k]=minmod(qp[0],qp[1],qp[2],ISp[0],ISp[1],ISp[2])+minmod(qn[0],qn[1],qn[2],ISn[0],ISn[1],ISn[2]);//j点相当于半点（j+1/2）
      }
}

//空间离散算子
/*x方向
入口：U
出口：Lf*/
void LLx(double U[Nx+5][Ny+5][4],double Lf[Nx+5][Ny+5][4])
{
     int i,j,k;
     static double FP[Nx+5][Ny+5][4],FN[Nx+5][Ny+5][4],F1[Nx+5][Ny+5][4];

     S_Wx(U,FP,FN);
     ENOX(F1,FP,FN);
     for(i=2;i<=Nx+2;i++)
     for(j=2;j<=Ny+2;j++)
     for(k=0;k<=3;k++)
     {
      Lf[i][j][k]=F1[i][j][k]-F1[i-1][j][k];
     }
}

/*y方向
入口：U
出口：Lg*/
void LLy(double U[Nx+5][Ny+5][4],double Lg[Nx+5][Ny+5][4])
{
     int i,j,k;
     static double GP[Nx+5][Ny+5][4],GN[Nx+5][Ny+5][4],G1[Nx+5][Ny+5][4];

     S_Wy(U,GP,GN);
     ENOY(G1,GP,GN);
     for(i=2;i<=Nx+2;i++)
     for(j=2;j<=Ny+2;j++)
     for(k=0;k<=3;k++)
     {
      Lg[i][j][k]=G1[i][j][k]-G1[i][j-1][k];
     }
}

/*ENO差分求解器
入口：U 上一时刻U矢量
        dx,dy,dt
出口：当前时刻的U矢量*/
void ENO_solver(double U[Nx+5][Ny+5][4],double dx,double dy,double dt)
{
     int i,j,k;
     double rrx,rry;
     static double Lf[Nx+5][Ny+5][4],Lg[Nx+5][Ny+5][4],U_mid[Nx+5][Ny+5][4];

     rrx=dt/dx;
     rry=dt/dy;

     LLx(U,Lf);
     for(i=2;i<=Nx+2;i++)
     for(j=2;j<=Ny+2;j++)
     for(k=0;k<=3;k++)
     {
      U_mid[i][j][k]=U[i][j][k]-rrx*Lf[i][j][k];
      }
     bound(U_mid);

     LLy(U_mid,Lg);
     for(i=2;i<=Nx+2;i++)
     for(j=2;j<=Ny+2;j++)
     for(k=0;k<=3;k++)
     {
      U[i][j][k]=U_mid[i][j][k]-rry*Lg[i][j][k];
     }
     bound(U);
}


/*输出函数*/ 

void Output(double U[Nx+5][Ny+5][4],double dx,double dy,double T)
{
 int i,j;
 FILE *fp;
 double rou,u,v,p;
 char filename[30];

 sprintf(filename,"%1.4f.plt",T); 
 fp=fopen(filename,"w");
 fprintf(fp,"TITLE     = \"Dataset\"\nVARIABLES = \"x\" \"y\" \"rou\" \"u\" \"v\" \"p\" \"E\"");
 fprintf(fp,"ZONE T=\"Zone 1\"\nI=%d J=%d K=%d ZONETYPE=Ordered\n",Ny+1,Nx+1,1);
 fprintf(fp,"DATAPACKING=POINT\n");


 for(i=2;i<=Nx+2;i++)
  for(j=2;j<=Ny+2;j++) 
  {
      if((i>Nx/6+2&&i<Nx/3+2&&j<Ny/4+2)||(i>Nx/2+2&&i<Nx*2/3+2&&j>Ny*3/4+2))
      {
         rou=0;
         u=0;
         v=0;
         p=0;
      }
      else
      {                      
         rou=U[i][j][0];
         u=U[i][j][1]/rou;
         v=U[i][j][2]/rou;
         p=(GAMA-1)*(U[i][j][3]-0.5*U[i][j][0]*(u*u+v*v));
      }           


fprintf(fp,"%20f%20f%20.10e%20.10e%20.10e%20.10e%20.10e\n",(i-2)*dx,(j-2)*dy,rou,u,v,p,U[i][j][3]); 

  }

 fclose(fp);

}

/*主函数*/
int main()
{
 double T=0,dx,dy,dt;
 int m=1;

 Init(U); 
 bound(U);

 dx=Lx/Nx;
 dy=Ly/Ny;

 Output(U,dx,dy,T); 

 while(T<TT)
 {
  dt=CFL(U,dx,dy);
  T+=dt; 
  printf("T=%10g    dt=%10g\n",T,dt);
  ENO_solver(U,dx,dy,dt);
  if(m*T0<=T&&m*T0>T-dt)
  {
   Output(U,dx,dy,T);
   m++;
  }
 }
 Output(U,dx,dy,T);
}

