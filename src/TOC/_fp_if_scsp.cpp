#include<stdio.h>
#include<math.h>
#include <string>

using namespace std;

#define ly 301// 
#define lx 601// 
#define Re 100//雷诺数
#define dx 0.01
#define dy 0.01
#define Dt 1e-3//时间步长
#define eps 1e-4//收敛限
#define C 1.8//人工压缩算法中的稳定性条件（不定值影响收敛！！！！！）
//全局变量
double U[lx+1][ly+1],V[lx+1][ly+1],P[lx+1][ly+1];
double Fux1[lx+1][ly+1],Fux2[lx+1][ly+1],Fuy1[lx+1][ly+1],Fuy2[lx+1][ly+1],Fvx1[lx+1][ly+1],Fvx2[lx+1][ly+1],
Fvy1[lx+1][ly+1],Fvy2[lx+1][ly+1],Sux[lx+1][ly+1],Suy[lx+1][ly+1],Svx[lx+1][ly+1],Svy[lx+1][ly+1];//1代表正2代表负
double dm; //判断收敛的值
/*-------------------------------------
进出口速度边界条件及上下速度边界条件（由于计算过程进出口的压力始终为1没变，初始条件设置了这里就没必要设了）
--------------------------------------*/
void BoundInUV()
{
	int i,j;
    for(j=(ly-1)/6+1;j<=(ly-1)*5/6;j++)
	{
		  if(j<=(ly-1)/6+0.1/dy)
			U[0][j]=10.0*((j-(ly-1)/6)-0.5)*dy;
		else if(j<=(ly-1)/6+1.9/dy)
			U[0][j]=1.0;
		else
			U[0][j]=10.0*(2.0-(j-(ly-1)/6-0.5)*dy);
	    V[0][j]=0;
		U[lx][j]=U[lx-1][j];
		V[lx][j]=V[lx-1][j]; //入口出口		
	} 
   		for(i=0;i<lx+1;i++)
		{
			if((i<=(lx-1)/6-1)|(i>=(lx-1)/3+1))
				U[i][(ly-1)/6]=-U[i][(ly-1)/6+1];
			else
           U[i][0]=-U[i][1];
		}
     for(j=0;j<(ly-1)/6;j++)
	 {
			V[(lx-1)/6][j]=-V[(lx-1)/6+1][j];
	    	U[(lx-1)/3+1][j]=-U[(lx-1)/3-2][j];
			V[(lx-1)/3+1][j]=-V[(lx-1)/3][j];
	}
	
	for(i=0;i<lx+1;i++)
	if((i<=(lx-1)/2-1)|(i>=(lx-1)*2/3+2))
			{
				U[i][(ly-1)*5/6+1]=-U[i][(ly-1)*5/6];
                V[i][(ly-1)*5/6+1]=-V[i][(ly-1)*5/6-1];
			}
			else
			{
				U[i][ly]=-U[i][ly-1];
		     	V[i][ly]=-V[i][ly-2];
			}
    for(j=(ly-1)*5/6+2;j<ly+1;j++)
	{
		V[(lx-1)/2][j]=-V[(lx-1)/2+1][j];
		U[(lx-1)*2/3+1][j]=-U[(lx-1)*2/3-1][j];
        V[(lx-1)*2/3+1][j]=-V[(lx-1)*2/3-2][j];
	}			                               //上下边界
	U[(lx-1)/3+1][(ly-1)/6]=U[(lx-1)/3*2+1][(ly-1)/6*5+1]=0;
	V[(lx-1)/6][(ly-1)/6-1]=V[(lx-1)/3+1][(ly-1)/6-1]=V[(lx-1)/2+1][(ly-1)/5*6+1]=V[(lx-1)/3*2+1][(ly-1)/6*5+1]=0;
}
/*-------------------------------------
上下压力边界条件
--------------------------------------*/
void BoundInP()
{
	int i,j;

	for(j=0;j<(ly-1)/6+1;j++)
		for(i=0;i<lx+1;i++)
		{
			if((i<=(lx-1)/6)|(i>=(lx-1)/3+1))
			{
				P[i][j]=0;
				P[i][(ly-1)/6]=P[i][(ly-1)/6+1];
			}
			else
				P[i][0]=P[i][1];
    	P[(lx-1)/6][j]=P[(lx-1)/6+1][j];
		P[(lx-1)/3+1][j]=P[(lx-1)/3][j];
		}
    for(j=(ly-1)*5/6+1;j<ly+1;j++)
		for(i=0;i<lx+1;i++)
		{
			if((i<=(lx-1)/2)|(i>=(lx-1)*2/3+1))
			{
				P[i][j]=0;
	            P[i][(ly-1)*5/6+1]=P[i][(ly-1)*5/6];
			}
			else
				P[i][ly]=P[i][ly-1];
		P[(lx-1)/2][j]=P[(lx-1)/2+1][j];
		P[(lx-1)*2/3+1][j]=P[(lx-1)*2/3][j];
		}
}

/*----------------------------------------
初始条件 
------------------------------------------*/
void Init()
{
	int i,j;
	for(i=0;i<lx+1;i++)
		for(j=0;j<ly+1;j++)
		{
			U[i][j]=V[i][j]=Fvy1[i][j]=0.0;
			P[i][j]=1.0;
		}
}

/*-------------------------------------------------
X方向对流项算子分裂  U，V都可用这一个函数求       ！
---------------------------------------------------*/
void  ConvDivX(double Velx[lx+1][ly+1],double Fx1[lx+1][ly+1],double Fx2[lx+1][ly+1],int X1,int X2,int Y1,int Y2)//中间段v只需求到ly-1点
{
	int i,j;
	for(j=Y1;j<=Y2;j++)
	{
		Fx1[X1][j]=0.5*(3.0*(Velx[X1+1][j]-Velx[X1][j])-(Velx[X1+2][j]-Velx[X1+1][j]));
		Fx2[X2][j]=0.5*(3.0*(Velx[X2][j]-Velx[X2-1][j])-(Velx[X2-1][j]-Velx[X2-2][j]));
	}                                                                        //左右边界
    for(j=Y1;j<=Y2;j++)
		for(i=X1+1;i<X2;i++)
		{
			Fx1[i][j]=(Velx[i+1][j]+4.0*Velx[i][j]-5.0*Velx[i-1][j]-2.0*Fx1[i-1][j])/4.0;
		    Fx2[X2+X1-i][j]=(5.0*Velx[X2+X1-i+1][j]-4.0*Velx[X2+X1-i][j]-Velx[X2+X1-i-1][j]-2.0*Fx2[X2+X1-i+1][j])/4.0;
		}
}
/*-------------------------------------------------
Y方向对流项算子分裂  U，V都可用这一个函数求   
---------------------------------------------------*/
void  ConvDivY(double V1[lx+1][ly+1],double Fy1[lx+1][ly+1],double Fy2[lx+1][ly+1],int X1,int X2,int Y1,int Y2)//1起点，2为终点
{
	int i,j;
	for(i=X1;i<=X2;i++)
	{
		Fy1[i][Y1]=0.5*(3.0*(V1[i][Y1+1]-V1[i][Y1])-(V1[i][Y1+2]-V1[i][Y1+1]));
		Fy2[i][Y2]=0.5*(3.0*(V1[i][Y2]-V1[i][Y2-1])-(V1[i][Y2-1]-V1[i][Y2-2]));
	} 
	for(i=X1;i<=X2;i++)
		for(j=Y1+1;j<Y2;j++)   
		{
		Fy1[i][j]=(V1[i][j+1]+4.0*V1[i][j]-5.0*V1[i][j-1]-2.0*Fy1[i][j-1])/4.0;
		Fy2[i][Y2+Y1-j]=(5.0*V1[i][Y2+Y1-j+1]-4*V1[i][Y2+Y1-j]-V1[i][Y2+Y1-j-1]-2.0*Fy2[i][Y2+Y1-j+1])/4.0;
		}
}
/*-------------------------------------------
求解X方向粘性项  追赶法    (端点出S为0)
--------------------------------------------*/
void ViscX(double U[lx+1][ly+1],double Sx[lx+1][ly+1],int X1,int X2,int Y1,int Y2)
{
	int i,j;
	double f[lx+1],B[lx+1],W[lx+1];  //临时变量 f：存储等号左边的数,B：存储追赶法第一步所求系数，W：存储追赶法LU分解L所对应的解
	for(j=Y1;j<=Y2;j++)
	{
		B[X1]=1.0/10.0;
		W[X1]=(U[X1+1][j]-2.0*U[X1][j]+U[X1-1][j])/5.0*6.0;
		for(i=X1+1;i<=X2;i++)
		{
			B[i]=1.0/(10.0-B[i-1]);
			f[i]=U[i+1][j]-2.0*U[i][j]+U[i-1][j];
			W[i]=(f[i]-1.0/12.0*W[i-1])/(5.0/6.0-1.0/12.0*B[i-1]);
		}
		Sx[X2][j]=W[X2];
		for(i=X2-1;i>=X1;i--)
			Sx[i][j]=W[i]-B[i]*Sx[i+1][j];
	}
}
/*-------------------------------------------
求解Y方向粘性项  追赶法   端点出S为0)
--------------------------------------------*/
void ViscY(double U[lx+1][ly+1],double Sy[lx+1][ly+1],int X1,int X2,int Y1,int Y2)
{
	int i,j;
	double f[ly+1],B[ly+1],W[ly+1];  //临时变量 f：存储等号左边的数,B：存储追赶法第一步所求系数，W：存储追赶法LU分解L所对应的解
	for(i=X1;i<=X2;i++)
	{
		B[Y1]=1.0/10.0;
		W[Y1]=(U[i][Y1+1]-2.0*U[i][Y1]+U[i][Y1-1])/5.0*6.0;
		for(j=Y1+1;j<=Y2;j++)
		{
			B[j]=1.0/(10.0-B[j-1]);
			f[j]=U[i][j+1]-2.0*U[i][j]+U[i][j-1];
			W[j]=(f[j]-1.0/12.0*W[j-1])/(5.0/6.0-1.0/12.0*B[j-1]);
		}
		Sy[i][Y2]=W[Y2];
		for(j=Y2-1;j>=Y1;j--)
			Sy[i][j]=W[j]-B[j]*Sy[i][j+1];
	}
}
/*-------------------------------------------
对流项 粘性项求解器（X方向分三块区域，Y方向分四块区域
--------------------------------------------*/
void SolveFS()
{	
    ConvDivX(U,Fux1,Fux2,0,lx,(ly-1)/6+1,(ly-1)*5/6);
	ConvDivX(V,Fvx1,Fvx2,0,lx,(ly-1)/6+1,(ly-1)*5/6-1);
    ConvDivX(U,Fux1,Fux2,(lx-1)/6,(lx-1)/3,1,(ly-1)/6);
	ConvDivX(V,Fvx1,Fvx2,(lx-1)/6,(lx-1)/3+1,1,(ly-1)/6);
    ConvDivX(U,Fux1,Fux2,(lx-1)/2,(lx-1)*2/3,(ly-1)*5/6+1,ly-1);
	ConvDivX(V,Fvx1,Fvx2,(lx-1)/2,(lx-1)*2/3+1,(ly-1)*5/6,ly-2); //X方向分三块
	ConvDivY(U,Fuy1,Fuy2,0,(lx-1)/6,(ly-1)/6,(ly-1)*5/6+1);
	ConvDivY(V,Fvy1,Fvy2,0,(lx-1)/6,(ly-1)/6,(ly-1)*5/6);
	ConvDivY(U,Fuy1,Fuy2,(lx-1)/6+1,(lx-1)/3-1,0,(ly-1)*5/6+1);
	ConvDivY(V,Fvy1,Fvy2,(lx-1)/6+1,(lx-1)/3,0,(ly-1)*5/6);
	ConvDivY(U,Fuy1,Fuy2,(lx-1)/3,(lx-1)/2,(ly-1)/6,(ly-1)*5/6+1);
	ConvDivY(V,Fvy1,Fvy2,(lx-1)/3+1,(lx-1)/2,(ly-1)/6,(ly-1)*5/6);
	ConvDivY(U,Fuy1,Fuy2,(lx-1)/2+1,(lx-1)*2/3-1,(ly-1)/6,ly);
	ConvDivY(V,Fvy1,Fvy2,(lx-1)/2+1,(lx-1)*2/3,(ly-1)/6,ly-1);
	ConvDivY(U,Fuy1,Fuy2,(lx-1)*2/3,lx,(ly-1)/6,(ly-1)*5/6+1);
	ConvDivY(V,Fvy1,Fvy2,(lx-1)*2/3+1,lx,(ly-1)/6,(ly-1)*5/6);//Y方向分五块
	ViscX(U,Sux,1,lx-1,(ly-1)/6+1,(ly-1)*5/6);
	ViscX(V,Svx,1,lx-1,(ly-1)/6+1,(ly-1)*5/6-1);
	ViscX(U,Sux,(lx-1)/6+1,(lx-1)/3-1,1,(ly-1)/6);
	ViscX(V,Svx,(lx-1)/6+1,(lx-1)/3,1,(ly-1)/6);
    ViscX(U,Sux,(lx-1)/2+1,(lx-1)*2/3-1,(ly-1)*5/6+1,ly-1);
	ViscX(V,Svx,(lx-1)/2+1,(lx-1)*2/3,(ly-1)*5/6,ly-2);//X方向分3块
    ViscY(U,Suy,1,(lx-1)/6,(ly-1)/6+1,(ly-1)*5/6);
	ViscY(V,Svy,1,(lx-1)/6,(ly-1)/6+1,(ly-1)*5/6-1);
	ViscY(U,Suy,(lx-1)/6+1,(lx-1)/3-1,1,(ly-1)*5/6);
	ViscY(V,Svy,(lx-1)/6+1,(lx-1)/3,1,(ly-1)*5/6-1);
	ViscY(U,Suy,(lx-1)/3,(lx-1)/2,(ly-1)/6+1,(ly-1)*5/6);
	ViscY(V,Svy,(lx-1)/3+1,(lx-1)/2,(ly-1)/6+1,(ly-1)*5/6-1);
	ViscY(U,Suy,(lx-1)/2+1,(lx-1)*2/3-1,(ly-1)/6+1,ly-1);
	ViscY(V,Svy,(lx-1)/2+1,(lx-1)*2/3,(ly-1)/6+1,ly-2);
	ViscY(U,Suy,(lx-1)*2/3,lx-1,(ly-1)/6+1,(ly-1)*5/6);
	ViscY(V,Svy,(lx-1)*2/3+1,lx-1,(ly-1)/6+1,(ly-1)*5/6-1);//Y方向分五块*/
}
/*-------------------------------------------
速度求解器(先进行算子分裂和粘性项求解再进行这一步）为什么u和v中的suy，svy形式不一样
---------------------------------------------*/
	
void SolveUV(int UX1,int UX2,int UY1,int UY2,int VX1,int VX2,int VY1,int VY2)//1为U，2为V
{
	int i,j;
    double UTemp[lx+1][ly+1],VTemp[lx+1][ly+1];
double uav,vav,adv,adv1,adv2,adv3,adv4,prs,vis;
		for(i=UX1;i<=UX2;i++)
		for(j=UY1;j<=UY2;j++)
		{
			vav=(V[i+1][j-1]+V[i+1][j+1]+V[i-1][j-1]+V[i-1][j+1])/4.0;
			adv1=(U[i][j]+fabs(U[i][j]))/2.0*Fux1[i][j]/dx;
			adv2=(U[i][j]-fabs(U[i][j]))/2.0*Fux2[i][j]/dx;
			adv3=(vav+fabs(vav))/2.0*Fuy1[i][j]/dy;
			adv4=(vav-fabs(vav))/2.0*Fuy2[i][j]/dy;
			adv=-adv1-adv2-adv3-adv4;
			vis=(Sux[i][j]/dx/dx+Suy[i][j]/dy/dy)/Re;
            prs=-(P[i+1][j]-P[i][j])/dx;
            UTemp[i][j]=(adv+vis+prs)*Dt+U[i][j];
		}
	for(i=VX1;i<=VX2;i++)
		for(j=VY1;j<=VY2;j++)
		{
			uav=(U[i-1][j+1]+U[i-1][j-1]+U[i+1][j+1]+U[i+1][j-1])/4.0;
		     adv1=(uav+fabs(uav))/2.0*Fvx1[i][j]/dx;
			adv2=(uav-fabs(uav))/2.0*Fvx2[i][j]/dx;
             adv3=(V[i][j]+fabs(V[i][j]))/2.0*Fvy1[i][j]/dy;
			adv4=(V[i][j]-fabs(V[i][j]))/2.0*Fvy2[i][j]/dy;
			adv=-adv1-adv2-adv3-adv4;
			vis=(Svx[i][j]/dx/dx+Svy[i][j]/dy/dy)/Re;
			prs=-(P[i][j+1]-P[i][j])/dy;
			VTemp[i][j]=(adv+vis+prs)*Dt+V[i][j];
		}
	for(i=UX1;i<=UX2;i++)
		for(j=UY1;j<=UY2;j++)
			U[i][j]=UTemp[i][j];
    for(i=VX1;i<=VX2;i++)
		for(j=VY1;j<=VY2;j++)
			V[i][j]=VTemp[i][j];
}
/*-------------------------------------------
压力求解器(先进行算子分裂和粘性项求解再进行这一步）
---------------------------------------------*/
void SolveP(int X1,int X2,int Y1,int Y2)
{
	int i,j;
	for(i=X1;i<=X2;i++)
		for(j=Y1;j<=Y2;j++)
			P[i][j]=P[i][j]-Dt*C*C*((U[i][j]-U[i-1][j])/dx+(V[i][j]-V[i][j-1])/dy);
}
/*-------------------------------------------------------
判断收敛
--------------------------------------------------------*/
void Judge(int X1,int X2,int Y1,int Y2)
{
	int i,j;
	double er;
	dm=0;
    for(i=X1;i<=X2;i++)
		for(j=Y1;j<=Y2;j++)	
		{
			er=fabs((U[i][j]-U[i-1][j])/dx+(V[i][j]-V[i][j-1])/dy);
			if(dm<er)dm=er; 
		}
}
/*--------------------------------------------------
输出函数
----------------------------------------------------*/
void output(const char *sn,double U[lx+1][ly+1],double V[lx+1][ly+1],double P[lx+1][ly+1])
{
	int i,j;
	FILE *fp;
	fp=fopen(sn,"w");
	fprintf(fp,"Title=\"computation result\"\n");
	fprintf(fp,"vars=\"x\",\"y\",\"u\",\"v\",\"p\"\n");
	fprintf(fp,"zone t=\"zone1\",I=%d,J=%d,F=POINT\n",ly,lx);
	for(i=1;i<=lx;i++)
		for(j=1;j<=ly;j++)
			fprintf(fp,"%16f%16f%20e%20e%20e\n",(i-1)*dx,(j-1)*dy,(U[i-1][j]+U[i][j])/2,(V[i][j-1]+V[i][j])/2,P[i][j]);
	fclose(fp);
}
/*--------------------------------------------------
主函数
----------------------------------------------------*/
int main()
{
	int n;
	Init();
    BoundInP();
	BoundInUV();	
		n=1;
	do
		{
			SolveFS(); 
			SolveUV(1,lx-1,(ly-1)/6+1,(ly-1)*5/6,1,lx-1,(ly-1)/6+1,(ly-1)*5/6-1);
			SolveUV((lx-1)/6+1,(lx-1)/3-1,1,(ly-1)/6,(lx-1)/6+1,(lx-1)/3,1,(ly-1)/6);
			SolveUV((lx-1)/2+1,(lx-1)*2/3-1,(ly-1)*5/6+1,ly-1,(lx-1)/2+1,(lx-1)*2/3,(ly-1)*5/6,ly-2);
            BoundInUV();
			SolveP(1,lx-1,(ly-1)/6+1,(ly-1)*5/6);
			SolveP((lx-1)/6+1,(lx-1)/3,1,(ly-1)/6);
			SolveP((lx-1)/2+1,(lx-1)*2/3,(ly-1)*5/6+1,ly-1);
			BoundInP();
            Judge((lx-1)/6+1,(lx-1)/3,1,(ly-1)/6);
			Judge((lx-1)/2+1,(lx-1)*2/3,(ly-1)*5/6+1,ly-1);
            Judge(1,lx-1,(ly-1)/6+1,(ly-1)*5/6);
			printf("%5d   %f\n",n,dm);//输出循环步数和每一步误差
			n++;
		}
	while(dm>eps);
    string output_str = string("result.plt");
    output(output_str.c_str(),U,V,P);
}