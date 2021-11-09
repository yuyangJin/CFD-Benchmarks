#include <cstdlib>
#include<iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

using namespace std;

#define I 200    //x方向网格数
#define J 200   //y方向网格数
#define GAMA 1.4  //气体常数
#define Lx 1.0    //计算区域x方向长度
#define Ly 1.0    //计算区域y方向长度
#define TT 0.31    //计算总时间
#define CFL 0.1
#define MIN(x,y)(((x)<(y))?(x):(y))
#define MAX(x,y)(((x)>(y))?(x):(y))
#define OPDT 5
#define Step 50

//全局变量 
double U[I+2][J+2][4],EA[I+2][J+2][4],EB[I+2][J+2][4],FA[I+2][J+2][4],FB[I+2][J+2][4],UX[I+2][J+2][4]; 
double k=-1, b=1;   
double UL1[I+2][J+2][4],UL2[I+2][J+2][4],UR1[I+2][J+2][4], UR2[I+2][J+2][4];
double UL3[I+2][J+2][4],UR3[I+2][J+2][4],UR4[I+2][J+2][4], UL4[I+2][J+2][4];

using std::cin;
using std::cout;
using std::endl;
using std::setw;
using std::sqrt;
using std::pow;
using std::fabs;

/*定义初始条t和边界条件*/
void Init()            
{
	int i,j;
	double rou1=1.0,u1=0,v1=0,a1=1,p1=0.71429,rou2;
	double Ms=5.0,p2,u2,v2=0;
	p2=p1*(1+2*GAMA/(GAMA+1)*(Ms*Ms-1));
	rou2=rou1*Ms*Ms/(1+(GAMA-1)/(GAMA+1)*(Ms*Ms-1));
	u2=u1+a1*2/(GAMA+1)*(Ms-1/Ms);

	for(i=0;i<=I+1;i++)
		for(j=0;j<=J+1;j++)  //设定区域
		{
			if((i<=2))     			{
				U[i][j][0]=rou2;
                U[i][j][1]=rou2*u2;
                U[i][j][2]=rou2*v2;
	            U[i][j][3]=p2/(GAMA-1)+rou2*(u2*u2+v2*v2)/2;

			}
			else
			{
				U[i][j][0]=rou1;
                U[i][j][1]=rou1*u1;
                U[i][j][2]=rou1*v1;
	            U[i][j][3]=p1/(GAMA-1)+rou1*(u1*u1+v1*v1)/2;

			}
		}

		return ;

}

double sign(double va)
{
	if(va>0)
	{return 1;}
	if(va<0)
	{return -1;}
	else
	{return 0;}
}

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

void bound(double U[I+2][J+2][4])
{
	int i,j,k;
	double den1=1,u1=0,v1=0,a1=1,p1=0.71429;
	double Ms=5.0,den2,p2,u2,v2=0;
	p2=p1*(1+2*GAMA/(GAMA+1)*(Ms*Ms-1));
	den2=den1*Ms*Ms/(1+(GAMA-1)/(GAMA+1)*(Ms*Ms-1));
	u2=u1+a1*2/(GAMA+1)*(Ms-1/Ms);


	for(i=0;i<2;i++)//左
	{
		for(j=0;j<J+1;j++)
		{
			U[i][j][0]=den2;
			U[i][j][1]=den2*u2;
			U[i][j][2]=den2*v2;
			U[i][j][3]=p2/(GAMA-1)+0.5*den2*(u2*u2+v2*v2);
		}
	}
	for(i=I;i<=I+1;i++)//右
	{
		for(j=173;j<=J+1;j++)
		{
			U[i][j][0]=U[i-2][j][0];
			U[i][j][1]=U[i-2][j][1];
			U[i][j][2]=U[i-2][j][2];
			U[i][j][3]=U[i-2][j][3];
		}
	}

	for(i=0;i<=29;i++)  	{
		for(j=0;j<=1;j++)
		{
			U[i][j][0]=U[i][4-j][0];
			U[i][j][1]=U[i][4-j][1];
			U[i][j][2]=(-1)*U[i][4-j][2];
			U[i][j][3]=U[i][4-j][3];
		}
	}
	//设定下边界的刚性壁面条件
	for(i=0;i<=29;i++)
		U[i][2][2]=0;

		for(i=2;i<=I-1;i++)
	{
		for(j=J;j<=J+1;j++)
		{
			U[i][j][0]=U[i][2*J-2-j][0];
			U[i][j][1]=U[i][2*J-2-j][1];
			U[i][j][2]=(-1)*U[i][2*J-2-j][2];
			U[i][j][3]=U[i][2*J-2-j][3];
		}	
				
	}
		//设定上边界的刚性壁面条件
	for(i=2;i<=I-1;i++)
		U[i][J-1][2]=0;
	

		//设定楔形面上的边界条件	
	for(i=30;i<=I;i++)   //此处可以改动
	{
		for(j=0;j<i-28;j++)
			if((j>=i-30))
			{ 
			       U[i][j][0]=U[28+j][i-28][0];
			       U[i][j][1]=U[28+j][i-28][2];
			       U[i][j][2]=U[28+j][i-28][1];
			       U[i][j][3]=U[28+j][i-28][3];
		        
			}
		
	}
	
	for(i=30;i<=I-1;i++)
	{
		U[i][i-28][0]=(U[i-1][i-27][0]+U[i+1][i-29][0])*0.5;
	   U[i][i-28][1]=0.5*(U[i-1][i-27][1]+U[i+1][i-29][1]);
	U[i][i-28][2]=U[i][i-28][1];
	U[i][i-28][3]=(U[i-1][i-27][3]+U[i+1][i-29][3])*0.5;
	}
}

//时间步长
double TIME(double U[I+2][J+2][4], double dx, double dy)
{
	int i,j;
	double u,v,rou,a,vel,p,maxvel;
	maxvel=0.00000001*0.00000001;
	for(i=0;i<=I+2;i++)    
	{
		for(j=2;j<=J-1;j++)
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



void U2EA(double UL3[I+2][J+2][4],double EA[I+2][J+2][4],int i,int j)
{

	double rou,u,v,p,a,h,f1,f2,f3;
	double DD,pL,pR,hL,hR;
	int s=0;

    rou=UL3[i][j][0];
	u=UL3[i][j][1]/UL3[i][j][0];
	v=UL3[i][j][2]/UL3[i][j][0];
	p=(GAMA-1)*(UL3[i][j][3]-0.5*(u*u+v*v)*rou);  	

	a=sqrt(GAMA*p/rou);
    h=a*a/(GAMA-1)+0.5*(u*u+v*v);

    f1=0.5*(u+sqrt(u*u+0.00001*0.00001));
	f2=0.5*(u-a+sqrt((u-a)*(u-a)+0.00001*0.00001));
	f3=0.5*(u+a+sqrt((u+a)*(u+a)+0.00001*0.00001));
	
	EA[i][j][0]=rou/2/GAMA*(2*(GAMA-1)*f1+f2+f3);
	EA[i][j][1]=rou/2/GAMA*(2*(GAMA-1)*u*f1+(u-a)*f2+(u+a)*f3);
	EA[i][j][2]=rou/2/GAMA*(2*(GAMA-1)*f1+f2+f3)*v;
	EA[i][j][3]=rou/2/GAMA*((GAMA-1)*(u*u+v*v)*f1+(h-a*u)*f2+(h+a*u)*f3);
	return ;
}

/*计算f-*/
void U2EB(double UR3[I+2][J+2][4],double EB[I+2][J+2][4],int i,int j)
{

	double rou,u,v,p,a,h,f1,f2,f3;
	double DD,pL,pR,hL,hR;
	int s=0;

	rou=UR3[i][j][0];
	u=UR3[i][j][1]/UR3[i][j][0];
	v=UR3[i][j][2]/UR3[i][j][0];
	p=(GAMA-1)*(UR3[i][j][3]-0.5*rou*(u*u+v*v));
	
	a=sqrt(GAMA*p/rou);
	h=a*a/(GAMA-1)+0.5*(u*u+v*v);

    f1=0.5*(u-sqrt(u*u+0.00001*0.00001));
	f2=0.5*(u-a-sqrt((u-a)*(u-a)+0.00001*0.00001));
	f3=0.5*(u+a-sqrt((u+a)*(u+a)+0.00001*0.00001));
			
	EB[i][j][0]=rou/2/GAMA*(2*(GAMA-1)*f1+f2+f3);	
	EB[i][j][1]=rou/2/GAMA*(2*(GAMA-1)*u*f1+(u-a)*f2+(u+a)*f3);	
	EB[i][j][2]=rou/2/GAMA*(2*(GAMA-1)*f1+f2+f3)*v;	
	EB[i][j][3]=rou/2/GAMA*((GAMA-1)*(u*u+v*v)*f1+(h-a*u)*f2+(h+a*u)*f2);
	return;
}

/*计算g+*/
void U2FA(double UL4[I+2][J+2][4],double FA[I+2][J+2][4],int i,int j)
{

	double rou,u,v,p,a,h,g1,g2,g3;
	double DD,pL,pR,hL,hR;
	int s=0;

	rou=UL4[i][j][0];
	u=UL4[i][j][1]/UL4[i][j][0];
	v=UL4[i][j][2]/UL4[i][j][0];
	p=(GAMA-1)*(UL4[i][j][3]-0.5*(u*u+v*v)*rou);
	
      
	a=sqrt(GAMA*p/rou);
	h=a*a/(GAMA-1)+0.5*(u*u+v*v);

	g1=0.5*(v+sqrt(v*v+0.00001*0.00001));
	g2=0.5*(v-a+sqrt((v-a)*(v-a)+0.00001*0.00001));
	g3=0.5*(v+a+sqrt((v+a)*(v+a)+0.00001*0.00001));
	
    FA[i][j][0]=rou/2/GAMA*(2*(GAMA-1)*g1+g2+g3);
    FA[i][j][1]=rou/2/GAMA*(2*(GAMA-1)*g1+g2+g3)*u;
    FA[i][j][2]=rou/2/GAMA*(2*(GAMA-1)*g1*v+(v-a)*g2+(v+a)*g3);
	FA[i][j][3]=rou/2/GAMA*((GAMA-1)*(u*u+v*v)*g1+(h-a*v)*g2+(h+a*v)*g3);
	return;
}

/*计算g-*/
void U2FB(double UR4[I+2][J+2][4],double FB[I+2][J+2][4],int i,int j)
{
	double rou,u,v,p,a,h,g1,g2,g3;
	double DD,pL,pR,hL,hR;
	int s=0;

	rou=UR4[i][j][0];
	u=UR4[i][j][1]/UR4[i][j][0];
	v=UR4[i][j][2]/UR4[i][j][0];
	p=(GAMA-1)*(UR4[i][j][3]-0.5*(u*u+v*v)*rou);
	
	a=sqrt(GAMA*p/rou);
	h=a*a/(GAMA-1)+0.5*(u*u+v*v);

	
    g1=0.5*(v-sqrt(v*v+0.00001*0.00001));
	g2=0.5*(v-a-sqrt((v-a)*(v-a)+0.00001*0.00001));
	g3=0.5*(v+a-sqrt((v+a)*(v+a)+0.00001*0.00001));	
	
    FB[i][j][0]=rou/2/GAMA*(2*(GAMA-1)*g1+g2+g3);
    FB[i][j][1]=rou/2/GAMA*(2*(GAMA-1)*g1+g2+g3)*u;
    FB[i][j][2]=rou/2/GAMA*(2*(GAMA-1)*g1*v+(v-a)*g2+(v+a)*g3);
	FB[i][j][3]=rou/2/GAMA*((GAMA-1)*(u*u+v*v)*g1+(h-a*v)*g2+(h+a*v)*g3);

    return;
}

int Output(double U[I+2][I+2][4],int m,int n)
{
	int i,j,k;
	double den,u,v,p;
	char Name[]={"TVD000.plt"};
	FILE *fp;

	if(n%Step==0)
	{
		Name[5]=char(m%10+48);
		Name[4]=char((m/10)%10+48);
		Name[3]=char(m/100+48);
		fp=fopen(Name,"w");
		fprintf(fp,"ZONE  I=%d J=%d  ",I+2,J+2);
		fprintf(fp,"DATAPACKING=POINT\n");
	for(i=0;i<I+2;i++)
		for(j=0;j<I+2;j++)
		{
			if(U[i][j][0]!=0)
		{
			den=U[i][j][0];
	u=U[i][j][1]/U[i][j][0];
	v=U[i][j][2]/U[i][j][0];
	p=(GAMA-1)*(U[i][j][3]-0.5*den*(u*u+v*v));
	
			}
		else
		{
			den=0;
	u=0;
	v=0;
	p=0;
		}
		fprintf(fp,"%f  %f  %f  %f   %f    %f\n",(i+0.)/100,(j+0.)/100,den,u,v,p);
		}
	fclose(fp);
		return 1;
	}
	else
		return 0;
}

int main()
{
	int i,j,s,p=0;
	double m=0,n=0;
	int m1=0,n1=0;
	double dx,dy;
	dx=Lx*1./I;
	dy=Ly*1./J;
	Init();
	
	for(p=0;p<=14000;p++)	{
	
       n=TIME(U,dx,dy);
	   m=n+m;


		for(i=0;i<=I+1;i++)
		{
			for(j=0;(j<=J+1);j++)
			{
				if(j>=i-(int(0.15*I)))
				if(U[i][j][0]!=0)
				{
					U2EA(U,EA,i,j); 
				    U2EB(U,EB,i,j); 
				    U2FA(U,FA,i,j);
				    U2FB(U,FB,i,j);
				
				}
            }
		}
		
		for(i=2;i<=I-1;i++)
		{
			for(j=2;(j<=J-1);j++)
			{
				if(j>=i-(int(0.15*I)-2))
				{
					if(U[i][j][0]!=0)
				  {
				
					for(s=0;s<=3;s++)
				    {   
					   U[i][j][s]=U[i][j][s]-n*I*(EA[i][j][s]+0.5*(minmod((EA[i][j][s]-EA[i-1][j][s]),(EA[i+1][j][s]-EA[i][j][s])))
						   +EB[i+1][j][s]-0.5*(minmod((EB[i+1][j][s]-EB[i][j][s]),(EB[i+2][j][s]-EB[i+1][j][s])))
						   -(EA[i-1][j][s]+0.5*(minmod((EA[i-1][j][s]-EA[i-2][j][s]),(EA[i][j][s]-EA[i-1][j][s])))
						   +EB[i][j][s]-0.5*(minmod((EB[i][j][s]-EB[i-1][j][s]),(EB[i+1][j][s]-EB[i][j][s])))))
						   -n*J*(FA[i][j][s]+0.5*(minmod((FA[i][j][s]-FA[i][j-1][s]),(FA[i][j+1][s]-FA[i][j][s])))
						   +FB[i][j+1][s]-0.5*(minmod((FB[i][j+1][s]-FB[i][j][s]),(FB[i][j+2][s]-FB[i][j+1][s])))
						   -(FA[i][j-1][s]+0.5*(minmod((FA[i][j-1][s]-FA[i][j-2][s]),(FA[i][j][s]-FA[i][j-1][s])))
						   +FB[i][j][s]-0.5*(minmod((FB[i][j][s]-FB[i][j-1][s]),(FB[i][j+1][s]-FB[i][j][s])))));
				    }
				  }
				}
			 }
		 }
		
		bound(U);
		cout<<p<<endl;
		cout<<m<<endl;		
        n1++;
        m1=m1+Output(U,m1,n1);
    } 
	return 0;
}
