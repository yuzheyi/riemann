#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define GAMA 1.4//气体常数
#define PI 3.141592654
#define L 2.0//计算区域
#define TT 0.4//总时间
#define Sf 0.8//时间步长因子//（与时间步长相关的量，用在函数CFL中）
#define J 1000//网格数

//全局变量
double U[J+2][3],Uf[J+2][3],Ef[J+2][3];
//U即A.1中的u，U[][0]=ρ，U[][1]=ρu，U[][0]=E。
//Uf和Ef为临时变量，在计算中会改变含义

/*
*计算时间步长
*入口: U，当前物理量，dx，网格宽度；
*返回: 时间步长。
*/
double CFL(double U[J+2][3],double dx)
{
    int i;
    double maxvel,p,u,vel;
    maxvel=1e-100;
    for(i=1;i<=J;i++)
    {
        u=U[i][1]/U[i][0];
        p=(GAMA-1)*(U[i][2]-0.5*U[i][0]*u*u);//（即A.3式）
        vel=sqrt(GAMA*p/U[i][0])+fabs(u);//fabs()求绝对值。？？时间步长的公式从哪来的
        if(vel>maxvel)maxvel=vel;
    }
    return Sf*dx/maxvel;
}

/*
*初始化
*入口: 无；
*出口: U， 已经给定的初始值，
*dx, 网格宽度。
*/
void Init(double U[J+2][3],double & dx)
{
    int i;
    double rou1=1.0 ,u1=0.0,p1=1.0; //初始条件
    double rou2=0.125,u2=0.0,p2=0.1;
    dx=L/J;//网格宽度赋值
    for(i=0;i<=J/2;i++)//气体1赋初值，0-500
    {
        U[i][0]=rou1;
        U[i][1]=rou1*u1;
        U[i][2]=p1/(GAMA-1)+rou1*u1*u1/2;//由A.3算出初始E值
    }
    for(i=J/2+1;i<=J+1;i++)//气体2赋初值，501-1001
    {
        U[i][0]=rou2;
        U[i][1]=rou2*u2;
        U[i][2]=p2/(GAMA-1)+rou2*u2*u2/2;
    }
}

/*
*边界条件
*入口: dx，网格宽度；
*出口: U， 已经给定的边界。
*/
void bound(double U[J+2][3],double dx)
{
    int k;
    //左边界
    for(k=0;k<3;k++)U[0][k]=U[1][k];
    //右边界
    for(k=0;k<3;k++)U[J+1][k]=U[J][k];
}

/*
*根据U计算E
*入口: U， 当前U矢量；
*出口: E， 计算得到的E矢量，
*U、E的定义见Euler方程组。即A.2中的u和f矢量
*/
void U2E(double U[3],double E[3])
{
    double u,p;
    u=U[1]/U[0];
    p=(GAMA-1)*(U[2]-0.5*U[1]*U[1]/U[0]);//即A.3式
    E[0]=U[1];
    E[1]=U[0]*u*u+p;
    E[2]=(U[2]+p)*u;
}



/*
*一维差分格式求解器
*入口: U， 上一时刻的U矢量，Uf、Ef，临时变量，
*dx，网格宽度，dt, 时间步长；
*出口: U， 计算得到的当前时刻U矢量。
*/
/*
思路：用当前的u，加入人工黏性滤波获得ū和对应的f，用ū和f计算出u(n+1/2)和f(n+1/2)，最后算出u(n+1)，时间向前推进一步。
*/

void MacCormack_1DSolver(double U[J+2][3],double Uf[J+2][3],double Ef[J+2][3],double dx,double dt)
{
    int i,k;
    double r,nu,q;
    r=dt/dx;//式A.4中的r
    nu=0.25;//即A.8中的η
    for(i=1;i<=J;i++)
    {
        q=fabs(fabs(U[i+1][0]-U[i][0])-fabs(U[i][0]-U[i-1][0]))
        /(fabs(U[i+1][0]-U[i][0])+fabs(U[i][0]-U[i-1][0])+1e-100); //开关函数A.6
        for(k=0;k<3;k++)
        {
            Ef[i][k]=U[i][k]+0.5*nu*q*(U[i+1][k]-2*U[i][k]+U[i-1][k]);//人工黏性项A.7，此时的Ef是A.7中的ū
        }
    }
    for(k=0;k<3;k++)
    for(i=1;i<=J;i++)
    {
        U[i][k]=Ef[i][k];
    }
    for(i=0;i<=J+1;i++)
    {
        U2E(U[i],Ef[i]);//此时的Ef是A.2中的f矢量
    }
    for(i=0;i<=J;i++)
    for(k=0;k<3;k++)
    {
        Uf[i][k]=U[i][k]-r*(Ef[i+1][k]-Ef[i][k]); //U(n+1/2)(i+1/2)，即A.4的第一式，此时Uf即u(n+1/2)？？后半部分Ef的i下标是否有误
    }
    for(i=0;i<=J;i++)
    {
        U2E(Uf[i],Ef[i]); //E(n+1/2)(i+1/2),用u(n+1/2)算出f(n+1/2)
    }
    for(i=1;i<=J;i++)
    for(k=0;k<3;k++)
    {
        U[i][k]=0.5*(U[i][k]+Uf[i][k])-0.5*r*(Ef[i][k]-Ef[i-1][k]); //U(n+1)(i)，即A.4的第二式，？？后半部分Ef的i下标是否有误
    }
}



/*
*输出结果, 用数据格式画图
*入口: U， 当前时刻U矢量，dx， 网格宽度；
*出口: 无。
*/
void Output(double U[J+2][3],double dx)
{
    int i;
    FILE *fp;
    double rou,u,p;
    fp=fopen("result.txt","w");
    for(i=0;i<=J+1;i++)
    {
        rou=U[i][0];
        u=U[i][1]/rou;
        p=(GAMA-1)*(U[i][2]-0.5*U[i][0]*u*u);
        fprintf(fp,"%20f%20.10e%20.10e%20.10e%20.10e\n",i*dx,rou,u,p,U[i][2]);
    }
    fclose(fp);
}


/*
*主函数
*入口: 无；
*出口: 无。
*/
int main()
{
    double T,dx,dt;
    Init(U,dx);
    T=0;
    while(T<TT)
    {
        dt=CFL(U,dx);
        T+=dt;
        printf("T=%10g dt=%10g\n",T,dt);
        MacCormack_1DSolver(U,Uf,Ef,dx,dt);
        bound(U,dx);
    }
    Output(U,dx);
}