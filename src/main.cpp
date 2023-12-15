#include <iostream>
#include <vector>
#include "Fvm/Fvm.h"
// #include "GaussSolver/gauss.h"
#include <cmath>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/Dense>
#include <Eigen/LU>
// #include <Eigen/Eigenvalues> 

#include <fstream>
#include <omp.h>

/**
 * @brief The main function generates grid points, constructs faces and cells.
 * 
 * @return int 
 */
Eigen::Matrix<double, 3, 3> jacobian(double u, double c,double h,double gamma) //雅可比矩阵绝对值
{

    Eigen::Matrix<double, 3, 3> A;
    A <<0 , 1 , 0 ,
    0.5*(gamma-3)*u*u , (3-gamma)*u , gamma-1 ,
    u*((gamma-2)*0.5*u*u-c*c/(gamma-1)) , c*c/(gamma-1)+0.5*(3-2*gamma)*u*u , gamma*u;
    return A;

}

Eigen::Matrix<double, 3, 3> jacobian_abs(double u, double c,double h,double gamma) //雅可比矩阵绝对值
{
    Eigen::Matrix<double, 3, 3> lambda_abs;//特征值
    lambda_abs <<
    abs(u), 0 , 0,
    0,abs(u-c),0,
    0,0,abs(u+c);

    Eigen::Matrix<double, 3, 1> R1;
    R1 <<
    1,u,0.5*u*u;

    Eigen::Matrix<double, 3, 1> R2;
    R2 <<
    1,u-c,h-u*c;

    Eigen::Matrix<double, 3, 1> R3;
    R3 <<
    1,u+c,h+u*c;

    Eigen::Matrix<double, 3, 3> R;  
    R.col(0) = R1;
    R.col(1) = R2;
    R.col(2) = R3;

    Eigen::Matrix<double, 3, 3> A_abs;
    A_abs = R*lambda_abs*R.inverse();
    return A_abs;
}


Eigen::Matrix<double, 3, 1> W2U(double rho, double velocity, double P, double gamma) //基本量转守恒量
{
    Eigen::Matrix<double, 3, 1> U;
    U(0) = rho;
    U(1) = rho*velocity;
    U(2) = P/(gamma-1)+0.5*rho*velocity*velocity;
    return U;
}

// Eigen::Matrix<double, 3, 1> U2W(double U[3], double gamma) //守恒量转基本量
// {
//     Eigen::Matrix<double, 3, 1> W;
//     W(0) = U[0];
//     W(1) = U[1]/U[0];
//     W(2) = (gamma-1)*(U[2]-0.5*U[1]*U[1]/U[0]);
//     return W;
// }

Eigen::Matrix<double, 3, 1> W2F(double rho, double velocity, double P, double gamma) //基本量转通量
{
    Eigen::Matrix<double, 3, 1> F;
    F(0) = rho*velocity;
    F(1) = rho*velocity*velocity+P;
    // F[2] = velocity*(0.5*rho*velocity*velocity+P*gamma/(gamma-1));
    F(2) = velocity*(P*gamma/(gamma - 1)+0.5*rho*velocity*velocity);
    return F;
}

Eigen::Matrix<double, 3, 1> Lax_Wendroff(double rho_l, double rho_r, double u_l, double u_r, double P_l, double P_r, double gamma, double deltaT, double deltaX) //Lax-Wendroff格式
{   


            //L-W格式
            //从基本量构造通量
			Eigen::Matrix<double, 3,1> F0 = W2F(rho_l, u_l, P_l, gamma);
            Eigen::Matrix<double, 3,1> F1 = W2F(rho_r, u_r, P_r, gamma);
            //从基本量构造守恒量,1
            Eigen::Matrix<double, 3,1> U0 = W2U(rho_l, u_l, P_l, gamma);
            Eigen::Matrix<double, 3,1> U1 = W2U(rho_r, u_r, P_r, gamma);
			

            Eigen::Matrix<double, 3, 1> U_half;//截面处对流项半时间步长
            U_half=0.5*(U0+U1) - deltaT/(2*deltaX)*(F1-F0);

            double velocity = U_half(1)/U_half(0);
            double rho = U_half(0);
            double P = (gamma-1)*(U_half(2) -0.5*U_half(1)*U_half(1)/U_half(0));
            Eigen::Matrix<double, 3, 1> F_half = W2F(rho, velocity, P, gamma);
            return F_half;
}


Eigen::Matrix<double, 3, 1> Roe(double rho_l, double rho_r, double u_l, double u_r, double P_l, double P_r, double gamma, double deltaT, double deltaX) //Roe格式
{   
            //Roe格式
            //从基本量构造通量
            Eigen::Matrix<double, 3,1> F0 = W2F(rho_l, u_l, P_l, gamma);
            Eigen::Matrix<double, 3,1> F1 = W2F(rho_r, u_r, P_r, gamma);
            //从基本量构造守恒量,1
            Eigen::Matrix<double, 3,1> U0 = W2U(rho_l, u_l, P_l, gamma);
            Eigen::Matrix<double, 3,1> U1 = W2U(rho_r, u_r, P_r, gamma);
            double c_l = sqrt(gamma*P_l/rho_l);
            double c_r = sqrt(gamma*P_r/rho_r);
            
            double H_r =  (U1(2)+P_r)/U1(0);
            double H_l = (U0(2)+P_l)/U0(0);
            double rho = pow(((sqrt(rho_l)+sqrt(rho_r))/2),2);
            double u = (sqrt(rho_l)*u_l+sqrt(rho_r)*u_r)/(sqrt(rho_l)+sqrt(rho_r));
            double H = (sqrt(rho_l)*H_l+sqrt(rho_r)*H_r)/(sqrt(rho_l)+sqrt(rho_r));
            double c = sqrt((gamma-1)*(H-0.5*u*u));
            double P = rho*(gamma-1)/gamma*(H-0.5*u*u);
            // Eigen::Matrix<double, 3, 1> W_r ;
            // W_r << rho_r, u_r, P_r;
            // Eigen::Matrix<double, 3, 1> W_l ;
            // W_l << rho_l, u_l, P_l;
            Eigen::Matrix<double, 3, 3> A = jacobian_abs(u, c, H , gamma);
            // Eigen::Matrix<double, 3, 3> A_l = jacobian(u_l, c_l, H_l , gamma);
            // Eigen::Matrix<double, 3, 3> A_r = jacobian(u_r, c_r, H_r , gamma);
            // Eigen::Matrix<double, 3,1> U = W2U(rho, u, P, gamma);

            Eigen::Matrix<double, 3, 1> F_half;

            F_half = 0.5*(F0+F1)-0.5*A*(U1-U0);//L-W

            // Eigen::Matrix<double, 3, 1> F_half = 0.5*A*(U0+U1)-0.5*deltaT/deltaX*A*A*(U1-U0);//L-W
            return F_half;
}

Eigen::Matrix<double, 3, 1> L_F(double rho_l, double rho_r, double u_l, double u_r, double P_l, double P_r, double gamma, double deltaT, double deltaX) //LF格式
{
            //从基本量构造通量
            Eigen::Matrix<double, 3,1> F0 = W2F(rho_l, u_l, P_l, gamma);
            Eigen::Matrix<double, 3,1> F1 = W2F(rho_r, u_r, P_r, gamma);
            //从基本量构造守恒量,1
            Eigen::Matrix<double, 3,1> U0 = W2U(rho_l, u_l, P_l, gamma);
            Eigen::Matrix<double, 3,1> U1 = W2U(rho_r, u_r, P_r, gamma);
            double c_l = sqrt(gamma*P_l/rho_l);
            double c_r = sqrt(gamma*P_r/rho_r);
            
            double H_r =  (U1(2)+P_r)/U1(0);
            double H_l = (U0(2)+P_l)/U0(0);
            double rho = pow(((sqrt(rho_l)+sqrt(rho_r))/2),2);
            double u = (sqrt(rho_l)*u_l+sqrt(rho_r)*u_r)/(sqrt(rho_l)+sqrt(rho_r));
            double H = (sqrt(rho_l)*H_l+sqrt(rho_r)*H_r)/(sqrt(rho_l)+sqrt(rho_r));
            double c = sqrt((gamma-1)*(H-0.5*u*u));
            Eigen::Matrix<double, 3, 3> A = jacobian(u, c, H , gamma);

            Eigen::Matrix<double, 3, 1> F_half;
            
            F_half = 0.5*A*(U0+U1)-0.5*deltaX/deltaT*(U1-U0);//L-f

            // if ((u-c)*u*(u+c) < 0)
            // {
            //     F_half = A*U1;
            // }
            // else if ((u-c)*u*(u+c) == 0)
            // {
            //     F_half = 0.5*A*(U0+U1)-0.5*deltaX/deltaT*(U1-U0);  
            // }
            // else
            // {
            //     F_half = A*U0;
            // }
            // F_half = 0.5*(F0+F1)-0.5*A*(U1-U0);//L-W
            return F_half;
}

int main(){
    //生成一维网格点
    const int number=400;
    // std::cout << "输入节点数: ";
    // std::cin >> number;
    std::vector<Point> mypoints(number);
    int count = 0;
    //定义
    double positionRange=1; //定义实际跨越范围
    double deltaT = 0.0002 ;//定义时间步长
    double totalT = 0.5; //定义总时间
    double num = 500;
    double Rg = 287; //气体常数
    double gamma = 1.4;
    // std::cout << "输入实际跨越范围: ";
    // std::cin >> positionRange;
    while (count < number) {
        mypoints[count].position = positionRange/(number-1)*count;
        // std::cout << mypoints[count].position << std::endl;
        count++;
    }
    
    
    //生成网格点
    std::vector<Face> faces(number);//定义面
    int test = faces.size();
    // Face faces[number]; 
    pointToFace(mypoints,faces);//构造面
    std::vector<Cell> cells(number-1);
    // Cell cells[number-1];//定义体
    pointToCell(mypoints,faces,cells);//构造体
    std::cout << "网格生成完成"<<std::endl;

    /*初始化*/
    for (int i = 0; i < number; i++) {
        faces[i].velocity = 0;
        cells[i].P = 0.1;
        cells[i].velocity = 0;
        cells[i].rho =0.125;
        Eigen::Matrix<double, 3, 1> temp = W2U(cells[i].rho,cells[i].velocity,cells[i].P,gamma);;
        cells[i].U[0] = temp(0);
        cells[i].U[1] = temp(1);
        cells[i].U[2] = temp(2);
        std::cout <<","<<i<<":"<<cells[i].P ;
    }
    // Define boundary conditions (first type boundary)



    for (int i = 0; i < number/2; i++) {


        cells[i].P = 3;
        cells[i].velocity = 0;
        cells[i].rho =3;
        Eigen::Matrix<double, 3,1> temp = W2U(cells[i].rho,cells[i].velocity,cells[i].P,gamma);;
        cells[i].U[0] = temp(0);
        cells[i].U[1] = temp(1);
        cells[i].U[2] = temp(2);
        // cells[i].P = 201325;
    }
    
    // for (int i = number; i > number/2; i--) {

    //     cells[i].velocity = 1;
    // }
    // for (int i = 0; i < number-1; i++) {

    //     cells[i].velocity = sin(2*3.1415926*cells[i].position);
    // }
    std::cout <<"初始化完成"<<std::endl;
   


    for(int i = 0; i < num; i++)//总时间迭代
    {
        std::cout <<"当前时间步："<<i*deltaT<<std::endl;

        // cells[0].P = 101325;
        // cells[0].T = 300;
        // cells[0].velocity = 0.01;
        // cells[0].velocity = cells[1].velocity;
        // cells[0].P = cells[1].P;
        // cells[0].rho = cells[1].rho ;
        // cells[0].rho = cells[0].P/Rg/cells[0].T ;
        // cells[0].velocity = 1 ;
        // cells[number-2].P = 101325;
        // cells[number-2].T = 300;
        // cells[number-2].P = cells[number-3].P;
        // cells[number-2].rho = cells[number-3].rho;
        // cells[number-2].velocity = cells[number-3].velocity;
        // faces[0].velocity = 1;
        // Eigen::Matrix<double,number-1 , 1>  temp;
        for(int j = 1; j < number-2; j++)//面迭代
        {
            double deltaX = abs(cells[j-1].position-cells[j].position);
            double div[2] = {(faces[j].position - cells[j-1].position), (faces[j].position - cells[j].position)};
            double normal[2] = {div[0] / fabs(div[0]),
                                div[1] / fabs(div[1])};
                                                            
            double P_r = cells[j].P;
            double P_l = cells[j-1].P;
            double rho_r = cells[j].rho;
            double rho_l = cells[j-1].rho;
            double u_r = cells[j].velocity;
            double u_l = cells[j-1].velocity;

            // Eigen::Matrix<double, 3, 1> F_half = Lax_Wendroff (rho_l, rho_r, u_l, u_r, P_l, P_r, gamma, deltaT, deltaX); 
            // Eigen::Matrix<double, 3, 1> F_half = Roe(rho_l, rho_r, u_l, u_r, P_l, P_r, gamma, deltaT, deltaX);
            Eigen::Matrix<double, 3, 1> F_half = L_F(rho_l, rho_r, u_l, u_r, P_l, P_r, gamma, deltaT, deltaX);  
            faces[j].F[0] =F_half(0);
            faces[j].F[1] =F_half(1);
            faces[j].F[2] =F_half(2);
        }


        int j =0;
        std::cout <<"单元速度："<<cells[j].velocity<<",";
        for( j = 1; j < number-3; j++)//单元迭代
        {

            double div[2]= {(cells[j].position - cells[j].faces[0]->position), (cells[j].position - cells[j].faces[1]->position)};
            double normal[2]={div[0] / fabs(div[0])*cells[j].faces[0]->area,
                                div[1] / fabs(div[1])*cells[j].faces[1]->area};
            Eigen::Matrix<double, 3, 1> F0;
            F0<<cells[j].faces[0]->F[0],cells[j].faces[0]->F[1],cells[j].faces[0]->F[2];
            Eigen::Matrix<double, 3, 1> F1;
            F1<<cells[j].faces[1]->F[0],cells[j].faces[1]->F[1],cells[j].faces[1]->F[2];

            // 获取当前单元的物理量
            // 计算守恒量 U
            Eigen::Matrix<double, 3, 1> U;
            U << cells[j].U[0],cells[j].U[1],cells[j].U[2];
            U  = U - deltaT/cells[j].volume*(F1-F0);
            // cells[j].velocity = cells[j].velocity + deltaT/cells[j].volume*(normal[0]*F0+normal[1]*F1);
            cells[j].U[0]=U(0);
            cells[j].U[1]=U(1);
            cells[j].U[2]=U(2);
            cells[j].rho = U(0);
            cells[j].velocity = U(1)/U(0);  
            cells[j].P = (gamma-1)*(U(2)-0.5*U(1)*U(1)/U(0));
            std::cout<<cells[j].velocity<<",";
        }
        j =number-2;
        std::cout<<cells[j].velocity<<std::endl;

        



        std::fstream f;
        f.open("cell_V.csv",std::ios::out|std::ios::app);
        // 遍历 vector，并将每个 Cell 的 position 和 physics 写入文件
        f << i*deltaT ;
        for (const auto& cell : cells) {
            f << ","<< cell.velocity ;
        }
        f <<std::endl;
        // 关闭文件
        f.close();
        



        f.open("cell_rho.csv",std::ios::out|std::ios::app);
        // 遍历 vector，并将每个 Cell 的 position 和 physics 写入文件
        f << i*deltaT ;
        for (const auto& cell : cells) {
            f << ","<< cell.rho ;
        }
        f <<std::endl;
        // 关闭文件
        f.close();



        f.open("cell_P.csv",std::ios::out|std::ios::app);
        // 遍历 vector，并将每个 Cell 的 position 和 physics 写入文件
        f << i*deltaT ;
        for (const auto& cell : cells) {
            f << ","<< cell.P ;
        }
        f <<std::endl;
        // 关闭文件
        f.close();

    }

    // 写入 position 和 physics 到文件


    return 0;
}
















