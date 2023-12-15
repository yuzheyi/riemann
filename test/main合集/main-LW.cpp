#include <iostream>
#include <vector>
#include "Fvm/Fvm.h"
// #include "GaussSolver/gauss.h"
#include <cmath>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <fstream>
#include <omp.h>

/**
 * @brief The main function generates grid points, constructs faces and cells.
 * 
 * @return int 
 */


Eigen::Matrix3d jacobian(double u, double rho,double P, double gamma) //雅可比矩阵
{
    double c = sqrt(gamma*P/rho);
    Eigen::Matrix3d A;
    A <<
    0 , 1 , 0,
    1/2*(gamma-3)*u*u , (3-gamma)*u , gamma-1,
    u*((gamma-2)/2*u*u-c*c/(gamma-1)) , c*c*(gamma-1)+(3-2*gamma)*u*u/2 , gamma*u;
    return A;
}


Eigen::Matrix<double, 3, 1> W2U(double rho, double velocity, double P, double gamma) //基本量转守恒量
{
    Eigen::Matrix<double, 3, 1> U;
    U(0) = rho;
    U(1) = rho*velocity;
    U(2) = P/(gamma-1)+0.5*rho*velocity*velocity;
    return U;
}

Eigen::Matrix<double, 3, 1> U2W(double U[3], double gamma) //守恒量转基本量
{
    Eigen::Matrix<double, 3, 1> W;
    W(0) = U[0];
    W(1) = U[1]/U[0];
    W(2) = (gamma-1)*(U[2]-0.5*U[1]*U[1]/U[0]);
    return W;
}

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
			Eigen::Martix<double, 3> F0 = W2F(rho_l, u_l, P_l, gamma);
            Eigen::Martix<double, 3> F1 = W2F(rho_r, u_r, P_r,, gamma);
            //从基本量构造守恒量
            Eigen::Martix<double, 3> U0 = W2U(rho_l, u_l, P_l, gamma);
            Eigen::Martix<double, 3> U1 = W2U(rho_r, u_r, P_r, gamma);
			

            Eigen::Matrix<double, 3, 1> U_half;//截面处对流项半时间步长
            U_half=0.5*(U0+U1) - deltaT/(2*deltaX)*(F1-F0);

            double velocity = U_half(1)/U_half(0);
            double rho = U_half(0);
            double P = (gamma-1)*(U_half(2) -0.5*U_half(1)*U_half(1)/U_half(0));
            Eigen::Matrix<double, 3, 1> F_half = W2F(rho, velocity, P, gamma);
            return F_half;



}

int main(){
    //生成一维网格点
    const int number=500;
    // std::cout << "输入节点数: ";
    // std::cin >> number;
    std::vector<Point> mypoints(number);
    int count = 0;
    //定义
    double positionRange=1; //定义实际跨越范围
    double deltaT = 0.001 ;//定义时间步长
    double totalT = 0.5; //定义总时间
    double num = 1000;
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
        std::array<double, 3> temp = W2U(cells[i].rho,cells[i].velocity,cells[i].P,gamma);;
        cells[i].U[0] = temp[0];
        cells[i].U[1] = temp[1];
        cells[i].U[2] = temp[2];
        std::cout <<","<<i<<":"<<cells[i].P ;
    }
    // Define boundary conditions (first type boundary)



    for (int i = 0; i < number/2; i++) {


        cells[i].P = 3;
        cells[i].velocity = 0;
        cells[i].rho =3;
        std::array<double, 3> temp = W2U(cells[i].rho,cells[i].velocity,cells[i].P,gamma);;
        cells[i].U[0] = temp[0];
        cells[i].U[1] = temp[1];
        cells[i].U[2] = temp[2];
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

            // 定义守恒量和物理量的数组
            const int numVariables = 3;  // 数组的大小(只计算速度)
            const int numCells = 2;     //依赖的网格数量
            Eigen::Matrix<double, numVariables, numCells> U;   // 守恒量构造
            Eigen::Matrix<double, numVariables, numCells> W;   // 基本量构造
            // 初始化示例矩阵 M，这里假设 M 是一个 3x1 的矩阵


            for (int k = 0; k < numCells; ++k) {
                // 当前单元的索引
                int currentCellIndex = j - 1 + k;
                // 获取当前单元的物理量
                double rho = cells[currentCellIndex].rho;
                double velocity = cells[currentCellIndex].velocity;
                double P = cells[currentCellIndex].P;
                // 计算基本量W 
                W.col(k) << rho,velocity,P;
                // 计算守恒量 U
                U.col(k) << cells[currentCellIndex].U[0],cells[currentCellIndex].U[1],cells[currentCellIndex].U[2];
                
                }
            
            // double rho_l=cells[j-1].rho;
            // double rho_r=cells[j].rho;
            // double u_l=cells[j-1].velocity;
            // double u_r=cells[j].velocity;
            // double H_l=cells[j-1].P/(gamma-1)+0.5*rho_l*u_l*u_l;
            // double H_r=cells[j].P/(gamma-1)+0.5*rho_r*u_r*u_r;

            // double u = (sqrt(rho_l)*u_l+sqrt(rho_r)*u_r)/(sqrt(rho_l)+sqrt(rho_r));
            // double H = (sqrt(rho_l)*H_l+sqrt(rho_r)*H_r)/(sqrt(rho_l)+sqrt(rho_r));
            // faces[j].velocity = u;
            // double c = sqrt((gamma-1)*(H-0.5*u*u));




            //L-W格式
			std::array<double, 3> F00 = W2F(cells[j-1].rho, cells[j-1].velocity, cells[j-1].P, gamma);
            std::array<double, 3> F11 = W2F(cells[j].rho, cells[j].velocity, cells[j].P, gamma);
            Eigen::Matrix<double, 3, 1> F0;
            F0 << F00[0],F00[1],F00[2];
            Eigen::Matrix<double, 3, 1> F1;
            F1 << F11[0],F11[1],F11[2];
			
            Eigen::Matrix<double, 3, 1> U_half;//截面处对流项半时间步长
            U_half=0.5*(U.col(0)+U.col(1)) - deltaT/(2*deltaX)*(F1-F0);

            double velocity = U_half(1)/U_half(0);
            double rho = U_half(0);
            double P = (gamma-1)*(U_half(2) -0.5*U_half(1)*U_half(1)/U_half(0));
            std::array<double, 3> F_half = W2F(rho, velocity, P, gamma);


            // // //FTCS（Forward time central space）格式:
            // Eigen::Matrix<double, numCells, 1> M1;//对流项矩阵
            // M1 << 0.5,
            // 0.5;

            //LF格式
            // double A = 0.5*(cells[j-1].velocity+cells[j].velocity);
            // Eigen::Matrix<double, numCells, 1> M1;//对流项矩阵
            // M1 << 0.5*A+0.5*div[0]/deltaT,
            // 0.5*A+0.5*div[1]/deltaT;


            //迎风格式,P/(gamma-1)+0.5*rho*velocity*velocity

            // //Roe平均
            // double rho_l=cells[j-1].rho;
            // double rho_r=cells[j].rho;
            // double u_l=cells[j-1].velocity;
            // double u_r=cells[j].velocity;
            // double H_l=cells[j-1].P/(gamma-1)+0.5*rho_l*u_l*u_l;
            // double H_r=cells[j].P/(gamma-1)+0.5*rho_r*u_r*u_r;

            // double u = (sqrt(rho_l)*u_l+sqrt(rho_r)*u_r)/(sqrt(rho_l)+sqrt(rho_r));
            // double H = (sqrt(rho_l)*H_l+sqrt(rho_r)*H_r)/(sqrt(rho_l)+sqrt(rho_r));
            // faces[j].velocity = u;
            // double c = sqrt((gamma-1)*(H-0.5*u*u));
            // Eigen::Matrix<double, 3, 1> M1;//对流项矩阵
            // // Eigen::MatrixXd AA = A;
            // if (AA.determinant() < 0)
            // {
            //     M1 = A * U.col(1);
            // }
            // else // AA >= 0
            // {
            //     M1 = A * U.col(0);
            // }
            // //Lax-Wendroff 格式
            // Eigen::Matrix<double, numCells, 1> M0;//对流项矩阵
            //     M0 << 
            //     0.5+0.25*deltaT/div[0],
            //     0.5+0.25*deltaT/div[1];
            // faces[j].velocity=F * M0;
            // double A = 0.5*(cells[j-1].velocity+cells[j].velocity);
            // Eigen::Matrix<double, 3, 1> M1;//对流项矩阵
            //     M1 =0.5*A*(U.col(1)+U.col(0))+0.5*deltaT/div[0]*A*A*U.col(0)+0.5*A*A*deltaT/div[1]*U.col(1);     
            faces[j].F[0] =F_half[0];
            faces[j].F[1] =F_half[1];
            faces[j].F[2] =F_half[2];
            faces[j].velocity = U_half(1)/U_half(0);
            faces[j].P = (gamma-1)*(U_half(2)-0.5*U_half(1)*U_half(1)/U_half(0));
            faces[j].rho = U_half(0);

        }

        // cells[0].U[0] = cells[1].U[0];
        // cells[0].U[1] = cells[1].U[1];
        // cells[0].U[2] = cells[1].U[2];
        // cells[number-2].U[0] = cells[number-3].U[0];
        // cells[number-2].U[1] = cells[number-3].U[1];
        // /cells[number-2].U[2] = cells[number-3].U[2];

        // faces[0].velocity = 1;
        // faces[0].P = faces[1].P ;
        // faces[0].rho = 1.2 ;
        // faces[0].F[0] = faces[0].velocity*faces[0].rho;
        // faces[0].F[1] = faces[0].velocity*faces[0].velocity*faces[0].rho+faces[0].P;
        // faces[0].F[2] = faces[0].velocity*(faces[0].F[1]/faces[0].rho+faces[0].P/(gamma-1));
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
            // U  = U + deltaT/cells[j].volume*(normal[0]*F0+normal[1]*F1);
            U  = U - deltaT/cells[j].volume*(F1-F0);
            // cells[j].velocity = cells[j].velocity + deltaT/cells[j].volume*(normal[0]*F0+normal[1]*F1);
            cells[j].U[0]=U(0);
            cells[j].U[1]=U(1);
            cells[j].U[2]=U(2);
            cells[j].rho = U(0);
            cells[j].velocity = U(1)/U(0);  
            cells[j].P = (gamma-1)*(U(2)-0.5*U(1)*U(1)/U(0));
            // cells[j].P = 0.5*(cells[j].faces[0]->P+cells[j].faces[1]->P);
            // double epsilon = 1e-10;

            // cells[j].rho = std::abs(U(0)) < epsilon ? 0 : U(0);
            // cells[j].velocity = std::abs(U(1)/U(0)) < epsilon ? 0 : U(1)/U(0);
            // cells[j].P = std::abs((gamma-1)*(U(2)-0.5*U(1)*U(1)/U(0))) < epsilon ? 0 : (gamma-1)*(U(2)-0.5*U(1)*U(1)/U(0));
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
















