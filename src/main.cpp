#include <iostream>
#include <vector>
#include "Fvm.h"
#include "gauss.h"
#include <cmath>
#include <fstream>
#include <windows.h>

/**
 * @brief The main function generates grid points, constructs faces and cells.
 * 
 * @return int 
 */

int main(){
    SetConsoleOutputCP(65001);/*设置输出格式为utf-8*/
    //生成一维网格点
    int number = 10;
    std::vector<Point> mypoints(number);
    int count = 0;
    //生成点的坐标
    while (count < number) {
        mypoints[count].position = count/5.0;
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

    test = cells[cells.size()].position;

    /*初始化*/
    for (int i = 0; i < number; i++) {
        faces[i].physics = 0;
        cells[i].physics = 0;
        std::cout << i<<":"<<cells[i].physics<<std::endl;
    }

    // Define boundary conditions (first type boundary)
    faces[0].physics = 1;
    faces[number-1].physics = 1;
    
    std::cout <<"初始化完成"<<std::endl;
   
    // double K[number-1][number-1] = {}; //构造系数矩阵
    std::vector<std::vector<double>> K(number - 1, std::vector<double>(number - 1, 0.0));
    // double F[number-1]={};//构造右端项
    std::vector<double> F(number - 1, 0.0);
    


    /*初始化边界条件，通过修改右端项的第一行与最后一行实现,其中point[0]=2,point[number]=0*/
    {
        int i = 0;
        double div[2]={(faces[i].position - cells[i].position),(cells[i+1].position - cells[i].position)};
        double normal[2] = {div[0]/fabs(div[0]),div[1]/fabs(div[1])};
        K[i][i] = -1/div[0]*normal[0] + -1/div[1]*normal[1];//系数矩阵
        K[i][i+1] = 1/div[1]*normal[1];
        double f_In = 1/div[0] * normal[0] * faces[i].physics;
        F[i] = 2*cells[i].volume - f_In; // assign a value to F[i]
        std::cout <<"定义边界向量系数:"<<F[i]<<std::endl;
        
    }
    {
        int i = number - 2 ;
        double div[2]={(cells[i-1].position - cells[i].position),(faces[i+1].position - cells[i].position)};
        double normal[2] = {div[0]/fabs(div[0]),div[1]/fabs(div[1])};
        K[i][i] = -1/div[0]*normal[0] + -1/div[1]*normal[1];
        K[i][i-1] = 1/div[0] * normal[0];
        double f_Out = 1/div[1] * faces[i+1].physics;
        F[i] = 2*cells[i].volume - f_Out; // assign a value to F[i]
        std::cout <<"定义边界向量系数:"<<F[i]<<std::endl;
    }


    /*初始化内部单元系数矩阵*/
    for(int i = 1; i < number-2; i++)
    {   
        double div[2]={(cells[i-1].position - cells[i].position),(cells[i+1].position - cells[i].position)};
        double normal[2] = {div[0]/fabs(div[0]),div[1]/fabs(div[1])};
        K[i][i] =  -1/div[0]*normal[0] + -1/div[1]*normal[1];
        K[i][i-1] = 1/div[0] * normal[0];
        K[i][i+1] = 1/div[1] * normal[1];
        F[i] = 2*cells[i].volume; // assign a value to F[i]
        std::cout <<"方阵系数:"<< K[i][i]<<"向量系数:"<<F[i]<<std::endl;
    }
    std::cout << "构造系数矩阵完成"<<std::endl; 

    


    /*求解矩阵*/
    gaussSolver(K,F,cells);



    // 写入 position 和 physics 到文件
    std::ofstream outFile("cell_data.csv");
    if (!outFile.is_open()) {
        std::cerr << "Error opening file for writing!" << std::endl;
        return 1;
    }
    outFile << "Position,Physics" <<std::endl;
    // 遍历 vector，并将每个 Cell 的 position 和 physics 写入文件
    for (const auto& cell : cells) {
        outFile << cell.position << ","<< cell.physics <<std::endl;
    }

    // 关闭文件
    outFile.close();

    std::cout << "Data written to file successfully." << std::endl;

    return 0;
}
















