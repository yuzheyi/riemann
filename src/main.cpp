#include <iostream>
#include <vector>
#include "Fvm.h"
// #include "gauss.h"
#include <cmath>
/**
 * @brief The main function generates grid points, constructs faces and cells.
 * 
 * @return int 
 */

int main() {
    //生成一维网格点
    int number = 300;
    Point mypoints[number];
    int count = 0;
    //生成点的坐标
    while (count < number) {
        mypoints[count].position = count/30.0;
        count++;
    }


    //生成网格点
    Face faces[number]; //定义面
    pointToFace(mypoints,number,faces);//构造面
    Cell cells[number-1];//定义体
    pointToCell(mypoints,number,faces,cells);//构造体
    // pointToCell(mypoints,number);//构造单元


    /*初始化*/
    for (int i = 0; i < number; i++) {
        faces[i].physics = 0;
        cells[i].physics = 0;
    }

    // Define boundary conditions (first type boundary)
    faces[0].physics = 1;
    faces[number-1].physics = 1;
    

   
    double K[number-1][number-1] = {}; //构造系数矩阵
    double F[number-1]={};//构造右端项

    /*初始化边界条件，通过修改右端项的第一行与最后一行实现,其中point[0]=2,point[number]=0*/
    {
        int i = 0;
        double div[2]={(faces[i].position - cells[i].position),(cells[i+1].position - cells[i].position)};
        double normal[2] = {div[0]/fabs(div[0]),div[1]/fabs(div[1])};
        K[i][i] = -1/div[0]*normal[0] + -1/div[1]*normal[1];//系数矩阵
        K[i][i+1] = 1/div[1]*normal[1];
        double f_In = 1/div[0] * normal[0] * faces[i].physics;
        F[i] = 2*cells[i].volume - f_In; // assign a value to F[i]
    }
    {
        int i = number - 2 ;
        double div[2]={(cells[i-1].position - cells[i].position),(faces[i+1].position - cells[i].position)};
        double normal[2] = {div[0]/fabs(div[0]),div[1]/fabs(div[1])};
        K[i][i] = -1/div[0]*normal[0] + -1/div[1]*normal[1];
        K[i][i-1] = 1/div[0] * normal[0];
        double f_Out = 1/div[1] * faces[i+1].physics;
        F[i] = 2*cells[i].volume - f_Out; // assign a value to F[i]
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
    }


    


    // size_t width = (height > 0) ? K[0].size() : 0;
    /*求解矩阵*/
    int totalNums = sizeof(K) / sizeof(K[0][0]);//总个数
    int cols = sizeof(K) / sizeof(K[0]);//行数
    int raws = totalNums / cols;//列数
    gaussSolver(&K[0][0],&F[0],cells,cols,raws);
    // x.solve(K,F);


    return 0;
}


