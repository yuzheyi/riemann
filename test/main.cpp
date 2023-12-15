// // main.cpp
// #include <iostream>
// #include <vector>
// #include "FluidProperties/FluidProperties.h"

// int main() {
//     FluidProperties fluid(1.0, 2.0, 3.0, 100000.0, 1.2, 1.1, 1.1, 10.0, 300.0);
//     std::vector<double> position = fluid.PointPosition();
//     position[0] = position[0] + 30;
//     double *f = &position[0];
//     std::cout << "xPosition: " << *f << std::endl;

//     return 0;
// }


#include <iostream>
#include <vector>
#include "Fvm.h"
#include <cmath>


/**
 * @brief The main function generates grid points, constructs faces and cells.
 * 
 * @return int 
 */
int main() {
    //生成一维网格点
    int number = 100;
    Point mypoints[number];
    int count = 0;
    while (count < number) {
        mypoints[count].position = count;
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
    faces[number-1].physics = 2;
    

   
    double K[number-1][number-1]; //构造系数矩阵
    double F[number-1];//构造右端项

    /*初始化边界条件，通过修改右端项的第一行与最后一行实现,其中point[0]=2,point[number]=0*/
    {
        int i = 0;
        double normal[2] = {(cells[i].position-faces[i].position)/fabs(cells[i].position-faces[i].position), (cells[i].position-cells[i+1].position)/fabs(cells[i].position-cells[i+1].position)};//表示通量方向
        K[i][i] = 1/(cells[i].position-faces[i].position) * normal[0]+1/(cells[i].position-cells[i+1].position) *  normal[1];//系数矩阵
        K[i][i+1] = -1/(cells[i].position-cells[i+1].position) * normal[1];
        double f_In = -1/(cells[i].position-faces[i].position) * normal[0] * faces[i].physics;
        F[i] = 2*cells[i].volume - f_In; // assign a value to F[i]
    }
    {
        int i = number - 2 ;
        double normal[2] = {(cells[i].position-cells[i-1].position)/fabs(cells[i].position-cells[i-1].position), (cells[i].position-faces[i+1].position)/fabs(cells[i].position-faces[i+1].position)};
        K[i][i] = 1/(cells[i].position-cells[i-1].position) * normal[0] + 1/(cells[i].position-faces[i+1].position) *  normal[1];
        K[i][i-1] = -1/(cells[i].position-cells[i-1].position) * normal[0];
        double f_Out = -1/(cells[i].position-faces[i+1].position) * normal[1] * faces[i+1].physics;
        F[i] = 2*cells[i].volume - f_Out; // assign a value to F[i]
    }



    for(int i = 1; i < number-2; i++)
    {
        double normal[2] = {(cells[i].position-cells[i-1].position)/fabs(cells[i].position-cells[i-1].position), (cells[i].position-cells[i+1].position)/fabs(cells[i].position-cells[i+1].position)};
        K[i][i] = 1/(cells[i].position-cells[i-1].position) * normal[0] + 1/(cells[i].position-cells[i+1].position) *  normal[1];
        K[i][i-1] = -1/(cells[i].position-cells[i-1].position) * normal[0];
        K[i][i+1] = -1/(cells[i].position-cells[i+1].position) * normal[1];
        F[i] = 2*cells[i].volume; // assign a value to F[i]
    }

    // for(i)
    //     {
    //         double normal[0]= (cells[i].position-cells[i-1].position)/fabs(cells[i].position-cells[i-1].position);
    //         double normal[1]= (cells[i].position-cells[i+1].position)/fabs(cells[i].position-cells[i+1].position);


    //         -1/(cells[i].position-cells[i-1].position) * normal[0]//矩阵中i-1的系数cells[i-1].physics 

    //         K[][]=1/(cells[i].position-cells[i-1].position) * normal[0] + 1/(cells[i].position-cells[i+1].position) *  normal[1]//矩阵中i的系数cells[i].physics

    //         -1/(cells[i].position-cells[i+1].position) * normal[1]//矩阵中i+1的系数cells[i+1].physics



    //     }

    return 0;
}


