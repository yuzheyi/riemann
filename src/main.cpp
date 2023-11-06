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
    //初始化
    for (int i = 0; i < number; i++) {
        faces[i].physics = 0;
        cells[i].physics = 0;
    }

    // Define boundary conditions (first type boundary)
    faces[0].physics = 1;
    faces[number-1].physics = 2;

    // // Calculate physical quantities for faces
    // for (int i = 0; i < number; i++) {
    //     double grad_F = (cells[i+1].physics - cells[i].physics) / (cells[i+1].position - cells[i].position);
    //     faces[i].physics = grad_F;
    // }

    // // Calculate physical quantities for cells
    // for (int i = 0; i < number-1; i++) {
    //     double div_F = (faces[i+1].physics - faces[i].physics) / (cells[i+1].position - cells[i].position);
    //     double grad = 2 / (cells[i+1].position - cells[i].position);
    //     cells[i].physics = div_F * cells[i].volume + cells[i].physics;
    // }

    return 0;
}


