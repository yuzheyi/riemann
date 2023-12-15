
#ifndef FVM_H
#define FVM_H
#include <vector>
struct Face {
    double area;           // 表示面积的属性
    double position;       // 表示位置的属性
    double velocity; // 表示速度的属性
    double P; // 表示压力的属性
    double rho; // 表示密度的属性
    double F[3];
};

struct Point {
    double position;  // 表示位置的属性
};


struct Cell {
    double position;  // 表示位置的属性
    double volume;  // 表示体积的属性
    double P; // 表示压力的属性
    double rho; // 表示密度的属性
    double velocity; // 表示速度的属性
    double T;  // 表示温度的属性
    double U[3];
    Face* faces[2];//准备定义face0为入口面，face1为出口面
};


void pointToFace(const std::vector<Point>& points, std::vector<Face>& faces);
void pointToCell(const std::vector<Point>& points, const std::vector<Face>& faces, std::vector<Cell>& cells);



// void pointToFace(Point* points, int number,Face* faces);
// void pointToCell(Point* points,int number,Face* faces,Cell* cells);
//定义梯度和散度


#endif
