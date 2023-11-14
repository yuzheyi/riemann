#include <vector>
#ifndef FVM_H
#define FVM_H

struct Face {
    double area;  // 表示面积的属性
    double position;  // 表示位置的属性d
    double physics;  // 表示物理量的属性
};

struct Point {
    double position;  // 表示位置的属性
};


struct Cell {
    double position;  // 表示位置的属性
    double volume;  // 表示体积的属性
    double physics; // 表示物理量的属性
    Face* faces[2];//准备定义face0为入口面，face1为出口面
};


void pointToFace(const std::vector<Point>& points, std::vector<Face>& faces);
void pointToCell(const std::vector<Point>& points, const std::vector<Face>& faces, std::vector<Cell>& cells);



// void pointToFace(Point* points, int number,Face* faces);
// void pointToCell(Point* points,int number,Face* faces,Cell* cells);
//定义梯度和散度


#endif
