#include "Fvm.h"
#include <cmath>
#include<iostream>


// void pointToFace(const std::vector<Point>& points, std::vector<Face>& faces);
void pointToCell(const std::vector<Point>& points, const std::vector<Face>& faces, std::vector<Cell>& cells)

// void pointToCell(Point* points,int number,Face* faces, Cell* cells)
{
    std::cout << "构造体开始"<<std::endl;
    // int pointLength = sizeof(points) / sizeof(points[0]);
    int count = 0;
    while(count< cells.size())
    {
        
        cells[count].position=points[count].position/2+points[count+1].position/2;
        cells[count].volume=fabs(points[count].position-points[count+1].position);
        cells[count].faces[0] = const_cast<Face*>(&faces[count]);
        cells[count].faces[1] = const_cast<Face*>(&faces[count + 1]);
        std::cout << cells[count].position<<",";
        count++;

    }
    std::cout << "构造体完成"<<std::endl;
}


void pointToFace(const std::vector<Point>& points, std::vector<Face>& faces)
// void pointToFace(Point* points, int number,Face* faces)
{
    std::cout << "构造面开始"<<std::endl;
    int count = 0;
    while(count<faces.size())
    {
        faces[count].area=1;
        faces[count].position=points[count].position;
        std::cout << faces[count].position<<",";
        count++;
    }
    std::cout << "构造面完成"<<std::endl;

}



 