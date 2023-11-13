#include "Fvm.h"
#include <cmath>
#include<iostream>

void pointToCell(Point* points,int number,Face* faces, Cell* cells)
{
    
    // int pointLength = sizeof(points) / sizeof(points[0]);
    int pointLength = number;
    int count = 0;
    while(count<pointLength-1)
    {
        cells[count].position=points[count].position/2+points[count+1].position/2;
        cells[count].volume=fabs(points[count].position-points[count+1].position);
        cells[count].faces[0] = &faces[count];
        cells[count].faces[1] = &faces[count+1];
        count++;
    }

}



void pointToFace(Point* points, int number,Face* faces)
{
    int pointLength = number;
    int count = 0;
    while(count<pointLength)
    {
        faces[count].area=1;
        faces[count].position=points[count].position;
        count++;
    }
}



 