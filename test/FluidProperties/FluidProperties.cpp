// FluidProperties.cpp
#include "FluidProperties.h"
#include <vector>
#include <algorithm>

//第一个FluidProperties是表示类
//在使用过程中先定义类FluidProperties，后取个体名字fluid
//使用方法FluidProperties fluid
//定义了流体所有的信息
FluidProperties::FluidProperties(double xPos, double yPos, double zPos, double pressure, double density, double xVelocity,double yVelocity,double zVelocity, double temperature) {
    xPosition = xPos;
    yPosition = yPos;
    zPosition = zPos;
    this->pressure = pressure;
    this->density = density;
    this->xVelocity = xVelocity;
    this->yVelocity = yVelocity;
    this->zVelocity = zVelocity;
    this->temperature = temperature;
}

std::vector<double> FluidProperties::PointPosition() const {
    return {xPosition, yPosition, zPosition};
}

//速度向量
std::vector<double> FluidProperties::velocity() const {
    return {xVelocity, yVelocity, zVelocity};
}

//守恒方程偏t项
std::vector<double> FluidProperties::tVector() const {
 return {density,density*xVelocity, density*yVelocity, density*zVelocity, density*temperature};
 }
//三个方向
// std::vector<double> FluidProperties::getVelocityVector(double velocityComponent) const {
//     double u = velocityComponent;
//     return u * std::vector<double>{density, density * u, density * yVelocity, density * zVelocity, density * temperature + pressure};//差一个c_v
// }
std::vector<double> FluidProperties::getVelocityVector(double velocityComponent) const {
    double u = velocityComponent;
    std::vector<double> values = {density, density * xVelocity*u, density * yVelocity*u, density * zVelocity*u, density * temperature * 715.8*u + pressure*u};

    std::vector<double> result;
    result.reserve(values.size());

    std::transform(values.begin(), values.end(), std::back_inserter(result), [u](double value) { return u * value; });

    return result;
}

std::vector<double> FluidProperties::uVector() const {
    return getVelocityVector(xVelocity);
}

std::vector<double> FluidProperties::vVector() const {
    return getVelocityVector(yVelocity);
}

std::vector<double> FluidProperties::wVector() const {
    return getVelocityVector(zVelocity);
}
