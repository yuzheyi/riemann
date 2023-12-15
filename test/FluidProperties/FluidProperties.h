#ifndef FLUIDPROPERTIES_H
#define FLUIDPROPERTIES_H

#include <vector>

class FluidProperties {
public:
    FluidProperties(double xPos, double yPos, double zPos, double pressure, double density, double xVelocity, double yVelocity, double zVelocity, double temperature);

    std::vector<double> PointPosition() const;
    std::vector<double> velocity() const;
    std::vector<double> getVelocityVector(double velocityComponent) const;
    std::vector<double> tVector() const;
    std::vector<double> uVector() const;
    std::vector<double> vVector() const;
    std::vector<double> wVector() const;

private:
    double xPosition;
    double yPosition;
    double zPosition;
    double pressure;
    double density;
    double xVelocity;
    double yVelocity;
    double zVelocity;
    double temperature;
};

#endif // FLUIDPROPERTIES_H


