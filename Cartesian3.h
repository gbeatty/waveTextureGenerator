#pragma once
#include <complex>

class Cartesian3
{
public:
    Cartesian3() :
        x(0.0),
        y(0.0),
        z(0.0)
    {

    }

    Cartesian3(double x, double y, double z) :
        x(x),
        y(y),
        z(z)
    {
        
    }

    Cartesian3(const Cartesian3& value) :
        x(value.x),
        y(value.y),
        z(value.z)
    {
        
    }

    double Magnitude() const
    {
        return std::sqrt(std::pow(x, 2.0) + std::pow(y, 2.0) + std::pow(z, 2.0));
    }

    Cartesian3 Normalize() const
    {
        double mag = Magnitude();
        if(mag == 0.0)
        {
            return Cartesian3(0.0, 0.0, 0.0);
        }

        return Cartesian3(x / mag, y / mag, z / mag);
    }

    double Dot(const Cartesian3& rhs) const
    {
        return (x*rhs.x) + (y*rhs.y) + (z*rhs.z);
    }

    Cartesian3 operator+(const double& value) const
    {
        return Cartesian3(x + value, y + value, z + value);
    }

    Cartesian3 operator*(const double& value) const
    {
        return Cartesian3(x * value, y * value, z * value);
    }

    Cartesian3 operator-() const
    {
        return Cartesian3(-x, -y, -z);
    }

    void operator=(const Cartesian3& value)
    {
        x = value.x;
        y = value.y;
        z = value.z;
    }

    double x;
    double y;
    double z;
};
