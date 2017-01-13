#pragma once

#include <iostream>
#include <cmath>

struct Vector
{
    float x;
    float y;
    Vector() : x(0), y(0) {};
    Vector(const float& x, const float&y) : x(x), y(y) {};

    Vector& operator+=(const Vector& rhs)
    {
        x += rhs.x;
        y += rhs.y;
        return *this;
    }

    Vector& operator-=(const Vector rhs)
    {
        x -= rhs.x;
        y -= rhs.y;
        return *this;
    }

    Vector& operator*=(float rhs)
    {
        x *= rhs;
        y *= rhs;
        return *this;
    }

    Vector& operator/=(float rhs)
    {
        x /= rhs;
        y /= rhs;
        return *this;
    }
};

typedef Vector Point;

static
Vector operator+(const Vector& lhs, const Vector rhs)
{
    return Vector(lhs.x + rhs.x, lhs.y + rhs.y);
}

static
Vector operator-(const Vector& lhs, const Vector rhs)
{
    return Vector(lhs.x - rhs.x, lhs.y - rhs.y);
}

static
Vector operator*(const Vector& lhs, float rhs)
{
    return Vector(rhs*lhs.x, rhs*lhs.y);
}

static
Vector operator/(const Vector& lhs, float rhs)
{
    return Vector(lhs.x/rhs, lhs.y/rhs);
}

static
float dot(const Vector& lhs, const Vector& rhs)
{
    return lhs.x * rhs.x + lhs.y * rhs.y;
}

static
float normSquared(const Vector& lhs)
{
    return lhs.x * lhs.x + lhs.y * lhs.y;
}

static
float norm(const Vector& lhs)
{
    return std::sqrt(normSquared(lhs));
}

struct Line
{
    Point pt0;
    Point pt1;
};

static
std::ostream& operator<<(std::ostream& os, const Point& pt)
{
    os << pt.x << ", " << pt.y;
    return os;
}
