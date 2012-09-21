/* 
 * File:   MDVec.h
 * Author: tingnan
 *
 * Created on February 8, 2011, 8:30 PM
 */

#include "md_header.h"
#include <vector>

#ifndef _MDVEC_H
#define	_MDVEC_H

class IntVec
{
public:
    int x;
    int y;
    int z;
    IntVec();

    //    int ToDimOne()
    //    {
    //        return x + y * MD_CELL_X + z * MD_CELL_X * MD_CELL_Y;
    //    }

    //    void ToDim23(int i)
    //    {
    //        z = i / (cell_length_x_ * cell_length_y_);
    //        i = i - z * cell_length_x_ * cell_length_y_;
    //        y = i / cell_length_x_;
    //        x = i - y * cell_length_x_;
    //    }

    //    int end()
    //    {
    //        if (x == cell_length_x_ - 1 && z == cell_length_z_ - 1 && y
    //            == cell_length_y_ - 1)
    //            return 1;
    //        else
    //            return 0;
    //    }

    //    int out()
    //    {
    //        if (x < 0 || y < 0 || z < 0 || x >= MD_CELL_X || y >= MD_CELL_Y || z >= MD_CELL_Z)
    //            return 1;
    //        else
    //            return 0;
    //    }

    //    int operator++()
    //    {
    //
    //        if (x != cell_length_x_ - 1)
    //        {
    //            x++;
    //            return 1;
    //        }
    //        else
    //            x = 0;
    //
    //        if (y != cell_length_y_ - 1)
    //        {
    //            y++;
    //            return 1;
    //        }
    //        else
    //            y = 0;
    //
    //        if (z != cell_length_z_ - 1)
    //        {
    //            z++;
    //            return 1;
    //        }
    //        return 0;
    //    }

    IntVec operator+(const IntVec & vb) const
    {
        IntVec result;
        result.x = x + vb.x;
        result.y = y + vb.y;
        result.z = z + vb.z;
        return result;
    }

    IntVec & operator=(const IntVec & vb)
    {
        x = vb.x;
        y = vb.y;
        z = vb.z;
        return *this;
    }

    IntVec & operator=(int b)
    {
        x = b;
        y = b;
        z = b;
        return *this;
    }

    friend bool equal(const IntVec & a, const IntVec & b)
    {
        if (a.x == b.x && a.z == b.z && a.y == b.y)
            return true;
        return false;
    }
};

class RealVec
{
public:
    double x;
    double z;
    double y;
    RealVec();
    RealVec(double, double);
    RealVec(double, double, double);

    RealVec operator+(const RealVec & vb) const
    {
        RealVec result;
        result.x = x + vb.x;
        result.y = y + vb.y;
        result.z = z + vb.z;
        return result;
    }

    RealVec operator-(const RealVec & vb) const
    {
        RealVec result;
        result.x = x - vb.x;
        result.y = y - vb.y;
        result.z = z - vb.z;
        return result;
    }

    double operator*(const RealVec & vb) const
    {
        double result = 0;
        result += x * vb.x;
        result += y * vb.y;
        result += z * vb.z;
        return result;
    }

    RealVec operator/(double vb) const
    {
        RealVec result;
        result.x = x / vb;
        result.y = y / vb;
        result.z = z / vb;
        return result;
    }

    friend RealVec operator*(double scalar, const RealVec & vb)
    {
        RealVec result;
        result.x = scalar * vb.x;
        result.y = scalar * vb.y;
        result.z = scalar * vb.z;
        return result;
    }

    friend RealVec operator*(const RealVec & vb, double scalar)
    {
        RealVec result;

        result.x = scalar * vb.x;
        result.y = scalar * vb.y;
        result.z = scalar * vb.z;
        return result;
    }

    RealVec & operator=(const RealVec & vb)
    {
        x = vb.x;
        y = vb.y;
        z = vb.z;
        return *this;
    }

    RealVec & operator=(double a)
    {
        x = a;
        y = a;
        z = a;
        return *this;
    }

    friend RealVec CrossProduct(const RealVec &, const RealVec &);

};

typedef std::vector <std::vector < double > > rmat;

#endif
