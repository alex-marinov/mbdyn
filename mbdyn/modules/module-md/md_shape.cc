/* 
 * File:   MDShape.cc
 * Author: tingnan
 * 
 * Created on February 8, 2011, 10:12 PM
 */

#include "md_shape.h"
#include "md_utility.h"
#include <cassert>
using namespace std;


MDShapeBase::MDShapeBase ()
{
    mass_ = 0.0;
    orientation_mat_.resize (3);
    int i;

    for (i = 0; i < 3; i++)
    {
        orientation_mat_[i].resize (3);
    }
}

MDShapeBase::~MDShapeBase ()
{
    orientation_mat_.clear ();
}

void
MDShapeBase::GetNodeAngularVel (RealVec & omega)
{
    omega = angular_velocity_;
}

void
MDShapeBase::GetNodeForce (RealVec & f)
{
    f = force_;
}

void
MDShapeBase::GetNodeMass (double &m)
{
    m = mass_;
}

void
MDShapeBase::GetNodeMoment (RealVec & moment)
{
    moment = moment_;
}

void
MDShapeBase::GetNodeOrientation (rmat & ori)
{
    ori = orientation_mat_;
}

void
MDShapeBase::GetNodePosition (RealVec & pos)
{
    pos = position_;
}

void
MDShapeBase::GetNodeVelocity (RealVec & vel)
{
    vel = velocity_;
}

bool MDShapeBase::GetSurfGeo (const RealVec &, double, RealVec &, RealVec &, double &, bool &)
{
    // no op
    return false;
}

void
MDShapeBase::NodesInit (const RealVec & pos, const RealVec & vel, const rmat & ori)
{
    position_ = pos;
    velocity_ = vel;
    orientation_mat_ = ori;
}

void
MDShapeBase::SetNodeAngularVel (const RealVec & omega)
{
    angular_velocity_ = omega;
}

void
MDShapeBase::SetNodeForce (const RealVec & f)
{
    force_ = f;
}

void
MDShapeBase::SetNodeMass (double m)
{
    mass_ = m;
}

void
MDShapeBase::SetNodeMoment (const RealVec & moment)
{
    moment_ = moment;
}

void
MDShapeBase::SetNodeOrientation (const rmat & ori)
{
    orientation_mat_ = ori;
}

void
MDShapeBase::SetNodePosition (const RealVec & pos)
{
    position_ = pos;
}

void
MDShapeBase::SetNodeVelocity (const RealVec & vel)
{
    velocity_ = vel;
}

void
MDSphere::SphereInit (double r, double m)
{
    radius_ = r;
    mass_ = m;
}

bool MDSphere::GetSurfGeo (const RealVec & sand_pos, double sand_radius, RealVec & surf_pos,
                           RealVec & surf_vel, double &virtual_radius, bool & if_inside)
{
    if ((sand_pos - position_) * (sand_pos - position_) < radius_ * radius_)
    {
        if_inside = 1;
        return false;
    }
    surf_pos = position_;
    surf_vel = velocity_;
    virtual_radius = radius_;
    return true;
}

bool MDCylinder::GetSurfGeo (const RealVec & sand_pos, double sand_radius, RealVec & surf_pos,
                             RealVec & surf_vel, double &virtual_radius, bool & if_inside)
{
    RealVec
        dr = sand_pos - position_;
    RealVec
        surf_drv,
        perp_drv;
    double
        proj_dr,
        perp_dr;

    proj_dr = dr * orientation_vec_;
    perp_drv = dr - orientation_vec_ * proj_dr;
    perp_dr = sqrt (perp_drv * perp_drv);
    if (abs (proj_dr) > 0.5 * length_)
    {
        // could interact with top or bottom
        if (abs (proj_dr) > sand_radius + 0.5 * length_)
        {
            return false;
        }
        else
        {
            if (perp_dr < radius_)
            {
                surf_drv = 0.5 * sign (proj_dr) * length_ * orientation_vec_ + perp_drv;
                surf_pos = position_ + surf_drv;
                surf_vel = velocity_ + CrossProduct (angular_velocity_, surf_drv);
                virtual_radius = 0;
                return true;
            }
            // now it is possible to collide with edge
            else
            {
                if (perp_dr < radius_ + sand_radius)
                {
                    surf_drv =
                        0.5 * sign (proj_dr) * length_ * orientation_vec_ +
                        perp_drv * (radius_ / perp_dr);
                    surf_pos = position_ + surf_drv;
                    surf_vel = velocity_ + CrossProduct (angular_velocity_, surf_drv);
                    virtual_radius = 0;
                    return true;
                }
                else
                    return false;
            }
        }
    }
    else
    {
        if (perp_dr > radius_ + sand_radius)
            return false;
        else
        {
            if (perp_dr < radius_)
            {
                if_inside = 1;
                return false;
            }
            surf_drv = perp_drv * radius_ / perp_dr + orientation_vec_ * proj_dr;
            surf_pos = position_ + surf_drv;
            surf_vel = velocity_ + CrossProduct (angular_velocity_, surf_drv);
            virtual_radius = 0;
            return true;
        }
    }
}


void
MDCylinder::CylinderInit (double r, double l, double m)
{
    radius_ = r;
    length_ = l;
    mass_ = m;
}

void
MDCylinder::SetNodeOrientation (const rmat & orientationMatrix)
{
    RealVec initVec (1, 0, 0);

    orientation_mat_ = orientationMatrix;
    orientation_vec_.x =
        orientation_mat_[0][0] * initVec.x + orientation_mat_[0][1] * initVec.y +
        orientation_mat_[0][2] * initVec.z;
    orientation_vec_.y =
        orientation_mat_[1][0] * initVec.x + orientation_mat_[1][1] * initVec.y +
        orientation_mat_[1][2] * initVec.z;
    orientation_vec_.z =
        orientation_mat_[2][0] * initVec.x + orientation_mat_[2][1] * initVec.y +
        orientation_mat_[2][2] * initVec.z;
}

void
MDCylinder::NodesInitCyl (const RealVec & pos, const RealVec & vel, const RealVec & ori_vec)
{
    position_ = pos;
    velocity_ = vel;
    orientation_vec_ = ori_vec;

}

bool MDRoundCylinder::GetSurfGeo (const RealVec & sand_pos, double sand_radius, RealVec & surf_pos,
                                  RealVec & surf_vel, double &virtual_radius, bool & if_inside)
{
    RealVec
        dr,
        tmpdr = sand_pos;

    dr = tmpdr - position_;
    RealVec
        surf_drv,
        perp_drv;
    double
        proj_dr,
        perp_dr;

    proj_dr = dr * orientation_vec_;
    perp_drv = dr - orientation_vec_ * proj_dr;
    perp_dr = sqrt (perp_drv * perp_drv);
    if (abs (proj_dr) > 0.5 * length_)
    {
        // could interact with top or bottom
        // now place a half ball at the top or bottom center
        surf_drv = 0.5 * sign (proj_dr) * length_ * orientation_vec_;
        surf_pos = position_ + surf_drv;
        tmpdr = sand_pos;
        if ((tmpdr - surf_pos) * (tmpdr - surf_pos) < radius_ * radius_)
        {
            if_inside = 1;
            return true;
        }
        surf_vel = velocity_ + CrossProduct (angular_velocity_, surf_drv);
        virtual_radius = radius_;
        return true;
    }
    else
    {
        if (perp_dr > radius_ + sand_radius)
            return false;
        else
        {
            if (perp_dr < radius_)
            {
                if_inside = 1;
                return false;
            }
            surf_drv = perp_drv * radius_ / perp_dr + orientation_vec_ * proj_dr;
            surf_pos = position_ + surf_drv;
            surf_vel = velocity_ + CrossProduct (angular_velocity_, surf_drv);
            virtual_radius = 0;
            return true;
        }
    }
}


void
MDRoundCylinder::CylinderInit (double r, double l, double m)
{
    radius_ = r;
    length_ = l;
    mass_ = m;
}

void
MDRoundCylinder::SetNodeOrientation (const rmat & orientationMatrix)
{
    RealVec initVec (1, 0, 0);

    orientation_mat_ = orientationMatrix;
    orientation_vec_.x =
        orientation_mat_[0][0] * initVec.x + orientation_mat_[0][1] * initVec.y +
        orientation_mat_[0][2] * initVec.z;
    orientation_vec_.y =
        orientation_mat_[1][0] * initVec.x + orientation_mat_[1][1] * initVec.y +
        orientation_mat_[1][2] * initVec.z;
    orientation_vec_.z =
        orientation_mat_[2][0] * initVec.x + orientation_mat_[2][1] * initVec.y +
        orientation_mat_[2][2] * initVec.z;
}

void
MDRoundCylinder::NodesInitCyl (const RealVec & pos, const RealVec & vel, const RealVec & ori_vec)
{
    position_ = pos;
    velocity_ = vel;
    orientation_vec_ = ori_vec;

}

MDSquare::MDSquare (void)
{
    RealVec vecx (1, 0, 0);
    RealVec vecy (0, 1, 0);

    direction_x_ = vecx;
    direction_y_ = vecy;
    direction_z_ = CrossProduct (direction_x_, direction_y_);
}

bool MDSquare::GetSurfGeo (const RealVec & sand_pos, double sand_radius, RealVec & surf_pos,
                           RealVec & surf_vel, double &virtual_radius, bool & if_inside)
{
    RealVec
        dr = sand_pos - position_;
    double
        drx = dr * direction_x_;
    double
        dry = dr * direction_y_;
    double
        drz = dr * direction_z_;
    RealVec
        surf_drv;

    if (abs (drx) < 0.5 * length_x_ && abs (dry) < 0.5 * length_y_ && abs (drz) < 0.5 * length_z_)
    {
        if_inside = 1;
        return false;
    }
    if (abs (drx) > 0.5 * length_x_ + sand_radius || abs (dry) > 0.5 * length_y_ + sand_radius
        || abs (drz) > 0.5 * length_z_ + sand_radius)
    {
        return false;
    }
    // collide with a vertex
    if (abs (drx) > 0.5 * length_x_ && abs (dry) > 0.5 * length_y_ && abs (drz) > 0.5 * length_z_)
    {
        surf_drv =
            0.5 * sign (drx) * length_x_ * direction_x_ +
            0.5 * sign (dry) * length_y_ * direction_y_ +
            0.5 * sign (drz) * length_z_ * direction_z_;
        surf_pos = position_ + surf_drv;
        surf_vel = velocity_ + CrossProduct (angular_velocity_, surf_drv);
        virtual_radius = 0;
        return true;
    }
    // collide with edge
    if (abs (drx) > 0.5 * length_x_ && abs (dry) > 0.5 * length_y_ && abs (drz) < 0.5 * length_z_)
    {
        // collide with z edge
        surf_drv =
            0.5 * sign (drx) * length_x_ * direction_x_ +
            0.5 * sign (dry) * length_y_ * direction_y_ + drz * direction_z_;
        surf_pos = position_ + surf_drv;
        surf_vel = velocity_ + CrossProduct (angular_velocity_, surf_drv);
        virtual_radius = 0;
        return true;
    }
    if (abs (drx) > 0.5 * length_x_ && abs (dry) < 0.5 * length_y_ && abs (drz) > 0.5 * length_z_)
    {
        // collide with y edge
        surf_drv =
            0.5 * sign (drx) * length_x_ * direction_x_ + dry * direction_y_ +
            0.5 * sign (drz) * length_z_ * direction_z_;
        surf_pos = position_ + surf_drv;
        surf_vel = velocity_ + CrossProduct (angular_velocity_, surf_drv);
        virtual_radius = 0;
        return true;
    }
    if (abs (drx) < 0.5 * length_x_ && abs (dry) > 0.5 * length_y_ && abs (drz) > 0.5 * length_z_)
    {
        // collide with x edge
        surf_drv =
            drx * direction_x_ + 0.5 * sign (dry) * length_y_ * direction_y_ +
            0.5 * sign (drz) * length_z_ * direction_z_;
        surf_pos = position_ + surf_drv;
        surf_vel = velocity_ + CrossProduct (angular_velocity_, surf_drv);
        virtual_radius = 0;
        return true;
    }
    // collide with a surface
    if (abs (drx) > 0.5 * length_x_ && abs (dry) < 0.5 * length_y_ && abs (drz) < 0.5 * length_z_)
    {
        // collide with yz surface
        surf_drv =
            0.5 * sign (drx) * length_x_ * direction_x_ + dry * direction_y_ + drz * direction_z_;
        surf_pos = position_ + surf_drv;
        surf_vel = velocity_ + CrossProduct (angular_velocity_, surf_drv);
        virtual_radius = 0;
        return true;
    }
    if (abs (drx) < 0.5 * length_x_ && abs (dry) > 0.5 * length_y_ && abs (drz) < 0.5 * length_z_)
    {
        // collide with xz surface
        surf_drv =
            drx * direction_x_ + 0.5 * sign (dry) * length_y_ * direction_y_ + drz * direction_z_;
        surf_pos = position_ + surf_drv;
        surf_vel = velocity_ + CrossProduct (angular_velocity_, surf_drv);
        virtual_radius = 0;
        return true;
    }
    if (abs (drx) < 0.5 * length_x_ && abs (dry) < 0.5 * length_y_ && abs (drz) > 0.5 * length_z_)
    {
        // collide with xy surface
        surf_drv =
            drx * direction_x_ + dry * direction_y_ + 0.5 * sign (drz) * length_z_ * direction_z_;
        surf_pos = position_ + surf_drv;
        surf_vel = velocity_ + CrossProduct (angular_velocity_, surf_drv);
        virtual_radius = 0;
        return true;
    }
    assert (0);
}

void
MDSquare::SquareInit (double lx, double ly, double lz, double m)
{
    length_x_ = lx;
    length_y_ = ly;
    length_z_ = lz;
    mass_ = m;
}

void
MDSquare::SetNodeOrientation (const rmat & orimat)
{
    RealVec vecx (1, 0, 0);
    RealVec vecy (0, 1, 0);

    orientation_mat_ = orimat;
    direction_x_.x =
        orientation_mat_[0][0] * vecx.x + orientation_mat_[0][1] * vecx.y +
        orientation_mat_[0][2] * vecx.z;
    direction_x_.y =
        orientation_mat_[1][0] * vecx.x + orientation_mat_[1][1] * vecx.y +
        orientation_mat_[1][2] * vecx.z;
    direction_x_.z =
        orientation_mat_[2][0] * vecx.x + orientation_mat_[2][1] * vecx.y +
        orientation_mat_[2][2] * vecx.z;
    direction_y_.x =
        orientation_mat_[0][0] * vecy.x + orientation_mat_[0][1] * vecy.y +
        orientation_mat_[0][2] * vecy.z;
    direction_y_.y =
        orientation_mat_[1][0] * vecy.x + orientation_mat_[1][1] * vecy.y +
        orientation_mat_[1][2] * vecy.z;
    direction_y_.z =
        orientation_mat_[2][0] * vecy.x + orientation_mat_[2][1] * vecy.y +
        orientation_mat_[2][2] * vecy.z;
    direction_z_ = CrossProduct (direction_x_, direction_y_);
}


MDArch::MDArch (void)
{
    RealVec vecx (1, 0, 0);
    RealVec vecy (0, 1, 0);

    direction_x_ = vecx;
    direction_y_ = vecy;
    direction_z_ = CrossProduct (direction_x_, direction_y_);
}

bool MDArch::GetSurfGeo (const RealVec & sand_pos, double sand_radius, RealVec & surf_pos,
                         RealVec & surf_vel, double &virtual_radius, bool & if_inside)
{
    RealVec
        dr = sand_pos - position_;
    double
        drx = dr * direction_x_;
    double
        dry = dr * direction_y_;
    double
        drz = dr * direction_z_;
    double
        theta;
    double
        PI = 3.1415926535897932384626433832795;

    if_inside = false;
    if (abs (drz) > 10000 * abs (dry))
        theta = PI / 2.0;
    else
    {
        if (dry > 0)
            theta = abs (atan (drz / dry));
        if (dry < 0)
            theta = PI - abs (atan (drz / dry));
    }
    RealVec
        surf_drv;

    if (abs (drx) > 0.5 * width_ + sand_radius)
    {
        // non-collision
        return false;
    }
    if ((dry * dry + drz * drz) < (radius_ - sand_radius) * (radius_ - sand_radius)
        || (dry * dry + drz * drz) > (radius_ + sand_radius) * (radius_ + sand_radius))
    {
        return false;
    }
    if (abs (drx) <= 0.5 * width_ + sand_radius && abs (drx) >= 0.5 * width_)
    {
        // collision with top(bottom)-edge and top(bottom)-corner
        if (theta < theta_)
        {
            // edge
            surf_drv =
                0.5 * sign (drx) * width_ * direction_x_ + radius_ * cos (theta) * direction_y_ +
                radius_ * sin (theta) * sign (drz) * direction_z_;
            surf_pos = position_ + surf_drv;
            surf_vel = velocity_ + CrossProduct (angular_velocity_, surf_drv);
            virtual_radius = 0;
            return true;
        }
        else
        {
            // corner
            surf_drv =
                0.5 * sign (drx) * width_ * direction_x_ + radius_ * cos (theta_) * direction_y_ +
                radius_ * sin (theta_) * sign (drz) * direction_z_;
            surf_pos = position_ + surf_drv;
            surf_vel = velocity_ + CrossProduct (angular_velocity_, surf_drv);
            virtual_radius = 0;
            return true;
        }
    }
    if (abs (drx) < 0.5 * width_)
    {
        // collision with side or side edge
        if (theta < theta_)
        {
            surf_drv =
                drx * direction_x_ + radius_ * cos (theta) * direction_y_ +
                radius_ * sin (theta) * sign (drz) * direction_z_;
            surf_pos = position_ + surf_drv;
            surf_vel = velocity_ + CrossProduct (angular_velocity_, surf_drv);
            virtual_radius = 0;
            return true;
        }
        else
        {
            surf_drv =
                drx * direction_x_ + radius_ * cos (theta_) * direction_y_ +
                radius_ * sin (theta_) * sign (drz) * direction_z_;
            surf_pos = position_ + surf_drv;
            surf_vel = velocity_ + CrossProduct (angular_velocity_, surf_drv);
            virtual_radius = 0;
            return true;
        }
    }
    assert (0);
}

void
MDArch::ArchInit (double r, double theta, double thickness, double width, double m)
{
    radius_ = r;
    theta_ = theta;
    thickness_ = thickness;     // currently not implied in the collision volume
    width_ = width;
    mass_ = m;
}

void
MDArch::SetNodeOrientation (const rmat & orimat)
{
    RealVec vecx (1, 0, 0);
    RealVec vecy (0, 1, 0);
    RealVec tmpvec;

    orientation_mat_ = orimat;
    direction_x_.x =
        orientation_mat_[0][0] * vecx.x + orientation_mat_[0][1] * vecx.y +
        orientation_mat_[0][2] * vecx.z;
    direction_x_.y =
        orientation_mat_[1][0] * vecx.x + orientation_mat_[1][1] * vecx.y +
        orientation_mat_[1][2] * vecx.z;
    direction_x_.z =
        orientation_mat_[2][0] * vecx.x + orientation_mat_[2][1] * vecx.y +
        orientation_mat_[2][2] * vecx.z;
    direction_y_.x =
        orientation_mat_[0][0] * vecy.x + orientation_mat_[0][1] * vecy.y +
        orientation_mat_[0][2] * vecy.z;
    direction_y_.y =
        orientation_mat_[1][0] * vecy.x + orientation_mat_[1][1] * vecy.y +
        orientation_mat_[1][2] * vecy.z;
    direction_y_.z =
        orientation_mat_[2][0] * vecy.x + orientation_mat_[2][1] * vecy.y +
        orientation_mat_[2][2] * vecy.z;
    direction_z_ = CrossProduct (direction_x_, direction_y_);
    tmpvec = direction_x_;
    direction_x_ = direction_y_;
    direction_y_ = -1.0 * direction_z_;
    direction_z_ = tmpvec;
    // cout << direction_x_.x << "\t" << direction_x_.y << "\t" << direction_x_.z << "\t" << endl;
    // cout << direction_y_.x << "\t" << direction_y_.y << "\t" << direction_y_.z << "\t" << endl;
    // cout << direction_z_.x << "\t" << direction_z_.y << "\t" << direction_z_.z << "\t" << endl;
}
