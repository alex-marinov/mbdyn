/* 
 * File:   MDShape.h
 * Author: tingnan
 *
 * Created on February 8, 2011, 10:12 PM
 * 
 * The MDShape class is defined here. basically a surface is composed of
 * 
 */

#include "md_vec.h"

#ifndef _MDSHAPE_H
#define	_MDSHAPE_H

class MDShape
{
  public:
    virtual bool GetSurfGeo (const RealVec &, double, RealVec &, RealVec &, double &, bool &) = 0;
    virtual void SetNodeMass (double) = 0;
    virtual void GetNodeMass (double &) = 0;
    virtual void GetNodeForce (RealVec &) = 0;
    virtual void GetNodePosition (RealVec &) = 0;
    virtual void GetNodeVelocity (RealVec &) = 0;
    virtual void GetNodeAngularVel (RealVec &) = 0;
    virtual void GetNodeOrientation (rmat &) = 0;
    virtual void GetNodeMoment (RealVec &) = 0;
    virtual void SetNodePosition (const RealVec &) = 0;
    virtual void SetNodeVelocity (const RealVec &) = 0;
    virtual void SetNodeAngularVel (const RealVec &) = 0;
    virtual void SetNodeForce (const RealVec &) = 0;
    virtual void SetNodeOrientation (const rmat &) = 0;
    virtual void SetNodeMoment (const RealVec &) = 0;
    virtual void NodesInit (const RealVec &, const RealVec &, const rmat &) = 0;
};

class MDShapeBase:public MDShape
{
  protected:
    double mass_;
    rmat orientation_mat_;
    RealVec velocity_;
    RealVec force_;
    RealVec position_;
    RealVec angular_velocity_;
    RealVec moment_;
    RealVec orientation_euler_;
  public:
            MDShapeBase ();
           ~MDShapeBase ();
    bool GetSurfGeo (const RealVec &, double, RealVec &, RealVec &, double &, bool &);
    // all following method should be common among objects, no need to redefine
    void GetNodeMass (double &);
    void GetNodeForce (RealVec &);
    void GetNodePosition (RealVec &);
    void GetNodeVelocity (RealVec &);
    void GetNodeAngularVel (RealVec &);
    void GetNodeOrientation (rmat &);
    void GetNodeMoment (RealVec &);
    void SetNodeMass (double);
    void SetNodePosition (const RealVec &);
    void SetNodeVelocity (const RealVec &);
    void SetNodeAngularVel (const RealVec &);
    void SetNodeForce (const RealVec &);
    void SetNodeOrientation (const rmat &);
    void SetNodeMoment (const RealVec &);
    void NodesInit (const RealVec &, const RealVec &, const rmat &);
};

class MDSphere:public MDShapeBase
{
  protected:
    double radius_;
  public:
    void NodesInitSph (const RealVec &, const RealVec &);
    void SphereInit (double, double);
    bool GetSurfGeo (const RealVec &, double, RealVec &, RealVec &, double &, bool &);
};

class MDCylinder:public MDShapeBase
{
  protected:
    double radius_;
    double length_;
    RealVec orientation_vec_;
  public:
            bool GetSurfGeo (const RealVec &, double, RealVec &, RealVec &, double &, bool &);
    void NodesInitCyl (const RealVec &, const RealVec &, const RealVec &);
    void SetNodeOrientation (const rmat &);
    void CylinderInit (double, double, double);
};

class MDRoundCylinder:public MDShapeBase
{
  protected:
    double radius_;
    double length_;
    RealVec orientation_vec_;
  public:
            bool GetSurfGeo (const RealVec &, double, RealVec &, RealVec &, double &, bool &);
    void NodesInitCyl (const RealVec &, const RealVec &, const RealVec &);
    void SetNodeOrientation (const rmat &);
    void CylinderInit (double, double, double);
};

class MDSquare:public MDShapeBase
{
  protected:
    RealVec direction_x_;
    RealVec direction_y_;
    RealVec direction_z_;
    double length_x_;
    double length_y_;
    double length_z_;
  public:
           MDSquare (void);
    bool GetSurfGeo (const RealVec &, double, RealVec &, RealVec &, double &, bool &);
    void SetNodeOrientation (const rmat &);
    void SquareInit (double, double, double, double);
};

class MDArch:public MDShapeBase
{
  protected:
    RealVec direction_x_;
    RealVec direction_y_;
    RealVec direction_z_;
    double radius_;
    double theta_;
    double thickness_;
    double width_;
  public:
           MDArch (void);
    bool GetSurfGeo (const RealVec &, double, RealVec &, RealVec &, double &, bool &);
    void SetNodeOrientation (const rmat &);
    void ArchInit (double, double, double, double, double);
};


#endif /* MDSHAPE_H */
