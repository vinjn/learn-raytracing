#ifndef __SMALLPT_SPHERE_H__
#define __SMALLPT_SPHERE_H__

#include "Utils.h"

// for now the only geometry is sphere

enum Refl_t 
{ 
    DIFFUSE,    // ¬˛∑¥…‰
    SPECULAR,   // ∏ﬂπ‚
    REFRACT,    // ’€…‰
};

struct Geometry
{
    double radius; 
    Vec pos, emission, color;
    Refl_t refl;
    Geometry(double rad_, Vec p_, Vec e_, Vec c_, Refl_t refl_):
    radius(rad_), pos(p_), emission(e_), color(c_), refl(refl_) {}
    virtual double intersect(const Ray &r) const = 0;
};

struct Sphere : public Geometry
{
    Sphere(double rad_, Vec p_, Vec e_, Vec c_, Refl_t refl_):Geometry(rad_, p_, e_, c_, refl_)
    {}

    virtual double intersect(const Ray &r) const 
    {   
        // returns distance, 0 if no hit
        Vec op = pos-r.orig; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
        double eps=1e-4, b=op.dot(r.dir), det=b*b-op.dot(op)+radius*radius;
        if (det<0) 
            return 0; 
        else 
            det=sqrt(det);
        double t = 0;
        if ((t = b - det) > eps)
            return t;
        else if ((t = b + det) > eps)
            return t;
        return 0;
    }
};


#endif __SMALLPT_SPHERE_H__