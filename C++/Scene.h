#ifndef __SMALLPT_SCENE_H__
#define __SMALLPT_SCENE_H__

// a scene composed of geometries

#include "Geometry.h"
#include <vector>

struct Scene
{
    std::vector<Geometry*> mGeomtries;

    virtual ~Scene()
    {
        size_t nGeom = mGeomtries.size();
        for (size_t i=0;i<nGeom;i++)
            delete mGeomtries[i];
    }

    inline bool intersect(const Ray &r, double &t, int &id) const
    {
        const double INF = t = 1e20;
        size_t nGeom = mGeomtries.size();
        for(size_t i = 0;i < nGeom;i++) 
        {
            double d = mGeomtries[i]->intersect(r);
            if (d && d < t)
            {
                t = d;
                id = i;
            }
        }
        return t < INF;
    }
};

Sphere spheres[] = 
{
    //   radius     position                emission        color               material
    Sphere(1e5,     Vec( 1e5+1,40.8,81.6),  Vec(),          Vec(.75,.25,.25),   DIFFUSE),// left wall
    Sphere(1e5,     Vec(-1e5+99,40.8,81.6), Vec(),          Vec(.25,.25,.75),   DIFFUSE),// right wall
    Sphere(1e5,     Vec(50,40.8, 1e5),      Vec(),          Vec(.75,.75,.75),   DIFFUSE),// back wall
    Sphere(1e5,     Vec(50,40.8,-1e5+170),  Vec(),          Vec(),              DIFFUSE),// front wall
    Sphere(1e5,     Vec(50, 1e5, 81.6),     Vec(),          Vec(.75,.75,.75),   DIFFUSE),// bottom wall
    Sphere(1e5,     Vec(50,-1e5+81.6,81.6), Vec(),          Vec(.75,.75,.75),   DIFFUSE),// top wall

    Sphere(16.5,    Vec(27,16.5,47),        Vec(),          Vec(1,1,1)*.999,    SPECULAR),// mirror 1
    Sphere(20.5,    Vec(73,16.5,78),        Vec(),          Vec(1,1,1)*.999,    SPECULAR),// mirror 2
    Sphere(10.5,    Vec(27,16.5,100),        Vec(),          Vec(1,1,1)*.999,    SPECULAR),// mirror 3

    Sphere(600,     Vec(50,681.6-.27,81.6), Vec(12,12,12),  Vec(),              DIFFUSE) // the white light
};

const int N_SPHERES = sizeof(spheres)/sizeof(Sphere);

#endif __SMALLPT_SCENE_H__