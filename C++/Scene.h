#ifndef __SMALLPT_SCENE_H__
#define __SMALLPT_SCENE_H__

// a scene composed of geometries

#include "Geometry.h"
#include "irrXML.h"
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

#endif __SMALLPT_SCENE_H__