#pragma once

// a scene composed of geometries

#include "Geometry.h"
#include "irrXML.h"
#include <vector>
#include <string>

struct IScene
{
    virtual ~IScene(){}

    static void printBuiltInSceneNames();
    static bool createBuiltInScene(IScene* pScene, const std::string& sceneName);

    virtual Geometry* getGeometry(int id) const = 0;  

    void addGeometries(Geometry* geoms[], size_t count)
    {
        for (size_t i = 0; i < count; i++)
        {
            addGeometry(geoms[i]);
        }
    }

    virtual bool intersect(const Ray& ray, double& hitDist, int& hitId) const = 0;

protected:
    virtual void addGeometry(Geometry* geom) = 0;
};

struct SimpleScene : public IScene
{
    Geometry* getGeometry(int id) const
    {
        return mGeomtries[id];
    }

    void addGeometry(Geometry* geom)
    {
        mGeomtries.push_back(geom);
    }

    virtual ~SimpleScene()
    {
        size_t nGeom = mGeomtries.size();
        for (size_t i=0; i < nGeom; i++)
            delete mGeomtries[i];
    }

    // brute-force
    bool intersect(const Ray& ray, double& hitDist, int& hitId) const
    {
        const double kInfinity = hitDist = 1e20;
        size_t nGeom = mGeomtries.size();
        for (size_t i = 0; i < nGeom; i++) 
        {
            double dist = mGeomtries[i]->intersect(ray);
            if (dist && dist < hitDist)
            {
                hitDist = dist;
                hitId = i;
            }
        }
        return hitDist < kInfinity;
    }

protected:
    std::vector<Geometry*> mGeomtries;
};
