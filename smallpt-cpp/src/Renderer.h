#pragma once

#include "Scene.h"

struct IRenderer
{
    virtual Vec render(const Scene& scene, const Ray& ray, int depth, RandomLCG& Xi) = 0;
};

struct SimpleRenderer : public IRenderer
{
    Vec render(const Scene& scene, const Ray& ray, int depth, RandomLCG& Xi)
    {
        double t;                               // distance to intersection
        int id = 0;                               // id of intersected object
        if (!scene.intersect(ray, t, id)) 
            return Vec(); // if miss, return black
        const Geometry& hitObj = *scene.getGeometry(id); 
        Vec hitPt = ray.orig + ray.dir*t;
        Vec n=(hitPt - hitObj.pos).norm();
        Vec nl = n.dot(ray.dir)<0 ? n:n*-1;
        Vec f = hitObj.color;
        double p = f.x>f.y && f.x>f.z ? f.x : f.y>f.z ? f.y : f.z; // p = max(f.x, f.y, f.z)

        if (++depth > 5 || !p) 
        {
            if (Xi() < p) 
                f = f * (1/p); 
            else 
                return hitObj.emission; // Russian roulette
        }
        if (depth > 100) 
            return hitObj.emission; // MILO

        if (hitObj.refl == DIFFUSE)
        {
            double r1 = 2*M_PI*Xi();
            double r2 = Xi();
            double r2s = sqrt(r2);
            Vec w = nl;
            Vec u=(fabs(w.x)>.1 ? Vec(0,1) : Vec(1)).cross(w);
            u.norm();
            Vec v = w.cross(u);
            Vec d = u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1 - r2);
            d.norm();

            return hitObj.emission + f.mult(render(scene, Ray(hitPt,d),depth,Xi));
        }
        else
        {
            Ray reflRay(hitPt, ray.dir - n*2*n.dot(ray.dir)); 
            return hitObj.emission + f.mult(render(scene, reflRay,depth,Xi));
        }
    }
};
