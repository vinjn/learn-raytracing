#pragma once

#include "Scene.h"

struct IRenderer
{
    // compute the radiance estimate along ray
    virtual Vec radiance(const IScene& scene, const Ray& ray, int depth, RandomLCG& Xi) = 0;
};


// diffuse version
struct DiffuseOnlyRenderer : public IRenderer
{
    Vec radiance(const IScene& scene, const Ray& ray, int depth, RandomLCG& Xi)
    {
        double hitDist;
        int hitId = 0;
        if (!scene.intersect(ray, hitDist, hitId)) 
            return Vec(); // if miss, return black

        const Geometry& hitObj = *scene.getGeometry(hitId); 
        Vec hitColor = hitObj.color;

        // use maximum reflectivity amount of Russian roulette
        double maxRefl = max(hitColor.r, hitColor.g, hitColor.b);
        if (++depth > 5 || !maxRefl) 
        {
            if (Xi() < maxRefl) 
                hitColor = hitColor / maxRefl; 
            else 
                return hitObj.emission;
        }

        Vec hitPt = ray.orig + ray.dir * hitDist;
        Vec hitNorm = (hitPt - hitObj.pos).norm();

        if (hitObj.refl == DIFFUSE)
        {
            // axis: newX/newY/newZ
            Vec newZ = hitNorm.dot(ray.dir) < 0 ? hitNorm : hitNorm * -1; // properly oriented surface normal;
            Vec newX = fabs(newZ.x) > .1 ? Vec::axisY() : Vec::axisX();
            newX = (cross(newX, newZ)).norm();
            Vec newY = cross(newZ, newX);

            // uniformly distributed random direction
            double angle = 2*M_PI*Xi(); // angle around
            double radius = Xi();   // distance from center
            double radius_sqrt = sqrt(radius);
            Vec sampDir = newX*cos(angle)*radius_sqrt + newY*sin(angle)*radius_sqrt + newZ*sqrt(1 - radius);
            sampDir.norm();

            return hitObj.emission + hitColor * radiance(scene, Ray(hitPt, sampDir), depth, Xi);
        }
        else
        {
            Ray reflRay(hitPt, reflect(ray.dir, hitNorm));
            return hitObj.emission + hitColor * radiance(scene, reflRay, depth, Xi);
        }
    }
};

// recursive version
struct SimpleRenderer : public IRenderer
{
    Vec radiance(const IScene& scene, const Ray& ray, int depth, RandomLCG& Xi)
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
        double p = max(f.r, f.g, f.b);

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
            Vec u = cross(fabs(w.x)>.1 ? Vec(0,1) : Vec(1), w);
            u.norm();
            Vec v = cross(w, u);
            Vec d = u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1 - r2);
            d.norm();

            return hitObj.emission + f * radiance(scene, Ray(hitPt,d),depth,Xi);
        }

        Ray reflRay(hitPt, reflect(ray.dir, n));

        if (hitObj.refl == SPECULAR)
        {
            return hitObj.emission + f * radiance(scene, reflRay,depth,Xi);
        }
        // REFRACT
        // otherwise we have a dielectric (glass) surface
        bool into = n.dot(nl)>0;                // Ray from outside going in?
        double nc=1, nt=1.5;
        double nnt=into?nc/nt:nt/nc;
        double ddn=ray.dir.dot(nl);
        double cos2t;
        // if total internal reflection, reflect
        if ((cos2t=1-nnt*nnt*(1-ddn*ddn))<0)    // Total internal reflection
            return hitObj.emission + f * radiance(scene, reflRay,depth,Xi);
        // otherwise, choose reflection or refraction
        Vec tdir = (ray.dir*nnt - n*((into?1:-1)*(ddn*nnt+sqrt(cos2t)))).norm();
        double a=nt-nc, b=nt+nc, R0=a*a/(b*b), c = 1-(into?-ddn:tdir.dot(n));
        double Re=R0+(1-R0)*c*c*c*c*c,Tr=1-Re,P=.25+.5*Re,RP=Re/P,TP=Tr/(1-P);
        return hitObj.emission + f * (depth>2 ? (Xi()<P ?   // Russian roulette
            radiance(scene, reflRay,depth,Xi)*RP:radiance(scene, Ray(hitPt,tdir),depth,Xi)*TP) :
        radiance(scene, reflRay,depth,Xi)*Re+radiance(scene, Ray(hitPt,tdir),depth,Xi)*Tr);
    }
};
