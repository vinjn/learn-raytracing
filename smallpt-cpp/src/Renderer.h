#pragma once

#include "Scene.h"

struct IRenderer
{
    // compute the radiance estimate along ray
    virtual Vec radiance(const IScene& scene, const Ray& ray, int depth, RandomLCG& Xi) = 0;
};

struct SimpleRenderer : public IRenderer
{
    Vec radiance(const IScene& scene, const Ray& ray, int depth, RandomLCG& Xi)
    {
        double t;                               // distance to intersection
        int id = 0;                               // id of intersected object
        if (!scene.intersect(ray, t, id)) 
            return Vec(); // if miss, return black
        const Geometry& hitObj = scene.getGeometry(id); 
        Vec hitPt = ray.orig + ray.dir*t;
        Vec n = (hitPt - hitObj.pos).norm();
        Vec nl = dot(n, ray.dir) < 0 ? n : n*-1;
        Vec hitColor = hitObj.color;
        double p = max(hitColor.r, hitColor.g, hitColor.b);

        if (++depth > 5) 
        {
            if (Xi() < p) 
                hitColor = hitColor * (1/p); 
            else 
                return hitObj.emission; // Russian roulette
        }

        if (hitObj.refl == DIFFUSE)
        {
            double r1 = 2*M_PI*Xi();
            double r2 = Xi();
            double r2s = sqrt(r2);

            Vec u, v, w = nl;
            buildUVFromW(u, v, w);
            Vec cosineV(cos(r1)*r2s, sin(r1)*r2s, sqrt(1 - r2));
            Vec d = transformVec(u, v, w, cosineV).norm();

            return hitObj.emission + hitColor * radiance(scene, Ray(hitPt,d), depth, Xi);
        }

        Ray reflRay(hitPt, reflect(ray.dir, n));

        if (hitObj.refl == SPECULAR)
        {
            return hitObj.emission + hitColor * radiance(scene, reflRay,depth,Xi);
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
            return hitObj.emission + hitColor * radiance(scene, reflRay,depth,Xi);
        // otherwise, choose reflection or refraction
        Vec tdir = (ray.dir*nnt - n*((into?1:-1)*(ddn*nnt+sqrt(cos2t)))).norm();
        double a=nt-nc, b=nt+nc, R0=a*a/(b*b), c = 1-(into?-ddn:tdir.dot(n));
        double Re=R0+(1-R0)*c*c*c*c*c,Tr=1-Re,P=.25+.5*Re,RP=Re/P,TP=Tr/(1-P);
        return hitObj.emission + hitColor * (depth>2 ? (Xi()<P ?   // Russian roulette
            radiance(scene, reflRay,depth,Xi)*RP:radiance(scene, Ray(hitPt,tdir),depth,Xi)*TP) :
        radiance(scene, reflRay,depth,Xi)*Re+radiance(scene, Ray(hitPt,tdir),depth,Xi)*Tr);
    }
};

struct ExplicitRenderer : public IRenderer
{
    Vec radiance(const IScene& scene, const Ray& ray, int depth, RandomLCG& Xi)
    {
        return radiance(scene, ray, depth, Xi, 1);
    }

private:
    Vec radiance(const IScene& scene, const Ray& ray, int depth, RandomLCG& Xi, int E)
    {
        double t;                               // distance to intersection
        int id = 0;                               // id of intersected object
        if (!scene.intersect(ray, t, id)) 
            return Vec(); // if miss, return black
        const Geometry& hitObj = scene.getGeometry(id); 
        Vec hitPt = ray.orig + ray.dir*t;
        Vec n=(hitPt - hitObj.pos).norm();
        Vec nl = n.dot(ray.dir)<0 ? n:n*-1;
        Vec hitColor = hitObj.color;
        double p = max(hitColor.r, hitColor.g, hitColor.b);

        if (++depth > 5 || !p) 
        {
            if (Xi() < p) 
                hitColor = hitColor * (1/p); 
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

            Vec u, v, w = nl;
            buildUVFromW(u, v, w);
            Vec cosineV(cos(r1)*r2s, sin(r1)*r2s, sqrt(1 - r2));
            Vec d = transformVec(u, v, w, cosineV).norm();

            // Loop over any lights
            Vec e;
            for (int i = 0; i < scene.getGeometryCount(); i++){
                const Geometry& geom = scene.getGeometry(i);
                if (geom.emission.x <= 0 && geom.emission.y <= 0 && geom.emission.z <= 0) 
                    continue; // skip non-lights

                Vec su, sv, sw = geom.pos - hitPt;
                buildUVFromW(su, sv, sw);
                double cos_a_max = sqrt(1 - geom.radius * geom.radius / dot(hitPt-geom.pos, hitPt-geom.pos));
                double eps1 = Xi(), eps2 = Xi();
                double cos_a = 1 - eps1 + eps1 * cos_a_max;
                double sin_a = sqrt(1 - cos_a * cos_a);
                double phi = 2 * M_PI * eps2;
                Vec randV(cos(phi) * sin_a, sin(phi) * sin_a, cos_a);
                Vec l = transformVec(su, sv, sw, randV).norm();
                if (scene.intersect(Ray(hitPt, l), t, id) && id == i)
                {  
                    // shadow ray
                    double omega = 2*M_PI*(1-cos_a_max);
                    e = e + hitColor * geom.emission * dot(l, nl) * omega * M_1_PI;  // 1/pi for brdf
                }
            }

            return hitObj.emission*E + e + hitColor * radiance(scene, Ray(hitPt, d), depth, Xi, 0);

        }

        Ray reflRay(hitPt, reflect(ray.dir, n));

        if (hitObj.refl == SPECULAR)
        {
            return hitObj.emission + hitColor * radiance(scene, reflRay, depth, Xi, 1);
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
            return hitObj.emission + hitColor * radiance(scene, reflRay,depth,Xi);
        // otherwise, choose reflection or refraction
        Vec tdir = (ray.dir*nnt - n*((into?1:-1)*(ddn*nnt+sqrt(cos2t)))).norm();
        double a=nt-nc, b=nt+nc, R0=a*a/(b*b), c = 1-(into?-ddn:tdir.dot(n));
        double Re=R0+(1-R0)*c*c*c*c*c,Tr=1-Re,P=.25+.5*Re,RP=Re/P,TP=Tr/(1-P);
        return hitObj.emission + hitColor * (depth>2 ? (Xi()<P ?   // Russian roulette
            radiance(scene, reflRay,depth,Xi)*RP:radiance(scene, Ray(hitPt,tdir),depth,Xi)*TP) :
        radiance(scene, reflRay,depth,Xi)*Re+radiance(scene, Ray(hitPt,tdir),depth,Xi)*Tr);
    }
};