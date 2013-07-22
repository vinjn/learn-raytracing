#pragma once

#include "Scene.h"

struct IRenderer
{
    virtual Vec radiance(const IScene& scene, const Ray& ray, int depth, RandomLCG& Xi) = 0;
};


// diffuse version
struct DiffuseOnlyRenderer : public IRenderer
{
    Vec radiance(const IScene& scene, const Ray& ray, int depth, RandomLCG& Xi)
    {
        double t;                               // distance to intersection
        int id = 0;                               // id of intersected object
        if (!scene.intersect(ray, t, id)) 
            return Vec(); // if miss, return black
        const Geometry& hitObj = *scene.getGeometry(id); 
        Vec hitPt = ray.orig + ray.dir*t;
        Vec hitNorm=(hitPt - hitObj.pos).norm();
        Vec nl = hitNorm.dot(ray.dir)<0 ? hitNorm:hitNorm*-1;
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

        //if (hitObj.refl == DIFFUSE)
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

            return hitObj.emission + f.mult(radiance(scene, Ray(hitPt,d),depth,Xi));
        }
    }
};

// non-recursive version
struct ForwardRenderer : public IRenderer
{
    Vec radiance(const IScene& scene, const Ray& r_, int depth_, RandomLCG& Xi)
    {
        double t;                               // distance to intersection
        int id = 0;                               // id of intersected object
        Ray r = r_;
        int depth = depth_;
        // L0 = Le0 + f0*(L1)
        //    = Le0 + f0*(Le1 + f1*L2)
        //    = Le0 + f0*(Le1 + f1*(Le2 + f2*(L3))
        //    = Le0 + f0*(Le1 + f1*(Le2 + f2*(Le3 + f3*(L4)))
        //    = ...
        //    = Le0 + f0*Le1 + f0*f1*Le2 + f0*f1*f2*Le3 + f0*f1*f2*f3*Le4 + ...
        // 
        // So:
        // F = 1
        // while (1){
        //   L + =  F*Lei
        //   F * =  fi
        // }
        Vec cl(0,0,0);   // accumulated color
        Vec cf(1,1,1);  // accumulated reflectance
        while (1)
        {
            if (!scene.intersect(r, t, id)) 
                return cl; // if miss, return black

            const Geometry &obj = *scene.getGeometry(id);        // the hit object
            Vec x = r.orig+r.dir*t;
            Vec n = (x-obj.pos).norm();
            Vec nl = n.dot(r.dir)<0?n:n*-1;
            Vec f = obj.color;
            double p = max(f.r, f.g, f.b);

            cl = cl + cf.mult(obj.emission);
            if (++depth>5) 
            {
                if (Xi()<p) 
                    f = f*(1/p); 
                return cl; //R.R.
            }
            cf = cf.mult(f);
            if (obj.refl == DIFFUSE)
            {
                // Ideal DIFFUSE reflection
                double r1 = 2*M_PI*Xi(), r2 = Xi(), r2s = sqrt(r2);
                Vec w = nl, u = ((fabs(w.x)>.1?Vec(0,1):Vec(1)).cross(w)).norm(), v = w.cross(u);
                Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).norm();
                //return obj.e + f.mult(radiance(Ray(x,d),depth,Xi));
                r = Ray(x,d);
                continue;
            } 

            Ray reflRay(x, Vec::reflect(r.dir, n));

            if (obj.refl == SPECULAR)
            {       
                // Ideal SPECULAR reflection
                //return obj.e + f.mult(radiance(Ray(x,r.dir-n*2*n.dot(r.dir)),depth,Xi));
                r = reflRay;
                continue;
            }
            
            // Ideal dielectric REFRACTION
            bool into = n.dot(nl)>0;                // Ray from outside going in?
            double nc = 1, nt = 1.5, nnt = into?nc/nt:nt/nc, ddn = r.dir.dot(nl), cos2t;
            if ((cos2t = 1-nnt*nnt*(1-ddn*ddn))<0)
            {  
                // Total internal reflection
                //return obj.e + f.mult(radiance(reflRay,depth,Xi));
                r = reflRay;
                continue;
            }
            Vec tdir = (r.dir*nnt - n*((into?1:-1)*(ddn*nnt+sqrt(cos2t)))).norm();
            double a = nt-nc, b = nt+nc, R0 = a*a/(b*b), c = 1-(into?-ddn:tdir.dot(n));
            double Re = R0+(1-R0)*c*c*c*c*c,Tr = 1-Re,P = .25+.5*Re,RP = Re/P,TP = Tr/(1-P);
            // return obj.e + f.mult(Xi()<P ?
            //                       radiance(reflRay,    depth,Xi)*RP:
            //                       radiance(Ray(x,tdir),depth,Xi)*TP);
            if (Xi()<P)
            {
                cf = cf*RP;
                r = reflRay;
            }
            else
            {
                cf = cf*TP;
                r = Ray(x,tdir);
            }
            continue;
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
            Vec u=(fabs(w.x)>.1 ? Vec(0,1) : Vec(1)).cross(w);
            u.norm();
            Vec v = w.cross(u);
            Vec d = u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1 - r2);
            d.norm();

            return hitObj.emission + f.mult(radiance(scene, Ray(hitPt,d),depth,Xi));
        }

        Ray reflRay(hitPt, Vec::reflect(ray.dir, n));

        if (hitObj.refl == SPECULAR)
        {
            return hitObj.emission + f.mult(radiance(scene, reflRay,depth,Xi));
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
            return hitObj.emission + f.mult(radiance(scene, reflRay,depth,Xi));
        // otherwise, choose reflection or refraction
        Vec tdir = (ray.dir*nnt - n*((into?1:-1)*(ddn*nnt+sqrt(cos2t)))).norm();
        double a=nt-nc, b=nt+nc, R0=a*a/(b*b), c = 1-(into?-ddn:tdir.dot(n));
        double Re=R0+(1-R0)*c*c*c*c*c,Tr=1-Re,P=.25+.5*Re,RP=Re/P,TP=Tr/(1-P);
        return hitObj.emission + f.mult(depth>2 ? (Xi()<P ?   // Russian roulette
            radiance(scene, reflRay,depth,Xi)*RP:radiance(scene, Ray(hitPt,tdir),depth,Xi)*TP) :
        radiance(scene, reflRay,depth,Xi)*Re+radiance(scene, Ray(hitPt,tdir),depth,Xi)*Tr);
    }
};
