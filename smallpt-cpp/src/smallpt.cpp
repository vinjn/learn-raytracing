// smallpt, a Path Tracer by Kevin Beason, 2008
// Usage: time ./smallpt 5000 && xv image.ppm
// refactored by vinjn.z@gmail.com

#include <stdlib.h> // Make : g++ -O3 -fopenmp smallpt.cpp -o smallpt
#include <stdio.h>  //        Remove "-fopenmp" for g++ version < 4.2
#include <time.h>       // MILO

#include "Utils.h"
#include "Scene.h"

Vec Cen(50,40.8,-860);

Sphere spheres[] = {//Scene: radius, position, emission, color, material
    // center 50 40.8 62
    // floor 0
    // back  0

    Sphere(1600, Vec(1,0,2)*3000, Vec(1,.9,.8)*1.2e1*1.56*2,Vec(), DIFFUSE), // sun
    Sphere(1560, Vec(1,0,2)*3500,Vec(1,.5,.05)*4.8e1*1.56*2, Vec(),  DIFFUSE), // horizon sun2
    //   Sphere(10000,Cen+Vec(0,0,-200), Vec(0.0627, 0.188, 0.569)*6e-2*8, Vec(.7,.7,1)*.25,  DIFF), // sky
    Sphere(10000,Cen+Vec(0,0,-200), Vec(0.00063842, 0.02001478, 0.28923243)*6e-2*8, Vec(.7,.7,1)*.25,  DIFFUSE), // sky

    Sphere(100000, Vec(50, -100000, 0),  Vec(),Vec(.3,.3,.3),DIFFUSE), // grnd
    Sphere(110000, Vec(50, -110048.5, 0),  Vec(.9,.5,.05)*4,Vec(),DIFFUSE),// horizon brightener
    Sphere(4e4, Vec(50, -4e4-30, -3000),  Vec(),Vec(.2,.2,.2),DIFFUSE),// mountains
    //  Sphere(3.99e4, Vec(50, -3.99e4+20.045, -3000),  Vec(),Vec(.7,.7,.7),DIFF),// mountains snow

    Sphere(26.5,Vec(22,26.5,42),   Vec(),Vec(1,1,1)*.596, SPECULAR), // white Mirr
    Sphere(13,Vec(75,13,82),   Vec(),Vec(.96,.96,.96)*.96, REFRACT),// Glas
    Sphere(22,Vec(87,22,24),   Vec(),Vec(.6,.6,.6)*.696, REFRACT)    // Glas2
};

Vec radiance(const Scene& scene, const Ray& ray, int depth, RandomLCG& Xi)
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

        return hitObj.emission + f.mult(radiance(scene, Ray(hitPt,d),depth,Xi));
    }
    else
    {
        Ray reflRay(hitPt, ray.dir - n*2*n.dot(ray.dir)); 
        return hitObj.emission + f.mult(radiance(scene, reflRay,depth,Xi));
    }
}


int main(int argc, char *argv[]){
    clock_t start = clock(); // MILO
    int w = 256, h = 256;
    int samps = argc==2 ? atoi(argv[1])/4 : 25; // # samples
    Ray cam(Vec(50,52,295.6), Vec(0,-0.042612,-1).norm()); // cam pos, dir
    Vec cx = Vec(w*.5135/h);
    Vec cy=(cx.cross(cam.dir)).norm()*.5135;
    Vec *c = new Vec[w*h];
    Vec r;

    // setup a simple scene             radius     position                emission        color               material
    Scene scene;

    scene.addGeometries(spheres, _countof(spheres));

    scene.addGeometry(new Sphere(1e5,     Vec( 1e5+1,40.8,81.6),   Vec(),          Vec(.75,.25,.25),   DIFFUSE));// left wall 
    scene.addGeometry(new Sphere(1e5,     Vec(-1e5+99,40.8,81.6),  Vec(),          Vec(.25,.25,.75),   DIFFUSE));// right wall
    scene.addGeometry(new Sphere(1e5,     Vec(50,40.8, 1e5),       Vec(),          Vec(.75,.75,.75),   DIFFUSE));// back wall
    scene.addGeometry(new Sphere(1e5,     Vec(50,40.8,-1e5+170),   Vec(),          Vec(),              DIFFUSE));// front wall
    scene.addGeometry(new Sphere(1e5,     Vec(50, 1e5, 81.6),      Vec(),          Vec(.75,.75,.75),   DIFFUSE));// bottom wall
    scene.addGeometry(new Sphere(1e5,     Vec(50,-1e5+81.6,81.6),  Vec(),          Vec(.75,.75,.75),   DIFFUSE));// top wall

    scene.addGeometry(new Sphere(16.5,    Vec(27,16.5,47),         Vec(),          Vec(1,1,1)*.999,    SPECULAR));// mirror 1
    scene.addGeometry(new Sphere(20.5,    Vec(73,16.5,78),         Vec(),          Vec(1,1,1)*.999,    SPECULAR));// mirror 2
    scene.addGeometry(new Sphere(10.5,    Vec(27,16.5,100),        Vec(),          Vec(1,1,1)*.999,    SPECULAR));// mirror 3

    scene.addGeometry(new Sphere(600,     Vec(50,681.6-.27,81.6),  Vec(12,12,12),  Vec(),              DIFFUSE)); // the white light

#pragma omp parallel for schedule(dynamic, 1) private(r)       // OpenMP
    for (int y = 0; y < h; y++){                       // Loop over image rows
        fprintf(stderr,"\rRendering (%d spp) %5.2f%%",samps*4,100.*y/(h - 1));
        RandomLCG Xi(y*y*y); // MILO
        for (unsigned short x = 0; x < w; x++)   // Loop cols
            for (int sy = 0, i=(h - y-1)*w+x; sy < 2; sy++)     // 2x2 subpixel rows
                for (int sx = 0; sx < 2; sx++, r = Vec()){        // 2x2 subpixel cols
                    for (int s = 0; s < samps; s++){
                        double dx = tentFilterRandom(Xi);
                        double dy = tentFilterRandom(Xi);
#if 1
                        Vec d = cx*( ( (sx+.5 + dx)/2 + x)/w - .5) +
                            cy*( ( (sy+.5 + dy)/2 + y)/h - .5) + cam.dir;
#else
                        Vec d = cx*(x/(double)w - .5 ) + cy*(y/(double)h - .5) + cam.dir;
#endif
                        r = r + radiance(scene, Ray(cam.orig+d*140,d.norm()),0,Xi)*(1./samps);
                    } // Camera rays are pushed ^^^^^ forward to start in interior
                    c[i] = c[i] + Vec(clamp(r.x),clamp(r.y),clamp(r.z))*.25;
                }
    }

    printf("\n%f sec\n", (float)(clock() - start)/CLOCKS_PER_SEC); // MILO
    FILE *f = fopen("image.ppm", "w");         // Write image to PPM file.
    fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
    for (int i = 0; i < w*h; i++)
        fprintf(f,"%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
}