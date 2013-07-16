// smallpt, a Path Tracer by Kevin Beason, 2008
// Usage: time ./smallpt 5000 && xv image.ppm
// refactored by vinjn.z@gmail.com

#include <stdlib.h> // Make : g++ -O3 -fopenmp smallpt.cpp -o smallpt
#include <stdio.h>  //        Remove "-fopenmp" for g++ version < 4.2
#include <time.h>       // MILO

#include "Utils.h"
#include "Scene.h"
#include "../getopt/getopt.h"

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


int main(int argc, char *argv[])
{
    clock_t start = clock(); // MILO
    int w = 256, h = 256;
    std::string sceneName = "cornell";
    int samps = 25; // # samples

    option long_options[] =
    {
        /* These options set a flag. */
        {"help",    no_argument,        0, 'h'},
        //{"brief",   no_argument,       &verbose_flag, 0},
        /* These options don't set a flag.
        We distinguish them by their indices. */
        {"viewport",required_argument,  0, 'v'},
        {"input",   required_argument,  0, 'i'},
        {"samples", required_argument,  0, 's'},
        {0, 0, 0, 0}
    };

    int opt;
    int option_index = 0;
    while ((opt = getopt_long (argc, argv, "hv:i:s:", long_options, &option_index)) != -1)
    {
        switch (opt)
        {
        case 'h':
            {
                fprintf(stdout, "smallpt --size 256 --scene cornell --samples 100\n\n");
                Scene::printBuiltInSceneNames();
                exit(0);
            }
        case 'v':{w = h = atoi(optarg); break;}
        case 'i':{sceneName = optarg; break;}
        case 's':{samps = atoi(optarg)/4; break;}
        default: exit(0);
        }
    }
    if (optind < argc) {
        printf ("non-option ARGV-elements: ");
        while (optind < argc)
            printf ("%s ", argv[optind++]);
        printf ("\n");
        exit(0);
    }

    Ray cam(Vec(50,52,295.6), Vec(0,-0.042612,-1).norm()); // cam pos, dir
    Vec cx = Vec(w*.5135/h);
    Vec cy=(cx.cross(cam.dir)).norm()*.5135;
    Vec *c = new Vec[w*h];
    Vec r;

    Scene scene;
    if (!Scene::createBuiltInScene(&scene, sceneName))
    {
        fprintf(stderr, "Invalid scene name: %s\n", sceneName.c_str());
        Scene::printBuiltInSceneNames();
        exit(1);
    }

    fprintf(stdout, "Rendering scene: %s\n", sceneName.c_str());

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
    char imageName[100];
    sprintf(imageName, "%s_%dX%d.ppm", sceneName.c_str(), w, h);
    FILE *f = fopen(imageName, "w");         // Write image to PPM file.
    fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
    for (int i = 0; i < w*h; i++)
        fprintf(f,"%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
    fclose(f);
}