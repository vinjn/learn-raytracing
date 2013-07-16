// smallpt, a Path Tracer by Kevin Beason, 2008
// Usage: time ./smallpt 5000 && xv image.ppm
// refactored by vinjn.z@gmail.com

#include <stdlib.h> // Make : g++ -O3 -fopenmp smallpt.cpp -o smallpt
#include <stdio.h>  //        Remove "-fopenmp" for g++ version < 4.2
#include <time.h>       // MILO

#include "Utils.h"
#include "Scene.h"
#include "../getopt/getopt.h"
#include "Renderer.h"

int main(int argc, char *argv[])
{
    clock_t start = clock(); // MILO
    int w = 256, h = 256;
    std::string sceneName = "cornell";
    int samps = 25; // # samples

    IRenderer* renderer = new SimpleRenderer();

    option long_options[] =
    {
        {"help",    no_argument,        0, 'h'},
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
                        r = r + renderer->render(scene, Ray(cam.orig+d*140,d.norm()),0,Xi)*(1./samps);
                    } // Camera rays are pushed ^^^^^ forward to start in interior
                    c[i] = c[i] + Vec(clamp(r.x),clamp(r.y),clamp(r.z))*.25;
                }
    }

    delete renderer;

    printf("\n%f sec\n", (float)(clock() - start)/CLOCKS_PER_SEC); // MILO

    {
        char imageName[100];
        sprintf(imageName, "%s_%dX%d.ppm", sceneName.c_str(), w, h);
        FILE *f = fopen(imageName, "w");         // Write image to PPM file.
        fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
        for (int i = 0; i < w*h; i++)
            fprintf(f,"%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
        fclose(f);
    }
}