// smallpt, a Path Tracer by Kevin Beason, 2008
// Usage: time ./smallpt 5000 && xv image.ppm
// enchanced by vinjn.z@gmail.com

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
    int width = 256, height = 256;
    std::string sceneName = "cornell";
    std::string rendererName = "simple";
    int nSamples = 25; // # samples

    option long_options[] =
    {
        {"help",    no_argument,        0, 'h'},
        {"viewport",required_argument,  0, 'v'},
        {"input",   required_argument,  0, 'i'},
        {"samples", required_argument,  0, 's'},
        {"renderer",required_argument,  0, 'r'},
        {0, 0, 0, 0}
    };

    int opt;
    int option_index = 0;
    while ((opt = getopt_long (argc, argv, "hv:i:s:r:", long_options, &option_index)) != -1)
    {
        switch (opt)
        {
        case 'h':
            {
                fprintf(stdout, "smallpt --renderer simple --viewport 256 --input cornell --samples 100\n\n");
                IScene::printBuiltInSceneNames();
                exit(0);
            }
        case 'v':{width = height = atoi(optarg); break;}
        case 'i':{sceneName = optarg; break;}
        case 's':{nSamples = atoi(optarg)/4; break;}
        case 'r':{rendererName = optarg; break;}
        default: exit(0);
        }
    }
    if (optind < argc) 
    {
        printf ("non-option ARGV-elements: ");
        while (optind < argc)
            printf ("%s ", argv[optind++]);
        printf ("\n");
        exit(0);
    }

    IRenderer* renderer = NULL;
    {
        if (rendererName == "diffuse")
            renderer = new DiffuseOnlyRenderer();
        else
            renderer = new SimpleRenderer();
    }

    const float CAMERA_SCALE = .5135;
    Ray cam(Vec(50,52,295.6), Vec(0,-0.042612,-1).norm()); // cam pos, dir
    Vec cx = Vec(width * CAMERA_SCALE / height);
    Vec cy = (cross(cx, cam.dir)).norm() * CAMERA_SCALE;
    Vec* colorBuf = new Vec[width*height];
    Vec r;

    IScene* scene = new SimpleScene();
    if (!IScene::createBuiltInScene(scene, sceneName))
    {
        fprintf(stderr, "Invalid scene name: %s\n", sceneName.c_str());
        IScene::printBuiltInSceneNames();
        exit(1);
    }

    char imageName[100];
    sprintf(imageName, "%s_%s_%dX%dX%d.ppm", 
        sceneName.c_str(), rendererName.c_str(), 
        width, height, nSamples*4);
    fprintf(stdout, "Rendering to %s\n", imageName);

#pragma omp parallel for schedule(dynamic, 1) private(r)       // OpenMP
    for (int y = 0; y < height; y++)
    {
        // Loop over image rows
        fprintf(stderr,"\rRendering (%d spp) %5.2f%%",nSamples*4,100.*y/(height - 1));
        RandomLCG Xi(y*y*y); // MILO
        for (unsigned short x = 0; x < width; x++)
        {
            // Loop cols
            int loc = (height - y - 1) * width + x;
            for (int sy = 0; sy < 2; sy++) 
            {
                // 2x2 subpixel rows
                for (int sx = 0; sx < 2; sx++, r = Vec())
                {
                    // 2x2 subpixel cols
                    for (int s = 0; s < nSamples; s++)
                    {
                        double dx = tentFilterRandom(Xi);
                        double dy = tentFilterRandom(Xi);
#if 1
                        Vec d = cx*( ( (sx+.5 + dx)/2 + x)/width - .5) +
                            cy*( ( (sy+.5 + dy)/2 + y)/height - .5) + cam.dir;
#else
                        Vec d = cx*(x/(double)width - .5 ) + cy*(y/(double)height - .5) + cam.dir;
#endif
                        r = r + renderer->radiance(*scene, Ray(cam.orig+d*140,d.norm()),0,Xi)*(1./nSamples);
                    } // Camera rays are pushed ^^^^^ forward to start in interior
                    colorBuf[loc] = colorBuf[loc] + Vec(clamp(r.x),clamp(r.y),clamp(r.z))*.25;
                }
            }
        }
    }

    delete renderer;
    delete scene;

    printf("\n%f sec\n", (float)(clock() - start)/CLOCKS_PER_SEC); // MILO

    {
        FILE *f = fopen(imageName, "w");         // Write image to PPM file.
        fprintf(f, "P3\n%d %d\n%d\n", width, height, 255);
        for (int i = 0; i < width*height; i++)
            fprintf(f,"%d %d %d ", toInt(colorBuf[i].x), toInt(colorBuf[i].y), toInt(colorBuf[i].z));
        fclose(f);
    }
}