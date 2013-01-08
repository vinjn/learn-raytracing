#define _USE_MATH_DEFINES
#include <math.h>   // smallpt, a Path Tracer by Kevin Beason, 2008
#include <stdlib.h> // Make : g++ -O3 -fopenmp smallpt.cpp -o smallpt
#include <stdio.h>  //        Remove "-fopenmp" for g++ version < 4.2
#include <time.h>       // MILO

struct RandomLCG {
    unsigned mSeed;
    RandomLCG(unsigned seed = 0) : mSeed(seed) {}
    float operator()() { mSeed = 214013 * mSeed + 2531011; return mSeed * (1.0f / 4294967296); }
};

struct Vec {        // Usage: time ./smallpt 5000 && xv image.ppm
    double x, y, z;                  // position, also color (r,g,b)
    Vec(double x_=0, double y_=0, double z_=0){ x=x_; y=y_; z=z_; }
    Vec operator+(const Vec &b) const { return Vec(x+b.x,y+b.y,z+b.z); }
    Vec operator-(const Vec &b) const { return Vec(x-b.x,y-b.y,z-b.z); }
    Vec operator*(double b) const { return Vec(x*b,y*b,z*b); }
    Vec mult(const Vec &b) const { return Vec(x*b.x,y*b.y,z*b.z); }
    Vec& norm(){ return *this = *this * (1/sqrt(x*x+y*y+z*z)); }
    double dot(const Vec &b) const { return x*b.x+y*b.y+z*b.z; } // cross:
    Vec operator%(const Vec &b){return Vec(y*b.z-z*b.y,z*b.x-x*b.z,x*b.y-y*b.x);}
};

struct Ray { Vec o, d; Ray(const Vec &o_, const Vec &d_) : o(o_), d(d_) {} };

enum Refl_t { 
    DIFFUSE, 
    SPECULAR, 
    REFRACTIVE
};  // material types, used in radiance()

struct Sphere {
    double rad;       // radius
    Vec p, e, c;      // position, emission, color
    Refl_t refl;      // reflection type (DIFFuse, SPECular, REFRactive)

    Sphere(double rad_, Vec p_, Vec e_, Vec c_, Refl_t refl_):
    rad(rad_), p(p_), e(e_), c(c_), refl(refl_) {}

    double intersect(const Ray &r) const { // returns distance, 0 if nohit
        Vec op = p-r.o; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
        double t, eps=1e-4, b=op.dot(r.d), det=b*b-op.dot(op)+rad*rad;
        if (det<0) 
            return 0; 
        else 
            det=sqrt(det);
        return (t=b-det)>eps ? t : ((t=b+det)>eps ? t : 0);
    }
};

Sphere spheres[] = {
    //Scene: radius, position,              emission,       color,              material
    Sphere(1e5,     Vec( 1e5+1,40.8,81.6),  Vec(),          Vec(.75,.25,.25),   DIFFUSE),//Left
    Sphere(1e5,     Vec(-1e5+99,40.8,81.6), Vec(),          Vec(.25,.25,.75),   DIFFUSE),//Rght
    Sphere(1e5,     Vec(50,40.8, 1e5),      Vec(),          Vec(.75,.75,.75),   DIFFUSE),//Back
    Sphere(1e5,     Vec(50,40.8,-1e5+170),  Vec(),          Vec(),              DIFFUSE),//Frnt
    Sphere(1e5,     Vec(50, 1e5, 81.6),     Vec(),          Vec(.75,.75,.75),   DIFFUSE),//Botm
    Sphere(1e5,     Vec(50,-1e5+81.6,81.6), Vec(),          Vec(.75,.75,.75),   DIFFUSE),//Top
    Sphere(16.5,    Vec(27,16.5,47),        Vec(),          Vec(1,1,1)*.999,    SPECULAR),//Mirr
    Sphere(16.5,    Vec(73,16.5,78),        Vec(),          Vec(1,1,1)*.999,    REFRACTIVE),//Glas
    Sphere(1.5,     Vec(50,81.6-16.5,81.6), Vec(4,4,4)*100,  Vec(),             DIFFUSE) //Lite
};

const int N_SPHERES = sizeof(spheres)/sizeof(Sphere);

inline double clamp(double x){ return x<0 ? 0 : x>1 ? 1 : x; }

inline int toInt(double x){ return int(pow(clamp(x),1/2.2)*255+.5); }

inline bool intersect(const Ray &r, double &t, int &id){
    const double INF = t = 1e20;
    for(int i = 0;i < N_SPHERES;i++) 
    {
        double d=spheres[i].intersect(r);
        if (d && d<t)
        {
            t=d;
            id=i;
        }
    }
    return t<INF;
}

Vec radiance(const Ray &r, int depth, RandomLCG &Xi, int E=1){
    double t;                               // distance to intersection
    int id=0;                               // id of intersected object
    if (!intersect(r, t, id)) 
        return Vec(); // if miss, return black
    const Sphere &obj = spheres[id];        // the hit object
    Vec x = r.o+r.d*t;
    Vec n=(x-obj.p).norm();
    Vec nl=n.dot(r.d)<0?n:n*-1;
    Vec f=obj.c;
    double p = f.x>f.y && f.x>f.z ? f.x : f.y>f.z ? f.y : f.z; // max refl

    if (++depth > 5 || !p) 
    {
        if (Xi()<p) 
            f = f*(1/p); 
        else 
            return obj.e; //R.R. Russian roulette
    }
    if (depth > 100) 
        return obj.e; // MILO

    if (obj.refl == DIFFUSE){                  // Ideal DIFFUSE reflection
        double r1=2*M_PI*Xi();
        double r2=Xi();
        double r2s=sqrt(r2);
        Vec w=nl;
        Vec u=(fabs(w.x)>.1 ? Vec(0,1) : Vec(1))%w;
        u.norm();
        Vec v=w%u;
        Vec d = u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2);
        d.norm();

        // Loop over any lights
        Vec e;
        for (int i=0; i<N_SPHERES; i++){
            const Sphere &s = spheres[i];
            if (s.e.x<=0 && s.e.y<=0 && s.e.z<=0) 
                continue; // skip non-lights

            Vec sw=s.p-x, su=((fabs(sw.x)>.1?Vec(0,1):Vec(1))%sw).norm(), sv=sw%su;
            double cos_a_max = sqrt(1-s.rad*s.rad/(x-s.p).dot(x-s.p));
            double eps1 = Xi(), eps2 = Xi();
            double cos_a = 1-eps1+eps1*cos_a_max;
            double sin_a = sqrt(1-cos_a*cos_a);
            double phi = 2*M_PI*eps2;
            Vec l = su*cos(phi)*sin_a + sv*sin(phi)*sin_a + sw*cos_a;
            l.norm();
            if (intersect(Ray(x,l), t, id) && id==i){  // shadow ray
                double omega = 2*M_PI*(1-cos_a_max);
                e = e + f.mult(s.e*l.dot(nl)*omega)*M_1_PI;  // 1/pi for brdf
            }
        }

        return obj.e*E + e + f.mult(radiance(Ray(x,d),depth,Xi,0));
    } 
    else
    {
        Ray reflRay(x, r.d-n*2*n.dot(r.d));     // Ideal dielectric REFRACTION
        if (obj.refl == SPECULAR){            // Ideal SPECULAR reflection
            return obj.e + f.mult(radiance(reflRay,depth,Xi));
        }

        bool into = n.dot(nl)>0;                // Ray from outside going in?
        double nc=1, nt=1.5, nnt=into?nc/nt:nt/nc, ddn=r.d.dot(nl), cos2t;
        if ((cos2t=1-nnt*nnt*(1-ddn*ddn))<0)    // Total internal reflection
            return obj.e + f.mult(radiance(reflRay,depth,Xi));
        Vec tdir = r.d*nnt - n*((into?1:-1)*(ddn*nnt+sqrt(cos2t)));
        tdir.norm();
        double a=nt-nc, b=nt+nc, R0=a*a/(b*b), c = 1-(into?-ddn:tdir.dot(n));
        double Re=R0+(1-R0)*c*c*c*c*c,Tr=1-Re,P=.25+.5*Re,RP=Re/P,TP=Tr/(1-P);
        return obj.e + f.mult(depth>2 ? (Xi()<P ?   // Russian roulette
            radiance(reflRay,depth,Xi)*RP:radiance(Ray(x,tdir),depth,Xi)*TP) :
        radiance(reflRay,depth,Xi)*Re+radiance(Ray(x,tdir),depth,Xi)*Tr);
    }

}

double centerRandom(RandomLCG& Xi)
{
    double r = 2*Xi();
    return r < 1 ? sqrt(r)-1: 1-sqrt(2-r);
}

int main(int argc, char *argv[]){
    clock_t start = clock(); // MILO
    int w=256, h=256;
    int samps = argc==2 ? atoi(argv[1])/4 : 25; // # samples
    Ray cam(Vec(50,52,295.6), Vec(0,-0.042612,-1).norm()); // cam pos, dir
    Vec cx=Vec(w*.5135/h);
    Vec cy=(cx%cam.d).norm()*.5135;
    Vec *c=new Vec[w*h];
    Vec r;

#pragma omp parallel for schedule(dynamic, 1) private(r)       // OpenMP
    for (int y=0; y<h; y++){                       // Loop over image rows
        fprintf(stderr,"\rRendering (%d spp) %5.2f%%",samps*4,100.*y/(h-1));
        RandomLCG Xi(y*y*y); // MILO
        for (unsigned short x=0; x<w; x++)   // Loop cols
            for (int sy=0, i=(h-y-1)*w+x; sy<2; sy++)     // 2x2 subpixel rows
                for (int sx=0; sx<2; sx++, r = Vec()){        // 2x2 subpixel cols
                    for (int s=0; s<samps; s++){
                        double dx = centerRandom(Xi);
                        double dy = centerRandom(Xi);
                        Vec d = cx*( ( (sx+.5 + dx)/2 + x)/w - .5) +
                            cy*( ( (sy+.5 + dy)/2 + y)/h - .5) + cam.d;
                        r = r + radiance(Ray(cam.o+d*140,d.norm()),0,Xi)*(1./samps);
                    } // Camera rays are pushed ^^^^^ forward to start in interior
                    c[i] = c[i] + Vec(clamp(r.x),clamp(r.y),clamp(r.z))*.25;
                }
    }

    printf("\n%f sec\n", (float)(clock() - start)/CLOCKS_PER_SEC); // MILO
    FILE *f = fopen("image.ppm", "w");         // Write image to PPM file.
    fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
    for (int i=0; i<w*h; i++)
        fprintf(f,"%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
}