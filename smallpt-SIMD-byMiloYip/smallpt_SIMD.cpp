#include <math.h>   // smallpt, a Path Tracer by Kevin Beason, 2008
#include <stdlib.h> // Make : g++ -O3 -fopenmp smallpt.cpp -o smallpt
#include <stdio.h>  //        Remove "-fopenmp" for g++ version < 4.2
#include <time.h>		// MILO
#include <windows.h>
#include <xnamath.h>
#define M_PI 3.141592653589793238462643	// MILO
struct RandomLCG {
	unsigned mSeed;
	RandomLCG(unsigned seed = 0) : mSeed(seed) {}
	float operator()() { mSeed = 214013 * mSeed + 2531011; return mSeed * (1.0f / 4294967296); }
};
struct Ray { XMVECTOR o, d; Ray(const XMVECTOR o_, const XMVECTOR d_) : o(o_), d(d_) {} };
enum Refl_t { DIFF, SPEC, REFR };  // material types, used in radiance()
struct Sphere {
  XMVECTOR p, e, c;      // position, emission, color
  XMVECTOR sqrad;
  float rad;       // radius
  Refl_t refl;      // reflection type (DIFFuse, SPECular, REFRactive)
  float maxRefl;
  Sphere(float rad_, XMVECTOR p_, XMVECTOR e_, XMVECTOR c_, Refl_t refl_):
    rad(rad_), p(p_), e(e_), c(c_), refl(refl_) { sqrad=XMVectorSet(rad * rad,0,0,0); 
	maxRefl = XMVectorGetX(c)>XMVectorGetY(c) && XMVectorGetX(c)>XMVectorGetZ(c) ? XMVectorGetX(c) : XMVectorGetY(c)>XMVectorGetZ(c) ? XMVectorGetY(c) : XMVectorGetZ(c);}
  float intersect(const Ray &r) const { // returns distance, 0 if nohit
    XMVECTOR op = p-r.o; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
    float t, eps=1e-4, b=XMVectorGetX(XMVector3Dot(op, r.d)), det=b*b-XMVectorGetX(XMVector3Dot(op, op)-sqrad);
    if (det<0) return 0; else det=sqrtf(det);
    return (t=b-det)>eps ? t : ((t=b+det)>eps ? t : 0);
  }
};
#define X 1000
Sphere spheres[] = {//Scene: radius, position, emission, color, material
  Sphere(X, XMVectorSet( X+1,40.8f,81.6f,0), XMVectorZero(),XMVectorSet(.75f,.25f,.25f,0),DIFF),//Left
  Sphere(X, XMVectorSet(-X+99,40.8f,81.6f,0),XMVectorZero(),XMVectorSet(.25f,.25f,.75f,0),DIFF),//Rght
  Sphere(X, XMVectorSet(50,40.8f, X,0),     XMVectorZero(),XMVectorSet(.75f,.75f,.75f,0),DIFF),//Back
  Sphere(X, XMVectorSet(50,40.8f,-X+170,0), XMVectorZero(),XMVectorZero(),            DIFF),//Frnt
  Sphere(X, XMVectorSet(50, X, 81.6f,0),    XMVectorZero(),XMVectorSet(.75f,.75f,.75f,0),DIFF),//Botm
  Sphere(X, XMVectorSet(50,-X+81.6f,81.6f,0),XMVectorZero(),XMVectorSet(.75f,.75f,.75f,0),DIFF),//Top
  Sphere(16.5,XMVectorSet(27,16.5f,47,0),       XMVectorZero(),XMVectorSet(1,1,1,0)*.999f, SPEC),//Mirr
  Sphere(16.5,XMVectorSet(73,16.5f,78,0),       XMVectorZero(),XMVectorSet(1,1,1,0)*.999f, REFR),//Glas
  Sphere(600, XMVectorSet(50,681.6f-.27f,81.6f,0),XMVectorSet(12,12,12,0),  XMVectorZero(), DIFF) //Lite
};
#undef X
inline float clamp(float x){ return x<0 ? 0 : x>1 ? 1 : x; }
inline int toInt(float x){ return int(powf(clamp(x),1/2.2f)*255+.5f); }
inline bool intersect(const Ray &r, float &t, int &id){
  float n=sizeof(spheres)/sizeof(Sphere), d, inf=t=1e20f;
  for(int i=int(n);i--;) if((d=spheres[i].intersect(r))&&d<t){t=d;id=i;}
  return t<inf;
}
XMVECTOR radiance(const Ray &r, int depth, RandomLCG &Xi){
  float t;                               // distance to intersection
  int id=0;                               // id of intersected object
  if (!intersect(r, t, id)) return XMVECTOR(); // if miss, return black
  const Sphere &obj = spheres[id];        // the hit object
  XMVECTOR x=r.o+r.d*t, n=XMVector3Normalize(x-obj.p), nl=XMVectorGetX(XMVector3Dot(n, r.d))<0?n:-n, refl=obj.c;
  float p = obj.maxRefl; // max refl
  if (++depth>5) if (Xi()<p) refl=refl*(1/p); else return obj.e; //R.R.
  if (depth > 100) return obj.e; // MILO
  if (obj.refl == DIFF){                  // Ideal DIFFUSE reflection
    float r1=2*XM_PI*Xi(), r2=Xi(), r2s=sqrtf(r2);
    XMVECTOR w=nl, u=XMVector3Normalize(XMVector3Cross((fabs(XMVectorGetX(w))>.1?XMVectorSet(0,1,0,0):XMVectorSet(1,0,0,0)), w)), v=XMVector3Cross(w, u);
    XMVECTOR d = XMVector3Normalize((u*(cosf(r1)*r2s) + v*(sinf(r1)*r2s) + w*sqrtf(1-r2)));
    return obj.e + refl * radiance(Ray(x,d),depth,Xi);
  } else if (obj.refl == SPEC)            // Ideal SPECULAR reflection
    return obj.e + refl * radiance(Ray(x,r.d-n*2*XMVector3Dot(n,r.d)),depth,Xi);
  Ray reflRay(x, r.d-n*2*XMVector3Dot(n, r.d));     // Ideal dielectric REFRACTION
  bool into = XMVectorGetX(XMVector3Dot(n,nl))>0;                // Ray from outside going in?
  float nc=1, nt=1.5, nnt=into?nc/nt:nt/nc, ddn=XMVectorGetX(XMVector3Dot(r.d,nl)), cos2t;
  if ((cos2t=1-nnt*nnt*(1-ddn*ddn))<0)    // Total internal reflection
    return obj.e + refl * radiance(reflRay,depth,Xi);
  XMVECTOR tdir = XMVector3Normalize(r.d*nnt - n*((into?1:-1)*(ddn*nnt+sqrtf(cos2t))));
  float a=nt-nc, b=nt+nc, R0=a*a/(b*b), c = 1-(into?-ddn:XMVectorGetX(XMVector3Dot(tdir,n)));
  float Re=R0+(1-R0)*c*c*c*c*c,Tr=1-Re,P=.25f+.5f*Re,RP=Re/P,TP=Tr/(1-P);
  return obj.e + refl * (depth>2 ? (Xi()<P ?   // Russian roulette
    radiance(reflRay,depth,Xi)*RP:radiance(Ray(x,tdir),depth,Xi)*TP) :
    radiance(reflRay,depth,Xi)*Re+radiance(Ray(x,tdir),depth,Xi)*Tr);
}
int main(int argc, char *argv[]){
  clock_t start = clock(); // MILO
  int w=256, h=256, samps = argc==2 ? atoi(argv[1])/4 : 25; // # samples
  Ray cam(XMVectorSet(50,52,295.6f,0), XMVector3Normalize(XMVectorSet(0,-0.042612f,-1,0))); // cam pos, dir
  XMVECTOR cx=XMVectorSet(w*.5135f/h,0,0,0), cy=(XMVector3Normalize(XMVector3Cross(cx,cam.d)))*.5135f, r, *c=(XMVECTOR*)_aligned_malloc(w*h*sizeof(XMVECTOR), 16);
  XMVECTOR invSamp = XMVectorReplicate(1.0f/samps);
#pragma omp parallel for schedule(dynamic, 1) private(r)       // OpenMP
  for (int y=0; y<h; y++){                       // Loop over image rows
    fprintf(stderr,"\rRendering (%d spp) %5.2f%%",samps*4,100.*y/(h-1));
	RandomLCG Xi(y*y*y); // MILO
    for (unsigned short x=0; x<w; x++)   // Loop cols
	  for (int sy=0, i=(h-y-1)*w+x; sy<2; sy++)     // 2x2 subpixel rows
        for (int sx=0; sx<2; sx++){        // 2x2 subpixel cols
		  r=XMVectorZero();
          for (int s=0; s<samps; s++){
            float r1=2*Xi(), dx=r1<1 ? sqrt(r1)-1: 1-sqrt(2-r1);
            float r2=2*Xi(), dy=r2<1 ? sqrt(r2)-1: 1-sqrt(2-r2);
            XMVECTOR d = cx*( ( (sx+.5f + dx)/2 + x)/w - .5f) +
                    cy*( ( (sy+.5f + dy)/2 + y)/h - .5f) + cam.d;
            r = r + radiance(Ray(cam.o+d*XMVectorReplicate(140),XMVector3Normalize(d)),0,Xi)*invSamp;
          } // Camera rays are pushed ^^^^^ forward to start in interior
          c[i] = c[i] + XMVectorSet(clamp(XMVectorGetX(r)),clamp(XMVectorGetY(r)),clamp(XMVectorGetZ(r)),0)*.25f;
        }
  }
  printf("\n%f sec\n", (float)(clock() - start)/CLOCKS_PER_SEC); // MILO
  FILE *refl = fopen("image.ppm", "w");         // Write image to PPM file.
  fprintf(refl, "P3\n%d %d\n%d\n", w, h, 255);
  for (int i=0; i<w*h; i++)
    fprintf(refl,"%d %d %d ", toInt(XMVectorGetX(c[i])), toInt(XMVectorGetY(c[i])), toInt(XMVectorGetZ(c[i])));
}