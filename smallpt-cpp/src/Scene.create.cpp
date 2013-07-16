#include "Scene.h"

namespace
{
    void createScene_cornell(Scene* pScene)
    {
        Geometry* geoms[] = {//Scene: radius, position, emission, color, material
            new Sphere(1e5,     Vec( 1e5+1,40.8,81.6),   Vec(),          Vec(.75,.25,.25),   DIFFUSE),// left wall 
            new Sphere(1e5,     Vec(-1e5+99,40.8,81.6),  Vec(),          Vec(.25,.25,.75),   DIFFUSE),// right wall
            new Sphere(1e5,     Vec(50,40.8, 1e5),       Vec(),          Vec(.75,.75,.75),   DIFFUSE),// back wall
            new Sphere(1e5,     Vec(50,40.8,-1e5+170),   Vec(),          Vec(),              DIFFUSE),// front wall
            new Sphere(1e5,     Vec(50, 1e5, 81.6),      Vec(),          Vec(.75,.75,.75),   DIFFUSE),// bottom wall
            new Sphere(1e5,     Vec(50,-1e5+81.6,81.6),  Vec(),          Vec(.75,.75,.75),   DIFFUSE),// top wall

            new Sphere(16.5,    Vec(27,16.5,47),         Vec(),          Vec(1,1,1)*.999,    SPECULAR),// mirror 1
            new Sphere(20.5,    Vec(73,16.5,78),         Vec(),          Vec(1,1,1)*.999,    SPECULAR),// mirror 2
            new Sphere(10.5,    Vec(27,16.5,100),        Vec(),          Vec(1,1,1)*.999,    SPECULAR),// mirror 3

            new Sphere(600,     Vec(50,681.6-.27,81.6),  Vec(12,12,12),  Vec(),              DIFFUSE), // the white light
        };
        pScene->addGeometries(geoms, _countof(geoms));
    }

    void createScene_sky(Scene* pScene)
    {
        // Idea stolen from Picogen http://picogen.org/ by phresnel/greenhybrid
        Vec Cen(50,40.8,-860);
        Geometry* geoms[] = {//Scene: radius, position, emission, color, material
            //center 50 40.8 62
            //floor 0
            //back  0

            new Sphere(1600, Vec(1,0,2)*3000, Vec(1,.9,.8)*1.2e1*1.56*2,Vec(), DIFFUSE), // sun
            new Sphere(1560, Vec(1,0,2)*3500,Vec(1,.5,.05)*4.8e1*1.56*2, Vec(),  DIFFUSE), // horizon sun2
            new Sphere(10000,Cen+Vec(0,0,-200), Vec(0.0627, 0.188, 0.569)*6e-2*8, Vec(.7,.7,1)*.25,  DIFFUSE), // sky
            new Sphere(10000,Cen+Vec(0,0,-200), Vec(0.00063842, 0.02001478, 0.28923243)*6e-2*8, Vec(.7,.7,1)*.25,  DIFFUSE), // sky

            new Sphere(100000, Vec(50, -100000, 0),  Vec(),Vec(.3,.3,.3),DIFFUSE), // grnd
            new Sphere(110000, Vec(50, -110048.5, 0),  Vec(.9,.5,.05)*4,Vec(),DIFFUSE),// horizon brightener
            new Sphere(4e4, Vec(50, -4e4-30, -3000),  Vec(),Vec(.2,.2,.2),DIFFUSE),// mountains
            new Sphere(3.99e4, Vec(50, -3.99e4+20.045, -3000),  Vec(),Vec(.7,.7,.7),DIFFUSE),// mountains snow

            new Sphere(26.5,Vec(22,26.5,42),   Vec(),Vec(1,1,1)*.596, SPECULAR), // white Mirr
            new Sphere(13,Vec(75,13,82),   Vec(),Vec(.96,.96,.96)*.96, REFRACT),// Glas
            new Sphere(22,Vec(87,22,24),   Vec(),Vec(.6,.6,.6)*.696, REFRACT)    // Glas2
        };
        pScene->addGeometries(geoms, _countof(geoms));
    }

    void createScene_nightsky(Scene* pScene)
    {
        Geometry* geoms[] = {//Scene: radius, position, emission, color, material
            //center 50 40.8 62
            //floor 0
            //back  0
            //rad       pos                   emis           col     refl
            new Sphere(1e3,   Vec(1,1,-2)*1e4,    Vec(1,1,1)*5e2,     Vec(), DIFFUSE), // moon
            new Sphere(3e2,   Vec(.6,.2,-2)*1e4,    Vec(1,1,1)*5e3,     Vec(), DIFFUSE), // moon

            new Sphere(2.5e3,   Vec(.82,.92,-2)*1e4,    Vec(1,1,1)*.8e2,     Vec(), DIFFUSE), // moon

            new Sphere(2.5e4, Vec(50, 0, 0),     Vec(1,1,1)*1e-3,    Vec(.2,.2,1)*0.0075, DIFFUSE), // sky
            new Sphere(2.5e4, Vec(50, 0, 0),  Vec(0.114, 0.133, 0.212)*1e-2,  Vec(.216,.384,1)*0.0007, DIFFUSE), // sky

            new Sphere(2.5e4, Vec(50, 0, 0),  Vec(0.114, 0.133, 0.212)*1e-2,  Vec(.216,.384,1)*0.003, DIFFUSE), // sky

            new Sphere(5e0,   Vec(-.2,0.16,-1)*1e4, Vec(1.00, 0.843, 0.698)*1e2,   Vec(), DIFFUSE),  // star
            new Sphere(5e0,   Vec(0,  0.18,-1)*1e4, Vec(1.00, 0.851, 0.710)*1e2,  Vec(), DIFFUSE),  // star
            new Sphere(5e0,   Vec(.3, 0.15,-1)*1e4, Vec(0.671, 0.780, 1.00)*1e2,   Vec(), DIFFUSE),  // star
            new Sphere(3.5e4,   Vec(600,-3.5e4+1, 300), Vec(),   Vec(.6,.8,1)*.01,  REFRACT),   //pool
            new Sphere(5e4,   Vec(-500,-5e4+0, 0),   Vec(),      Vec(1,1,1)*.35,  DIFFUSE),    //hill
            new Sphere(16.5,  Vec(27,0,47),         Vec(),              Vec(1,1,1)*.33, DIFFUSE), //hut
            new Sphere(7,     Vec(27+8*sqrt(2.0),0,47+8*sqrt(2.0)),Vec(),  Vec(1,1,1)*.33,  DIFFUSE), //door
            new Sphere(500,   Vec(-1e3,-300,-3e3), Vec(),  Vec(1,1,1)*.351,    DIFFUSE),  //mnt
            new Sphere(830,   Vec(0,   -500,-3e3), Vec(),  Vec(1,1,1)*.354,    DIFFUSE),  //mnt
            new Sphere(490,  Vec(1e3,  -300,-3e3), Vec(),  Vec(1,1,1)*.352,    DIFFUSE),  //mnt
        };
        pScene->addGeometries(geoms, _countof(geoms));
    }

    void createScene_island(Scene* pScene)
    {
        // Inspired by cover of "Time Planet Earth: An Illustrated History"
        Vec Cen(50,-20,-860);
        Geometry* geoms[] = {//Scene: radius, position, emission, color, material
            // center 50 40.8 62
            // floor 0
            // back  0
            //     rad       pos                   emis           col     refl

            new Sphere(160,  Cen+Vec(0, 600, -500),Vec(1,1,1)*2e2, Vec(),  DIFFUSE), // sun
            new Sphere(800, Cen+Vec(0,-880,-9120),Vec(1,1,1)*2e1, Vec(),  DIFFUSE), // horizon
            new Sphere(10000,Cen+Vec(0,0,-200), Vec(0.0627, 0.188, 0.569)*1e0, Vec(1,1,1)*.4,  DIFFUSE), // sky

            //  new Sphere(1000, Cen+Vec(0,-1080,-8020),Vec(1,1,1)*2e1, Vec(),  DIFFUSE), // horizon
            //  new Sphere(10000,Cen+Vec(0,0,-200), Vec(0.0627, 0.188, 0.569)*1e0, Vec(1,1,1)*.3,  DIFFUSE), // sky

            //  new Sphere(800, Cen+Vec(0,-720,-200),Vec(),  Vec(0, 0.588, 0.8),  REFRACT), // water
            //  new Sphere(800, Cen+Vec(0,-720,-200),Vec(),  Vec(0.106, 0.725, 0.949),  REFRACT), // water
            //  new Sphere(800, Cen+Vec(0,-720,-200),Vec(),  Vec(0.110, 0.988, 0.945),  REFRACT), // water
            new Sphere(800, Cen+Vec(0,-720,-200),Vec(),  Vec(0.110, 0.898, 1.00)*.996,  REFRACT), // water
            new Sphere(790, Cen+Vec(0,-720,-200),Vec(),  Vec(.4,.3,.04)*.6,    DIFFUSE), // earth
            new Sphere(325, Cen+Vec(0,-255,-50), Vec(),  Vec(.4,.3,.04)*.8,       DIFFUSE), // island
            new Sphere(275, Cen+Vec(0,-205,-33), Vec(),  Vec(.02,.3,.02)*.75,      DIFFUSE), // grass
        };
        pScene->addGeometries(geoms, _countof(geoms));
    }

    void createScene_vista(Scene* pScene)
    {
        Vec Cen(50,-20,-860);
        Geometry* geoms[] = {//Scene: radius, position, emission, color, material
            // center 50 40.8 62
            // floor 0
            // back  0
            //     rad       pos                   emis           col     refl

            new Sphere(8000, Cen+Vec(0,-8000,-900),Vec(1,.4,.1)*5e-1, Vec(),  DIFFUSE), // sun
            new Sphere(1e4,  Cen+Vec(), Vec(0.631, 0.753, 1.00)*3e-1, Vec(1,1,1)*.5,  DIFFUSE), // sky

            new Sphere(150,  Cen+Vec(-350,0, -100),Vec(),  Vec(1,1,1)*.3,  DIFFUSE), // mnt
            new Sphere(200,  Cen+Vec(-210,0,-100), Vec(),  Vec(1,1,1)*.3,  DIFFUSE), // mnt
            new Sphere(145,  Cen+Vec(-210,85,-100),Vec(),  Vec(1,1,1)*.8,  DIFFUSE), // snow
            new Sphere(150,  Cen+Vec(-50,0,-100),  Vec(),  Vec(1,1,1)*.3,  DIFFUSE), // mnt
            new Sphere(150,  Cen+Vec(100,0,-100),  Vec(),  Vec(1,1,1)*.3,  DIFFUSE), // mnt
            new Sphere(125,  Cen+Vec(250,0,-100),  Vec(),  Vec(1,1,1)*.3,  DIFFUSE), // mnt
            new Sphere(150,  Cen+Vec(375,0,-100),  Vec(),  Vec(1,1,1)*.3,  DIFFUSE), // mnt

            new Sphere(2500, Cen+Vec(0,-2400,-500),Vec(),  Vec(1,1,1)*.1,  DIFFUSE), // mnt base

            new Sphere(8000, Cen+Vec(0,-8000,200), Vec(),  Vec(.2,.2,1),    REFRACT), // water
            new Sphere(8000, Cen+Vec(0,-8000,1100),Vec(),  Vec(0,.3,0),     DIFFUSE), // grass
            new Sphere(8   , Cen+Vec(-75, -5, 850),Vec(),  Vec(0,.3,0),     DIFFUSE), // bush
            new Sphere(30,   Cen+Vec(0,   23, 825),Vec(),  Vec(1,1,1)*.996, REFRACT), // ball

            new Sphere(30,  Cen+Vec(200,280,-400),  Vec(),  Vec(1,1,1)*.8,  DIFFUSE),   // clouds
            new Sphere(37,  Cen+Vec(237,280,-400),  Vec(),  Vec(1,1,1)*.8,  DIFFUSE),   // clouds
            new Sphere(28,  Cen+Vec(267,280,-400),  Vec(),  Vec(1,1,1)*.8,  DIFFUSE),   // clouds

            new Sphere(40,  Cen+Vec(150,280,-1000),  Vec(),  Vec(1,1,1)*.8,  DIFFUSE),  // clouds
            new Sphere(37,  Cen+Vec(187,280,-1000),  Vec(),  Vec(1,1,1)*.8,  DIFFUSE),  // clouds

            new Sphere(40,  Cen+Vec(600,280,-1100),  Vec(),  Vec(1,1,1)*.8,  DIFFUSE),  // clouds
            new Sphere(37,  Cen+Vec(637,280,-1100),  Vec(),  Vec(1,1,1)*.8,  DIFFUSE),  // clouds

            new Sphere(37,  Cen+Vec(-800,280,-1400),  Vec(),  Vec(1,1,1)*.8,  DIFFUSE), // clouds
            new Sphere(37,  Cen+Vec(0,280,-1600),  Vec(),  Vec(1,1,1)*.8,  DIFFUSE),    // clouds
            new Sphere(37,  Cen+Vec(537,280,-1800),  Vec(),  Vec(1,1,1)*.8,  DIFFUSE),  // clouds

        };
        pScene->addGeometries(geoms, _countof(geoms));
    }

    void createScene_overlap(Scene* pScene)
    {
        double D=50;
        double R=40;
        Geometry* geoms[] = {//Scene: radius, position, emission, color, material
            new Sphere(150, Vec(50+75,28,62), Vec(1,1,1)*0e-3, Vec(1,.9,.8)*.93, REFRACT),
            new Sphere(28,  Vec(50+5,-28,62), Vec(1,1,1)*1e1, Vec(1,1,1)*0, DIFFUSE),
            new Sphere(300, Vec(50,28,62), Vec(1,1,1)*0e-3, Vec(1,1,1)*.93, SPECULAR)
        };

        pScene->addGeometries(geoms, _countof(geoms));
    }

    void createScene_wada(Scene* pScene)
    {
        double R=60;
        //double R=120;
        double T=30*M_PI/180.;
        double D=R/cos(T);
        double Z=60;
        Geometry* geoms[] = {//Scene: radius, position, emission, color, material
            // center 50 40.8 62
            // floor 0
            // back  0
            new Sphere(1e5, Vec(50, 100, 0),      Vec(1,1,1)*3e0, Vec(), DIFFUSE), // sky
            new Sphere(1e5, Vec(50, -1e5-D-R, 0), Vec(),     Vec(.1,.1,.1),DIFFUSE),           //grnd

            new Sphere(R, Vec(50,40.8,62)+Vec( cos(T),sin(T),0)*D, Vec(), Vec(1,.3,.3)*.999, SPECULAR), //red
            new Sphere(R, Vec(50,40.8,62)+Vec(-cos(T),sin(T),0)*D, Vec(), Vec(.3,1,.3)*.999, SPECULAR), //grn
            new Sphere(R, Vec(50,40.8,62)+Vec(0,-1,0)*D,         Vec(), Vec(.3,.3,1)*.999, SPECULAR), //blue
            new Sphere(R, Vec(50,40.8,62)+Vec(0,0,-1)*D,       Vec(), Vec(.53,.53,.53)*.999, SPECULAR), //back
            new Sphere(R, Vec(50,40.8,62)+Vec(0,0,1)*D,      Vec(), Vec(1,1,1)*.999, REFRACT), //front

            //   new Sphere(R, Vec(50,35,Z)+Vec( cos(T),sin(T),0)*D, Vec(1,1,1)*1e-1, Vec(1,1,1)*.999, SPECULAR), //red
            //   new Sphere(R, Vec(50,35,Z)+Vec(-cos(T),sin(T),0)*D, Vec(1,1,1)*1e-1, Vec(1,1,1)*.999, SPECULAR), //grn
            //   new Sphere(R, Vec(50,35,Z)+Vec(0,-1,0)*D,           Vec(1,1,1)*1e-1, Vec(1,1,1)*.999, SPECULAR), //blue
            //   new Sphere(R, Vec(50,35,Z)+Vec(0,0,-1)*D*1.6,       Vec(1,1,1)*0e-1, Vec(0.275, 0.612, 0.949)*.999, SPECULAR), //back
            //  new Sphere(R, Vec(50,40.8,62)+Vec(0,0,1)*D*.2877,          Vec(1,1,1)*0e-1, Vec(1,1,1)*.999, REFRACT), //front

        };
        pScene->addGeometries(geoms, _countof(geoms));
    }

    void createScene_wada2(Scene* pScene)
    {
        //double R=60;
        double R=120;     // radius
        double T=30*M_PI/180.;
        double D=R/cos(T);     //distance
        // double D=60;     //distance
        // double R=D*sqrt(2);
        double Z=62;
        Vec C=Vec(0.275, 0.612, 0.949);
        Geometry* geoms[] = {//Scene: radius, position, emission, color, material

            new Sphere(R, Vec(50,28,Z)+Vec( cos(T),sin(T),0)*D,    C*6e-2,Vec(1,1,1)*.996, SPECULAR), //red
            new Sphere(R, Vec(50,28,Z)+Vec(-cos(T),sin(T),0)*D,    C*6e-2,Vec(1,1,1)*.996, SPECULAR), //grn
            new Sphere(R, Vec(50,28,Z)+Vec(0,-1,0)*D,              C*6e-2,Vec(1,1,1)*.996, SPECULAR), //blue
            new Sphere(R, Vec(50,28,Z)+Vec(0,0,-1)*R*2*sqrt(2./3.),C*0e-2,Vec(1,1,1)*.996, SPECULAR), //back
            //  new Sphere(1e5, Vec(50,28,Z)+Vec(0,0,1e5+170),   Vec(1,1,1)*0,Vec(1,1,1)*.996, SPECULAR), //front
            //  new Sphere(2*R*2*sqrt(2./3.)-R*2*sqrt(2./3.)/3., Vec(50,28,Z)+Vec(0,0,-R*2*sqrt(2./3.)/3.),   Vec(1,1,1)*0,Vec(1,1,1)*.3333, SPECULAR), //front
            new Sphere(2*2*R*2*sqrt(2./3.)-R*2*sqrt(2./3.)/3., Vec(50,28,Z)+Vec(0,0,-R*2*sqrt(2./3.)/3.),   Vec(1,1,1)*0,Vec(1,1,1)*.5, SPECULAR), //front
        };

        pScene->addGeometries(geoms, _countof(geoms));
    }

    void createScene_forest(Scene* pScene)
    {
        Vec tc(0.0588, 0.361, 0.0941);
        Vec sc = Vec(1,1,1)*.7;
        Geometry* geoms[] = {//Scene: radius, position, emission, color, material
            // center 50 40.8 62
            // floor 0
            // back  0
            //  new Sphere(1e5, Vec(50, 1e5+100, 0),  Vec(1,1,1)*1,Vec(),DIFFUSE), //lite
            //  new Sphere(1e5, Vec(50, -1e5, 0),  Vec(),Vec(.3,.3,.1),DIFFUSE), //grnd
            //  new Sphere(1e5, Vec(50, 1e5+100, 0),  Vec(0.761, 0.875, 1.00)*1.3,Vec(),DIFFUSE),
            //  //lite
            new Sphere(1e5, Vec(50, 1e5+130, 0),  Vec(1,1,1)*1.3,Vec(),DIFFUSE), //lite
            new Sphere(1e2, Vec(50, -1e2+2, 47),  Vec(),Vec(1,1,1)*.7,DIFFUSE), //grnd

            new Sphere(1e4, Vec(50, -30, 300)+Vec(-sin(50*M_PI/180),0,cos(50*M_PI/180))*1e4, Vec(), Vec(1,1,1)*.99,SPECULAR),// mirr L
            new Sphere(1e4, Vec(50, -30, 300)+Vec(sin(50*M_PI/180),0,cos(50*M_PI/180))*1e4,  Vec(), Vec(1,1,1)*.99,SPECULAR),// mirr R
            new Sphere(1e4, Vec(50, -30, -50)+Vec(-sin(30*M_PI/180),0,-cos(30*M_PI/180))*1e4,Vec(), Vec(1,1,1)*.99,SPECULAR),// mirr FL
            new Sphere(1e4, Vec(50, -30, -50)+Vec(sin(30*M_PI/180),0,-cos(30*M_PI/180))*1e4, Vec(), Vec(1,1,1)*.99,SPECULAR),// mirr


            new Sphere(4, Vec(50,6*.6,47),   Vec(),Vec(.13,.066,.033), DIFFUSE),//"tree"
            new Sphere(16,Vec(50,6*2+16*.6,47),   Vec(), tc,  DIFFUSE),//"tree"
            new Sphere(11,Vec(50,6*2+16*.6*2+11*.6,47),   Vec(), tc,  DIFFUSE),//"tree"
            new Sphere(7, Vec(50,6*2+16*.6*2+11*.6*2+7*.6,47),   Vec(), tc,  DIFFUSE),//"tree"

            new Sphere(15.5,Vec(50,1.8+6*2+16*.6,47),   Vec(), sc,  DIFFUSE),//"tree"
            new Sphere(10.5,Vec(50,1.8+6*2+16*.6*2+11*.6,47),   Vec(), sc,  DIFFUSE),//"tree"
            new Sphere(6.5, Vec(50,1.8+6*2+16*.6*2+11*.6*2+7*.6,47),   Vec(), sc,  DIFFUSE),//"tree"
        };

        pScene->addGeometries(geoms, _countof(geoms));
    }
}

#define SCENE_ENTRY_LIST \
    ONE_ENTRY(cornell);\
    ONE_ENTRY(sky);\
    ONE_ENTRY(nightsky);\
    ONE_ENTRY(island);\
    ONE_ENTRY(vista);\
    ONE_ENTRY(overlap);\
    ONE_ENTRY(wada);\
    ONE_ENTRY(wada2);\
    ONE_ENTRY(forest);

bool Scene::createBuiltInScene( Scene* pScene, const std::string& sceneName )
{
#define ONE_ENTRY(tag) if (sceneName == #tag) { createScene_##tag(pScene); return true; }
    SCENE_ENTRY_LIST;
#undef ONE_ENTRY
    return false;
}

void Scene::printBuiltInSceneNames()
{
    fprintf(stdout, "Built-in scene names: \n");
#define ONE_ENTRY(tag) fprintf(stdout, "\t%s\n", #tag);
    SCENE_ENTRY_LIST;
#undef ONE_ENTRY
}
