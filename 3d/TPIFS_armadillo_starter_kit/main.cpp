#include <stdlib.h>
#include <cmath>
#include <GL/glut.h>
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <ctime>
#include <cstdlib> 
#include <chrono>
#include <fstream>

#include "forme.h"
#include "config.h"

using namespace std;



struct IndiceCube{
    bool x;
    bool y;
    bool z;
    string toString(){
        return to_string(x)+to_string(y)+to_string(z);
    }
    int toInt(){
        int v = 0;
        if(x)
            v+=4;
        if(y)
            v+=2;
        if(z)
            v++;
        return v;
    }
    IndiceCube add(IndiceCube b)
    {
        IndiceCube c = {x,y,z};
        if(b.x)
            c.x=!x;
        if(b.y)
            c.y=!y;
        if(b.z)
            c.z=!z;
        return c;
    }
};

struct EdgeCube{
    IndiceCube P0;
    IndiceCube P1;
    IndiceCube coupe;
    Point inter;
    string toString(){
        string ret = P0.toString()+":"+P1.toString()+":";
        if(coupe.x)
            ret+="x";
        else
        {
            if(coupe.y)
                ret+="y";
            else
                ret+="z";
        }
        return ret;
    }
    void order()
    {
        if(P0.toInt()>P1.toInt())
        {
            IndiceCube t = P0;
            P0 = P1;
            P1 = t;
        }
    }
};

void affichage(void);

void clavier(unsigned char touche,int x,int y);
void affiche_repere(void);
void dessinZone(void);

void mouse(int, int, int, int);
void mouseMotion(int, int);
//void reshape(int,int);

double val(Point P);
double func(Point P);
double dist(Point P);
double dist2(Point P1, Point P2);
Point hsv2rgb(double h, double s, double v);
void initMap();

map<string,Forme> guide;
map<string,bool> marque;
vector<Point> pile;
map<string,EdgeCube> CubePlane;

Point A = {-1,-1,-1};
Point B = {-1,-1,1};
Point C = {-1,1,1};
Point D = {-1,1,-1};
Point E = {1,-1,-1};
Point F = {1,-1,1};
Point G = {1,1,1};
Point H = {1,1,-1};

IndiceCube c000 = {0,0,0};
IndiceCube c001 = {0,0,1};
IndiceCube c010 = {0,1,0};
IndiceCube c011 = {0,1,1};
IndiceCube c100 = {1,0,0};
IndiceCube c101 = {1,0,1};
IndiceCube c110 = {1,1,0};
IndiceCube c111 = {1,1,1};

IndiceCube cX = {1,0,0};
IndiceCube cY = {0,1,0};
IndiceCube cZ = {0,0,1};

double va = -1;
double vb = -1;
double vc = -1;
double vd = -1;
double ve = -1;
double vf = -1;
double vg = -1;
double vh = -1;

int nb = 26;


// variables globales pour OpenGL
bool mouseLeftDown;
bool mouseRightDown;
bool mouseMiddleDown;
float mouseX, mouseY;
float cameraAngleX;
float cameraAngleY;
float cameraDistance=0.;

//------------------------------------------------------
int main(int argc,char **argv)
{
    initMap();
    srand(time(0));
  /* initialisation de glut et creation
     de la fenetre */
  glutInit(&argc,argv);
  glutInitDisplayMode(GLUT_RGB);
  glutInitWindowPosition(200,200);
  glutInitWindowSize(600,600);
  glutCreateWindow("ifs");

  /* Initialisation d'OpenGL */
  glClearColor(0.0,0.0,0.0,0.0);
  glColor3f(1.0,1.0,1.0);
  glPointSize(1.0);
	
  /* enregistrement des fonctions de rappel */
  glutDisplayFunc(affichage);
  glutKeyboardFunc(clavier);
  glutMouseFunc(mouse);
  glutMotionFunc(mouseMotion);

// initialisation des tranformations

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	 gluPerspective(60.0f,(GLfloat)200/(GLfloat)200,0.1f,100.0f);
	glMatrixMode(GL_MODELVIEW);
	//glLoadIdentity();	
//	glScalef(.7,.7,.7);
 gluLookAt(0.,0.,2., 0.,0.,0., 0.,1.,0.);
  //  glTranslatef(0.0,0.0,-5.0);

/* Entree dans la boucle principale glut */
  glutMainLoop();
  return 0;
}
//------------------------------------------------------
void affiche_repere(void)
{
  glBegin(GL_LINES);
  glColor3f(1.0,0.0,0.0);
  glVertex2f(0.,0.);
  glVertex2f(1.,0.);
  glEnd(); 

	 glBegin(GL_LINES);
  glColor3f(0.0,1.0,0.0);
  glVertex2f(0.,0.);
  glVertex2f(0.,1.);
  glEnd(); 
   glBegin(GL_LINES);
  glColor3f(0.0,0.0,1.0);
  glVertex3f(0.,0.,0.);
  glVertex3f(0.,0.,1.);
  glEnd(); 
}

void initMap()
{
    guide["00000000"] = Forme(0,-1,-1,-1,-1);
    guide["00000001"] = Forme(1,7,-1,-1,-1);
    guide["00000010"] = Forme(1,6,-1,-1,-1);
    guide["00000011"] = Forme(2,6,7,-1,-1);
    guide["00000100"] = Forme(1,5,-1,-1,-1);
    guide["00000101"] = Forme(2,5,7,-1,-1);
    guide["00000110"] = Forme(3,5,6,-1,-1);
    guide["00000111"] = Forme(4,5,6,7,-1);
    guide["00001000"] = Forme(1,4,-1,-1,-1);
    guide["00001001"] = Forme(3,4,7,-1,-1);
    guide["00001010"] = Forme(2,4,6,-1,-1);
    guide["00001011"] = Forme(4,4,6,7,-1);
    guide["00001100"] = Forme(2,4,5,-1,-1);
    guide["00001101"] = Forme(4,4,5,7,-1);
    guide["00001110"] = Forme(4,4,5,6,-1);
    guide["00001111"] = Forme(5,4,5,6,7);
    guide["00010000"] = Forme(1,3,-1,-1,-1);
    guide["00010001"] = Forme(2,3,7,-1,-1);
    guide["00010010"] = Forme(3,3,6,-1,-1);
    guide["00010011"] = Forme(4,3,6,7,-1);
    guide["00010100"] = Forme(3,3,5,-1,-1);
    guide["00010101"] = Forme(4,3,5,7,-1);
    guide["00010110"] = Forme(12,3,5,6,-1);
    guide["00010111"] = Forme(8,3,5,6,7);
    guide["00011000"] = Forme(10,3,4,-1,-1);
    guide["00011001"] = Forme(11,4,3,7,-1);
    guide["00011010"] = Forme(11,3,4,6,-1);
    guide["00011011"] = Forme(9,3,4,6,7);
    guide["00011100"] = Forme(11,3,4,5,-1);
    guide["00011101"] = Forme(9,3,4,5,7);
    guide["00011110"] = Forme(6,3,4,5,6);
    guide["00011111"] = Forme(4,0,1,2,-1);
    guide["00100000"] = Forme(1,2,-1,-1,-1);
    guide["00100001"] = Forme(3,2,7,-1,-1);
    guide["00100010"] = Forme(2,2,6,-1,-1);
    guide["00100011"] = Forme(4,2,6,7,-1);
    guide["00100100"] = Forme(10,2,5,-1,-1);
    guide["00100101"] = Forme(11,2,5,7,-1);
    guide["00100110"] = Forme(11,5,2,6,-1);
    guide["00100111"] = Forme(9,2,5,6,7);
    guide["00101000"] = Forme(3,2,4,-1,-1);
    guide["00101001"] = Forme(12,2,4,7,-1);
    guide["00101010"] = Forme(4,2,4,6,-1);
    guide["00101011"] = Forme(8,2,4,6,7);
    guide["00101100"] = Forme(11,2,4,5,-1);
    guide["00101101"] = Forme(6,2,4,5,7);
    guide["00101110"] = Forme(9,2,4,5,6);
    guide["00101111"] = Forme(4,0,1,3,-1);
    guide["00110000"] = Forme(2,2,3,-1,-1);
    guide["00110001"] = Forme(4,2,3,7,-1);
    guide["00110010"] = Forme(4,2,3,6,-1);
    guide["00110011"] = Forme(5,2,3,6,7);
    guide["00110100"] = Forme(11,5,2,3,-1);
    guide["00110101"] = Forme(9,2,3,5,7);
    guide["00110110"] = Forme(6,6,2,3,5);
    guide["00110111"] = Forme(4,0,1,4,-1);
    guide["00111000"] = Forme(11,4,2,3,-1);
    guide["00111001"] = Forme(6,4,2,3,7);
    guide["00111010"] = Forme(9,2,3,4,6);
    guide["00111011"] = Forme(4,0,1,5,-1);
    guide["00111100"] = Forme(13,2,3,4,5);
    guide["00111101"] = Forme(11,6,0,1,-1);
    guide["00111110"] = Forme(11,7,0,1,-1);
    guide["00111111"] = Forme(2,0,1,-1,-1);
    guide["01000000"] = Forme(1,1,-1,-1,-1);
    guide["01000001"] = Forme(3,1,7,-1,-1);
    guide["01000010"] = Forme(10,1,6,-1,-1);
    guide["01000011"] = Forme(11,1,6,7,-1);
    guide["01000100"] = Forme(2,1,5,-1,-1);
    guide["01000101"] = Forme(4,1,5,7,-1);
    guide["01000110"] = Forme(11,6,1,5,-1);
    guide["01000111"] = Forme(9,1,5,6,7);
    guide["01001000"] = Forme(3,1,4,-1,-1);
    guide["01001001"] = Forme(12,1,4,7,-1);
    guide["01001010"] = Forme(11,1,4,6,-1);
    guide["01001011"] = Forme(6,1,4,6,7);
    guide["01001100"] = Forme(4,1,4,5,-1);
    guide["01001101"] = Forme(8,1,4,5,7);
    guide["01001110"] = Forme(9,1,4,5,6);
    guide["01001111"] = Forme(4,0,2,3,-1);
    guide["01010000"] = Forme(2,1,3,-1,-1);
    guide["01010001"] = Forme(4,1,3,7,-1);
    guide["01010010"] = Forme(11,6,1,3,-1);
    guide["01010011"] = Forme(9,1,3,6,7);
    guide["01010100"] = Forme(4,1,3,5,-1);
    guide["01010101"] = Forme(5,1,3,5,7);
    guide["01010110"] = Forme(9,1,3,5,6);
    guide["01010111"] = Forme(4,0,2,4,-1);
    guide["01011000"] = Forme(11,4,1,3,-1);
    guide["01011001"] = Forme(6,4,1,3,7);
    guide["01011010"] = Forme(13,1,3,4,6);
    guide["01011011"] = Forme(11,5,0,2,-1);
    guide["01011100"] = Forme(9,1,3,4,5);
    guide["01011101"] = Forme(4,0,2,6,-1);
    guide["01011110"] = Forme(11,7,0,2,-1);
    guide["01011111"] = Forme(2,0,2,-1,-1);
    guide["01100000"] = Forme(3,1,2,-1,-1);
    guide["01100001"] = Forme(12,1,2,7,-1);
    guide["01100010"] = Forme(11,1,2,6,-1);
    guide["01100011"] = Forme(6,1,2,6,7);
    guide["01100100"] = Forme(11,2,1,5,-1);
    guide["01100101"] = Forme(6,2,1,5,7);
    guide["01100110"] = Forme(13,1,5,2,6);
    guide["01100111"] = Forme(11,3,0,4,-1);
    guide["01101000"] = Forme(12,1,2,4,-1);
    guide["01101001"] = Forme(7,1,2,4,7);
    guide["01101010"] = Forme(6,1,2,4,6);
    guide["01101011"] = Forme(12,0,3,5,-1);
    guide["01101100"] = Forme(6,2,1,4,5);
    guide["01101101"] = Forme(12,0,3,6,-1);
    guide["01101110"] = Forme(11,0,3,7,-1);
    guide["01101111"] = Forme(3,0,3,-1,-1);
    guide["01110000"] = Forme(4,1,2,3,-1);
    guide["01110001"] = Forme(8,1,2,3,7);
    guide["01110010"] = Forme(9,1,2,3,6);
    guide["01110011"] = Forme(4,0,4,5,-1);
    guide["01110100"] = Forme(9,1,2,3,5);
    guide["01110101"] = Forme(4,0,4,6,-1);
    guide["01110110"] = Forme(11,7,0,4,-1);
    guide["01110111"] = Forme(2,0,4,-1,-1);
    guide["01111000"] = Forme(6,4,1,2,3);
    guide["01111001"] = Forme(12,0,5,6,-1);
    guide["01111010"] = Forme(11,0,5,7,-1);
    guide["01111011"] = Forme(3,0,5,-1,-1);
    guide["01111100"] = Forme(11,0,6,7,-1);
    guide["01111101"] = Forme(3,0,6,-1,-1);
    guide["01111110"] = Forme(10,0,7,-1,-1);
    guide["01111111"] = Forme(1,0,-1,-1,-1);
}

double interpolation(double a, double b)
{
    double aa = abs(a);
    double ab = abs(b);
    if(aa+ab!=0)
        return ab/(aa+ab);
    return 0.5;
}

double dist2(Point P1, Point P2)
{
    double x = P2.x - P1.x;
    double y = P2.y - P1.y;
    double z = P2.z - P1.z;
    return sqrt(x*x+y*y+z*z);
}

double dist(Point P)
{
    return dist2(P,{0,0,0});
}

Point hsv2rgb(double h, double s, double v)
{
    Point C = {0,0,0};
    h = 360*h;
    while(h>=360)
    {
        h-=360;
    }
    while(h<0)
    {
        h+=360;
    }
    double ti = floor(h/60.0);
    double f = h/60.0 - ti;
    double l = v * ( 1 - s);
    double m = v * (1 - f * s);
    double n = v * (1 - (1-f) * s);
    switch((int)ti)
    {
    case 0:
        C = {v,n,l};
        break;
    case 1:
        C = {m,v,l};
        break;
    case 2:
        C = {l,v,n};
        break;
    case 3:
        C = {l,m,v};
        break;
    case 4:
        C = {n,l,v};
        break;
    case 5:
        C = {v,l,m};
        break;
    }
    return C;
}
//-----------------------------------------------------
double val(Point P)
{
    if(va<0 || vb<0 || vc<0 || vd<0 || ve<0 || vf<0 || vg<0 || vh<0)
    {
        va = (double)(rand()) / (double)(RAND_MAX);
        vb = (double)(rand()) / (double)(RAND_MAX);
        vc = (double)(rand()) / (double)(RAND_MAX);
        vd = (double)(rand()) / (double)(RAND_MAX);
        ve = (double)(rand()) / (double)(RAND_MAX);
        vf = (double)(rand()) / (double)(RAND_MAX);
        vg = (double)(rand()) / (double)(RAND_MAX);
        vh = (double)(rand()) / (double)(RAND_MAX);
    }
    double xd = (P.x - A.x) / (E.x - A.x);
    double yd = (P.y - A.y) / (D.y - A.y);
    double zd = (P.z - A.z) / (B.z - A.z);
    double v00 = xd*ve + (1-xd)*va;
    double v01 = xd*vf + (1-xd)*vb;
    double v11 = xd*vg + (1-xd)*vc;
    double v10 = xd*vh + (1-xd)*vd;
    double v0 = yd*v10 + (1-yd)*v00;
    double v1 = yd*v11 + (1-yd)*v01;
    double v = zd*v1 + (1-zd)*v0;
    return v;
}

double func(Point P)
{
    double v = val(P);
    if(abs(P.x)>1 || abs(P.y)>1 || abs(P.z)>1)
        return -1;
    return v;
}

void dessinZone()
{
    glPointSize(5.0f);
    glBegin(GL_POINTS);
    Point Co = hsv2rgb(va,1,1);
    glColor3f(Co.x,Co.y,Co.z);
    glVertex3f(A.x,A.y,A.z);
    Co = hsv2rgb(vb,1,1);
    glColor3f(Co.x,Co.y,Co.z);
    glVertex3f(B.x,B.y,B.z);
    Co = hsv2rgb(vc,1,1);
    glColor3f(Co.x,Co.y,Co.z);
    glVertex3f(C.x,C.y,C.z);
    Co = hsv2rgb(vd,1,1);
    glColor3f(Co.x,Co.y,Co.z);
    glVertex3f(D.x,D.y,D.z);
    Co = hsv2rgb(ve,1,1);
    glColor3f(Co.x,Co.y,Co.z);
    glVertex3f(E.x,E.y,E.z);
    Co = hsv2rgb(vf,1,1);
    glColor3f(Co.x,Co.y,Co.z);
    glVertex3f(F.x,F.y,F.z);
    Co = hsv2rgb(vg,1,1);
    glColor3f(Co.x,Co.y,Co.z);
    glVertex3f(G.x,G.y,G.z);
    Co = hsv2rgb(vh,1,1);
    glColor3f(Co.x,Co.y,Co.z);
    glVertex3f(H.x,H.y,H.z);
    glEnd();
    glPointSize(1.0f);
    double nb = 1;
    glBegin(GL_POINTS);
    for(int i=-nb; i<=nb; i++)
    {
        double x = (double)i/nb;
        for(int j=-nb; j<=nb; j++)
        {
            double y = (double)j/nb;
            for(int k=-nb; k<=nb; k++)
            {
                double z = (double)k/nb;
                Point P = {x,y,z};
                double v = func(P);
                Co = hsv2rgb(v,1,1);
                if(v==-1)
                    Co = {0,0,0};
                glColor3f(Co.x,Co.y,Co.z);
                glVertex3f(x,y,z);
            }
        }
    }
    glEnd();
}

double distFromSurface(Point P0, Point P1, double target)
{
    double s = 0;
    for(int i=0; i<=10; i++)
    {
        double t = (double)i/10.0;
        Point P = {0,0,0};
        P.x = t * P1.x + (1-t) * P0.x;
        P.y = t * P1.y + (1-t) * P0.y;
        P.z = t * P1.z + (1-t) * P0.z;
        s+=abs(target-val(P));
    }
    return s;
}

void marchingCube(Point origine, double pas, double target)
{
    marque.clear();
    pile.push_back({0,0,0});
    int st = 0;
    while(pile.size()!=0)
    {
        cout<<"Pile "<<pile.size()<<endl;
        Point Pi = pile[0];
        string k = to_string(Pi.x)+":"+to_string(Pi.y)+":"+to_string(Pi.z);
        if(marque.find(k) == marque.end())
        {
            marque[k]=true;
            Point P = {origine.x+pas*Pi.x, origine.y+pas*Pi.y, origine.z+pas*Pi.z};
            Point P000 = {P.x-pas/2, P.y-pas/2, P.z-pas/2};
            Point P001 = {P.x-pas/2, P.y-pas/2, P.z+pas/2};
            Point P011 = {P.x-pas/2, P.y+pas/2, P.z+pas/2};
            Point P010 = {P.x-pas/2, P.y+pas/2, P.z-pas/2};
            Point P100 = {P.x+pas/2, P.y-pas/2, P.z-pas/2};
            Point P101 = {P.x+pas/2, P.y-pas/2, P.z+pas/2};
            Point P111 = {P.x+pas/2, P.y+pas/2, P.z+pas/2};
            Point P110 = {P.x+pas/2, P.y+pas/2, P.z-pas/2};

            double v000 = target - func(P000);
            double v001 = target - func(P001);
            double v011 = target - func(P011);
            double v010 = target - func(P010);
            double v100 = target - func(P100);
            double v101 = target - func(P101);
            double v111 = target - func(P111);
            double v110 = target - func(P110);
            {
                if(v000 == 0)
                {
                    if(v001>0 || v010>0 || v100>0)
                        v000 = 0.00001;
                }
                if(v001 == 0)
                {
                    if(v000>0 || v011>0 || v101>0)
                        v001 = 0.00001;
                }
                if(v011 == 0)
                {
                    if(v010>0 || v001>0 || v111>0)
                        v011 = 0.00001;
                }
                if(v010 == 0)
                {
                    if(v011>0 || v000>0 || v110>0)
                        v010 = 0.00001;
                }
                if(v100 == 0)
                {
                    if(v101>0 || v110>0 || v000>0)
                        v100 = 0.00001;
                }
                if(v101 == 0)
                {
                    if(v100>0 || v111>0 || v001>0)
                        v101 = 0.00001;
                }
                if(v111 == 0)
                {
                    if(v110>0 || v101>0 || v011>0)
                        v111 = 0.00001;
                }
                if(v110 == 0)
                {
                    if(v111>0 || v100>0 || v010>0)
                        v110 = 0.00001;
                }
            }
            double lim = target+1;
            bool stop = false;
            {
                if(v000==lim)
                {
                    double sum = 0;
                    int nb = 0;
                    if(v001!=lim)
                    {
                        sum+=v001;
                        nb++;
                    }
                    if(v010!=lim)
                    {
                        sum+=v010;
                        nb++;
                    }
                    if(v100!=lim)
                    {
                        sum+=v100;
                        nb++;
                    }
                    if(nb<2)
                    {
                        stop = true;
                    }
                    else
                        v000 = sum/(double)nb;
                }
                if(v001==lim)
                {
                    double sum = 0;
                    int nb = 0;
                    if(v000!=lim)
                    {
                        sum+=v000;
                        nb++;
                    }
                    if(v011!=lim)
                    {
                        sum+=v011;
                        nb++;
                    }
                    if(v101!=lim)
                    {
                        sum+=v101;
                        nb++;
                    }
                    if(nb<2)
                    {
                        stop = true;
                    }
                    else
                        v001 = sum/(double)nb;
                }
                if(v011==lim)
                {
                    double sum = 0;
                    int nb = 0;
                    if(v010!=lim)
                    {
                        sum+=v010;
                        nb++;
                    }
                    if(v001!=lim)
                    {
                        sum+=v001;
                        nb++;
                    }
                    if(v111!=lim)
                    {
                        sum+=v111;
                        nb++;
                    }
                    if(nb<2)
                    {
                        stop = true;
                    }
                    else
                        v011 = sum/(double)nb;
                }
                if(v010==lim)
                {
                    double sum = 0;
                    int nb = 0;
                    if(v011!=lim)
                    {
                        sum+=v011;
                        nb++;
                    }
                    if(v000!=lim)
                    {
                        sum+=v000;
                        nb++;
                    }
                    if(v110!=lim)
                    {
                        sum+=v110;
                        nb++;
                    }
                    if(nb<2)
                    {
                        stop = true;
                    }
                    else
                        v010 = sum/(double)nb;
                }
                if(v100==lim)
                {
                    double sum = 0;
                    int nb = 0;
                    if(v101!=lim)
                    {
                        sum+=v101;
                        nb++;
                    }
                    if(v110!=lim)
                    {
                        sum+=v110;
                        nb++;
                    }
                    if(v000!=lim)
                    {
                        sum+=v000;
                        nb++;
                    }
                    if(nb<2)
                    {
                        stop = true;
                    }
                    else
                        v100 = sum/(double)nb;
                }
                if(v101==lim)
                {
                    double sum = 0;
                    int nb = 0;
                    if(v100!=lim)
                    {
                        sum+=v100;
                        nb++;
                    }
                    if(v111!=lim)
                    {
                        sum+=v111;
                        nb++;
                    }
                    if(v001!=lim)
                    {
                        sum+=v001;
                        nb++;
                    }
                    if(nb<2)
                    {
                        stop = true;
                    }
                    else
                        v101 = sum/(double)nb;
                }
                if(v111==lim)
                {
                    double sum = 0;
                    int nb = 0;
                    if(v110!=lim)
                    {
                        sum+=v110;
                        nb++;
                    }
                    if(v101!=lim)
                    {
                        sum+=v101;
                        nb++;
                    }
                    if(v011!=lim)
                    {
                        sum+=v011;
                        nb++;
                    }
                    if(nb<2)
                    {
                        stop = true;
                    }
                    else
                        v111 = sum/(double)nb;
                }
                if(v110==lim)
                {
                    double sum = 0;
                    int nb = 0;
                    if(v111!=lim)
                    {
                        sum+=v111;
                        nb++;
                    }
                    if(v100!=lim)
                    {
                        sum+=v100;
                        nb++;
                    }
                    if(v010!=lim)
                    {
                        sum+=v010;
                        nb++;
                    }
                    if(nb<2)
                    {
                        stop = true;
                    }
                    else
                        v110 = sum/(double)nb;
                }
            }
            if(!stop)
            {
                float c = (double)st/100.0;
                Point C = hsv2rgb(c,1,1);
                glColor3f(C.x,C.y,C.z);
                /*
                glBegin(GL_LINES);
                    glVertex3f(P000.x,P000.y,P000.z);
                    glVertex3f(P001.x,P001.y,P001.z);
                    glVertex3f(P000.x,P000.y,P000.z);
                    glVertex3f(P100.x,P100.y,P100.z);
                    glVertex3f(P000.x,P000.y,P000.z);
                    glVertex3f(P010.x,P010.y,P010.z);
                    glVertex3f(P001.x,P001.y,P001.z);
                    glVertex3f(P011.x,P011.y,P011.z);
                    glVertex3f(P001.x,P001.y,P001.z);
                    glVertex3f(P101.x,P101.y,P101.z);
                    glVertex3f(P100.x,P100.y,P100.z);
                    glVertex3f(P101.x,P101.y,P101.z);
                    glVertex3f(P100.x,P100.y,P100.z);
                    glVertex3f(P110.x,P110.y,P110.z);
                    glVertex3f(P010.x,P010.y,P010.z);
                    glVertex3f(P011.x,P011.y,P011.z);
                    glVertex3f(P010.x,P010.y,P010.z);
                    glVertex3f(P110.x,P110.y,P110.z);
                    glVertex3f(P011.x,P011.y,P011.z);
                    glVertex3f(P111.x,P111.y,P111.z);
                    glVertex3f(P101.x,P101.y,P101.z);
                    glVertex3f(P111.x,P111.y,P111.z);
                    glVertex3f(P111.x,P111.y,P111.z);
                    glVertex3f(P110.x,P110.y,P110.z);
                glEnd();
                */
                glPointSize(2.0f);
                glBegin(GL_POINTS);
                    glVertex3f(P.x,P.y,P.z);
                glEnd();
                CubePlane.clear();
                int n000 = 0;
                int n001 = 0;
                int n010 = 0;
                int n011 = 0;
                int n100 = 0;
                int n101 = 0;
                int n110 = 0;
                int n111 = 0;
                int nx = 0;
                int ny = 0;
                int nz = 0;

                
                {
                    if((v000<=0 && v001>0) || (v000>0 && v001<=0))
                    {
                        Point inter = {0,0,0};
                        double vint = interpolation(v000,v001);
                        inter.x = P000.x*vint + (1-vint)*P001.x;
                        inter.y = P000.y*vint + (1-vint)*P001.y;
                        inter.z = P000.z*vint + (1-vint)*P001.z;
                        EdgeCube ec = {c000,c001,cZ,inter};
                        n000++;
                        n001++;
                        nz++;
                        CubePlane[ec.toString()]=ec;
                    }
                    if((v000<=0 && v010>0) || (v000>0 && v010<=0))
                    {
                        Point inter = {0,0,0};
                        double vint = interpolation(v000,v010);
                        inter.x = P000.x*vint + (1-vint)*P010.x;
                        inter.y = P000.y*vint + (1-vint)*P010.y;
                        inter.z = P000.z*vint + (1-vint)*P010.z;
                        EdgeCube ec = {c000,c010,cY,inter};
                        n000++;
                        n010++;
                        ny++;
                        CubePlane[ec.toString()]=ec;
                    }
                    if((v000<=0 && v100>0) || (v000>0 && v100<=0))
                    {
                        Point inter = {0,0,0};
                        double vint = interpolation(v000,v100);
                        inter.x = P000.x*vint + (1-vint)*P100.x;
                        inter.y = P000.y*vint + (1-vint)*P100.y;
                        inter.z = P000.z*vint + (1-vint)*P100.z;
                        EdgeCube ec = {c000,c100,cX,inter};
                        n000++;
                        n100++;
                        nx++;
                        CubePlane[ec.toString()]=ec;
                    }

                    if((v111<=0 && v110>0) || (v111>0 && v110<=0))
                    {
                        Point inter = {0,0,0};
                        double vint = interpolation(v111,v110);
                        inter.x = P111.x*vint + (1-vint)*P110.x;
                        inter.y = P111.y*vint + (1-vint)*P110.y;
                        inter.z = P111.z*vint + (1-vint)*P110.z;
                        EdgeCube ec = {c110,c111,cZ,inter};
                        n111++;
                        n110++;
                        nz++;
                        CubePlane[ec.toString()]=ec;
                    }
                    if((v111<=0 && v101>0) || (v111>0 && v101<=0))
                    {
                        Point inter = {0,0,0};
                        double vint = interpolation(v111,v101);
                        inter.x = P111.x*vint + (1-vint)*P101.x;
                        inter.y = P111.y*vint + (1-vint)*P101.y;
                        inter.z = P111.z*vint + (1-vint)*P101.z;
                        EdgeCube ec = {c101,c111,cY,inter};
                        n111++;
                        n101++;
                        ny++;
                        CubePlane[ec.toString()]=ec;
                    }
                    if((v111<=0 && v011>0) || (v111>0 && v011<=0))
                    {
                        Point inter = {0,0,0};
                        double vint = interpolation(v111,v011);
                        inter.x = P111.x*vint + (1-vint)*P011.x;
                        inter.y = P111.y*vint + (1-vint)*P011.y;
                        inter.z = P111.z*vint + (1-vint)*P011.z;
                        EdgeCube ec = {c011,c111,cX,inter};
                        n111++;
                        n011++;
                        nx++;
                        CubePlane[ec.toString()]=ec;
                    }

                    if((v011<=0 && v001>0) || (v011>0 && v001<=0))
                    {
                        Point inter = {0,0,0};
                        double vint = interpolation(v011,v001);
                        inter.x = P011.x*vint + (1-vint)*P001.x;
                        inter.y = P011.y*vint + (1-vint)*P001.y;
                        inter.z = P011.z*vint + (1-vint)*P001.z;
                        EdgeCube ec = {c001,c011,cY,inter};
                        n011++;
                        n001++;
                        ny++;
                        CubePlane[ec.toString()]=ec;
                    }
                    if((v101<=0 && v001>0) || (v101>0 && v001<=0))
                    {
                        Point inter = {0,0,0};
                        double vint = interpolation(v101,v001);
                        inter.x = P101.x*vint + (1-vint)*P001.x;
                        inter.y = P101.y*vint + (1-vint)*P001.y;
                        inter.z = P101.z*vint + (1-vint)*P001.z;
                        EdgeCube ec = {c001,c101,cX,inter};
                        n101++;
                        n001++;
                        nx++;
                        CubePlane[ec.toString()]=ec;
                    }

                    if((v010<=0 && v011>0) || (v010>0 && v011<=0))
                    {
                        Point inter = {0,0,0};
                        double vint = interpolation(v010,v011);
                        inter.x = P010.x*vint + (1-vint)*P011.x;
                        inter.y = P010.y*vint + (1-vint)*P011.y;
                        inter.z = P010.z*vint + (1-vint)*P011.z;
                        EdgeCube ec = {c010,c011,cZ,inter};
                        n010++;
                        n011++;
                        nz++;
                        CubePlane[ec.toString()]=ec;
                    }
                    if((v010<=0 && v110>0) || (v010>0 && v110<=0))
                    {
                        Point inter = {0,0,0};
                        double vint = interpolation(v010,v110);
                        inter.x = P010.x*vint + (1-vint)*P110.x;
                        inter.y = P010.y*vint + (1-vint)*P110.y;
                        inter.z = P010.z*vint + (1-vint)*P110.z;
                        EdgeCube ec = {c010,c110,cX,inter};
                        n010++;
                        n110++;
                        nx++;
                        CubePlane[ec.toString()]=ec;
                    }

                    if((v101<=0 && v100>0) || (v101>0 && v100<=0))
                    {
                        Point inter = {0,0,0};
                        double vint = interpolation(v101,v100);
                        inter.x = P101.x*vint + (1-vint)*P100.x;
                        inter.y = P101.y*vint + (1-vint)*P100.y;
                        inter.z = P101.z*vint + (1-vint)*P100.z;
                        EdgeCube ec = {c100,c101,cZ,inter};
                        n101++;
                        n100++;
                        nz++;
                        CubePlane[ec.toString()]=ec;
                    }
                    if((v110<=0 && v100>0) || (v110>0 && v100<=0))
                    {
                        Point inter = {0,0,0};
                        double vint = interpolation(v110,v100);
                        inter.x = P110.x*vint + (1-vint)*P100.x;
                        inter.y = P110.y*vint + (1-vint)*P100.y;
                        inter.z = P110.z*vint + (1-vint)*P100.z;
                        EdgeCube ec = {c100,c110,cY,inter};
                        n100++;
                        n110++;
                        ny++;
                        CubePlane[ec.toString()]=ec;
                    }
                }

                int s = n000 + n001 + n010 + n011 + n100 + n101 + n110 + n111 + nz + ny + nx;
                if(s==0)
                {
                    double v0 = abs(v000+v001+v010+v011)/4.0;
                    double v1 = abs(v100+v101+v110+v111)/4.0;
                    double v2 = abs(v000+v010+v100+v110)/4.0;
                    double v3 = abs(v001+v011+v101+v111)/4.0;
                    double v4 = abs(v000+v001+v100+v101)/4.0;
                    double v5 = abs(v010+v011+v110+v111)/4.0;
                    double m = min(min(min(v0,v1),min(v2,v3)),min(v4,v5));
                    int ra = rand()%25;
                    if((v0==m && ra>=6) || ra==0 || ra==24)
                    {
                        cout<<"Dir v0 "<<m<<" "<<v0<<endl;
                        pile.push_back({Pi.x-1,Pi.y,Pi.z});
                    }
                    if((v1==m && ra>=6) || ra==1 || ra==24)
                    {
                        cout<<"Dir v1 "<<m<<" "<<v1<<endl;
                        pile.push_back({Pi.x+1,Pi.y,Pi.z});
                    }
                    if((v2==m && ra>=6) || ra==2 || ra==24)
                    {
                        cout<<"Dir v2 "<<m<<" "<<v2<<endl;
                        pile.push_back({Pi.x,Pi.y,Pi.z-1});
                    }
                    if((v3==m && ra>=6) || ra==3 || ra==24)
                    {
                        cout<<"Dir v3 "<<m<<" "<<v3<<endl;
                        pile.push_back({Pi.x,Pi.y,Pi.z+1});
                    }
                    if((v4==m && ra>=6) || ra==4 || ra==24)
                    {
                        cout<<"Dir v4 "<<m<<" "<<v4<<endl;
                        pile.push_back({Pi.x,Pi.y-1,Pi.z});
                    }
                    if((v5==m && ra>=6) || ra==5 || ra==24)
                    {
                        cout<<"Dir v5 "<<m<<" "<<v5<<endl;
                        pile.push_back({Pi.x,Pi.y+1,Pi.z});
                    }
                }
                else
                {
                    cout<<n000<<"\t"<<n001<<"\t"<<n010<<"\t"<<n011<<"\t"<<endl;
                    cout<<n100<<"\t"<<n101<<"\t"<<n110<<"\t"<<n111<<"\t"<<endl;

                    if(n000==3)
                    {
                        n000=0;
                        IndiceCube i = c000;
                        IndiceCube i0 = i.add(cX);
                        IndiceCube i1 = i.add(cY);
                        IndiceCube i2 = i.add(cZ);
                        int ip0 = i0.toInt();
                        int ip1 = i1.toInt();
                        int ip2 = i2.toInt();
                        cout<<i.toString()<<"\t"<<i0.toString()<<"\t"<<i1.toString()<<"\t"<<i2.toString()<<endl;
                        switch(ip0)
                        {
                            case 0:
                                n000--;
                                break;
                            case 1:
                                n001--;
                                break;
                            case 2:
                                n010--;
                                break;
                            case 3:
                                n011--;
                                break;
                            case 4:
                                n100--;
                                break;
                            case 5:
                                n101--;
                                break;
                            case 6:
                                n110--;
                                break;
                            case 7:
                                n111--;
                        }
                        switch(ip1)
                        {
                            case 0:
                                n000--;
                                break;
                            case 1:
                                n001--;
                                break;
                            case 2:
                                n010--;
                                break;
                            case 3:
                                n011--;
                                break;
                            case 4:
                                n100--;
                                break;
                            case 5:
                                n101--;
                                break;
                            case 6:
                                n110--;
                                break;
                            case 7:
                                n111--;
                        }
                        switch(ip2)
                        {
                            case 0:
                                n000--;
                                break;
                            case 1:
                                n001--;
                                break;
                            case 2:
                                n010--;
                                break;
                            case 3:
                                n011--;
                                break;
                            case 4:
                                n100--;
                                break;
                            case 5:
                                n101--;
                                break;
                            case 6:
                                n110--;
                                break;
                            case 7:
                                n111--;
                        }
                        Point P0 = {0,0,0};
                        Point P1 = {0,0,0};
                        Point P2 = {0,0,0};
                        if(i.toInt()<i0.toInt())
                        {
                            P0 = CubePlane[i.toString()+":"+i0.toString()+":x"].inter;
                        }
                        else
                        {
                            P0 = CubePlane[i0.toString()+":"+i.toString()+":x"].inter;
                        }
                        if(i.toInt()<i1.toInt())
                        {
                            P1 = CubePlane[i.toString()+":"+i1.toString()+":y"].inter;
                        }
                        else
                        {
                            P1 = CubePlane[i1.toString()+":"+i.toString()+":y"].inter;
                        }
                        if(i.toInt()<i2.toInt())
                        {
                            P2 = CubePlane[i.toString()+":"+i2.toString()+":z"].inter;
                        }
                        else
                        {
                            P2 = CubePlane[i2.toString()+":"+i.toString()+":z"].inter;
                        }
                        glBegin(GL_TRIANGLES);
                        glVertex3f(P0.x,P0.y,P0.z);
                        glVertex3f(P1.x,P1.y,P1.z);
                        glVertex3f(P2.x,P2.y,P2.z);
                        glEnd();
                    }

                    if(n001==3)
                    {
                        n001=0;
                        IndiceCube i = c001;
                        IndiceCube i0 = i.add(cX);
                        IndiceCube i1 = i.add(cY);
                        IndiceCube i2 = i.add(cZ);
                        int ip0 = i0.toInt();
                        int ip1 = i1.toInt();
                        int ip2 = i2.toInt();
                        cout<<i.toString()<<"\t"<<i0.toString()<<"\t"<<i1.toString()<<"\t"<<i2.toString()<<endl;
                        switch(ip0)
                        {
                            case 0:
                                n000--;
                                break;
                            case 1:
                                n001--;
                                break;
                            case 2:
                                n010--;
                                break;
                            case 3:
                                n011--;
                                break;
                            case 4:
                                n100--;
                                break;
                            case 5:
                                n101--;
                                break;
                            case 6:
                                n110--;
                                break;
                            case 7:
                                n111--;
                        }
                        switch(ip1)
                        {
                            case 0:
                                n000--;
                                break;
                            case 1:
                                n001--;
                                break;
                            case 2:
                                n010--;
                                break;
                            case 3:
                                n011--;
                                break;
                            case 4:
                                n100--;
                                break;
                            case 5:
                                n101--;
                                break;
                            case 6:
                                n110--;
                                break;
                            case 7:
                                n111--;
                        }
                        switch(ip2)
                        {
                            case 0:
                                n000--;
                                break;
                            case 1:
                                n001--;
                                break;
                            case 2:
                                n010--;
                                break;
                            case 3:
                                n011--;
                                break;
                            case 4:
                                n100--;
                                break;
                            case 5:
                                n101--;
                                break;
                            case 6:
                                n110--;
                                break;
                            case 7:
                                n111--;
                        }
                        Point P0 = {0,0,0};
                        Point P1 = {0,0,0};
                        Point P2 = {0,0,0};
                        if(i.toInt()<i0.toInt())
                        {
                            P0 = CubePlane[i.toString()+":"+i0.toString()+":x"].inter;
                        }
                        else
                        {
                            P0 = CubePlane[i0.toString()+":"+i.toString()+":x"].inter;
                        }
                        if(i.toInt()<i1.toInt())
                        {
                            P1 = CubePlane[i.toString()+":"+i1.toString()+":y"].inter;
                        }
                        else
                        {
                            P1 = CubePlane[i1.toString()+":"+i.toString()+":y"].inter;
                        }
                        if(i.toInt()<i2.toInt())
                        {
                            P2 = CubePlane[i.toString()+":"+i2.toString()+":z"].inter;
                        }
                        else
                        {
                            P2 = CubePlane[i2.toString()+":"+i.toString()+":z"].inter;
                        }
                        glBegin(GL_TRIANGLES);
                        glVertex3f(P0.x,P0.y,P0.z);
                        glVertex3f(P1.x,P1.y,P1.z);
                        glVertex3f(P2.x,P2.y,P2.z);
                        glEnd();
                    }

                    if(n010==3)
                    {
                        n010=0;
                        IndiceCube i = c010;
                        IndiceCube i0 = i.add(cX);
                        IndiceCube i1 = i.add(cY);
                        IndiceCube i2 = i.add(cZ);
                        int ip0 = i0.toInt();
                        int ip1 = i1.toInt();
                        int ip2 = i2.toInt();
                        cout<<i.toString()<<"\t"<<i0.toString()<<"\t"<<i1.toString()<<"\t"<<i2.toString()<<endl;
                        switch(ip0)
                        {
                            case 0:
                                n000--;
                                break;
                            case 1:
                                n001--;
                                break;
                            case 2:
                                n010--;
                                break;
                            case 3:
                                n011--;
                                break;
                            case 4:
                                n100--;
                                break;
                            case 5:
                                n101--;
                                break;
                            case 6:
                                n110--;
                                break;
                            case 7:
                                n111--;
                        }
                        switch(ip1)
                        {
                            case 0:
                                n000--;
                                break;
                            case 1:
                                n001--;
                                break;
                            case 2:
                                n010--;
                                break;
                            case 3:
                                n011--;
                                break;
                            case 4:
                                n100--;
                                break;
                            case 5:
                                n101--;
                                break;
                            case 6:
                                n110--;
                                break;
                            case 7:
                                n111--;
                        }
                        switch(ip2)
                        {
                            case 0:
                                n000--;
                                break;
                            case 1:
                                n001--;
                                break;
                            case 2:
                                n010--;
                                break;
                            case 3:
                                n011--;
                                break;
                            case 4:
                                n100--;
                                break;
                            case 5:
                                n101--;
                                break;
                            case 6:
                                n110--;
                                break;
                            case 7:
                                n111--;
                        }
                        Point P0 = {0,0,0};
                        Point P1 = {0,0,0};
                        Point P2 = {0,0,0};
                        if(i.toInt()<i0.toInt())
                        {
                            P0 = CubePlane[i.toString()+":"+i0.toString()+":x"].inter;
                        }
                        else
                        {
                            P0 = CubePlane[i0.toString()+":"+i.toString()+":x"].inter;
                        }
                        if(i.toInt()<i1.toInt())
                        {
                            P1 = CubePlane[i.toString()+":"+i1.toString()+":y"].inter;
                        }
                        else
                        {
                            P1 = CubePlane[i1.toString()+":"+i.toString()+":y"].inter;
                        }
                        if(i.toInt()<i2.toInt())
                        {
                            P2 = CubePlane[i.toString()+":"+i2.toString()+":z"].inter;
                        }
                        else
                        {
                            P2 = CubePlane[i2.toString()+":"+i.toString()+":z"].inter;
                        }
                        glBegin(GL_TRIANGLES);
                        glVertex3f(P0.x,P0.y,P0.z);
                        glVertex3f(P1.x,P1.y,P1.z);
                        glVertex3f(P2.x,P2.y,P2.z);
                        glEnd();
                    }

                    if(n100==3)
                    {
                        n100=0;
                        IndiceCube i = c100;
                        IndiceCube i0 = i.add(cX);
                        IndiceCube i1 = i.add(cY);
                        IndiceCube i2 = i.add(cZ);
                        int ip0 = i0.toInt();
                        int ip1 = i1.toInt();
                        int ip2 = i2.toInt();
                        cout<<i.toString()<<"\t"<<i0.toString()<<"\t"<<i1.toString()<<"\t"<<i2.toString()<<endl;
                        switch(ip0)
                        {
                            case 0:
                                n000--;
                                break;
                            case 1:
                                n001--;
                                break;
                            case 2:
                                n010--;
                                break;
                            case 3:
                                n011--;
                                break;
                            case 4:
                                n100--;
                                break;
                            case 5:
                                n101--;
                                break;
                            case 6:
                                n110--;
                                break;
                            case 7:
                                n111--;
                        }
                        switch(ip1)
                        {
                            case 0:
                                n000--;
                                break;
                            case 1:
                                n001--;
                                break;
                            case 2:
                                n010--;
                                break;
                            case 3:
                                n011--;
                                break;
                            case 4:
                                n100--;
                                break;
                            case 5:
                                n101--;
                                break;
                            case 6:
                                n110--;
                                break;
                            case 7:
                                n111--;
                        }
                        switch(ip2)
                        {
                            case 0:
                                n000--;
                                break;
                            case 1:
                                n001--;
                                break;
                            case 2:
                                n010--;
                                break;
                            case 3:
                                n011--;
                                break;
                            case 4:
                                n100--;
                                break;
                            case 5:
                                n101--;
                                break;
                            case 6:
                                n110--;
                                break;
                            case 7:
                                n111--;
                        }
                        Point P0 = {0,0,0};
                        Point P1 = {0,0,0};
                        Point P2 = {0,0,0};
                        if(i.toInt()<i0.toInt())
                        {
                            P0 = CubePlane[i.toString()+":"+i0.toString()+":x"].inter;
                        }
                        else
                        {
                            P0 = CubePlane[i0.toString()+":"+i.toString()+":x"].inter;
                        }
                        if(i.toInt()<i1.toInt())
                        {
                            P1 = CubePlane[i.toString()+":"+i1.toString()+":y"].inter;
                        }
                        else
                        {
                            P1 = CubePlane[i1.toString()+":"+i.toString()+":y"].inter;
                        }
                        if(i.toInt()<i2.toInt())
                        {
                            P2 = CubePlane[i.toString()+":"+i2.toString()+":z"].inter;
                        }
                        else
                        {
                            P2 = CubePlane[i2.toString()+":"+i.toString()+":z"].inter;
                        }
                        glBegin(GL_TRIANGLES);
                        glVertex3f(P0.x,P0.y,P0.z);
                        glVertex3f(P1.x,P1.y,P1.z);
                        glVertex3f(P2.x,P2.y,P2.z);
                        glEnd();
                    }

                    if(n011==3)
                    {
                        n011=0;
                        IndiceCube i = c011;
                        IndiceCube i0 = i.add(cX);
                        IndiceCube i1 = i.add(cY);
                        IndiceCube i2 = i.add(cZ);
                        int ip0 = i0.toInt();
                        int ip1 = i1.toInt();
                        int ip2 = i2.toInt();
                        cout<<i.toString()<<"\t"<<i0.toString()<<"\t"<<i1.toString()<<"\t"<<i2.toString()<<endl;
                        switch(ip0)
                        {
                            case 0:
                                n000--;
                                break;
                            case 1:
                                n001--;
                                break;
                            case 2:
                                n010--;
                                break;
                            case 3:
                                n011--;
                                break;
                            case 4:
                                n100--;
                                break;
                            case 5:
                                n101--;
                                break;
                            case 6:
                                n110--;
                                break;
                            case 7:
                                n111--;
                        }
                        switch(ip1)
                        {
                            case 0:
                                n000--;
                                break;
                            case 1:
                                n001--;
                                break;
                            case 2:
                                n010--;
                                break;
                            case 3:
                                n011--;
                                break;
                            case 4:
                                n100--;
                                break;
                            case 5:
                                n101--;
                                break;
                            case 6:
                                n110--;
                                break;
                            case 7:
                                n111--;
                        }
                        switch(ip2)
                        {
                            case 0:
                                n000--;
                                break;
                            case 1:
                                n001--;
                                break;
                            case 2:
                                n010--;
                                break;
                            case 3:
                                n011--;
                                break;
                            case 4:
                                n100--;
                                break;
                            case 5:
                                n101--;
                                break;
                            case 6:
                                n110--;
                                break;
                            case 7:
                                n111--;
                        }
                        Point P0 = {0,0,0};
                        Point P1 = {0,0,0};
                        Point P2 = {0,0,0};
                        if(i.toInt()<i0.toInt())
                        {
                            P0 = CubePlane[i.toString()+":"+i0.toString()+":x"].inter;
                        }
                        else
                        {
                            P0 = CubePlane[i0.toString()+":"+i.toString()+":x"].inter;
                        }
                        if(i.toInt()<i1.toInt())
                        {
                            P1 = CubePlane[i.toString()+":"+i1.toString()+":y"].inter;
                        }
                        else
                        {
                            P1 = CubePlane[i1.toString()+":"+i.toString()+":y"].inter;
                        }
                        if(i.toInt()<i2.toInt())
                        {
                            P2 = CubePlane[i.toString()+":"+i2.toString()+":z"].inter;
                        }
                        else
                        {
                            P2 = CubePlane[i2.toString()+":"+i.toString()+":z"].inter;
                        }
                        glBegin(GL_TRIANGLES);
                        glVertex3f(P0.x,P0.y,P0.z);
                        glVertex3f(P1.x,P1.y,P1.z);
                        glVertex3f(P2.x,P2.y,P2.z);
                        glEnd();
                    }

                    if(n101==3)
                    {
                        n101=0;
                        IndiceCube i = c101;
                        IndiceCube i0 = i.add(cX);
                        IndiceCube i1 = i.add(cY);
                        IndiceCube i2 = i.add(cZ);
                        int ip0 = i0.toInt();
                        int ip1 = i1.toInt();
                        int ip2 = i2.toInt();
                        cout<<i.toString()<<"\t"<<i0.toString()<<"\t"<<i1.toString()<<"\t"<<i2.toString()<<endl;
                        switch(ip0)
                        {
                            case 0:
                                n000--;
                                break;
                            case 1:
                                n001--;
                                break;
                            case 2:
                                n010--;
                                break;
                            case 3:
                                n011--;
                                break;
                            case 4:
                                n100--;
                                break;
                            case 5:
                                n101--;
                                break;
                            case 6:
                                n110--;
                                break;
                            case 7:
                                n111--;
                        }
                        switch(ip1)
                        {
                            case 0:
                                n000--;
                                break;
                            case 1:
                                n001--;
                                break;
                            case 2:
                                n010--;
                                break;
                            case 3:
                                n011--;
                                break;
                            case 4:
                                n100--;
                                break;
                            case 5:
                                n101--;
                                break;
                            case 6:
                                n110--;
                                break;
                            case 7:
                                n111--;
                        }
                        switch(ip2)
                        {
                            case 0:
                                n000--;
                                break;
                            case 1:
                                n001--;
                                break;
                            case 2:
                                n010--;
                                break;
                            case 3:
                                n011--;
                                break;
                            case 4:
                                n100--;
                                break;
                            case 5:
                                n101--;
                                break;
                            case 6:
                                n110--;
                                break;
                            case 7:
                                n111--;
                        }
                        Point P0 = {0,0,0};
                        Point P1 = {0,0,0};
                        Point P2 = {0,0,0};
                        if(i.toInt()<i0.toInt())
                        {
                            P0 = CubePlane[i.toString()+":"+i0.toString()+":x"].inter;
                        }
                        else
                        {
                            P0 = CubePlane[i0.toString()+":"+i.toString()+":x"].inter;
                        }
                        if(i.toInt()<i1.toInt())
                        {
                            P1 = CubePlane[i.toString()+":"+i1.toString()+":y"].inter;
                        }
                        else
                        {
                            P1 = CubePlane[i1.toString()+":"+i.toString()+":y"].inter;
                        }
                        if(i.toInt()<i2.toInt())
                        {
                            P2 = CubePlane[i.toString()+":"+i2.toString()+":z"].inter;
                        }
                        else
                        {
                            P2 = CubePlane[i2.toString()+":"+i.toString()+":z"].inter;
                        }
                        glBegin(GL_TRIANGLES);
                        glVertex3f(P0.x,P0.y,P0.z);
                        glVertex3f(P1.x,P1.y,P1.z);
                        glVertex3f(P2.x,P2.y,P2.z);
                        glEnd();
                    }

                    if(n110==3)
                    {
                        n110=0;
                        IndiceCube i = c110;
                        IndiceCube i0 = i.add(cX);
                        IndiceCube i1 = i.add(cY);
                        IndiceCube i2 = i.add(cZ);
                        int ip0 = i0.toInt();
                        int ip1 = i1.toInt();
                        int ip2 = i2.toInt();
                        cout<<i.toString()<<"\t"<<i0.toString()<<"\t"<<i1.toString()<<"\t"<<i2.toString()<<endl;
                        switch(ip0)
                        {
                            case 0:
                                n000--;
                                break;
                            case 1:
                                n001--;
                                break;
                            case 2:
                                n010--;
                                break;
                            case 3:
                                n011--;
                                break;
                            case 4:
                                n100--;
                                break;
                            case 5:
                                n101--;
                                break;
                            case 6:
                                n110--;
                                break;
                            case 7:
                                n111--;
                        }
                        switch(ip1)
                        {
                            case 0:
                                n000--;
                                break;
                            case 1:
                                n001--;
                                break;
                            case 2:
                                n010--;
                                break;
                            case 3:
                                n011--;
                                break;
                            case 4:
                                n100--;
                                break;
                            case 5:
                                n101--;
                                break;
                            case 6:
                                n110--;
                                break;
                            case 7:
                                n111--;
                        }
                        switch(ip2)
                        {
                            case 0:
                                n000--;
                                break;
                            case 1:
                                n001--;
                                break;
                            case 2:
                                n010--;
                                break;
                            case 3:
                                n011--;
                                break;
                            case 4:
                                n100--;
                                break;
                            case 5:
                                n101--;
                                break;
                            case 6:
                                n110--;
                                break;
                            case 7:
                                n111--;
                        }
                        Point P0 = {0,0,0};
                        Point P1 = {0,0,0};
                        Point P2 = {0,0,0};
                        if(i.toInt()<i0.toInt())
                        {
                            P0 = CubePlane[i.toString()+":"+i0.toString()+":x"].inter;
                        }
                        else
                        {
                            P0 = CubePlane[i0.toString()+":"+i.toString()+":x"].inter;
                        }
                        if(i.toInt()<i1.toInt())
                        {
                            P1 = CubePlane[i.toString()+":"+i1.toString()+":y"].inter;
                        }
                        else
                        {
                            P1 = CubePlane[i1.toString()+":"+i.toString()+":y"].inter;
                        }
                        if(i.toInt()<i2.toInt())
                        {
                            P2 = CubePlane[i.toString()+":"+i2.toString()+":z"].inter;
                        }
                        else
                        {
                            P2 = CubePlane[i2.toString()+":"+i.toString()+":z"].inter;
                        }
                        glBegin(GL_TRIANGLES);
                        glVertex3f(P0.x,P0.y,P0.z);
                        glVertex3f(P1.x,P1.y,P1.z);
                        glVertex3f(P2.x,P2.y,P2.z);
                        glEnd();
                    }

                    if(n111==3)
                    {
                        n111=0;
                        IndiceCube i = c111;
                        IndiceCube i0 = i.add(cX);
                        IndiceCube i1 = i.add(cY);
                        IndiceCube i2 = i.add(cZ);
                        int ip0 = i0.toInt();
                        int ip1 = i1.toInt();
                        int ip2 = i2.toInt();
                        cout<<i.toString()<<"\t"<<i0.toString()<<"\t"<<i1.toString()<<"\t"<<i2.toString()<<endl;
                        switch(ip0)
                        {
                            case 0:
                                n000--;
                                break;
                            case 1:
                                n001--;
                                break;
                            case 2:
                                n010--;
                                break;
                            case 3:
                                n011--;
                                break;
                            case 4:
                                n100--;
                                break;
                            case 5:
                                n101--;
                                break;
                            case 6:
                                n110--;
                                break;
                            case 7:
                                n111--;
                        }
                        switch(ip1)
                        {
                            case 0:
                                n000--;
                                break;
                            case 1:
                                n001--;
                                break;
                            case 2:
                                n010--;
                                break;
                            case 3:
                                n011--;
                                break;
                            case 4:
                                n100--;
                                break;
                            case 5:
                                n101--;
                                break;
                            case 6:
                                n110--;
                                break;
                            case 7:
                                n111--;
                        }
                        switch(ip2)
                        {
                            case 0:
                                n000--;
                                break;
                            case 1:
                                n001--;
                                break;
                            case 2:
                                n010--;
                                break;
                            case 3:
                                n011--;
                                break;
                            case 4:
                                n100--;
                                break;
                            case 5:
                                n101--;
                                break;
                            case 6:
                                n110--;
                                break;
                            case 7:
                                n111--;
                        }
                        Point P0 = {0,0,0};
                        Point P1 = {0,0,0};
                        Point P2 = {0,0,0};
                        if(i.toInt()<i0.toInt())
                        {
                            P0 = CubePlane[i.toString()+":"+i0.toString()+":x"].inter;
                        }
                        else
                        {
                            P0 = CubePlane[i0.toString()+":"+i.toString()+":x"].inter;
                        }
                        if(i.toInt()<i1.toInt())
                        {
                            P1 = CubePlane[i.toString()+":"+i1.toString()+":y"].inter;
                        }
                        else
                        {
                            P1 = CubePlane[i1.toString()+":"+i.toString()+":y"].inter;
                        }
                        if(i.toInt()<i2.toInt())
                        {
                            P2 = CubePlane[i.toString()+":"+i2.toString()+":z"].inter;
                        }
                        else
                        {
                            P2 = CubePlane[i2.toString()+":"+i.toString()+":z"].inter;
                        }
                        glBegin(GL_TRIANGLES);
                        glVertex3f(P0.x,P0.y,P0.z);
                        glVertex3f(P1.x,P1.y,P1.z);
                        glVertex3f(P2.x,P2.y,P2.z);
                        glEnd();
                    }

                    if(n000>0)
                    {
                        IndiceCube i = c000;
                        int* targ = &n000;
                        IndiceCube i0 = i.add(cX);
                        IndiceCube i1 = i.add(cY);
                        IndiceCube i2 = i.add(cZ);
                        EdgeCube e0 = {i,i0,cX,{0,0,0}};
                        EdgeCube e1 = {i,i1,cY,{0,0,0}};
                        EdgeCube e2 = {i,i2,cZ,{0,0,0}};
                        e0.order();
                        e1.order();
                        e2.order();
                        string chx = "";
                        if(CubePlane.find(e0.toString())!=CubePlane.end())
                        {
                            chx = e0.toString();
                            if(*targ!=2)
                            {
                                targ = &nx;
                            }
                        }
                        else
                        {
                            if(CubePlane.find(e1.toString())!=CubePlane.end())
                            {
                                chx = e1.toString();
                                if(*targ!=2)
                                {
                                    targ = &ny;
                                }
                            }
                            else
                            {
                                if(CubePlane.find(e2.toString())!=CubePlane.end())
                                {
                                    chx = e2.toString();
                                    if(*targ!=2)
                                    {
                                        targ = &nz;
                                    }
                                }
                            }
                        }
                        vector<Point> face;
                        int st = 0;
                        while(*targ!=0 && st<3)
                        {
                            st++;
                            cout<<n000<<"\t"<<n001<<"\t"<<n010<<"\t"<<n011<<"\t"<<endl;
                            cout<<n100<<"\t"<<n101<<"\t"<<n110<<"\t"<<n111<<"\t"<<endl;
                            EdgeCube ec = CubePlane[chx];
                            cout<<ec.toString()<<endl;
                            CubePlane.erase(ec.toString());
                            face.push_back(ec.inter);
                            int ip0 = ec.P0.toInt();
                            int ip1 = ec.P1.toInt();
                            int ipA = ec.coupe.toInt();
                            vector<EdgeCube> liste;
                            switch(ip0)
                            {
                                case 0:
                                    n000--;
                                    break;
                                case 1:
                                    n001--;
                                    break;
                                case 2:
                                    n010--;
                                    break;
                                case 3:
                                    n011--;
                                    break;
                                case 4:
                                    n100--;
                                    break;
                                case 5:
                                    n101--;
                                    break;
                                case 6:
                                    n110--;
                                    break;
                                case 7:
                                    n111--;
                            }
                            switch(ip1)
                            {
                                case 0:
                                    n000--;
                                    break;
                                case 1:
                                    n001--;
                                    break;
                                case 2:
                                    n010--;
                                    break;
                                case 3:
                                    n011--;
                                    break;
                                case 4:
                                    n100--;
                                    break;
                                case 5:
                                    n101--;
                                    break;
                                case 6:
                                    n110--;
                                    break;
                                case 7:
                                    n111--;
                            }
                            switch(ipA)
                            {
                                case 1:
                                    nz--;
                                    liste.push_back({ec.P0.add(cY),ec.P1.add(cY),cX,{0,0,0}});
                                    liste.push_back({ec.P0.add(cZ),ec.P1.add(cZ),cX,{0,0,0}});
                                    break;
                                case 2:
                                    liste.push_back({ec.P0.add(cX),ec.P1.add(cX),cY,{0,0,0}});
                                    liste.push_back({ec.P0.add(cZ),ec.P1.add(cZ),cY,{0,0,0}});
                                    ny--;
                                    break;
                                case 4:
                                    liste.push_back({ec.P0.add(cX),ec.P1.add(cX),cZ,{0,0,0}});
                                    liste.push_back({ec.P0.add(cY),ec.P1.add(cY),cZ,{0,0,0}});
                                    nx--;
                                    break;
                            }
                            if(ec.P0.toInt()!=i.toInt())
                                i = ec.P0;
                            else
                                i = ec.P1;
                            i0 = i.add(cX);
                            i1 = i.add(cY);
                            i2 = i.add(cZ);
                            e0 = {i,i0,cX,{0,0,0}};
                            e1 = {i,i1,cX,{0,0,0}};
                            e2 = {i,i2,cX,{0,0,0}};
                            e0.order();
                            e1.order();
                            e2.order();
                            liste.push_back(e0);
                            liste.push_back(e1);
                            liste.push_back(e2);
                            EdgeCube e = ec;
                            double min = -1;
                            for(int i=0; i<liste.size(); i++)
                            {
                                cout<<e.toString()<<endl;
                                if(CubePlane.find(liste[i].toString())!=CubePlane.end())
                                {
                                    EdgeCube et = CubePlane[liste[i].toString()];
                                    double dist = distFromSurface(ec.inter,et.inter,target);
                                    if(min==-1 || dist<min)
                                    {
                                        e = et;
                                        min = dist;
                                    }
                                    
                                }
                            }
                            chx = e.toString();
                        }
                        glBegin(GL_LINE_LOOP);
                        for(int i=0; i<face.size(); i++)
                        {
                            Point P = face[i];
                            glVertex3f(P.x,P.y,P.z);
                        }
                        glEnd();
                    }

                    if(n001>0)
                    {
                        IndiceCube i = c001;
                        int* targ = &n001;
                        IndiceCube i0 = i.add(cX);
                        IndiceCube i1 = i.add(cY);
                        IndiceCube i2 = i.add(cZ);
                        EdgeCube e0 = {i,i0,cX,{0,0,0}};
                        EdgeCube e1 = {i,i1,cY,{0,0,0}};
                        EdgeCube e2 = {i,i2,cZ,{0,0,0}};
                        e0.order();
                        e1.order();
                        e2.order();
                        string chx = "";
                        if(CubePlane.find(e0.toString())!=CubePlane.end())
                        {
                            chx = e0.toString();
                            if(*targ!=2)
                            {
                                targ = &nx;
                            }
                        }
                        else
                        {
                            if(CubePlane.find(e1.toString())!=CubePlane.end())
                            {
                                chx = e1.toString();
                                if(*targ!=2)
                                {
                                    targ = &ny;
                                }
                            }
                            else
                            {
                                if(CubePlane.find(e2.toString())!=CubePlane.end())
                                {
                                    chx = e2.toString();
                                    if(*targ!=2)
                                    {
                                        targ = &nz;
                                    }
                                }
                            }
                        }
                        vector<Point> face;
                        int st = 0;
                        while(*targ!=0 && st<3)
                        {
                            st++;
                            cout<<n000<<"\t"<<n001<<"\t"<<n010<<"\t"<<n011<<"\t"<<endl;
                            cout<<n100<<"\t"<<n101<<"\t"<<n110<<"\t"<<n111<<"\t"<<endl;
                            EdgeCube ec = CubePlane[chx];
                            cout<<ec.toString()<<endl;
                            CubePlane.erase(ec.toString());
                            face.push_back(ec.inter);
                            int ip0 = ec.P0.toInt();
                            int ip1 = ec.P1.toInt();
                            int ipA = ec.coupe.toInt();
                            vector<EdgeCube> liste;
                            switch(ip0)
                            {
                                case 0:
                                    n000--;
                                    break;
                                case 1:
                                    n001--;
                                    break;
                                case 2:
                                    n010--;
                                    break;
                                case 3:
                                    n011--;
                                    break;
                                case 4:
                                    n100--;
                                    break;
                                case 5:
                                    n101--;
                                    break;
                                case 6:
                                    n110--;
                                    break;
                                case 7:
                                    n111--;
                            }
                            switch(ip1)
                            {
                                case 0:
                                    n000--;
                                    break;
                                case 1:
                                    n001--;
                                    break;
                                case 2:
                                    n010--;
                                    break;
                                case 3:
                                    n011--;
                                    break;
                                case 4:
                                    n100--;
                                    break;
                                case 5:
                                    n101--;
                                    break;
                                case 6:
                                    n110--;
                                    break;
                                case 7:
                                    n111--;
                            }
                            switch(ipA)
                            {
                                case 1:
                                    nz--;
                                    liste.push_back({ec.P0.add(cY),ec.P1.add(cY),cX,{0,0,0}});
                                    liste.push_back({ec.P0.add(cZ),ec.P1.add(cZ),cX,{0,0,0}});
                                    break;
                                case 2:
                                    liste.push_back({ec.P0.add(cX),ec.P1.add(cX),cY,{0,0,0}});
                                    liste.push_back({ec.P0.add(cZ),ec.P1.add(cZ),cY,{0,0,0}});
                                    ny--;
                                    break;
                                case 4:
                                    liste.push_back({ec.P0.add(cX),ec.P1.add(cX),cZ,{0,0,0}});
                                    liste.push_back({ec.P0.add(cY),ec.P1.add(cY),cZ,{0,0,0}});
                                    nx--;
                                    break;
                            }
                            if(ec.P0.toInt()!=i.toInt())
                                i = ec.P0;
                            else
                                i = ec.P1;
                            i0 = i.add(cX);
                            i1 = i.add(cY);
                            i2 = i.add(cZ);
                            e0 = {i,i0,cX,{0,0,0}};
                            e1 = {i,i1,cX,{0,0,0}};
                            e2 = {i,i2,cX,{0,0,0}};
                            e0.order();
                            e1.order();
                            e2.order();
                            liste.push_back(e0);
                            liste.push_back(e1);
                            liste.push_back(e2);
                            EdgeCube e = ec;
                            double min = -1;
                            for(int i=0; i<liste.size(); i++)
                            {
                                if(CubePlane.find(liste[i].toString())!=CubePlane.end())
                                {
                                    EdgeCube et = CubePlane[liste[i].toString()];
                                    double dist = distFromSurface(ec.inter,et.inter,target);
                                    if(min==-1 || dist<min)
                                    {
                                        e = et;
                                        min = dist;
                                    }
                                    
                                }
                            }
                            chx = e.toString();
                        }
                        glBegin(GL_LINE_LOOP);
                        for(int i=0; i<face.size(); i++)
                        {
                            Point P = face[i];
                            glVertex3f(P.x,P.y,P.z);
                        }
                        glEnd();
                    }

                    if(n010>0)
                    {
                        IndiceCube i = c010;
                        int* targ = &n010;
                        IndiceCube i0 = i.add(cX);
                        IndiceCube i1 = i.add(cY);
                        IndiceCube i2 = i.add(cZ);
                        EdgeCube e0 = {i,i0,cX,{0,0,0}};
                        EdgeCube e1 = {i,i1,cY,{0,0,0}};
                        EdgeCube e2 = {i,i2,cZ,{0,0,0}};
                        e0.order();
                        e1.order();
                        e2.order();
                        string chx = "";
                        if(CubePlane.find(e0.toString())!=CubePlane.end())
                        {
                            chx = e0.toString();
                            if(*targ!=2)
                            {
                                targ = &nx;
                            }
                        }
                        else
                        {
                            if(CubePlane.find(e1.toString())!=CubePlane.end())
                            {
                                chx = e1.toString();
                                if(*targ!=2)
                                {
                                    targ = &ny;
                                }
                            }
                            else
                            {
                                if(CubePlane.find(e2.toString())!=CubePlane.end())
                                {
                                    chx = e2.toString();
                                    if(*targ!=2)
                                    {
                                        targ = &nz;
                                    }
                                }
                            }
                        }
                        vector<Point> face;
                        int st = 0;
                        while(*targ!=0 && st<3)
                        {
                            st++;
                            cout<<n000<<"\t"<<n001<<"\t"<<n010<<"\t"<<n011<<"\t"<<endl;
                            cout<<n100<<"\t"<<n101<<"\t"<<n110<<"\t"<<n111<<"\t"<<endl;
                            EdgeCube ec = CubePlane[chx];
                            cout<<ec.toString()<<endl;
                            CubePlane.erase(ec.toString());
                            face.push_back(ec.inter);
                            int ip0 = ec.P0.toInt();
                            int ip1 = ec.P1.toInt();
                            int ipA = ec.coupe.toInt();
                            vector<EdgeCube> liste;
                            switch(ip0)
                            {
                                case 0:
                                    n000--;
                                    break;
                                case 1:
                                    n001--;
                                    break;
                                case 2:
                                    n010--;
                                    break;
                                case 3:
                                    n011--;
                                    break;
                                case 4:
                                    n100--;
                                    break;
                                case 5:
                                    n101--;
                                    break;
                                case 6:
                                    n110--;
                                    break;
                                case 7:
                                    n111--;
                            }
                            switch(ip1)
                            {
                                case 0:
                                    n000--;
                                    break;
                                case 1:
                                    n001--;
                                    break;
                                case 2:
                                    n010--;
                                    break;
                                case 3:
                                    n011--;
                                    break;
                                case 4:
                                    n100--;
                                    break;
                                case 5:
                                    n101--;
                                    break;
                                case 6:
                                    n110--;
                                    break;
                                case 7:
                                    n111--;
                            }
                            switch(ipA)
                            {
                                case 1:
                                    nz--;
                                    liste.push_back({ec.P0.add(cY),ec.P1.add(cY),cX,{0,0,0}});
                                    liste.push_back({ec.P0.add(cZ),ec.P1.add(cZ),cX,{0,0,0}});
                                    break;
                                case 2:
                                    liste.push_back({ec.P0.add(cX),ec.P1.add(cX),cY,{0,0,0}});
                                    liste.push_back({ec.P0.add(cZ),ec.P1.add(cZ),cY,{0,0,0}});
                                    ny--;
                                    break;
                                case 4:
                                    liste.push_back({ec.P0.add(cX),ec.P1.add(cX),cZ,{0,0,0}});
                                    liste.push_back({ec.P0.add(cY),ec.P1.add(cY),cZ,{0,0,0}});
                                    nx--;
                                    break;
                            }
                            if(ec.P0.toInt()!=i.toInt())
                                i = ec.P0;
                            else
                                i = ec.P1;
                            i0 = i.add(cX);
                            i1 = i.add(cY);
                            i2 = i.add(cZ);
                            e0 = {i,i0,cX,{0,0,0}};
                            e1 = {i,i1,cX,{0,0,0}};
                            e2 = {i,i2,cX,{0,0,0}};
                            e0.order();
                            e1.order();
                            e2.order();
                            liste.push_back(e0);
                            liste.push_back(e1);
                            liste.push_back(e2);
                            EdgeCube e = ec;
                            double min = -1;
                            for(int i=0; i<liste.size(); i++)
                            {
                                if(CubePlane.find(liste[i].toString())!=CubePlane.end())
                                {
                                    EdgeCube et = CubePlane[liste[i].toString()];
                                    double dist = distFromSurface(ec.inter,et.inter,target);
                                    if(min==-1 || dist<min)
                                    {
                                        e = et;
                                        min = dist;
                                    }
                                    
                                }
                            }
                            chx = e.toString();
                        }
                        glBegin(GL_LINE_LOOP);
                        for(int i=0; i<face.size(); i++)
                        {
                            Point P = face[i];
                            glVertex3f(P.x,P.y,P.z);
                        }
                        glEnd();
                    }

                    if(n100>0)
                    {
                        IndiceCube i = c100;
                        int* targ = &n100;
                        IndiceCube i0 = i.add(cX);
                        IndiceCube i1 = i.add(cY);
                        IndiceCube i2 = i.add(cZ);
                        EdgeCube e0 = {i,i0,cX,{0,0,0}};
                        EdgeCube e1 = {i,i1,cY,{0,0,0}};
                        EdgeCube e2 = {i,i2,cZ,{0,0,0}};
                        e0.order();
                        e1.order();
                        e2.order();
                        string chx = "";
                        if(CubePlane.find(e0.toString())!=CubePlane.end())
                        {
                            chx = e0.toString();
                            if(*targ!=2)
                            {
                                targ = &nx;
                            }
                        }
                        else
                        {
                            if(CubePlane.find(e1.toString())!=CubePlane.end())
                            {
                                chx = e1.toString();
                                if(*targ!=2)
                                {
                                    targ = &ny;
                                }
                            }
                            else
                            {
                                if(CubePlane.find(e2.toString())!=CubePlane.end())
                                {
                                    chx = e2.toString();
                                    if(*targ!=2)
                                    {
                                        targ = &nz;
                                    }
                                }
                            }
                        }
                        vector<Point> face;
                        int st = 0;
                        while(*targ!=0 && st<3)
                        {
                            st++;
                            cout<<n000<<"\t"<<n001<<"\t"<<n010<<"\t"<<n011<<"\t"<<endl;
                            cout<<n100<<"\t"<<n101<<"\t"<<n110<<"\t"<<n111<<"\t"<<endl;
                            EdgeCube ec = CubePlane[chx];
                            cout<<ec.toString()<<endl;
                            CubePlane.erase(ec.toString());
                            face.push_back(ec.inter);
                            int ip0 = ec.P0.toInt();
                            int ip1 = ec.P1.toInt();
                            int ipA = ec.coupe.toInt();
                            vector<EdgeCube> liste;
                            switch(ip0)
                            {
                                case 0:
                                    n000--;
                                    break;
                                case 1:
                                    n001--;
                                    break;
                                case 2:
                                    n010--;
                                    break;
                                case 3:
                                    n011--;
                                    break;
                                case 4:
                                    n100--;
                                    break;
                                case 5:
                                    n101--;
                                    break;
                                case 6:
                                    n110--;
                                    break;
                                case 7:
                                    n111--;
                            }
                            switch(ip1)
                            {
                                case 0:
                                    n000--;
                                    break;
                                case 1:
                                    n001--;
                                    break;
                                case 2:
                                    n010--;
                                    break;
                                case 3:
                                    n011--;
                                    break;
                                case 4:
                                    n100--;
                                    break;
                                case 5:
                                    n101--;
                                    break;
                                case 6:
                                    n110--;
                                    break;
                                case 7:
                                    n111--;
                            }
                            switch(ipA)
                            {
                                case 1:
                                    nz--;
                                    liste.push_back({ec.P0.add(cY),ec.P1.add(cY),cX,{0,0,0}});
                                    liste.push_back({ec.P0.add(cZ),ec.P1.add(cZ),cX,{0,0,0}});
                                    break;
                                case 2:
                                    liste.push_back({ec.P0.add(cX),ec.P1.add(cX),cY,{0,0,0}});
                                    liste.push_back({ec.P0.add(cZ),ec.P1.add(cZ),cY,{0,0,0}});
                                    ny--;
                                    break;
                                case 4:
                                    liste.push_back({ec.P0.add(cX),ec.P1.add(cX),cZ,{0,0,0}});
                                    liste.push_back({ec.P0.add(cY),ec.P1.add(cY),cZ,{0,0,0}});
                                    nx--;
                                    break;
                            }
                            if(ec.P0.toInt()!=i.toInt())
                                i = ec.P0;
                            else
                                i = ec.P1;
                            i0 = i.add(cX);
                            i1 = i.add(cY);
                            i2 = i.add(cZ);
                            e0 = {i,i0,cX,{0,0,0}};
                            e1 = {i,i1,cX,{0,0,0}};
                            e2 = {i,i2,cX,{0,0,0}};
                            e0.order();
                            e1.order();
                            e2.order();
                            liste.push_back(e0);
                            liste.push_back(e1);
                            liste.push_back(e2);
                            EdgeCube e = ec;
                            double min = -1;
                            for(int i=0; i<liste.size(); i++)
                            {
                                if(CubePlane.find(liste[i].toString())!=CubePlane.end())
                                {
                                    EdgeCube et = CubePlane[liste[i].toString()];
                                    double dist = distFromSurface(ec.inter,et.inter,target);
                                    if(min==-1 || dist<min)
                                    {
                                        e = et;
                                        min = dist;
                                    }
                                    
                                }
                            }
                            chx = e.toString();
                        }
                        glBegin(GL_LINE_LOOP);
                        for(int i=0; i<face.size(); i++)
                        {
                            Point P = face[i];
                            glVertex3f(P.x,P.y,P.z);
                        }
                        glEnd();
                    }

                    if(n011>0)
                    {
                        IndiceCube i = c011;
                        int* targ = &n011;
                        IndiceCube i0 = i.add(cX);
                        IndiceCube i1 = i.add(cY);
                        IndiceCube i2 = i.add(cZ);
                        EdgeCube e0 = {i,i0,cX,{0,0,0}};
                        EdgeCube e1 = {i,i1,cY,{0,0,0}};
                        EdgeCube e2 = {i,i2,cZ,{0,0,0}};
                        e0.order();
                        e1.order();
                        e2.order();
                        string chx = "";
                        if(CubePlane.find(e0.toString())!=CubePlane.end())
                        {
                            chx = e0.toString();
                            if(*targ!=2)
                            {
                                targ = &nx;
                            }
                        }
                        else
                        {
                            if(CubePlane.find(e1.toString())!=CubePlane.end())
                            {
                                chx = e1.toString();
                                if(*targ!=2)
                                {
                                    targ = &ny;
                                }
                            }
                            else
                            {
                                if(CubePlane.find(e2.toString())!=CubePlane.end())
                                {
                                    chx = e2.toString();
                                    if(*targ!=2)
                                    {
                                        targ = &nz;
                                    }
                                }
                            }
                        }
                        vector<Point> face;
                        int st = 0;
                        while(*targ!=0 && st<3)
                        {
                            st++;
                            cout<<n000<<"\t"<<n001<<"\t"<<n010<<"\t"<<n011<<"\t"<<endl;
                            cout<<n100<<"\t"<<n101<<"\t"<<n110<<"\t"<<n111<<"\t"<<endl;
                            EdgeCube ec = CubePlane[chx];
                            cout<<ec.toString()<<endl;
                            CubePlane.erase(ec.toString());
                            face.push_back(ec.inter);
                            int ip0 = ec.P0.toInt();
                            int ip1 = ec.P1.toInt();
                            int ipA = ec.coupe.toInt();
                            vector<EdgeCube> liste;
                            switch(ip0)
                            {
                                case 0:
                                    n000--;
                                    break;
                                case 1:
                                    n001--;
                                    break;
                                case 2:
                                    n010--;
                                    break;
                                case 3:
                                    n011--;
                                    break;
                                case 4:
                                    n100--;
                                    break;
                                case 5:
                                    n101--;
                                    break;
                                case 6:
                                    n110--;
                                    break;
                                case 7:
                                    n111--;
                            }
                            switch(ip1)
                            {
                                case 0:
                                    n000--;
                                    break;
                                case 1:
                                    n001--;
                                    break;
                                case 2:
                                    n010--;
                                    break;
                                case 3:
                                    n011--;
                                    break;
                                case 4:
                                    n100--;
                                    break;
                                case 5:
                                    n101--;
                                    break;
                                case 6:
                                    n110--;
                                    break;
                                case 7:
                                    n111--;
                            }
                            switch(ipA)
                            {
                                case 1:
                                    nz--;
                                    liste.push_back({ec.P0.add(cY),ec.P1.add(cY),cX,{0,0,0}});
                                    liste.push_back({ec.P0.add(cZ),ec.P1.add(cZ),cX,{0,0,0}});
                                    break;
                                case 2:
                                    liste.push_back({ec.P0.add(cX),ec.P1.add(cX),cY,{0,0,0}});
                                    liste.push_back({ec.P0.add(cZ),ec.P1.add(cZ),cY,{0,0,0}});
                                    ny--;
                                    break;
                                case 4:
                                    liste.push_back({ec.P0.add(cX),ec.P1.add(cX),cZ,{0,0,0}});
                                    liste.push_back({ec.P0.add(cY),ec.P1.add(cY),cZ,{0,0,0}});
                                    nx--;
                                    break;
                            }
                            if(ec.P0.toInt()!=i.toInt())
                                i = ec.P0;
                            else
                                i = ec.P1;
                            i0 = i.add(cX);
                            i1 = i.add(cY);
                            i2 = i.add(cZ);
                            e0 = {i,i0,cX,{0,0,0}};
                            e1 = {i,i1,cX,{0,0,0}};
                            e2 = {i,i2,cX,{0,0,0}};
                            e0.order();
                            e1.order();
                            e2.order();
                            liste.push_back(e0);
                            liste.push_back(e1);
                            liste.push_back(e2);
                            EdgeCube e = ec;
                            double min = -1;
                            for(int i=0; i<liste.size(); i++)
                            {
                                if(CubePlane.find(liste[i].toString())!=CubePlane.end())
                                {
                                    EdgeCube et = CubePlane[liste[i].toString()];
                                    double dist = distFromSurface(ec.inter,et.inter,target);
                                    if(min==-1 || dist<min)
                                    {
                                        e = et;
                                        min = dist;
                                    }
                                    
                                }
                            }
                            chx = e.toString();
                        }
                        glBegin(GL_LINE_LOOP);
                        for(int i=0; i<face.size(); i++)
                        {
                            Point P = face[i];
                            glVertex3f(P.x,P.y,P.z);
                        }
                        glEnd();
                    }

                    if(n101>0)
                    {
                        IndiceCube i = c101;
                        int* targ = &n101;
                        IndiceCube i0 = i.add(cX);
                        IndiceCube i1 = i.add(cY);
                        IndiceCube i2 = i.add(cZ);
                        EdgeCube e0 = {i,i0,cX,{0,0,0}};
                        EdgeCube e1 = {i,i1,cY,{0,0,0}};
                        EdgeCube e2 = {i,i2,cZ,{0,0,0}};
                        e0.order();
                        e1.order();
                        e2.order();
                        string chx = "";
                        if(CubePlane.find(e0.toString())!=CubePlane.end())
                        {
                            chx = e0.toString();
                            if(*targ!=2)
                            {
                                targ = &nx;
                            }
                        }
                        else
                        {
                            if(CubePlane.find(e1.toString())!=CubePlane.end())
                            {
                                chx = e1.toString();
                                if(*targ!=2)
                                {
                                    targ = &ny;
                                }
                            }
                            else
                            {
                                if(CubePlane.find(e2.toString())!=CubePlane.end())
                                {
                                    chx = e2.toString();
                                    if(*targ!=2)
                                    {
                                        targ = &nz;
                                    }
                                }
                            }
                        }
                        vector<Point> face;
                        int st = 0;
                        while(*targ!=0 && st<3)
                        {
                            st++;
                            cout<<n000<<"\t"<<n001<<"\t"<<n010<<"\t"<<n011<<"\t"<<endl;
                            cout<<n100<<"\t"<<n101<<"\t"<<n110<<"\t"<<n111<<"\t"<<endl;
                            EdgeCube ec = CubePlane[chx];
                            cout<<ec.toString()<<endl;
                            CubePlane.erase(ec.toString());
                            face.push_back(ec.inter);
                            int ip0 = ec.P0.toInt();
                            int ip1 = ec.P1.toInt();
                            int ipA = ec.coupe.toInt();
                            vector<EdgeCube> liste;
                            switch(ip0)
                            {
                                case 0:
                                    n000--;
                                    break;
                                case 1:
                                    n001--;
                                    break;
                                case 2:
                                    n010--;
                                    break;
                                case 3:
                                    n011--;
                                    break;
                                case 4:
                                    n100--;
                                    break;
                                case 5:
                                    n101--;
                                    break;
                                case 6:
                                    n110--;
                                    break;
                                case 7:
                                    n111--;
                            }
                            switch(ip1)
                            {
                                case 0:
                                    n000--;
                                    break;
                                case 1:
                                    n001--;
                                    break;
                                case 2:
                                    n010--;
                                    break;
                                case 3:
                                    n011--;
                                    break;
                                case 4:
                                    n100--;
                                    break;
                                case 5:
                                    n101--;
                                    break;
                                case 6:
                                    n110--;
                                    break;
                                case 7:
                                    n111--;
                            }
                            switch(ipA)
                            {
                                case 1:
                                    nz--;
                                    liste.push_back({ec.P0.add(cY),ec.P1.add(cY),cX,{0,0,0}});
                                    liste.push_back({ec.P0.add(cZ),ec.P1.add(cZ),cX,{0,0,0}});
                                    break;
                                case 2:
                                    liste.push_back({ec.P0.add(cX),ec.P1.add(cX),cY,{0,0,0}});
                                    liste.push_back({ec.P0.add(cZ),ec.P1.add(cZ),cY,{0,0,0}});
                                    ny--;
                                    break;
                                case 4:
                                    liste.push_back({ec.P0.add(cX),ec.P1.add(cX),cZ,{0,0,0}});
                                    liste.push_back({ec.P0.add(cY),ec.P1.add(cY),cZ,{0,0,0}});
                                    nx--;
                                    break;
                            }
                            if(ec.P0.toInt()!=i.toInt())
                                i = ec.P0;
                            else
                                i = ec.P1;
                            i0 = i.add(cX);
                            i1 = i.add(cY);
                            i2 = i.add(cZ);
                            e0 = {i,i0,cX,{0,0,0}};
                            e1 = {i,i1,cX,{0,0,0}};
                            e2 = {i,i2,cX,{0,0,0}};
                            e0.order();
                            e1.order();
                            e2.order();
                            liste.push_back(e0);
                            liste.push_back(e1);
                            liste.push_back(e2);
                            EdgeCube e = ec;
                            double min = -1;
                            for(int i=0; i<liste.size(); i++)
                            {
                                if(CubePlane.find(liste[i].toString())!=CubePlane.end())
                                {
                                    EdgeCube et = CubePlane[liste[i].toString()];
                                    double dist = distFromSurface(ec.inter,et.inter,target);
                                    if(min==-1 || dist<min)
                                    {
                                        e = et;
                                        min = dist;
                                    }
                                    
                                }
                            }
                            chx = e.toString();
                        }
                        glBegin(GL_LINE_LOOP);
                        for(int i=0; i<face.size(); i++)
                        {
                            Point P = face[i];
                            glVertex3f(P.x,P.y,P.z);
                        }
                        glEnd();
                    }

                    if(n110>0)
                    {
                        IndiceCube i = c110;
                        int* targ = &n110;
                        IndiceCube i0 = i.add(cX);
                        IndiceCube i1 = i.add(cY);
                        IndiceCube i2 = i.add(cZ);
                        EdgeCube e0 = {i,i0,cX,{0,0,0}};
                        EdgeCube e1 = {i,i1,cY,{0,0,0}};
                        EdgeCube e2 = {i,i2,cZ,{0,0,0}};
                        e0.order();
                        e1.order();
                        e2.order();
                        string chx = "";
                        if(CubePlane.find(e0.toString())!=CubePlane.end())
                        {
                            chx = e0.toString();
                            if(*targ!=2)
                            {
                                targ = &nx;
                            }
                        }
                        else
                        {
                            if(CubePlane.find(e1.toString())!=CubePlane.end())
                            {
                                chx = e1.toString();
                                if(*targ!=2)
                                {
                                    targ = &ny;
                                }
                            }
                            else
                            {
                                if(CubePlane.find(e2.toString())!=CubePlane.end())
                                {
                                    chx = e2.toString();
                                    if(*targ!=2)
                                    {
                                        targ = &nz;
                                    }
                                }
                            }
                        }
                        vector<Point> face;
                        int st = 0;
                        while(*targ!=0 && st<3)
                        {
                            st++;
                            cout<<n000<<"\t"<<n001<<"\t"<<n010<<"\t"<<n011<<"\t"<<endl;
                            cout<<n100<<"\t"<<n101<<"\t"<<n110<<"\t"<<n111<<"\t"<<endl;
                            EdgeCube ec = CubePlane[chx];
                            cout<<ec.toString()<<endl;
                            CubePlane.erase(ec.toString());
                            face.push_back(ec.inter);
                            int ip0 = ec.P0.toInt();
                            int ip1 = ec.P1.toInt();
                            int ipA = ec.coupe.toInt();
                            vector<EdgeCube> liste;
                            switch(ip0)
                            {
                                case 0:
                                    n000--;
                                    break;
                                case 1:
                                    n001--;
                                    break;
                                case 2:
                                    n010--;
                                    break;
                                case 3:
                                    n011--;
                                    break;
                                case 4:
                                    n100--;
                                    break;
                                case 5:
                                    n101--;
                                    break;
                                case 6:
                                    n110--;
                                    break;
                                case 7:
                                    n111--;
                            }
                            switch(ip1)
                            {
                                case 0:
                                    n000--;
                                    break;
                                case 1:
                                    n001--;
                                    break;
                                case 2:
                                    n010--;
                                    break;
                                case 3:
                                    n011--;
                                    break;
                                case 4:
                                    n100--;
                                    break;
                                case 5:
                                    n101--;
                                    break;
                                case 6:
                                    n110--;
                                    break;
                                case 7:
                                    n111--;
                            }
                            switch(ipA)
                            {
                                case 1:
                                    nz--;
                                    liste.push_back({ec.P0.add(cY),ec.P1.add(cY),cX,{0,0,0}});
                                    liste.push_back({ec.P0.add(cZ),ec.P1.add(cZ),cX,{0,0,0}});
                                    break;
                                case 2:
                                    liste.push_back({ec.P0.add(cX),ec.P1.add(cX),cY,{0,0,0}});
                                    liste.push_back({ec.P0.add(cZ),ec.P1.add(cZ),cY,{0,0,0}});
                                    ny--;
                                    break;
                                case 4:
                                    liste.push_back({ec.P0.add(cX),ec.P1.add(cX),cZ,{0,0,0}});
                                    liste.push_back({ec.P0.add(cY),ec.P1.add(cY),cZ,{0,0,0}});
                                    nx--;
                                    break;
                            }
                            if(ec.P0.toInt()!=i.toInt())
                                i = ec.P0;
                            else
                                i = ec.P1;
                            i0 = i.add(cX);
                            i1 = i.add(cY);
                            i2 = i.add(cZ);
                            e0 = {i,i0,cX,{0,0,0}};
                            e1 = {i,i1,cX,{0,0,0}};
                            e2 = {i,i2,cX,{0,0,0}};
                            e0.order();
                            e1.order();
                            e2.order();
                            liste.push_back(e0);
                            liste.push_back(e1);
                            liste.push_back(e2);
                            EdgeCube e = ec;
                            double min = -1;
                            for(int i=0; i<liste.size(); i++)
                            {
                                if(CubePlane.find(liste[i].toString())!=CubePlane.end())
                                {
                                    EdgeCube et = CubePlane[liste[i].toString()];
                                    double dist = distFromSurface(ec.inter,et.inter,target);
                                    if(min==-1 || dist<min)
                                    {
                                        e = et;
                                        min = dist;
                                    }
                                    
                                }
                            }
                            chx = e.toString();
                        }
                        glBegin(GL_LINE_LOOP);
                        for(int i=0; i<face.size(); i++)
                        {
                            Point P = face[i];
                            glVertex3f(P.x,P.y,P.z);
                        }
                        glEnd();
                    }

                    if(n111>0)
                    {
                        IndiceCube i = c111;
                        int* targ = &n111;
                        IndiceCube i0 = i.add(cX);
                        IndiceCube i1 = i.add(cY);
                        IndiceCube i2 = i.add(cZ);
                        EdgeCube e0 = {i,i0,cX,{0,0,0}};
                        EdgeCube e1 = {i,i1,cY,{0,0,0}};
                        EdgeCube e2 = {i,i2,cZ,{0,0,0}};
                        e0.order();
                        e1.order();
                        e2.order();
                        string chx = "";
                        if(CubePlane.find(e0.toString())!=CubePlane.end())
                        {
                            chx = e0.toString();
                            if(*targ!=2)
                            {
                                targ = &nx;
                            }
                        }
                        else
                        {
                            if(CubePlane.find(e1.toString())!=CubePlane.end())
                            {
                                chx = e1.toString();
                                if(*targ!=2)
                                {
                                    targ = &ny;
                                }
                            }
                            else
                            {
                                if(CubePlane.find(e2.toString())!=CubePlane.end())
                                {
                                    chx = e2.toString();
                                    if(*targ!=2)
                                    {
                                        targ = &nz;
                                    }
                                }
                            }
                        }
                        vector<Point> face;
                        int st = 0;
                        while(*targ!=0 && st<3)
                        {
                            st++;
                            cout<<n000<<"\t"<<n001<<"\t"<<n010<<"\t"<<n011<<"\t"<<endl;
                            cout<<n100<<"\t"<<n101<<"\t"<<n110<<"\t"<<n111<<"\t"<<endl;
                            EdgeCube ec = CubePlane[chx];
                            cout<<ec.toString()<<endl;
                            CubePlane.erase(ec.toString());
                            face.push_back(ec.inter);
                            int ip0 = ec.P0.toInt();
                            int ip1 = ec.P1.toInt();
                            int ipA = ec.coupe.toInt();
                            vector<EdgeCube> liste;
                            switch(ip0)
                            {
                                case 0:
                                    n000--;
                                    break;
                                case 1:
                                    n001--;
                                    break;
                                case 2:
                                    n010--;
                                    break;
                                case 3:
                                    n011--;
                                    break;
                                case 4:
                                    n100--;
                                    break;
                                case 5:
                                    n101--;
                                    break;
                                case 6:
                                    n110--;
                                    break;
                                case 7:
                                    n111--;
                            }
                            switch(ip1)
                            {
                                case 0:
                                    n000--;
                                    break;
                                case 1:
                                    n001--;
                                    break;
                                case 2:
                                    n010--;
                                    break;
                                case 3:
                                    n011--;
                                    break;
                                case 4:
                                    n100--;
                                    break;
                                case 5:
                                    n101--;
                                    break;
                                case 6:
                                    n110--;
                                    break;
                                case 7:
                                    n111--;
                            }
                            switch(ipA)
                            {
                                case 1:
                                    nz--;
                                    liste.push_back({ec.P0.add(cY),ec.P1.add(cY),cX,{0,0,0}});
                                    liste.push_back({ec.P0.add(cZ),ec.P1.add(cZ),cX,{0,0,0}});
                                    break;
                                case 2:
                                    liste.push_back({ec.P0.add(cX),ec.P1.add(cX),cY,{0,0,0}});
                                    liste.push_back({ec.P0.add(cZ),ec.P1.add(cZ),cY,{0,0,0}});
                                    ny--;
                                    break;
                                case 4:
                                    liste.push_back({ec.P0.add(cX),ec.P1.add(cX),cZ,{0,0,0}});
                                    liste.push_back({ec.P0.add(cY),ec.P1.add(cY),cZ,{0,0,0}});
                                    nx--;
                                    break;
                            }
                            if(ec.P0.toInt()!=i.toInt())
                                i = ec.P0;
                            else
                                i = ec.P1;
                            i0 = i.add(cX);
                            i1 = i.add(cY);
                            i2 = i.add(cZ);
                            e0 = {i,i0,cX,{0,0,0}};
                            e1 = {i,i1,cX,{0,0,0}};
                            e2 = {i,i2,cX,{0,0,0}};
                            e0.order();
                            e1.order();
                            e2.order();
                            liste.push_back(e0);
                            liste.push_back(e1);
                            liste.push_back(e2);
                            EdgeCube e = ec;
                            double min = -1;
                            for(int i=0; i<liste.size(); i++)
                            {
                                if(CubePlane.find(liste[i].toString())!=CubePlane.end())
                                {
                                    EdgeCube et = CubePlane[liste[i].toString()];
                                    double dist = distFromSurface(ec.inter,et.inter,target);
                                    if(min==-1 || dist<min)
                                    {
                                        e = et;
                                        min = dist;
                                    }
                                    
                                }
                            }
                            chx = e.toString();
                        }
                        glBegin(GL_LINE_LOOP);
                        for(int i=0; i<face.size(); i++)
                        {
                            Point P = face[i];
                            glVertex3f(P.x,P.y,P.z);
                        }
                        glEnd();
                    }

                    pile.push_back({Pi.x+1,Pi.y,Pi.z});
                    pile.push_back({Pi.x-1,Pi.y,Pi.z});
                    pile.push_back({Pi.x,Pi.y+1,Pi.z});
                    pile.push_back({Pi.x,Pi.y-1,Pi.z});
                    pile.push_back({Pi.x,Pi.y,Pi.z+1});
                    pile.push_back({Pi.x,Pi.y,Pi.z-1});
                    
                }
            }
            
        }
        pile.erase(pile.begin());
        st++;
    }
    cout<<"end"<<endl;
}

//------------------------------------------------------
void affichage(void)
{
    glMatrixMode(GL_MODELVIEW);
    /* effacement de l'image avec la couleur de fond */
    glClear(GL_COLOR_BUFFER_BIT);
    glPushMatrix();
    glTranslatef(0,0,cameraDistance);
    glRotatef(cameraAngleX,1.,0.,0.)	;
    glRotatef(cameraAngleY,0.,1.,0.);
    //affiche_repere();
    //dessinZone();
    //marchingCube({0,0,0},0.1,0.5);
    
    Config c;
    c.max = 8;
    c = c.add(nb);
    Point P[8];
    P[0]={-0.5,-0.5,-0.5};
    P[1]={-0.5,-0.5,0.5};
    P[2]={-0.5,0.5,-0.5};
    P[3]={-0.5,0.5,0.5};
    P[4]={0.5,-0.5,-0.5};
    P[5]={0.5,-0.5,0.5};
    P[6]={0.5,0.5,-0.5};
    P[7]={0.5,0.5,0.5};
    double v[8];
    for(int i=0; i<8; i++)
    {
        int x = 1;
        if(c.p.size()>7-i && c.p[7-i])
        {
            x=-1;
        }
        v[i]=0.5*x;
    }
    cout<<c.toString()<<endl;
    glPointSize(5.0f);
    glBegin(GL_POINTS);
    for(int i=0; i<8; i++)
    {
        if(v[i]>0)
        {
            glColor3f(0,1,0);
        }
        else
        {
            glColor3f(1,0,0);
        }
        glVertex3f(P[i].x,P[i].y,P[i].z);
    }
    glEnd();
    if(guide.find(c.toString())==guide.end())
        c.inv();
    cout<<c.toString()<<endl;
    Forme f = guide[c.toString()];
    cout<<f.toString()<<endl;
    vector<Zone> vz = f.dessin(P,v);
    for(int i=0; i<vz.size(); i++)
    {
        Zone z = vz[i];
        for(int j=0; j<z.faces.size(); j++)
        {
            Triangle t = z.faces[j];
            glColor3f(1,1,1);
            glBegin(GL_TRIANGLES);
            glVertex3f(t.P0.x,t.P0.y,t.P0.z);
            glVertex3f(t.P1.x,t.P1.y,t.P1.z);
            glVertex3f(t.P2.x,t.P2.y,t.P2.z);
            glEnd();
            glColor3f(0,0,0);
            glBegin(GL_LINE_LOOP);
            glVertex3f(t.P0.x,t.P0.y,t.P0.z);
            glVertex3f(t.P1.x,t.P1.y,t.P1.z);
            glVertex3f(t.P2.x,t.P2.y,t.P2.z);
            glEnd();
        }
    }
    /*
    do{
        cout<<c.toString()<<endl;
        c = c.next();
    }while(c.val()!=0);
    */
    
    glPopMatrix();
    /* on force l'affichage du resultat */
    glFlush();
}

//------------------------------------------------------


//------------------------------------------------------
void clavier(unsigned char touche,int x,int y)
{

  switch (touche)
    {
    case '+': //* affichage du carre plein 
        nb++;
      glutPostRedisplay();
      break;
    case '-': //* affichage du carre plein 
        nb--;
        if(nb<0) nb=0;
      glutPostRedisplay();
      break;
    case 'f': //* affichage en mode fil de fer 
      glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
      glutPostRedisplay();
      break;
      case 'p': //* affichage du carre plein 
      glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
      glutPostRedisplay();
      break;
  case 's' : //* Affichage en mode sommets seuls 
      glPolygonMode(GL_FRONT_AND_BACK,GL_POINT);
      glutPostRedisplay();
      break;
    case 'q' : //*la touche 'q' permet de quitter le programme 
      exit(0);
    }
    
}
void mouse(int button, int state, int x, int y)
{
    mouseX = x;
    mouseY = y;

    if(button == GLUT_LEFT_BUTTON)
    {
        if(state == GLUT_DOWN)
        {
            mouseLeftDown = true;
        }
        else if(state == GLUT_UP)
            mouseLeftDown = false;
    }

    else if(button == GLUT_RIGHT_BUTTON)
    {
        if(state == GLUT_DOWN)
        {
            mouseRightDown = true;
        }
        else if(state == GLUT_UP)
            mouseRightDown = false;
    }

    else if(button == GLUT_MIDDLE_BUTTON)
    {
        if(state == GLUT_DOWN)
        {
            mouseMiddleDown = true;
        }
        else if(state == GLUT_UP)
            mouseMiddleDown = false;
    }
}


void mouseMotion(int x, int y)
{
    if(mouseLeftDown)
    {
        cameraAngleY += (x - mouseX);
        cameraAngleX += (y - mouseY);
        mouseX = x;
        mouseY = y;
    }
    if(mouseRightDown)
    {
        cameraDistance += (y - mouseY) * 0.2f;
        mouseY = y;
    }

    glutPostRedisplay();
}

    
    
