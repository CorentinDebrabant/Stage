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
void setSize();

map<string,Forme> guide;
map<string,bool> marque;
vector<Point> pile;
map<string,EdgeCube> CubePlane;
vector<double> targets;

vector<vector<vector<double>>> champ;
int resoC = 10;

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
double taille = 2;
int fonction = 0;

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
    double v = 0;
    switch(fonction)
    {
        case 1:
            {
                if(champ.size()==0)
                {
                    for(int i=-resoC; i<=resoC; i++)
                    {
                        vector<vector<double>> cy;
                        for(int j=-resoC; j<=resoC; j++)
                        {
                            vector<double> cz;
                            for(int k=-resoC; k<=resoC; k++)
                            {
                                if(k==-resoC && j==-resoC && i==-resoC)
                                {
                                    double x = ((double)(rand()) / (double)(RAND_MAX));
                                    cz.push_back(x);
                                }
                                else
                                {
                                    vector<double> val;
                                    if(i!=-resoC)
                                    {
                                        val.push_back(champ[resoC+i-1][resoC+j][resoC+k]);
                                        if(j!=-resoC)
                                        {
                                            val.push_back(champ[resoC+i-1][resoC+j-1][resoC+k]);
                                            if(k!=-resoC)
                                            {
                                                val.push_back(champ[resoC+i-1][resoC+j-1][resoC+k-1]);
                                            }
                                        }
                                        if(k!=-resoC)
                                        {
                                            val.push_back(champ[resoC+i-1][resoC+j][resoC+k-1]);
                                        }
                                    }
                                    
                                    if(j!=-resoC)
                                    {
                                        val.push_back(cy[cy.size()-1][resoC+k]);
                                        if(k!=-resoC)
                                        {
                                            val.push_back(cy[cy.size()-1][resoC+k-1]);
                                        }
                                    }
                                    
                                    if(k!=-resoC)
                                    {
                                        val.push_back(cz[cz.size()-1]);
                                    }
                                    
                                    double avg = 0;
                                    int nb = 0;
                                    for(int b=0; b<val.size(); b++)
                                    {
                                        nb++;
                                        avg+=val[b];
                                    }
                                    avg/=(double)nb;
                                    double x = ((double)(rand()) / (double)(RAND_MAX))/5-0.1;
                                    cz.push_back(avg+x);
                                }
                            }
                            cy.push_back(cz);
                        }
                        champ.push_back(cy);
                    }
                }
                double p = (taille/(resoC*2));
                int x = resoC + (int)floor(P.x/p);
                int y = resoC + (int)floor(P.y/p);
                int z = resoC + (int)floor(P.z/p);
                if(x==resoC*2)
                    x=resoC*2-1;
                if(y==resoC*2)
                    y=resoC*2-1;
                if(z==resoC*2)
                    z=resoC*2-1;
                double v0 = champ[x][y][z];
                double v1 = champ[x][y][z+1];
                double v2 = champ[x][y+1][z];
                double v3 = champ[x][y+1][z+1];
                double v4 = champ[x+1][y][z];
                double v5 = champ[x+1][y][z+1];
                double v6 = champ[x+1][y+1][z];
                double v7 = champ[x+1][y+1][z+1];
                Point S = {(x-resoC)*p,(y-resoC)*p,(z-resoC)*p};
                double xd = (P.x - S.x) / (2*p);
                double yd = (P.y - S.y) / (2*p);
                double zd = (P.z - S.z) / (2*p);
                double v00 = xd*v4 + (1-xd)*v0;
                double v01 = xd*v5 + (1-xd)*v1;
                double v11 = xd*v7 + (1-xd)*v3;
                double v10 = xd*v6 + (1-xd)*v2;
                double vr0 = yd*v10 + (1-yd)*v00;
                double vr1 = yd*v11 + (1-yd)*v01;
                v = zd*vr1 + (1-zd)*vr0;
            }
            break;
        default:
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
                    cout<<va<<" "<<vb<<" "<<vc<<" "<<vd<<" "<<ve<<" "<<vf<<" "<<vg<<" "<<vh<<endl;
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
                v = zd*v1 + (1-zd)*v0;
            }
            break;
    }
    return v;
}

double func(Point P)
{
    if(abs(P.x)>1 || abs(P.y)>1 || abs(P.z)>1)
        return -1;
    return val(P);
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
    float v = 0;
    glPointSize(10.0f);
    glBegin(GL_POINTS);
    while(v<=1.01)
    {
        Co = hsv2rgb(v,1,1);
        glColor3f(Co.x,Co.y,Co.z);
        glVertex3f(-1.5,0,v);
        v+=0.1;
    }
    glEnd();
    glPointSize(1.0f);
    double nb = 3;
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

void marchingCube(Point origine, double pas, double target/*, vector<vector<Point>> *C*/)
{
    //cout<<"Start "<<pas<<" "<<target<<" Départ en ("<<origine.x<<", "<<origine.y<<", "<<origine.z<<") et valeur = "<<func(origine)<<endl;
    pile.push_back({0,0,0});
    int nbP = 5;
    int np = (taille/pas)/nbP;
    bool plop = false;
    int st = 0;
    while(pile.size()!=0)
    {
        if(pile.size()>10 && !plop)
        {
            plop = true;
            
        }
        //cout<<"Pile "<<pile.size()<<endl;
        Point Pi = pile[0];
        string k = to_string(Pi.x)+":"+to_string(Pi.y)+":"+to_string(Pi.z);
        if(marque.find(k) == marque.end())
        {
            marque[k]=true;
            Point P = {origine.x+pas*Pi.x, origine.y+pas*Pi.y, origine.z+pas*Pi.z};
            //cout<<func(P)<<endl;
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
                //cout<<"x"<<endl;
                float c = (double)st/100.0;
                Point C = hsv2rgb(c,1,1);
                //glColor3f(C.x,C.y,C.z);
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
                Config g;
                g.max = 8;
                if(v000>0)
                    g=g.inv(7);
                if(v001>0)
                    g=g.inv(6);
                if(v010>0)
                    g=g.inv(5);
                if(v011>0)
                    g=g.inv(4);
                if(v100>0)
                    g=g.inv(3);
                if(v101>0)
                    g=g.inv(2);
                if(v110>0)
                    g=g.inv(1);
                if(v111>0)
                    g=g.inv(0);
                if(guide.find(g.toString())==guide.end())
                    g.inv();
                Forme f = guide[g.toString()];
                if(f.type==0)
                {
                    //cout<<"move"<<endl;
                    /*glColor3f(0,1,0);
                    glPointSize(3.0f);
                    glBegin(GL_POINTS);
                    glVertex3f(P.x,P.y,P.z);
                    glEnd();*/
                    double v0 = abs(v000);
                    double v1 = abs(v001);
                    double v2 = abs(v010);
                    double v3 = abs(v011);
                    double v4 = abs(v100);
                    double v5 = abs(v101);
                    double v6 = abs(v110);
                    double v7 = abs(v111);
                    double m = min(min(min(v0,v1),min(v2,v3)),min(min(v4,v5),min(v6,v7)));
                    int ra = rand()%15;
                    ra = 7;
                    if(v0==m || v1==m || v2==m || v3==m || ra==0 || ra>=13)
                    {
                        pile.push_back({Pi.x-1,Pi.y,Pi.z});
                    }
                    if(v4==m || v5==m || v6==m || v7==m || ra==1 || ra>=13)
                    {
                        pile.push_back({Pi.x+1,Pi.y,Pi.z});
                    }
                    if(v0==m || v2==m || v4==m || v6==m || ra==2 || ra>=13)
                    {
                        pile.push_back({Pi.x,Pi.y,Pi.z-1});
                    }
                    if(v1==m || v3==m || v5==m || v7==m || ra==3 || ra>=13)
                    {
                        pile.push_back({Pi.x,Pi.y,Pi.z+1});
                    }
                    if(v0==m || v1==m || v4==m || v6==m || ra==4 || ra>=13)
                    {
                        pile.push_back({Pi.x,Pi.y-1,Pi.z});
                    }
                    if(v2==m || v3==m || v5==m || v7==m || ra==5 || ra>=13)
                    {
                        pile.push_back({Pi.x,Pi.y+1,Pi.z});
                    }
                }
                else
                {
                    //cout<<"dessin"<<endl;
                    /*if((int)(Pi.x+Pi.y+Pi.z)%np==0)
                    {
                        glPointSize(10.0f);
                        glColor3f(1,1,1);
                        glBegin(GL_POINTS);
                        glVertex3f(P.x,P.y,P.z);
                        glEnd();
                    }*/
                    /*glColor3f(1,0,0);
                    glPointSize(3.0f);
                    glBegin(GL_POINTS);
                    glVertex3f(P.x,P.y,P.z);
                    glEnd();*/
                    Point P[8] = {P000,P001,P010,P011,P100,P101,P110,P111};
                    double v[8]= {v000,v001,v010,v011,v100,v101,v110,v111};
                    vector<Zone> vz = f.dessin(P,v);
                    for(int i=0; i<vz.size(); i++)
                    {
                        Zone z = vz[i];
                        Point Co = hsv2rgb(target,1,1);
                        glColor3f(Co.x,Co.y,Co.z);
                        for(int j=0; j<z.faces.size(); j++)
                        {
                            Triangle t = z.faces[j];
                            //glColor3f(1,1,1);
                            /*glBegin(GL_TRIANGLES);
                            glVertex3f(t.P0.x,t.P0.y,t.P0.z);
                            glVertex3f(t.P1.x,t.P1.y,t.P1.z);
                            glVertex3f(t.P2.x,t.P2.y,t.P2.z);
                            glEnd();*/
                            
                            //glColor3f(0,0,0);
                            glBegin(GL_LINE_LOOP);
                            glVertex3f(t.P0.x,t.P0.y,t.P0.z);
                            glVertex3f(t.P1.x,t.P1.y,t.P1.z);
                            glVertex3f(t.P2.x,t.P2.y,t.P2.z);
                            glEnd();
                        }
                    }
                    pile.push_back({Pi.x-1,Pi.y,Pi.z});
                    pile.push_back({Pi.x+1,Pi.y,Pi.z});
                    pile.push_back({Pi.x,Pi.y,Pi.z-1});
                    pile.push_back({Pi.x,Pi.y,Pi.z+1});
                    pile.push_back({Pi.x,Pi.y-1,Pi.z});
                    pile.push_back({Pi.x,Pi.y+1,Pi.z});
                }
            }
            
        }
        pile.erase(pile.begin());
        st++;
    }
    //cout<<"Fin"<<endl;
}

void lancerMarchingCube(double pas, double target)
{
    marque.clear();
    
    marchingCube({0,0,0},pas,target);
    bool stop = false;
    Point P[8] = {{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0}};
    bool bp[8] = {false,false,false,false,false,false,false,false};
    while(!stop)
    {
        int n = 0;
        for(int i=0; i<8; i++)
        {
            Point Pc = P[i];
            Config c;
            c.max = 3;
            c.inv();
            c.inv();
            c = c.add(i);
            if(func(Pc)!=-1 && !bp[i])
            {
                if(c.p[0])
                    Pc.x=pas;
                else
                    Pc.x-=pas;
                if(c.p[1])
                    Pc.y+=pas;
                else
                    Pc.y-=pas;
                if(c.p[2])
                    Pc.z+=pas;
                else
                    Pc.z-=pas;
            }
            if(func(Pc)==-1)
            {
                n++;
                if(!bp[i])
                {
                    bp[i]=true;
                    if(c.p[0])
                        Pc.x-=pas*2;
                    else
                        Pc.x+=pas*2;
                    if(c.p[1])
                        Pc.y-=pas*2;
                    else
                        Pc.y+=pas*2;
                    if(c.p[2])
                        Pc.z-=pas*2;
                    else
                        Pc.z+=pas*2;
                }
            }
            P[i]=Pc;        
        }
        if(n==8)
            stop=true;
    }
    for(int i=0; i<8; i++)
    {
        marchingCube(P[i],pas,target);
    }
}

//------------------------------------------------------
void affichage(void)
{
    targets.clear();
    glMatrixMode(GL_MODELVIEW);
    /* effacement de l'image avec la couleur de fond */
    glClear(GL_COLOR_BUFFER_BIT);
    glPushMatrix();
    glTranslatef(0,0,cameraDistance);
    glRotatef(cameraAngleX,1.,0.,0.)	;
    glRotatef(cameraAngleY,0.,1.,0.);
    //affiche_repere();
    dessinZone();
    /*targets.push_back(0.1);
    targets.push_back(0.2);
    targets.push_back(0.3);
    targets.push_back(0.4);
    targets.push_back(0.5);
    targets.push_back(0.6);
    targets.push_back(0.7);
    targets.push_back(0.8);
    targets.push_back(0.9);*/
    targets.push_back(0.5);
    targets.push_back(0.6);
    for(int i=0; i<targets.size(); i++)
    {
        
        double v = targets[i];
        lancerMarchingCube(0.1,v);
    }
    /*Config c;
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
    Config co;
    co.max = 3;
    co = co.add(0);
    for(int i=0; i<8; i++)
    {
        Point C = {0,0,0};
        switch(i)
        {
            case 1:
                C = {0,0,0.5};
                break;
            case 2:
                C = {0,0.5,0};
                break;
            case 3:
                C = {0,0.5,0.5};
                break;
            case 4:
                C = {0.5,0,0};
                break;
            case 5:
                C = {0.5,0,0.5};
                break;
            case 6:
                C = {0.5,0.5,0};
                break;
            case 7:
                C = {0.5,0.5,0.5};
                break;
        }
        if(v[i]>0)
        {
            glColor3f(C.x,C.y,C.z);
        }
        else
        {
            glColor3f(C.x*2,C.y*2,C.z*2);
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
    }*/
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

    
    
