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
#include <algorithm>

using namespace std;
using namespace chrono;

void affichage(void);
ofstream debug;
    
void clavier(unsigned char touche,int x,int y);
void affiche_repere(void);
void initMap();
void mouse(int, int, int, int);
void mouseMotion(int, int);

//void reshape(int,int);

double normeMin = .001;
double normeMax = .005;

float color=0;

int nbF = -1;
int nbE = 0;

struct Point{
    double x;
    double y;
    double z;
};

struct PointCourbe{
    Point P;
    float val;
    int isoC;
    int isoInd;
    int graC;
    int indGra;
};


struct PointReturn{
    Point P;
    bool isOk;
};

struct Repere{
    Point P;
    Point N;
    Point T;
};

struct Face{
    vector<PointCourbe> sommets;
};

struct half_edge{
    int id;
    Point incident;
    float densityI;
    Point controleI;
    float densityCI;
    Point controleO;
    float densityCO;
    Point origine;
    float densityO;
    int face;
    half_edge* next;
    half_edge* previous;
    half_edge* opposite;
};

struct heFace{
    int id;
    float density;
    half_edge* incidente;
};

struct RetourDensite{
    int nb;
    float somme;
};

Point addP(Point P1, Point P2, double mult);
double getAngle(Point V0, Point V1);

map<string,bool> marque;
vector<Point> pile;

half_edge nullHe = {-1,{0,0,0},0,{0,0,0},0,{0,0,0},0,{0,0,0},0,-1,nullptr,nullptr,nullptr};

map<string,half_edge*> hes;
vector<heFace> Faces;

vector<float> targets;
map<int,vector<vector<Point>>> courbes;
map<int,vector<vector<PointCourbe>>> courbeFinale;
map<string,int> msTable;
vector<vector<PointCourbe>> courbeGrad;
//vector<PointFunc> listeAlea;

vector<Face> listeFace;

float msMoy = 0;
int nbMs = 0;
float cbMoy = 0;
int nbCb = 0;

float moy = 0;
int ng=0;
float moy2=0;
int n2=0;

float pX = 0.2;
float pY = 0.2;

bool dessinne = true;

int nb = 0;
//vector<Point> P = {};

struct Square{
    Point P;
    int type;
    int courbe;
    bool fini = false;
    double pas;
};

map<string,bool> Passage;

float a1 = -1;
float b1 = -1;
float c1 = -1;
float d1 = -1;

static int nbi=5 ;
static int nb_transfo=0 ;

const float xMin = -1;
const float xMax = 1;
const float yMin = -1;
const float yMax = 1;

vector<Point> surface;
vector<int> pointsImportants;

const float limite = 0.0001;
const float x_size = xMax-xMin;
const float y_size = yMax-yMin;

// variables globales pour OpenGL
bool mouseLeftDown;
bool mouseRightDown;
bool mouseMiddleDown;
float mouseX, mouseY;
float cameraAngleX;
float cameraAngleY;
float cameraDistance=0;
float decX = 0;
float decY = 0;
int dir = -1;

void idle()
{
    cout<<pX<<" "<<pY<<endl;
    if(abs(pY)>0.0051)
        dir=-dir;
    pX+=(0.00001*dir);
    pY+=(0.00001*dir);
    glutPostRedisplay();
}

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
  //glutIdleFunc(idle);

// initialisation des tranformations

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	 gluPerspective(60.0f,(GLfloat)200/(GLfloat)200,0.0001f,100.0f);
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

//-----------------------------------------------------

Point hsv2rgb(int h,float s, float v)
{
    Point C = {0,0,0};
    /*h = h%360;
    float c = s*v;
    float x = c * ( 1 -abs((h/60)%2 - 1));
    float m = v-c;
    cout<<c<<" "<<x<<" "<<m<<" "<<h<<endl;
    if(h<60)
    {
        C = {c,x,0};
    }
    else
    {
        if(h<120)
        {
            C = {x,c,0};
        }
        else
        {
            if(h<180)
            {
                C = {0,c,x};
            }
            else
            {
                if(h<240)
                {
                    C = {0,x,c};
                }
                else
                {
                    if(h<300)
                    {
                        C = {x,0,c};
                    }
                    else
                    {
                        C = {c,0,x};
                    }
                }
            }
        }
    }*/
    int ti = (int)floor(h/60)%6;
    float f = (float)h/60 - ti;
    float l = v * ( 1 - s);
    float m = v * (1 - f * s);
    float n = v * (1 - (1-f) * s);
    switch(ti)
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

float produitScalaire(Point V1, Point V2)
{
    return V1.x*V2.x+V1.y*V2.y+V1.z*V2.z;
}

float interpolation(float a, float b)
{
    float aa = abs(a);
    float ab = abs(b);
    if(aa+ab!=0)
        return ab/(aa+ab);
    return 0.5;
}

float dist2(Point P1, Point P2)
{
    float x = P2.x - P1.x;
    float y = P2.y - P1.y;
    float z = P2.z - P1.z;
    return sqrt(x*x+y*y+z*z);
}

float dist(float x, float y, float z)
{
    return sqrt(x*x+y*y+z*z);
}

double val(double x, double y, double z)
{
    /*
    if(a1==-1)
    {
        a1 = (float)(rand()) / (float)(RAND_MAX);
    }
    if(b1==-1)
    {
        b1 = (float)(rand()) / (float)(RAND_MAX);
    }
    if(c1==-1)
    {
        c1 = (float)(rand()) / (float)(RAND_MAX);
    }
    if(d1==-1)
    {
        d1 = (float)(rand()) / (float)(RAND_MAX);
    }
    if(a1>=0.2 && b1>=0.2 && c1>=0.2 && d1>=0.2)
    {
        a1 = 0;
    }
    if(a1<=0.8 && b1<=0.8 && c1<=0.8 && d1<=0.8)
    {
        c1 = 1;
    }
    Point P = {x,y,z};
    Point A = {-1,-1,0};
    Point B = {1,-1,0};
    Point C = {1,1,0};
    Point D = {-1,1,0};
    float Ap = dist2(A,P);
    float Bp = dist2(B,P);
    float Cp = dist2(C,P);
    float Dp = dist2(D,P);
    float iac = interpolation(Ap,Cp);
    float ibd = interpolation(Bp,Dp);
    float v1 = a1*iac + (1-iac)*c1;
    float v2 = b1*ibd + (1-ibd)*d1;
    return (v1+v2)/2;
    */
    
    double v = cos(2*x+3*y+0.5)+cos(5*x+2*y-0.1);
    v/=4;
    v+=0.5;
    return v;
    
}

double func(float x, float y, float z)
{
    if(x>=-1 && y>=-1 && x<=1 && y<=1)
        return val(x,y,z);
    else
        return -1;
    //return (sqrt(2)-dist(x,y,z))/sqrt(2);
    /*
    if(abs(y)<1 && abs(x)<=1)
        return max(abs(x),0.01f);
    return 0;
    */
    /*
    if(listeAlea.size()==0)
    {
        bool isZero = false;
        bool isOne = false;
        Point A = {-1,-1,0};
        Point B = {1,-1,0};
        Point C = {1,1,0};
        Point D = {-1,1,0};
        a1 = (float)(rand()) / (float)(RAND_MAX);
        b1 = (float)(rand()) / (float)(RAND_MAX);
        c1 = (float)(rand()) / (float)(RAND_MAX);
        d1 = (float)(rand()) / (float)(RAND_MAX);
        listeAlea.push_back({A,a1});
        listeAlea.push_back({B,b1});
        listeAlea.push_back({C,c1});
        listeAlea.push_back({D,d1});
        if(a1<=0.2 || b1<=0.2 || c1<=0.2 || d1<=0.2)
            isZero = true;
        if(a1>=0.8 || b1>=0.8 || c1>=0.8 || d1>=0.8)
            isOne = true;
        while(!isZero || !isOne || listeAlea.size()<6)
        {
            float x = 2* (float)(rand()) / (float)(RAND_MAX) -1;
            float y = 2* (float)(rand()) / (float)(RAND_MAX) -1;
            Point P = {x,y,0};
            float v = (float)(rand()) / (float)(RAND_MAX);
            listeAlea.push_back({P,v});
            if(v<=0.2)
                isZero = true;
            if(v>=0.8)
                isOne = true;
        }
        cout<<listeAlea.size()<<endl;
    }
    Point P = {x,y,z};
    if(x>=-1 && y>=-1 && x<=1 && y<=1)
    {
        float maxDist = 0;
        for(int i=0; i<listeAlea.size(); i++)
        {
            float d = dist2(P,listeAlea[i].P);
            if(d>maxDist)
                maxDist =d;
        }
        float sumDistInv = 0;
        float sumVal = 0;
        for(int i=0; i<listeAlea.size(); i++)
        {
            PointFunc Pf = listeAlea[i];
            float d = maxDist-dist2(P,Pf.P);
            sumVal += d * Pf.val;
            sumDistInv+= d;
        }
        return sumVal/sumDistInv;
    }
    return 0;
    */
    
    
    /*
    Point P = {x,y,z};
    Point A = {-1,-1,0};
    Point B = {1,-1,0};
    Point C = {1,1,0};
    Point D = {-1,1,0};
    float Ap = dist2(A,P);
    float Bp = dist2(B,P);
    float Cp = dist2(C,P);
    float Dp = dist2(D,P);
    if(x>=-1 && y>=-1 && x<=1 && y<=1)
        return min(min(Ap,Bp),min(Cp,Dp));
    return 0;*/
}

void initMap()
{
    msTable["0000"]=0;
    msTable["0001"]=1;
    msTable["0010"]=2;
    msTable["0011"]=3;
    msTable["0100"]=4;
    msTable["0101"]=5;
    msTable["0110"]=6;
    msTable["0111"]=7;
    msTable["1000"]=7;
    msTable["1001"]=6;
    msTable["1010"]=5;
    msTable["1011"]=4;
    msTable["1100"]=3;
    msTable["1101"]=2;
    msTable["1110"]=1;
    msTable["1111"]=0;
}

PointReturn findPoint(Point P, float target, int ind, float pas)
{
    float v0 = target-func(P.x,P.y,P.z);
    if(abs(v0)<limite)
    {
        return {P,true};
    }
    if(ind>3000)
    {
        return {P,false};
    }
    float v1 = target-func(P.x-pas,P.y,P.z);
    float v2 = target-func(P.x,P.y-pas,P.z);
    float v3 = target-func(P.x+pas,P.y,P.z);
    float v4 = target-func(P.x,P.y+pas,P.z);
    int s0 = v0/abs(v0);
    int s1 = v1/abs(v1);
    int s2 = v2/abs(v2);
    int s3 = v3/abs(v3);
    int s4 = v4/abs(v4);
    if(s0==s1 && s1==s2 && s2==s3 && s3==s4)
    {
        int dir = -1;
        while(true)
        {
            dir = rand()%4;
            int r = rand()%20;
            switch(dir)
            {
                case 0:
                    if(abs(v1)<abs(v0) || r==0)
                    {
                        Point nP = {P.x-pas,P.y,P.z};
                        return findPoint(nP,target,ind+1,pas);
                    }
                    break;
                case 1:
                    if(abs(v2)<abs(v0) || r==0)
                    {
                        Point nP = {P.x,P.y-pas,P.z};
                        return findPoint(nP,target,ind+1,pas);
                    }
                    break;
                case 2:
                    if(abs(v3)<abs(v0) || r==0)
                    {
                        Point nP = {P.x+pas,P.y,P.z};
                        return findPoint(nP,target,ind+1,pas);
                    }
                    break;
                case 3:
                    if(abs(v4)<abs(v0) || r==0)
                    {
                        Point nP = {P.x,P.y+pas,P.z};
                        return findPoint(nP,target,ind+1,pas);
                    }
                    break;
            }
        }
    }
    else
    {
        Point sum = {0,0,0};
        float n = 0;
        if(s0!=s1)
        {
            float inter = interpolation(v0,v1);
            sum.x+=inter*P.x + (1-inter)*(P.x-pas);
            sum.y+=inter*P.y + (1-inter)*(P.y);
            sum.z+=inter*P.z + (1-inter)*(P.z);
            n++;
        }
        if(s0!=s2)
        {
            float inter = interpolation(v0,v2);
            sum.x+=inter*P.x + (1-inter)*(P.x);
            sum.y+=inter*P.y + (1-inter)*(P.y-pas);
            sum.z+=inter*P.z + (1-inter)*(P.z);
            n++;
        }
        if(s0!=s3)
        {
            float inter = interpolation(v0,v3);
            sum.x+=inter*P.x + (1-inter)*(P.x+pas);
            sum.y+=inter*P.y + (1-inter)*(P.y);
            sum.z+=inter*P.z + (1-inter)*(P.z);
            n++;
        }
        if(s0!=s4)
        {
            float inter = interpolation(v0,v4);
            sum.x+=inter*P.x + (1-inter)*(P.x);
            sum.y+=inter*P.y + (1-inter)*(P.y+pas);
            sum.z+=inter*P.z + (1-inter)*(P.z);
            n++;
        }
        sum.x/=n;
        sum.y/=n;
        sum.z/=n;
        return findPoint(sum, target, ind+1, pas/2.0);
    }
}

/*Point findGradient(Point P)
{
    auto tstart= system_clock::now();
    float nb = 3;
    float start = 0;
    float end = M_PI*2;
    float pas = (end-start)/nb;
    int step = 0;
    float ori = func(P.x,P.y,P.z);
    float retVal;
    float compMax = 0;
    float ind = 0;
    while(step<200)
    {
        float max = 0;
        int id = -1;
        for(int i=0; i<=nb; i++)
        {
            float angle = start + pas*i;
            Point cs = {0.0001*cos(angle),0.0001*sin(angle),0};
            float val = func(P.x+cs.x,P.y+cs.y,P.z);
            if(abs(ori-val)>max)
            {
                max=abs(ori-val);
                id=i;
            }
            ind++;
        }
        step++;
        retVal = start + pas * id;
        end = start + pas * (id+1);
        start = start + pas * (id-1);
        pas = (end-start)/nb;
        if(abs(compMax-max)<limite*0.0001)
        {
            break;
        }
        compMax=max;
    }
    Point R1 = {cos(retVal),sin(retVal),0};
    auto tend= system_clock::now();
    auto elapsed = duration_cast<nanoseconds>(tend - tstart);
    //cout<<"Gradient : "<<elapsed.count()<<endl;
    moy+=elapsed.count();
    ng++;
    return R1;
}*/

Repere getRepere(Point P)
{
    double dtO = 0.00000001;
    double dt = dtO;
    double target = val(P.x,P.y,P.z);
    /*
    float p1 = target - func(P.x-dt,P.y-dt,P.z);
    float p2 = target - func(P.x+dt,P.y-dt,P.z);
    float p3 = target - func(P.x+dt,P.y+dt,P.z);
    float p4 = target - func(P.x-dt,P.y+dt,P.z);
    if(p1==target || p2==target || p3==target || p4==target)
    {
        target = 0;
        dt*=2;
    }
    */
    double p1 = val(P.x+dt,P.y,P.z);
    double p2 = val(P.x,P.y+dt,P.z);
    bool p = true;
    while(abs(target-p1)<0.000000000001 && abs(target-p2)<0.000000000001)
    {
        if(p)
            dt/=2;
        else
            dt*=2;
        p1 = val(P.x+dt,P.y,P.z);
        p2 = val(P.x,P.y+dt,P.z);
        if(dt<0.000000000001)
        {
            dt=dtO;
            p=false;
        }
    }
    //cout<<"p0 "<<target<<" "<<p1<<" "<<p2<<endl;
    //cout<<"P "<<P.x<<" "<<P.x+dt<<" "<<P.y<<" "<<P.y+dt<<endl;
    p1 = target-p1;
    p2 = target-p2;
    //cout<<"p1 "<<p1<<" "<<p2<<" "<<p3<<" "<<p4<<endl;
    float o = -M_PI/2;
    Point N = {0,0,0};
    Point T = {0,0,0};
    N.x = p1/dt;
    N.y = p2/dt;
    float m = dist(N.x,N.y,N.z);
    N.x/=m;
    N.y/=m;
    N.z/=m;
    //cout<<N.x<<" "<<N.y<<endl;
    T.x = N.x * cos(o) - N.y * sin(o);
    T.y = N.x * sin(o) + N.y * cos(o);
    T.z = N.z;
    float m2 = dist(T.x,T.y,T.z);
    T.x/=m2;
    T.y/=m2;
    T.z/=m2;
    return {P,N,T};
    /*
    if(p1==-1 || p2==-1)
    {
        target = -1;
        dt*=2;
        cout<<"U"<<endl;
    }
    if(target!=-1)
    {
        /*
        if(p1==0 && p2>0 && p4>0)
        {
            p1 = 0.001;
        }
        if(p2==0 && p1>0 && p3>0)
        {
            p2 = 0.001;
        }
        if(p3==0 && p2>0 && p4>0)
        {
            p3 = 0.001;
        }
        if(p4==0 && p1>0 && p3>0)
        {
            p4 = 0.001;
        }
        debug<<P.x<<"\t"<<P.y<<"\t"<<target<<"\t"<<p1<<"\t"<<p2<<"\t"<<p3<<"\t"<<p4<<endl;
        Point R = {0,0,0};
        int chx = msTable[to_string(p1>0)+to_string(p2>0)+to_string(p3>0)+to_string(p4>0)];
        float i1 = 0;
        float i2 = 0;
        float x = 0;
        float y = 0;
        float o = -M_PI/2;
        float m = 0;
        float x1 = 0;
        float x2 = 0;
        float y1v = 0;
        float y2v = 0;
        Point T = {0,0,0};
        Point mid = {0,0,0};
        switch(chx)
        {
            case 1:
                i1 = interpolation(p4,p1); 
                y = P.y-dt+i1*2*dt;
                i2 = interpolation(p3,p4);
                x = P.x-dt+i2*2*dt;
                mid = {(P.x-dt+x)/2,(P.y-dt+y)/2,P.z};
                T = {x-P.x+dt,P.y+dt-y,0};
                m = dist(T.x,T.y,T.z);
                T.x/=m;
                T.y/=m;
                T.z/=m;
                R.x = T.x * cos(o) - T.y * sin(o);
                R.y = T.x * sin(o) + T.y * cos(o);
                R.z = T.z;
                break;
            case 2:
                i1 = interpolation(p3,p2); 
                y = P.y-dt+i1*2*dt;
                i2 = interpolation(p3,p4);
                x = P.x-dt+i2*2*dt;
                mid = {(P.x+dt+x)/2,(P.y+dt+y)/2,P.z};
                T = {P.x-dt-x,y-P.y-dt,0};
                m = dist(T.x,T.y,T.z);
                T.x/=m;
                T.y/=m;
                T.z/=m;
                R.x = T.x * cos(o) - T.y * sin(o);
                R.y = T.x * sin(o) + T.y * cos(o);
                R.z = T.z;
                break;
            case 3:
                i1 = interpolation(p3,p2); 
                y1v = P.y-dt+i1*2*dt;
                i2 = interpolation(p4,p1);
                y2v = P.y-dt+i2*2*dt;
                mid = {P.x,(y1v+y2v)/2,P.z};
                T = {2*dt,y1v-y2v,0};
                m = dist(T.x,T.y,T.z);
                T.x/=m;
                T.y/=m;
                T.z/=m;
                R.x = T.x * cos(o) - T.y * sin(o);
                R.y = T.x * sin(o) + T.y * cos(o);
                R.z = T.z;
                break;
            case 4:
                i1 = interpolation(p3,p2); 
                y = P.y-dt+i1*2*dt;
                i2 = interpolation(p2,p1);
                x = P.x-dt+i2*2*dt;
                mid = {(P.x+dt+x)/2,(P.y-dt+y)/2,P.z};
                T = {P.x+dt-x,y-P.y+dt,0};
                m = dist(T.x,T.y,T.z);
                T.x/=m;
                T.y/=m;
                T.z/=m;
                R.x = T.x * cos(o) - T.y * sin(o);
                R.y = T.x * sin(o) + T.y * cos(o);
                R.z = T.z;
                break;
            case 6:
                i1 = interpolation(p2,p1); 
                x1 = P.x-dt+i1*2*dt;
                i2 = interpolation(p3,p4);
                x2 = P.x-dt+i2*2*dt;
                mid = {(x1+x2)/2,P.y,P.z};
                T = {(x1-x2),-2*dt,0};
                m = dist(T.x,T.y,T.z);
                T.x/=m;
                T.y/=m;
                T.z/=m;
                R.x = T.x * cos(o) - T.y * sin(o);
                R.y = T.x * sin(o) + T.y * cos(o);
                R.z = T.z;
                break;
            case 7:
                i1 = interpolation(p4,p1); 
                y = P.y-dt+i1*2*dt;
                i2 = interpolation(p2,p1);
                x= P.x-dt+i2*2*dt;
                mid = {(P.x-dt+x)/2,(P.y-dt+y)/2,P.z};
                T = {x-P.x+dt,P.y-dt-y,0};
                m = dist(T.x,T.y,T.z);
                T.x/=m;
                T.y/=m;
                T.z/=m;
                R.x = T.x * cos(o) - T.y * sin(o);
                R.y = T.x * sin(o) + T.y * cos(o);
                R.z = T.z;
                break;
        }
        return {P,R,T};
        */  
    //}
    /*else
    {
        cout<<endl;
        float ap = M_PI/50;
        float angle = 0;
        Point P1 = {0,0,0};
        Point P2 = {0,0,0};
        int pi = 0;
        int mult = 1;
        while(pi<2)
        {
            pi = 0;
            Point f = {P.x+cos(angle)*dt*mult,P.y+sin(angle)*dt*mult,P.z};
            float vf = func(f.x,f.y,f.z);
            bool type = (vf==-1);
            for(int i=1; i<100; i++)
            {
                angle = i*ap;
                f = {P.x+cos(angle)*dt*mult,P.y+sin(angle)*dt*mult,P.z};
                vf = func(f.x,f.y,f.z);
                //cout<<type<<" "<<vf<<" "<<mult<<endl;
                if((vf==-1)!=type)
                {
                    type=!type;
                    if(pi==0)
                    {
                        P1=f;
                    }
                    else
                    {
                        P2=f;
                    }
                    pi++;
                }
            }
            mult++;
        }
        Point T = {0,0,0};
        T.x = abs(P2.x-P1.x);
        T.y = abs(P2.y-P1.y);
        T.z = abs(P2.z-P1.z);
        cout<<"P1 "<<P1.x<<" "<<P1.y<<" "<<P1.z<<" "<<func(P1.x,P1.y,P1.z)<<endl;
        cout<<"P2 "<<P2.x<<" "<<P2.y<<" "<<P2.z<<" "<<func(P2.x,P2.y,P2.z)<<endl;
        float m = dist(T.x,T.y,T.z);
        T.x/=m;
        T.y/=m;
        T.z/=m;
        Point R = {0,0,0};
        float o = -M_PI/2;
        R.x = T.x * cos(o) - T.y * sin(o);
        R.y = T.x * sin(o) + T.y * cos(o);
        R.z = T.z;
        P1 = {P.x+0.05f*R.x,P.y+0.05f*R.y,P.z+0.05f*R.z};
        cout<<p1<<" "<<p2<<endl;
        cout<<"P0 "<<P.x<<" "<<P.y<<" "<<P.z<<" "<<dist(P.x,P.y,P.z)<<endl;
        cout<<"nP1 "<<P1.x<<" "<<P1.y<<" "<<P1.z<<" "<<dist(P1.x,P1.y,P1.z)<<endl;
        cout<<"R "<<R.x<<" "<<R.y<<" "<<R.z<<endl;
        if(dist(P1.x,P1.y,P1.z)>dist(P.x,P.y,P.z))
        {
            P1 = {P.x-0.05f*R.x,P.y-0.05f*R.y,P.z-0.05f*R.z};
            cout<<"nP2 "<<P1.x<<" "<<P1.y<<" "<<P1.z<<" "<<dist(P1.x,P1.y,P1.z)<<endl;
            //return getRepere(P1);
        }
        cout<<endl;
        if(func(P1.x,P1.y,P1.z)==-1)
            exit(0);
        return getRepere(P1);
    }*/
}


Repere getRepere0(Repere R)
{
    float dt = 0.001;
    Point P = R.P;
    float ap = M_PI/50;
    float angle = 0;
    Point P1 = {0,0,0};
    Point P2 = {0,0,0};
    int pi = 0;
    int mult = 1;
    while(pi==0)
    {
        Point f = {P.x+cos(angle)*dt*mult,P.y+sin(angle)*dt*mult,P.z};
        float vf = func(f.x,f.y,f.z);
        bool type = (vf==0);
        for(int i=1; i<100; i++)
        {
            angle = i*ap;
            f = {P.x+cos(angle)*dt*mult,P.y+sin(angle)*dt*mult,P.z};
            vf = func(f.x,f.y,f.z);
            if((vf==0)!=type)
            {
                type=!type;
                if(pi==0)
                {
                    P1=f;
                }
                else
                {
                    P2=f;
                }
                pi++;
            }
        }
        mult++;
    }
    Point T = {0,0,0};
    if(P1.x<P2.x)
    {
        T.x = P2.x-P1.x;
        T.y = P2.y-P1.y;
        T.z = P2.z-P1.z;
    }
    else
    {
        T.x = P1.x-P2.x;
        T.y = P1.y-P2.y;
        T.z = P1.z-P2.z;
    }
    float m = dist(T.x,T.y,T.z);
    T.x/=m;
    T.y/=m;
    T.z/=m;
    return {P,R.N,T};
}

void courbe(Point P, Point PrevN, int dir, float target)
{
    int g = 0;
    float val = func(P.x,P.y,P.z);
    //cout<<P.x<<" "<<P.y<<" "<<P.z<<" "<<val<<" "<<dir<<endl;
    if(val==0)
        return;
    Point N = {0,0,0};
    float dt = 0.01;
    if(dir!=0)
    {
        float val = func(P.x,P.y,P.z);
        auto start= system_clock::now();
        int nh = 0;
        while(val!=0)
        {
            Repere R = getRepere(P);
            N = R.N;
            Point T = R.T;
            Point nP = {P.x+T.x*dt*dir,P.y+T.y*dt*dir,P.z+T.z*dt*dir};
            Repere R2 = getRepere(nP);
            T = R2.T;
            Point nP2 = {P.x+T.x*dt*dir,P.y+T.y*dt*dir,P.z+T.z*dt*dir};
            float val1 = target-func(nP.x,nP.y,nP.z);
            float val2 = target-func(nP2.x,nP2.y,nP2.z);
            float inter = interpolation(val1,val2);
            Point nP3 = {0,0,0};
            nP3.x = nP.x * inter + (1-inter) * nP2.x;
            nP3.y = nP.y * inter + (1-inter) * nP2.y;
            nP3.z = nP.z * inter + (1-inter) * nP2.z;
            glBegin(GL_LINES);
            glVertex3f(P.x,P.y,P.z);
            glVertex3f(nP3.x,nP3.y,nP3.z);
            glEnd();
            P = nP3;
            PrevN = N;
            val = func(P.x,P.y,P.z);
            nh++;
            auto end= system_clock::now();
            auto elapsed = duration_cast<nanoseconds>(end - start);
            //cout<<"1 Point : "<<elapsed.count()<<endl;
            start= system_clock::now();
            moy2+=elapsed.count();
            n2++;
        }
    }
    else
    {
        N = PrevN;
        float o = M_PI/2;
        Point T1 = {0,0,N.z};
        T1.x = N.x * cos(o) - N.y * sin(o);
        T1.y = N.x * sin(o) + N.y * cos(o);
        Point nP1 = {P.x+T1.x*dt,P.y+T1.y*dt,P.z+T1.z*dt};
        o = -M_PI/2;
        Point T2 = {0,0,N.z};
        T2.x = N.x * cos(o) - N.y * sin(o);
        T2.y = N.x * sin(o) + N.y * cos(o);
        Point nP2 = {P.x+T2.x*dt,P.y+T2.y*dt,P.z+T2.z*dt};
        glBegin(GL_LINES);
        glColor3f(1,0,1);
        glVertex3f(P.x,P.y,P.z);
        glVertex3f(nP1.x,nP1.y,nP1.z);
        glVertex3f(P.x,P.y,P.z);
        glVertex3f(nP2.x,nP2.y,nP2.z);
        glEnd();
        courbe(nP1,N,1,target);
        courbe(nP2,N,-1,target);
    }    
}

void marchingSquare(float x, float y, float pas, int dir, float target, int ind, vector<Point> *P, int ordre )
{
    //cout<<ind<<" "<<dir<<" "<<ordre<<" "<<x<<" "<<y<<endl;
    if(ind<((x_size/pas)*(y_size/pas))/4)
    {
        int ix = round(x*1000);
        int iy = round(y*1000);
        int ip = round(pas*1000);
        string k = to_string(ix)+to_string(iy)+to_string(ip);
        //cout<<x<<" "<<y<<" "<<dir<<" "<<endl;
        if(Passage.find(k) == Passage.end())
        {
            Passage[k]=true;
            float p1 = target-func(x,y,0);
            float p2 = target-func(x+pas,y,0);
            float p3 = target-func(x+pas,y+pas,0);
            float p4 = target-func(x,y+pas,0);
            /*
            glColor3f(1,1,1);
            glBegin(GL_LINE_LOOP);
            glVertex3f(x,y,0);
            glVertex3f(x+pas,y,0);
            glVertex3f(x+pas,y+pas,0);
            glVertex3f(x,y+pas,0);
            glEnd();
            */
            if(p1==target+1 || p2==target+1 || p3==target+1 || p4==target+1)
            {
                if(dir==0)
                {
                    float nx = x;
                    float ny = y;
                    if(x<0)
                    {
                        nx+=abs(pas);
                    }
                    else
                    {
                        nx-=abs(pas);
                    }
                    if(y<0)
                    {
                        ny+=abs(pas);
                    }
                    else
                    {
                        ny-=abs(pas);
                    }
                    marchingSquare(nx,ny,pas,0,target,ind+1,P,0);
                }
            }
            else
            {
                Point P1 = {x,y,0};
                Point P2 = {x+pas,y,0};
                Point P3 = {x+pas,y+pas,0};
                Point P4 = {x,y+pas,0};
                float ap12 = abs(p1)+abs(p2);
                float ap14 = abs(p1)+abs(p4);
                float ap32 = abs(p3)+abs(p2);
                float ap34 = abs(p3)+abs(p4);
                //cout<<"12: "<<ap12<<", 14: "<<ap14<<", 32: "<<ap32<<", 34: "<<ap34<<" "<<p1<<" "<<p2<<" "<<p3<<" "<<p4<<endl;
                
                
                
                if(p1<0)
                {
                    if(p2<0)
                    {
                        if(p3<0)
                        {
                            if(p4<0)
                            {
                                //cout<<"AIE"<<endl;
                                if(ap12 < ap14)
                                {
                                    if(ap12 < ap32)
                                    {
                                        if(ap12<ap34 && ap12!=2*target)
                                        {
                                            marchingSquare(x,y-pas,pas,2,target, ind+1,P,ordre);
                                        }
                                        else
                                        {
                                            marchingSquare(x,y+pas,pas,4,target, ind+1,P,ordre);
                                        }
                                    }
                                    else
                                    {
                                        if(ap32<ap34 && ap32!=2*target)
                                        {
                                            marchingSquare(x+pas,y,pas,3,target, ind+1,P,ordre);
                                        }
                                        else
                                        {
                                            marchingSquare(x,y+pas,pas,4,target, ind+1,P,ordre);
                                        }
                                    }
                                }
                                else
                                {
                                    if(ap14 < ap32)
                                    {
                                        if(ap14<ap34 && ap14!=2*target)
                                        {
                                            marchingSquare(x-pas,y,pas,1,target, ind+1,P,ordre);
                                        }
                                        else
                                        {
                                            marchingSquare(x,y+pas,pas,4,target, ind+1,P,ordre);
                                        }
                                    }
                                    else
                                    {
                                        if(ap32<ap34 && ap32!=2*target)
                                        {
                                            marchingSquare(x+pas,y,pas,3,target, ind+1,P,ordre);
                                        }
                                        else
                                        {
                                            marchingSquare(x,y+pas,pas,4,target, ind+1,P,ordre);
                                        }
                                    }
                                }
                            }
                            else
                            {
                                //Calcul point
                                if(dir!=2 && ap34!=2*target)
                                {
                                    float inte = interpolation(p3,p4);
                                    Point p = {0,0,0};
                                    p.x = inte * P3.x + (1-inte)*P4.x;
                                    p.y = inte * P3.y + (1-inte)*P4.y;
                                    p.z = inte * P3.z + (1-inte)*P4.z;
                                    float np = target - func(p.x,p.y,p.z);
                                    if(np>limite)
                                    {
                                        if(np<0)
                                        {
                                            inte = interpolation(p4,np);
                                            p.x = inte * P4.x + (1-inte)*p.x;
                                            p.y = inte * P4.y + (1-inte)*p.y;
                                            p.z = inte * P4.z + (1-inte)*p.z;
                                        }
                                        else
                                        {
                                            inte = interpolation(np,p3);
                                            p.x = inte * p.x + (1-inte)*P3.x;
                                            p.y = inte * p.y + (1-inte)*P3.y;
                                            p.z = inte * p.z + (1-inte)*P3.z;
                                        }
                                    }
                                    //cout<<"Dessin "<<(p1<0)<<(p2<0)<<(p3<0)<<(p4<0)<<ordre<<" "<<dir<<" "<<ind<<endl;
                                    if(ordre==0)
                                    {
                                        P->push_back(p);
                                        nb++;
                                        marchingSquare(x,y+pas,pas,4,target, ind+1,P,1);
                                    }
                                    else
                                    {
                                        if(ordre==1)
                                        {
                                            P->push_back(p);
                                            nb++;
                                        }
                                        else
                                        {
                                            P->insert(P->begin(),p);
                                            nb++;
                                        }
                                        marchingSquare(x,y+pas,pas,4,target, ind+1,P,ordre);
                                    }
                                    
                                }
                                if(dir!=3 && ap14!=2*target)
                                {
                                    float inte = interpolation(p1,p4);
                                    Point p = {0,0,0};
                                    p.x = inte * P1.x + (1-inte)*P4.x;
                                    p.y = inte * P1.y + (1-inte)*P4.y;
                                    p.z = inte * P1.z + (1-inte)*P4.z;
                                    float np = target - func(p.x,p.y,p.z);
                                    if(np>limite)
                                    {
                                        if(np<0)
                                        {
                                            inte = interpolation(p4,np);
                                            p.x = inte * P4.x + (1-inte)*p.x;
                                            p.y = inte * P4.y + (1-inte)*p.y;
                                            p.z = inte * P4.z + (1-inte)*p.z;
                                        }
                                        else
                                        {
                                            inte = interpolation(np,p1);
                                            p.x = inte * p.x + (1-inte)*P1.x;
                                            p.y = inte * p.y + (1-inte)*P1.y;
                                            p.z = inte * p.z + (1-inte)*P1.z;
                                        }
                                    }
                                    //cout<<"Dessin "<<(p1<0)<<(p2<0)<<(p3<0)<<(p4<0)<<ordre<<" "<<dir<<" "<<ind<<endl;
                                    if(ordre==0)
                                    {
                                        P->insert(P->begin(),p);
                                        nb++;
                                        marchingSquare(x-pas,y,pas,1,target, ind+1,P,-1);
                                    }
                                    else
                                    {
                                        if(ordre==1)
                                        {
                                            P->push_back(p);
                                            nb++;
                                        }
                                        else
                                        {
                                            P->insert(P->begin(),p);
                                            nb++;
                                        }
                                        marchingSquare(x-pas,y,pas,1,target, ind+1,P,ordre);
                                    }
                                }

                            }
                        }
                        else
                        {
                            if(p4<0)
                            {
                                //Calcul point
                                if(dir!=2 && ap34!=2*target)
                                {
                                    float inte = interpolation(p3,p4);
                                    Point p = {0,0,0};
                                    p.x = inte * P3.x + (1-inte)*P4.x;
                                    p.y = inte * P3.y + (1-inte)*P4.y;
                                    p.z = inte * P3.z + (1-inte)*P4.z;
                                    float np = target - func(p.x,p.y,p.z);
                                    if(np>limite)
                                    {
                                        if(np>0)
                                        {
                                            inte = interpolation(p4,np);
                                            p.x = inte * P4.x + (1-inte)*p.x;
                                            p.y = inte * P4.y + (1-inte)*p.y;
                                            p.z = inte * P4.z + (1-inte)*p.z;
                                        }
                                        else
                                        {
                                            inte = interpolation(np,p3);
                                            p.x = inte * p.x + (1-inte)*P3.x;
                                            p.y = inte * p.y + (1-inte)*P3.y;
                                            p.z = inte * p.z + (1-inte)*P3.z;
                                        }
                                    }
                                    //cout<<"Dessin "<<(p1<0)<<(p2<0)<<(p3<0)<<(p4<0)<<ordre<<" "<<dir<<" "<<ind<<endl;
                                    if(ordre==0)
                                    {
                                        P->push_back(p);
                                        nb++;
                                        marchingSquare(x,y+pas,pas,4,target, ind+1,P,1);
                                    }
                                    else
                                    {
                                        if(ordre==1)
                                        {
                                            P->push_back(p);
                                            nb++;
                                        }
                                        else
                                        {
                                            P->insert(P->begin(),p);
                                            nb++;
                                        }
                                        marchingSquare(x,y+pas,pas,4,target, ind+1,P,ordre);
                                    }
                                }
                                if(dir!=1 && ap32!=2*target)
                                {
                                    float inte = interpolation(p3,p2);
                                    Point p = {0,0,0};
                                    p.x = inte * P3.x + (1-inte)*P2.x;
                                    p.y = inte * P3.y + (1-inte)*P2.y;
                                    p.z = inte * P3.z + (1-inte)*P2.z;
                                    float np = target - func(p.x,p.y,p.z);
                                    if(np>limite)
                                    {
                                        if(np>0)
                                        {
                                            inte = interpolation(p2,np);
                                            p.x = inte * P2.x + (1-inte)*p.x;
                                            p.y = inte * P2.y + (1-inte)*p.y;
                                            p.z = inte * P2.z + (1-inte)*p.z;
                                        }
                                        else
                                        {
                                            inte = interpolation(np,p3);
                                            p.x = inte * p.x + (1-inte)*P3.x;
                                            p.y = inte * p.y + (1-inte)*P3.y;
                                            p.z = inte * p.z + (1-inte)*P3.z;
                                        }
                                    }
                                    //cout<<"Dessin "<<(p1<0)<<(p2<0)<<(p3<0)<<(p4<0)<<ordre<<" "<<dir<<" "<<ind<<endl;
                                    if(ordre==0)
                                    {
                                        P->insert(P->begin(),p);
                                        nb++;
                                        marchingSquare(x+pas,y,pas,3,target, ind+1,P,-1);
                                    }
                                    else
                                    {
                                        if(ordre==1)
                                        {
                                            P->push_back(p);
                                            nb++;
                                        }
                                        else
                                        {
                                            P->insert(P->begin(),p);
                                            nb++;
                                        }
                                        marchingSquare(x+pas,y,pas,3,target, ind+1,P,ordre);
                                    }
                                }
                            }
                            else
                            {
                                //Calcul point
                                if(dir!=1 && ap32!=2*target)
                                {
                                    float inte = interpolation(p3,p2);
                                    Point p = {0,0,0};
                                    p.x = inte * P3.x + (1-inte)*P2.x;
                                    p.y = inte * P3.y + (1-inte)*P2.y;
                                    p.z = inte * P3.z + (1-inte)*P2.z;
                                    float np = target - func(p.x,p.y,p.z);
                                    if(np>limite)
                                    {
                                        if(np>0)
                                        {
                                            inte = interpolation(p2,np);
                                            p.x = inte * P2.x + (1-inte)*p.x;
                                            p.y = inte * P2.y + (1-inte)*p.y;
                                            p.z = inte * P2.z + (1-inte)*p.z;
                                        }
                                        else
                                        {
                                            inte = interpolation(np,p3);
                                            p.x = inte * p.x + (1-inte)*P3.x;
                                            p.y = inte * p.y + (1-inte)*P3.y;
                                            p.z = inte * p.z + (1-inte)*P3.z;
                                        }
                                    }
                                    //cout<<"Dessin "<<(p1<0)<<(p2<0)<<(p3<0)<<(p4<0)<<ordre<<" "<<dir<<" "<<ind<<endl;
                                    if(ordre==0)
                                    {
                                        P->push_back(p);
                                        nb++;
                                        marchingSquare(x+pas,y,pas,3,target, ind+1,P,1);
                                    }
                                    else
                                    {
                                        if(ordre==1)
                                        {
                                            P->push_back(p);
                                            nb++;
                                        }
                                        else
                                        {
                                            P->insert(P->begin(),p);
                                            nb++;
                                        }
                                        marchingSquare(x+pas,y,pas,3,target, ind+1,P,ordre);
                                    }
                                }
                                if(dir!=3 && ap14!=2*target)
                                {
                                    float inte = interpolation(p1,p4);
                                    Point p = {0,0,0};
                                    p.x = inte * P1.x + (1-inte)*P4.x;
                                    p.y = inte * P1.y + (1-inte)*P4.y;
                                    p.z = inte * P1.z + (1-inte)*P4.z;
                                    float np = target - func(p.x,p.y,p.z);
                                    if(np>limite)
                                    {
                                        if(np<0)
                                        {
                                            inte = interpolation(p4,np);
                                            p.x = inte * P4.x + (1-inte)*p.x;
                                            p.y = inte * P4.y + (1-inte)*p.y;
                                            p.z = inte * P4.z + (1-inte)*p.z;
                                        }
                                        else
                                        {
                                            inte = interpolation(np,p1);
                                            p.x = inte * p.x + (1-inte)*P1.x;
                                            p.y = inte * p.y + (1-inte)*P1.y;
                                            p.z = inte * p.z + (1-inte)*P1.z;
                                        }
                                    }
                                    //cout<<"Dessin "<<(p1<0)<<(p2<0)<<(p3<0)<<(p4<0)<<ordre<<" "<<dir<<" "<<ind<<endl;
                                    if(ordre==0)
                                    {
                                        P->insert(P->begin(),p);
                                        nb++;
                                        marchingSquare(x-pas,y,pas,1,target, ind+1,P,-1);
                                    }
                                    else
                                    {
                                        if(ordre==1)
                                        {
                                            P->push_back(p);
                                            nb++;
                                        }
                                        else
                                        {
                                            P->insert(P->begin(),p);
                                            nb++;
                                        }
                                        marchingSquare(x-pas,y,pas,1,target, ind+1,P,ordre);
                                    }
                                }
                            }
                        }
                    }
                    else
                    {
                        if(p3<0)
                        {
                            if(p4<0)
                            {
                                //Calcul point
                                if(dir!=1 && ap32!=2*target)
                                {
                                    float inte = interpolation(p3,p2);
                                    Point p = {0,0,0};
                                    p.x = inte * P3.x + (1-inte)*P2.x;
                                    p.y = inte * P3.y + (1-inte)*P2.y;
                                    p.z = inte * P3.z + (1-inte)*P2.z;
                                    float np = target - func(p.x,p.y,p.z);
                                    if(np>limite)
                                    {
                                        if(np<0)
                                        {
                                            inte = interpolation(p2,np);
                                            p.x = inte * P2.x + (1-inte)*p.x;
                                            p.y = inte * P2.y + (1-inte)*p.y;
                                            p.z = inte * P2.z + (1-inte)*p.z;
                                        }
                                        else
                                        {
                                            inte = interpolation(np,p3);
                                            p.x = inte * p.x + (1-inte)*P3.x;
                                            p.y = inte * p.y + (1-inte)*P3.y;
                                            p.z = inte * p.z + (1-inte)*P3.z;
                                        }
                                    }
                                    //cout<<"Dessin "<<(p1<0)<<(p2<0)<<(p3<0)<<(p4<0)<<ordre<<" "<<dir<<" "<<ind<<endl;
                                    if(ordre==0)
                                    {
                                        P->push_back(p);
                                        nb++;
                                        marchingSquare(x+pas,y,pas,3,target, ind+1,P,1);
                                    }
                                    else
                                    {
                                        if(ordre==1)
                                        {
                                            P->push_back(p);
                                            nb++;
                                        }
                                        else
                                        {
                                            P->insert(P->begin(),p);
                                            nb++;
                                        }
                                        marchingSquare(x+pas,y,pas,3,target, ind+1,P,ordre);
                                    }
                                }
                                if(dir!=4 && ap12!=2*target)
                                {
                                    float inte = interpolation(p1,p2);
                                    Point p = {0,0,0};
                                    p.x = inte * P1.x + (1-inte)*P2.x;
                                    p.y = inte * P1.y + (1-inte)*P2.y;
                                    p.z = inte * P1.z + (1-inte)*P2.z;
                                    float np = target - func(p.x,p.y,p.z);
                                    if(np>limite)
                                    {
                                        if(np<0)
                                        {
                                            inte = interpolation(p2,np);
                                            p.x = inte * P2.x + (1-inte)*p.x;
                                            p.y = inte * P2.y + (1-inte)*p.y;
                                            p.z = inte * P2.z + (1-inte)*p.z;
                                        }
                                        else
                                        {
                                            inte = interpolation(np,p1);
                                            p.x = inte * p.x + (1-inte)*P1.x;
                                            p.y = inte * p.y + (1-inte)*P1.y;
                                            p.z = inte * p.z + (1-inte)*P1.z;
                                        }
                                    }
                                    //cout<<"Dessin "<<(p1<0)<<(p2<0)<<(p3<0)<<(p4<0)<<ordre<<" "<<dir<<" "<<ind<<endl;
                                    if(ordre==0)
                                    {
                                        P->insert(P->begin(),p);
                                        nb++;
                                        marchingSquare(x,y-pas,pas,2,target, ind+1,P,-1);
                                    }
                                    else
                                    {
                                        if(ordre==1)
                                        {
                                            P->push_back(p);
                                            nb++;
                                        }
                                        else
                                        {
                                            P->insert(P->begin(),p);
                                            nb++;
                                        }
                                        marchingSquare(x,y-pas,pas,2,target, ind+1,P,ordre);
                                    }
                                }
                            }
                            else
                            {
                                //Double Passage
                            }
                        }
                        else
                        {
                            if(p4<0)
                            {
                                //Calcul point
                                if(dir!=2 && ap34!=2*target)
                                {
                                    float inte = interpolation(p3,p4);
                                    Point p = {0,0,0};
                                    p.x = inte * P3.x + (1-inte)*P4.x;
                                    p.y = inte * P3.y + (1-inte)*P4.y;
                                    p.z = inte * P3.z + (1-inte)*P4.z;
                                    float np = target - func(p.x,p.y,p.z);
                                    if(np>limite)
                                    {
                                        if(np>0)
                                        {
                                            inte = interpolation(p4,np);
                                            p.x = inte * P4.x + (1-inte)*p.x;
                                            p.y = inte * P4.y + (1-inte)*p.y;
                                            p.z = inte * P4.z + (1-inte)*p.z;
                                        }
                                        else
                                        {
                                            inte = interpolation(np,p3);
                                            p.x = inte * p.x + (1-inte)*P3.x;
                                            p.y = inte * p.y + (1-inte)*P3.y;
                                            p.z = inte * p.z + (1-inte)*P3.z;
                                        }
                                    }
                                    //cout<<"Dessin "<<(p1<0)<<(p2<0)<<(p3<0)<<(p4<0)<<ordre<<" "<<dir<<" "<<ind<<endl;
                                    if(ordre==0)
                                    {
                                        P->push_back(p);
                                        nb++;
                                        marchingSquare(x,y+pas,pas,4,target, ind+1,P,1);
                                    }
                                    else
                                    {
                                        if(ordre==1)
                                        {
                                            P->push_back(p);
                                            nb++;
                                        }
                                        else
                                        {
                                            P->insert(P->begin(),p);
                                            nb++;
                                        }
                                        marchingSquare(x,y+pas,pas,4,target, ind+1,P,ordre);
                                    }
                                }
                                if(dir!=4 && ap12!=2*target)
                                {
                                    float inte = interpolation(p1,p2);
                                    Point p = {0,0,0};
                                    p.x = inte * P1.x + (1-inte)*P2.x;
                                    p.y = inte * P1.y + (1-inte)*P2.y;
                                    p.z = inte * P1.z + (1-inte)*P2.z;
                                    float np = target - func(p.x,p.y,p.z);
                                    if(np>limite)
                                    {
                                        if(np<0)
                                        {
                                            inte = interpolation(p2,np);
                                            p.x = inte * P2.x + (1-inte)*p.x;
                                            p.y = inte * P2.y + (1-inte)*p.y;
                                            p.z = inte * P2.z + (1-inte)*p.z;
                                        }
                                        else
                                        {
                                            inte = interpolation(np,p1);
                                            p.x = inte * p.x + (1-inte)*P1.x;
                                            p.y = inte * p.y + (1-inte)*P1.y;
                                            p.z = inte * p.z + (1-inte)*P1.z;
                                        }
                                    }
                                    //cout<<"Dessin "<<(p1<0)<<(p2<0)<<(p3<0)<<(p4<0)<<ordre<<" "<<dir<<" "<<ind<<endl;
                                    if(ordre==0)
                                    {
                                        P->insert(P->begin(),p);
                                        nb++;
                                        marchingSquare(x,y-pas,pas,2,target, ind+1,P,-1);
                                    }
                                    else
                                    {
                                        if(ordre==1)
                                        {
                                            P->push_back(p);
                                            nb++;
                                        }
                                        else
                                        {
                                            P->insert(P->begin(),p);
                                            nb++;
                                        }
                                        marchingSquare(x,y-pas,pas,2,target, ind+1,P,ordre);
                                    }
                                } 
                            }
                            else
                            {
                                //Calcul point
                                if(dir!=3 && ap14!=2*target)
                                {
                                    float inte = interpolation(p1,p4);
                                    Point p = {0,0,0};
                                    p.x = inte * P1.x + (1-inte)*P4.x;
                                    p.y = inte * P1.y + (1-inte)*P4.y;
                                    p.z = inte * P1.z + (1-inte)*P4.z;
                                    float np = target - func(p.x,p.y,p.z);
                                    if(np>limite)
                                    {
                                        if(np<0)
                                        {
                                            inte = interpolation(p4,np);
                                            p.x = inte * P4.x + (1-inte)*p.x;
                                            p.y = inte * P4.y + (1-inte)*p.y;
                                            p.z = inte * P4.z + (1-inte)*p.z;
                                        }
                                        else
                                        {
                                            inte = interpolation(np,p1);
                                            p.x = inte * p.x + (1-inte)*P1.x;
                                            p.y = inte * p.y + (1-inte)*P1.y;
                                            p.z = inte * p.z + (1-inte)*P1.z;
                                        }
                                    }
                                    //cout<<"Dessin "<<(p1<0)<<(p2<0)<<(p3<0)<<(p4<0)<<ordre<<" "<<dir<<" "<<ind<<endl;
                                    if(ordre==0)
                                    {
                                        P->push_back(p);
                                        nb++;
                                        marchingSquare(x-pas,y,pas,1,target, ind+1,P,1);
                                    }
                                    else
                                    {
                                        if(ordre==1)
                                        {
                                            P->push_back(p);
                                            nb++;
                                        }
                                        else
                                        {
                                            P->insert(P->begin(),p);
                                            nb++;
                                        }
                                        marchingSquare(x-pas,y,pas,1,target, ind+1,P,ordre);
                                    }
                                }
                                if(dir!=4 && ap12!=2*target)
                                {
                                    float inte = interpolation(p1,p2);
                                    Point p = {0,0,0};
                                    p.x = inte * P1.x + (1-inte)*P2.x;
                                    p.y = inte * P1.y + (1-inte)*P2.y;
                                    p.z = inte * P1.z + (1-inte)*P2.z;
                                    float np = target - func(p.x,p.y,p.z);
                                    if(np>limite)
                                    {
                                        if(np<0)
                                        {
                                            inte = interpolation(p2,np);
                                            p.x = inte * P2.x + (1-inte)*p.x;
                                            p.y = inte * P2.y + (1-inte)*p.y;
                                            p.z = inte * P2.z + (1-inte)*p.z;
                                        }
                                        else
                                        {
                                            inte = interpolation(np,p1);
                                            p.x = inte * p.x + (1-inte)*P1.x;
                                            p.y = inte * p.y + (1-inte)*P1.y;
                                            p.z = inte * p.z + (1-inte)*P1.z;
                                        }
                                    }
                                    //cout<<"Dessin "<<(p1<0)<<(p2<0)<<(p3<0)<<(p4<0)<<ordre<<" "<<dir<<" "<<ind<<endl;
                                    if(ordre==0)
                                    {
                                        P->insert(P->begin(),p);
                                        nb++;
                                        marchingSquare(x,y-pas,pas,2,target, ind+1,P,-1);
                                    }
                                    else
                                    {
                                        if(ordre==1)
                                        {
                                            P->push_back(p);
                                            nb++;
                                        }
                                        else
                                        {
                                            P->insert(P->begin(),p);
                                            nb++;
                                        }
                                        marchingSquare(x,y-pas,pas,2,target, ind+1,P,ordre);
                                    }
                                }
                            }
                        }
                    }
                }
                else
                {
                    if(p2<0)
                    {
                        if(p3<0)
                        {
                            if(p4<0)
                            {
                                //Calcul point
                                if(dir!=3 && ap14!=2*target)
                                {
                                    float inte = interpolation(p1,p4);
                                    Point p = {0,0,0};
                                    p.x = inte * P1.x + (1-inte)*P4.x;
                                    p.y = inte * P1.y + (1-inte)*P4.y;
                                    p.z = inte * P1.z + (1-inte)*P4.z;
                                    float np = target - func(p.x,p.y,p.z);
                                    if(np>limite)
                                    {
                                        if(np>0)
                                        {
                                            inte = interpolation(p4,np);
                                            p.x = inte * P4.x + (1-inte)*p.x;
                                            p.y = inte * P4.y + (1-inte)*p.y;
                                            p.z = inte * P4.z + (1-inte)*p.z;
                                        }
                                        else
                                        {
                                            inte = interpolation(np,p1);
                                            p.x = inte * p.x + (1-inte)*P1.x;
                                            p.y = inte * p.y + (1-inte)*P1.y;
                                            p.z = inte * p.z + (1-inte)*P1.z;
                                        }
                                    }
                                    //cout<<"Dessin "<<(p1<0)<<(p2<0)<<(p3<0)<<(p4<0)<<ordre<<" "<<dir<<" "<<ind<<endl;
                                    if(ordre==0)
                                    {
                                        P->push_back(p);
                                        nb++;
                                        marchingSquare(x-pas,y,pas,1,target, ind+1,P,1);
                                    }
                                    else
                                    {
                                        if(ordre==1)
                                        {
                                            P->push_back(p);
                                            nb++;
                                        }
                                        else
                                        {
                                            P->insert(P->begin(),p);
                                            nb++;
                                        }
                                        marchingSquare(x-pas,y,pas,1,target, ind+1,P,ordre);
                                    }
                                }
                                if(dir!=4 && ap12!=2*target)
                                {
                                    float inte = interpolation(p1,p2);
                                    Point p = {0,0,0};
                                    p.x = inte * P1.x + (1-inte)*P2.x;
                                    p.y = inte * P1.y + (1-inte)*P2.y;
                                    p.z = inte * P1.z + (1-inte)*P2.z;
                                    float np = target - func(p.x,p.y,p.z);
                                    if(np>limite)
                                    {
                                        if(np>0)
                                        {
                                            inte = interpolation(p2,np);
                                            p.x = inte * P2.x + (1-inte)*p.x;
                                            p.y = inte * P2.y + (1-inte)*p.y;
                                            p.z = inte * P2.z + (1-inte)*p.z;
                                        }
                                        else
                                        {
                                            inte = interpolation(np,p1);
                                            p.x = inte * p.x + (1-inte)*P1.x;
                                            p.y = inte * p.y + (1-inte)*P1.y;
                                            p.z = inte * p.z + (1-inte)*P1.z;
                                        }
                                    }
                                    //cout<<"Dessin "<<(p1<0)<<(p2<0)<<(p3<0)<<(p4<0)<<ordre<<" "<<dir<<" "<<ind<<endl;
                                    if(ordre==0)
                                    {
                                        P->insert(P->begin(),p);
                                        nb++;
                                        marchingSquare(x,y-pas,pas,2,target, ind+1,P,-1);
                                    }
                                    else
                                    {
                                        if(ordre==1)
                                        {
                                            P->push_back(p);
                                            nb++;
                                        }
                                        else
                                        {
                                            P->insert(P->begin(),p);
                                            nb++;
                                        }
                                        marchingSquare(x,y-pas,pas,2,target, ind+1,P,ordre);
                                    }
                                }
                            }
                            else
                            {
                                //Calcul point
                                if(dir!=2 && ap34!=2*target)
                                {
                                    float inte = interpolation(p3,p4);
                                    Point p = {0,0,0};
                                    p.x = inte * P3.x + (1-inte)*P4.x;
                                    p.y = inte * P3.y + (1-inte)*P4.y;
                                    p.z = inte * P3.z + (1-inte)*P4.z;
                                    float np = target - func(p.x,p.y,p.z);
                                    if(np>limite)
                                    {
                                        if(np<0)
                                        {
                                            inte = interpolation(p4,np);
                                            p.x = inte * P4.x + (1-inte)*p.x;
                                            p.y = inte * P4.y + (1-inte)*p.y;
                                            p.z = inte * P4.z + (1-inte)*p.z;
                                        }
                                        else
                                        {
                                            inte = interpolation(np,p3);
                                            p.x = inte * p.x + (1-inte)*P3.x;
                                            p.y = inte * p.y + (1-inte)*P3.y;
                                            p.z = inte * p.z + (1-inte)*P3.z;
                                        }
                                    }
                                    //cout<<"Dessin "<<(p1<0)<<(p2<0)<<(p3<0)<<(p4<0)<<ordre<<" "<<dir<<" "<<ind<<endl;
                                    if(ordre==0)
                                    {
                                        P->push_back(p);
                                        nb++;
                                        marchingSquare(x,y+pas,pas,4,target, ind+1,P,1);
                                    }
                                    else
                                    {
                                        if(ordre==1)
                                        {
                                            P->push_back(p);
                                            nb++;
                                        }
                                        else
                                        {
                                            P->insert(P->begin(),p);
                                            nb++;
                                        }
                                        marchingSquare(x,y+pas,pas,4,target, ind+1,P,ordre);
                                    }
                                }
                                if(dir!=4 && ap12!=2*target)
                                {
                                    float inte = interpolation(p1,p2);
                                    Point p = {0,0,0};
                                    p.x = inte * P1.x + (1-inte)*P2.x;
                                    p.y = inte * P1.y + (1-inte)*P2.y;
                                    p.z = inte * P1.z + (1-inte)*P2.z;
                                    float np = target - func(p.x,p.y,p.z);
                                    if(np>limite)
                                    {
                                        if(np>0)
                                        {
                                            inte = interpolation(p2,np);
                                            p.x = inte * P2.x + (1-inte)*p.x;
                                            p.y = inte * P2.y + (1-inte)*p.y;
                                            p.z = inte * P2.z + (1-inte)*p.z;
                                        }
                                        else
                                        {
                                            inte = interpolation(np,p1);
                                            p.x = inte * p.x + (1-inte)*P1.x;
                                            p.y = inte * p.y + (1-inte)*P1.y;
                                            p.z = inte * p.z + (1-inte)*P1.z;
                                        }
                                    }
                                    //cout<<"Dessin "<<(p1<0)<<(p2<0)<<(p3<0)<<(p4<0)<<ordre<<" "<<dir<<" "<<ind<<endl;
                                    if(ordre==0)
                                    {
                                        P->insert(P->begin(),p);
                                        nb++;
                                        marchingSquare(x,y-pas,pas,2,target, ind+1,P,-1);
                                    }
                                    else
                                    {
                                        if(ordre==1)
                                        {
                                            P->push_back(p);
                                            nb++;
                                        }
                                        else
                                        {
                                            P->insert(P->begin(),p);
                                            nb++;
                                        }
                                        marchingSquare(x,y-pas,pas,2,target, ind+1,P,ordre);
                                    }
                                }
                            }
                        }
                        else
                        {
                            if(p4<0)
                            {
                                //Double Passage
                            }
                            else
                            {
                                //Calcul point
                                if(dir!=1 && ap32!=2*target)
                                {
                                    float inte = interpolation(p3,p2);
                                    Point p = {0,0,0};
                                    p.x = inte * P3.x + (1-inte)*P2.x;
                                    p.y = inte * P3.y + (1-inte)*P2.y;
                                    p.z = inte * P3.z + (1-inte)*P2.z;
                                    float np = target - func(p.x,p.y,p.z);
                                    if(np>limite)
                                    {
                                        if(np>0)
                                        {
                                            inte = interpolation(p2,np);
                                            p.x = inte * P2.x + (1-inte)*p.x;
                                            p.y = inte * P2.y + (1-inte)*p.y;
                                            p.z = inte * P2.z + (1-inte)*p.z;
                                        }
                                        else
                                        {
                                            inte = interpolation(np,p3);
                                            p.x = inte * p.x + (1-inte)*P3.x;
                                            p.y = inte * p.y + (1-inte)*P3.y;
                                            p.z = inte * p.z + (1-inte)*P3.z;
                                        }
                                    }
                                    //cout<<"Dessin "<<(p1<0)<<(p2<0)<<(p3<0)<<(p4<0)<<ordre<<" "<<dir<<" "<<ind<<endl;
                                    if(ordre==0)
                                    {
                                        P->push_back(p);
                                        nb++;
                                        marchingSquare(x+pas,y,pas,3,target, ind+1,P,1);
                                    }
                                    else
                                    {
                                        if(ordre==1)
                                        {
                                            P->push_back(p);
                                            nb++;
                                        }
                                        else
                                        {
                                            P->insert(P->begin(),p);
                                            nb++;
                                        }
                                        marchingSquare(x+pas,y,pas,3,target, ind+1,P,ordre);
                                    }
                                }
                                if(dir!=4 && ap12!=2*target)
                                {
                                    float inte = interpolation(p1,p2);
                                    Point p = {0,0,0};
                                    p.x = inte * P1.x + (1-inte)*P2.x;
                                    p.y = inte * P1.y + (1-inte)*P2.y;
                                    p.z = inte * P1.z + (1-inte)*P2.z;
                                    float np = target - func(p.x,p.y,p.z);
                                    if(np>limite)
                                    {
                                        if(np>0)
                                        {
                                            inte = interpolation(p2,np);
                                            p.x = inte * P2.x + (1-inte)*p.x;
                                            p.y = inte * P2.y + (1-inte)*p.y;
                                            p.z = inte * P2.z + (1-inte)*p.z;
                                        }
                                        else
                                        {
                                            inte = interpolation(np,p1);
                                            p.x = inte * p.x + (1-inte)*P1.x;
                                            p.y = inte * p.y + (1-inte)*P1.y;
                                            p.z = inte * p.z + (1-inte)*P1.z;
                                        }
                                    }
                                    //cout<<"Dessin "<<(p1<0)<<(p2<0)<<(p3<0)<<(p4<0)<<ordre<<" "<<dir<<" "<<ind<<endl;
                                    if(ordre==0)
                                    {
                                        P->insert(P->begin(),p);
                                        nb++;
                                        marchingSquare(x,y-pas,pas,2,target, ind+1,P,-1);
                                    }
                                    else
                                    {
                                        if(ordre==1)
                                        {
                                            P->push_back(p);
                                            nb++;
                                        }
                                        else
                                        {
                                            P->insert(P->begin(),p);
                                            nb++;
                                        }
                                        marchingSquare(x,y-pas,pas,2,target, ind+1,P,ordre);
                                    }
                                }
                            }
                        }
                    }
                    else
                    {
                        if(p3<0)
                        {
                            if(p4<0)
                            {
                                //Calcul point
                                if(dir!=3 && ap14!=2*target)
                                {
                                    float inte = interpolation(p1,p4);
                                    Point p = {0,0,0};
                                    p.x = inte * P1.x + (1-inte)*P4.x;
                                    p.y = inte * P1.y + (1-inte)*P4.y;
                                    p.z = inte * P1.z + (1-inte)*P4.z;
                                    float np = target - func(p.x,p.y,p.z);
                                    if(np>limite)
                                    {
                                        if(np>0)
                                        {
                                            inte = interpolation(p4,np);
                                            p.x = inte * P4.x + (1-inte)*p.x;
                                            p.y = inte * P4.y + (1-inte)*p.y;
                                            p.z = inte * P4.z + (1-inte)*p.z;
                                        }
                                        else
                                        {
                                            inte = interpolation(np,p1);
                                            p.x = inte * p.x + (1-inte)*P1.x;
                                            p.y = inte * p.y + (1-inte)*P1.y;
                                            p.z = inte * p.z + (1-inte)*P1.z;
                                        }
                                    }
                                    //cout<<"Dessin "<<(p1<0)<<(p2<0)<<(p3<0)<<(p4<0)<<ordre<<" "<<dir<<" "<<ind<<endl;
                                    if(ordre==0)
                                    {
                                        P->push_back(p);
                                        nb++;
                                        marchingSquare(x-pas,y,pas,1,target, ind+1,P,1);
                                    }
                                    else
                                    {
                                        if(ordre==1)
                                        {
                                            P->push_back(p);
                                            nb++;
                                        }
                                        else
                                        {
                                            P->insert(P->begin(),p);
                                            nb++;
                                        }
                                        marchingSquare(x-pas,y,pas,1,target, ind+1,P,ordre);
                                    }
                                }
                                if(dir!=1 && ap32!=2*target)
                                {
                                    float inte = interpolation(p3,p2);
                                    Point p = {0,0,0};
                                    p.x = inte * P3.x + (1-inte)*P2.x;
                                    p.y = inte * P3.y + (1-inte)*P2.y;
                                    p.z = inte * P3.z + (1-inte)*P2.z;
                                    float np = target - func(p.x,p.y,p.z);
                                    if(np>limite)
                                    {
                                        if(np<0)
                                        {
                                            inte = interpolation(p2,np);
                                            p.x = inte * P2.x + (1-inte)*p.x;
                                            p.y = inte * P2.y + (1-inte)*p.y;
                                            p.z = inte * P2.z + (1-inte)*p.z;
                                        }
                                        else
                                        {
                                            inte = interpolation(np,p3);
                                            p.x = inte * p.x + (1-inte)*P3.x;
                                            p.y = inte * p.y + (1-inte)*P3.y;
                                            p.z = inte * p.z + (1-inte)*P3.z;
                                        }
                                    }
                                    //cout<<"Dessin "<<(p1<0)<<(p2<0)<<(p3<0)<<(p4<0)<<ordre<<" "<<dir<<" "<<ind<<endl;
                                    if(ordre==0)
                                    {
                                        P->insert(P->begin(),p);
                                        nb++;
                                        marchingSquare(x+pas,y,pas,3,target, ind+1,P,-1);
                                    }
                                    else
                                    {
                                        if(ordre==1)
                                        {
                                            P->push_back(p);
                                            nb++;
                                        }
                                        else
                                        {
                                            P->insert(P->begin(),p);
                                            nb++;
                                        }
                                        marchingSquare(x+pas,y,pas,3,target, ind+1,P,ordre);
                                    }
                                }
                            }
                            else
                            {
                                //Calcul Point
                                if(dir!=2 && ap34!=2*target)
                                {
                                    float inte = interpolation(p3,p4);
                                    Point p = {0,0,0};
                                    p.x = inte * P3.x + (1-inte)*P4.x;
                                    p.y = inte * P3.y + (1-inte)*P4.y;
                                    p.z = inte * P3.z + (1-inte)*P4.z;
                                    float np = target - func(p.x,p.y,p.z);
                                    if(np>limite)
                                    {
                                        if(np<0)
                                        {
                                            inte = interpolation(p4,np);
                                            p.x = inte * P4.x + (1-inte)*p.x;
                                            p.y = inte * P4.y + (1-inte)*p.y;
                                            p.z = inte * P4.z + (1-inte)*p.z;
                                        }
                                        else
                                        {
                                            inte = interpolation(np,p3);
                                            p.x = inte * p.x + (1-inte)*P3.x;
                                            p.y = inte * p.y + (1-inte)*P3.y;
                                            p.z = inte * p.z + (1-inte)*P3.z;
                                        }
                                    }
                                    //cout<<"Dessin "<<(p1<0)<<(p2<0)<<(p3<0)<<(p4<0)<<ordre<<" "<<dir<<" "<<ind<<endl;
                                    if(ordre==0)
                                    {
                                        P->push_back(p);
                                        nb++;
                                        marchingSquare(x,y+pas,pas,4,target, ind+1,P,1);
                                    }
                                    else
                                    {
                                        if(ordre==1)
                                        {
                                            P->push_back(p);
                                            nb++;
                                        }
                                        else
                                        {
                                            P->insert(P->begin(),p);
                                            nb++;
                                        }
                                        marchingSquare(x,y+pas,pas,4,target, ind+1,P,ordre);
                                    }
                                }
                                if(dir!=1 && ap32!=2*target)
                                {
                                    float inte = interpolation(p3,p2);
                                    Point p = {0,0,0};
                                    p.x = inte * P3.x + (1-inte)*P2.x;
                                    p.y = inte * P3.y + (1-inte)*P2.y;
                                    p.z = inte * P3.z + (1-inte)*P2.z;
                                    float np = target - func(p.x,p.y,p.z);
                                    if(np>limite)
                                    {
                                        if(np<0)
                                        {
                                            inte = interpolation(p2,np);
                                            p.x = inte * P2.x + (1-inte)*p.x;
                                            p.y = inte * P2.y + (1-inte)*p.y;
                                            p.z = inte * P2.z + (1-inte)*p.z;
                                        }
                                        else
                                        {
                                            inte = interpolation(np,p3);
                                            p.x = inte * p.x + (1-inte)*P3.x;
                                            p.y = inte * p.y + (1-inte)*P3.y;
                                            p.z = inte * p.z + (1-inte)*P3.z;
                                        }
                                    }
                                    //cout<<"Dessin "<<(p1<0)<<(p2<0)<<(p3<0)<<(p4<0)<<ordre<<" "<<dir<<" "<<ind<<endl;
                                    if(ordre==0)
                                    {
                                        P->insert(P->begin(),p);
                                        nb++;
                                        marchingSquare(x+pas,y,pas,3,target, ind+1,P,-1);
                                    }
                                    else
                                    {
                                        if(ordre==1)
                                        {
                                            P->push_back(p);
                                            nb++;
                                        }
                                        else
                                        {
                                            P->insert(P->begin(),p);
                                            nb++;
                                        }
                                        marchingSquare(x+pas,y,pas,3,target, ind+1,P,ordre);
                                    }
                                }
                            }
                        }
                        else
                        {
                            if(p4<0)
                            {
                                //Calcul point
                                if(dir!=2 && ap34!=2*target)
                                {
                                    float inte = interpolation(p3,p4);
                                    Point p = {0,0,0};
                                    p.x = inte * P3.x + (1-inte)*P4.x;
                                    p.y = inte * P3.y + (1-inte)*P4.y;
                                    p.z = inte * P3.z + (1-inte)*P4.z;
                                    float np = target - func(p.x,p.y,p.z);
                                    if(np>limite)
                                    {
                                        if(np>0)
                                        {
                                            inte = interpolation(p4,np);
                                            p.x = inte * P4.x + (1-inte)*p.x;
                                            p.y = inte * P4.y + (1-inte)*p.y;
                                            p.z = inte * P4.z + (1-inte)*p.z;
                                        }
                                        else
                                        {
                                            inte = interpolation(np,p3);
                                            p.x = inte * p.x + (1-inte)*P3.x;
                                            p.y = inte * p.y + (1-inte)*P3.y;
                                            p.z = inte * p.z + (1-inte)*P3.z;
                                        }
                                    }
                                    //cout<<"Dessin "<<(p1<0)<<(p2<0)<<(p3<0)<<(p4<0)<<ordre<<" "<<dir<<" "<<ind<<endl;
                                    if(ordre==0)
                                    {
                                        P->push_back(p);
                                        nb++;
                                        marchingSquare(x,y+pas,pas,4,target, ind+1,P,1);
                                    }
                                    else
                                    {
                                        if(ordre==1)
                                        {
                                            P->push_back(p);
                                            nb++;
                                        }
                                        else
                                        {
                                            P->insert(P->begin(),p);
                                            nb++;
                                        }
                                        marchingSquare(x,y+pas,pas,4,target, ind+1,P,ordre);
                                    }
                                }
                                if(dir!=3 && ap14!=2*target)
                                {
                                    float inte = interpolation(p1,p4);
                                    Point p = {0,0,0};
                                    p.x = inte * P1.x + (1-inte)*P4.x;
                                    p.y = inte * P1.y + (1-inte)*P4.y;
                                    p.z = inte * P1.z + (1-inte)*P4.z;
                                    float np = target - func(p.x,p.y,p.z);
                                    if(np>limite)
                                    {
                                        if(np>0)
                                        {
                                            inte = interpolation(p4,np);
                                            p.x = inte * P4.x + (1-inte)*p.x;
                                            p.y = inte * P4.y + (1-inte)*p.y;
                                            p.z = inte * P4.z + (1-inte)*p.z;
                                        }
                                        else
                                        {
                                            inte = interpolation(np,p1);
                                            p.x = inte * p.x + (1-inte)*P1.x;
                                            p.y = inte * p.y + (1-inte)*P1.y;
                                            p.z = inte * p.z + (1-inte)*P1.z;
                                        }
                                    }
                                    //cout<<"Dessin "<<(p1<0)<<(p2<0)<<(p3<0)<<(p4<0)<<ordre<<" "<<dir<<" "<<ind<<endl;
                                    if(ordre==0)
                                    {
                                        P->insert(P->begin(),p);
                                        nb++;
                                        marchingSquare(x-pas,y,pas,1,target, ind+1,P,-1);
                                    }
                                    else
                                    {
                                        if(ordre==1)
                                        {
                                            P->push_back(p);
                                            nb++;
                                        }
                                        else
                                        {
                                            P->insert(P->begin(),p);
                                            nb++;
                                        }
                                        marchingSquare(x-pas,y,pas,1,target, ind+1,P,ordre);
                                    }
                                }
                            }
                            else
                            {
                                //cout<<"AIE"<<endl;
                                if(ap12 < ap14)
                                {
                                    if(ap12 < ap32)
                                    {
                                        if(ap12<ap34 && ap12!=2*target)
                                        {
                                            marchingSquare(x,y-pas,pas,2,target, ind+1,P,ordre);
                                        }
                                        else
                                        {
                                            marchingSquare(x,y+pas,pas,4,target, ind+1,P,ordre);
                                        }
                                    }
                                    else
                                    {
                                        if(ap32<ap34 && ap32!=2*target)
                                        {
                                            marchingSquare(x+pas,y,pas,3,target, ind+1,P,ordre);
                                        }
                                        else
                                        {
                                            marchingSquare(x,y+pas,pas,4,target, ind+1,P,ordre);
                                        }
                                    }
                                }
                                else
                                {
                                    if(ap14 < ap32)
                                    {
                                        if(ap14<ap34 && ap14!=2*target)
                                        {
                                            marchingSquare(x-pas,y,pas,1,target, ind+1,P,ordre);
                                        }
                                        else
                                        {
                                            marchingSquare(x,y+pas,pas,4,target, ind+1,P,ordre);
                                        }
                                    }
                                    else
                                    {
                                        if(ap32<ap34 && ap32!=2*target)
                                        {
                                            marchingSquare(x+pas,y,pas,3,target, ind+1,P,ordre);
                                        }
                                        else
                                        {
                                            marchingSquare(x,y+pas,pas,4,target, ind+1,P,ordre);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void dessin(vector<Point> P, float r, float g, float b)
{
    glBegin(GL_LINE_STRIP);
    for(int h=0; h<P.size(); h++)
    {
        glColor3f((float)h/(float)P.size(),(float)h/(float)P.size(),(float)h/(float)P.size());
        glVertex3f(P[h].x,P[h].y,P[h].z);
    }
    glEnd();
    /*int ph = P.size()/10;
    if(ph<1)
        ph=1;
    for(int h=1; h<P.size(); h++)
    {
        if((h-1)%ph==0)
        {
            Point mid = {P[h].x,P[h].y,P[h].z};
            mid.x+=P[h-1].x;
            mid.y+=P[h-1].y;
            mid.z+=P[h-1].z;
            mid.x/=2;
            mid.y/=2;
            mid.z/=2;
            Point t = {0,0,0};
            t.x = P[h].x-P[h-1].x;
            t.y = P[h].y-P[h-1].y;
            t.z = P[h].z-P[h-1].z;
            float mt = dist(t.x,t.y,t.z);
            t.x/=mt;
            t.y/=mt;
            t.z/=mt;
            float o = M_PI/2;
            Point N = {0,0,t.z};
            N.x = t.x * cos(o) - t.y * sin(o);
            N.y = t.x * sin(o) + t.y * cos(o);
            N.x/=10.0;
            N.y/=10.0;
            N.z/=10.0;
            
            glBegin(GL_LINES);
            glColor3f(0,0,1);
            glVertex3f(mid.x-N.x,mid.y-N.y,mid.z-N.z);
            glColor3f(0,1,0);
            glVertex3f(mid.x+N.x,mid.y+N.y,mid.z+N.z);
            glEnd();
        }
    } */
}

void lancerMarchingSquare(float target, Point Color)
{
    Passage.clear();
    vector<Point> P;
    marchingSquare(xMin,yMin,0.01,0,target,0,&P,0);
    dessin(P,Color.x,Color.y,Color.z);
    if(P.size()>0)
        courbes[target*1000].push_back(P);
    vector<Point> P1;
    marchingSquare(xMax,yMin,0.01,0,target,0,&P1,0);
    dessin(P1,Color.x,Color.y,Color.z);
    if(P1.size()>0)
        courbes[target*1000].push_back(P1);
    vector<Point> P2;
    marchingSquare(xMax,yMax,0.01,0,target,0,&P2,0);
    dessin(P2,Color.x,Color.y,Color.z);
    if(P2.size()>0)
        courbes[target*1000].push_back(P2);
    vector<Point> P3;
    marchingSquare(xMin,yMax,0.01,0,target,0,&P3,0);
    dessin(P3,Color.x,Color.y,Color.z);
    if(P3.size()>0)
        courbes[target*1000].push_back(P3);
    vector<Point> P4;
    marchingSquare(0,0,0.01,0,target,0,&P4,0);
    dessin(P4,Color.x,Color.y,Color.z);
    if(P4.size()>0)
        courbes[target*1000].push_back(P4);
}

void newMarchingSquare(double target, double pas)
{
    map<string,Square*> marque;
    vector<Square*> squares;
    vector<Point> pile;
    pile.push_back({0,0,pas});
    int n = 0;
    double min = 0.0000000001;
    while(pile.size()!=0)
    {
        Point Pi = pile[0];
        string k = to_string(Pi.x)+":"+to_string(Pi.y);
        if(marque.find(k)==marque.end())
        {
            Point P = {Pi.x*Pi.z,Pi.y*Pi.z,0};
            if(Pi.z<pas)
                cout<<pile.size()<<" "<<P.x<<" "<<P.y<<" "<<Pi.z<<endl;
            Point P0 = {P.x-Pi.z,P.y-Pi.z,P.z};
            Point P1 = {P.x-Pi.z,P.y+Pi.z,P.z};
            Point P2 = {P.x+Pi.z,P.y+Pi.z,P.z};
            Point P3 = {P.x+Pi.z,P.y-Pi.z,P.z};
            double v0 = target - func(P0.x,P0.y,P0.z);
            double v1 = target - func(P1.x,P1.y,P1.z);
            double v2 = target - func(P2.x,P2.y,P2.z);
            double v3 = target - func(P3.x,P3.y,P3.z);
            if((v0==target+1 || v1==target+1 || v2==target+1 || v3==target+1))
            {
                if(v0!=target+1 || v1!=target+1 || v2!=target+1 || v3!=target+1)
                {
                    if(Pi.z>pas/4)
                    {
                        glPointSize(1.0f);
                        glBegin(GL_POINTS);
                        glColor3f(1,1,1);
                        glVertex3f(P.x,P.y,P.z);
                        glEnd();
                        int nb = squares.size();
                        squares.push_back(new Square);
                        squares[nb]->P = P;
                        squares[nb]->type = -1;
                        squares[nb]->courbe = -1;
                        squares[nb]->pas = Pi.z;
                        marque[k]=squares[nb];
                        pile.push_back({Pi.x-0.5,Pi.y-0.5,Pi.z/2});
                        pile.push_back({Pi.x-0.5,Pi.y+0.5,Pi.z/2});
                        pile.push_back({Pi.x+0.5,Pi.y+0.5,Pi.z/2});
                        pile.push_back({Pi.x+0.5,Pi.y-0.5,Pi.z/2});
                    }
                }
            }
            else
            {
                if(abs(Pi.z-pas)<min)
                {
                    glPointSize(1.0f);
                    glBegin(GL_POINTS);
                    glColor3f(0,1,0);
                    glVertex3f(P.x,P.y,P.z);
                    glEnd();
                    if(v0==0)
                    {
                        if((v1>0 && v3<0) || (v1<0 && v3>0))
                        {
                            if(abs(v1)<abs(v3))
                                v0 = copysign(min,v1);
                            else
                                v0 = copysign(min,v3);
                        }
                        else
                        {
                            if(abs(v1)>abs(v3))
                                v0 = copysign(min,-v1);
                            else
                                v0 = copysign(min,-v3);
                        }
                    }
                    if(v1==0)
                    {
                        if((v0>0 && v2<0) || (v0<0 && v2>0))
                        {
                            if(abs(v0)<abs(v2))
                                v1 = copysign(min,v0);
                            else
                                v1 = copysign(min,v2);
                        }
                        else
                        {
                            v1 = copysign(min,-v0);
                        }
                    }
                    if(v2==0)
                    {
                        if((v1>0 && v3<0) || (v1<0 && v3>0))
                        {
                            if(abs(v1)<abs(v3))
                                v2 = copysign(min,v1);
                            else
                                v2 = copysign(min,v3);
                        }
                        else
                        {
                            v2 = copysign(min,-v1);
                        }
                    }
                    if(v3==0)
                    {
                        if((v0>0 && v2<0) || (v0<0 && v2>0))
                        {
                            if(abs(v0)<abs(v2))
                                v3 = copysign(min,v0);
                            else
                                v3 = copysign(min,v2);
                        }
                        else
                        {
                            v3 = copysign(min,-v0);
                        }
                    }
                    string s = to_string(v0>0)+to_string(v1>0)+to_string(v2>0)+to_string(v3>0);
                    int ch = msTable[s];
                    int nb = squares.size();
                    squares.push_back(new Square);
                    squares[nb]->P = P;
                    squares[nb]->type = ch;
                    squares[nb]->pas = Pi.z;
                    if(ch!=0)
                    {
                        squares[nb]->courbe = n;
                        n++;
                    }
                    else
                        squares[nb]->courbe = -1;
                    marque[k]=squares[nb];
                    pile.push_back({Pi.x-1,Pi.y,Pi.z});
                    pile.push_back({Pi.x+1,Pi.y,Pi.z});
                    pile.push_back({Pi.x,Pi.y-1,Pi.z});
                    pile.push_back({Pi.x,Pi.y+1,Pi.z});
                }
            }
        }
        pile.erase(pile.begin());
    }
}

double getAngle(Point V0, Point V1)
{
    double d = dist(V0.x,V0.y,V0.z)*dist(V1.x,V1.y,V1.z);
    if(d==0)
        return 0;
    double c = (V0.x*V1.x+V0.y*V1.y+V0.z*V1.z);
    double o = acos(clamp(c/d,-1.0,1.0));
    return abs(o*180.0/M_PI);
}

Repere getVecteurGradient(Point P, Point PrevN, double dt, double off)
{
    double val = func(P.x,P.y,P.z);
    int d = 1;
    Repere R = getRepere(P);
    Point N = R.N;
    double ang1 = getAngle(PrevN,N);
    double ang2 = 0;
    bool plop = false;
    if(ang1>160)
    {
        if(ang1<172)
        {
            cout<<"0PPP "<<ang1<<endl;
            plop=true;
        }
        //cout<<"a0 "<<ang1<<endl;
        /*glPointSize(20.0f);
        glColor3f(0,1,1);
        glBegin(GL_POINTS);
        glVertex3f(P.x,P.y,P.z);
        glEnd();*/
        d=-1;
    }
    else
    {
        if(ang1>8)
        {
            cout<<"0PPP "<<ang1<<endl;
            plop=true;
            /*cout<<"N0 "<<PrevN.x<<" "<<PrevN.y<<endl;
            cout<<"N1 "<<N.x<<" "<<N.y<<endl;
            double p1 = func(P.x+0.0001,P.y,P.z);
            double p2 = func(P.x,P.y+0.0001,P.z);
            cout<<(1.0-val)*100000<<" "<<(1-p1)*100000<<" "<<(1-p2)*100000<<endl;*/
            glPointSize(10.0f);
            glColor3f(1,1,1);
            glBegin(GL_POINTS);
            glVertex3f(P.x,P.y,P.z);
            glEnd();
        }
    }
    Point nP3 = {0,0,0};
    if(isnan(R.N.x))
    {
        cout<<"E"<<endl;
        N = PrevN;
        double m = dist(N.x,N.y,N.z);
        double r = m * dt + off;
        nP3 = {P.x+N.x*r,P.y+N.y*r,P.z+N.z*r};
    }
    else
    {
        double m = dist(N.x,N.y,N.z);
        double r = m * dt + off;
        double m2 = 0;
        double r2 = 0;
        Point nP = {P.x+N.x*r*d,P.y+N.y*r*d,P.z+N.z*r*d};
        Repere R2 = getRepere(nP);
        double x0 = 0;
        int nb0 = 1;
        int prev = (int)floor(val*pow(10,5));
        while(x0<r)
        {
            Point temp = addP(P,N,x0);
            double vt = func(temp.x,temp.y,temp.z);
            int ivt = (int)floor(vt*pow(10,5));
            if(ivt!=prev)
            {
                nb0++;
                prev=ivt;
            }
            x0+=1/pow(10,3);
        }
        if(!isnan(R2.N.x))
        {
            int d2 = d;
            double ang2 = getAngle(R2.N,N);
            /*cout<<"a "<<ang1<<" "<<ang2<<endl;
            cout<<"N0 "<<PrevN.x<<" "<<PrevN.y<<endl;
            cout<<"N1 "<<N.x<<" "<<N.y<<endl;*/
            if(ang2>160)
            {
                //cout<<"a1 "<<ang2<<endl;
                d2=-d;
                if(ang2<172)
                {
                    cout<<"OH1 "<<ang2<<endl;
                    plop=true;
                }
            }
            else
            {
                if(ang2>8)
                {
                    cout<<"OH2 "<<ang2<<endl;
                    plop=true;
                    /*cout<<"1PPP "<<ang2<<" "<<ang1<<" "<<getAngle(R2.N,PrevN)<<endl;
                    cout<<"N0 "<<PrevN.x<<" "<<PrevN.y<<endl;
                    cout<<"N1 "<<N.x<<" "<<N.y<<endl;
                    
                    double p1 = func(P.x+0.0001,P.y,P.z);
                    double p2 = func(P.x,P.y+0.0001,P.z);
                    double p3 = func(nP.x+0.0001,nP.y,nP.z);
                    double p4 = func(nP.x,nP.y+0.0001,nP.z);
                    double v = func(nP.x,nP.y,nP.z);
                    cout<<(1.0-val)*100000<<" "<<(1-p1)*100000<<" "<<(1-p2)*100000<<" ,nP : "<<(1.0-v)*100000<<" "<<(1-p3)*100000<<" "<<(1-p4)*100000<<endl;*/
                    glPointSize(10.0f);
                    glColor3f(1,1,0);
                    glBegin(GL_POINTS);
                    glVertex3f(P.x,P.y,P.z);
                    glEnd();
                }
            }
            N = R2.N;
            m2 = dist(N.x,N.y,N.z);
            r2 = m2 * dt + off;
            Point nP2 = {P.x+N.x*r2*d2,P.y+N.y*r2*d2,P.z+N.z*r2*d2};
            x0 = 0;
            int nb2 = 1;
            prev = (int)floor(val*pow(10,5));
            while(x0<r)
            {
                Point temp = addP(P,N,x0);
                double vt = func(temp.x,temp.y,temp.z);
                int ivt = (int)floor(vt*pow(10,5));
                if(ivt!=prev)
                {
                    nb2++;
                    prev=ivt;
                }
                x0+=1/pow(10,3);
            }
            //float val1 = val-func(nP.x,nP.y,nP.z);
            //float val2 = val-func(nP2.x,nP2.y,nP2.z);
            float inter = interpolation((double)nb2,(double)nb0);
            
            if(plop)
                cout<<nb0<<" "<<nb2<<" "<<inter<<endl;
            /*cout<<"N2 "<<N.x<<" "<<N.y<<endl;
            Point N3 = {nP.x-P.x,nP.y-P.y,nP.z-P.z};
            cout<<"N3 "<<N3.x<<" "<<N3.y<<endl;
            cout<<"v "<<val1<<" "<<val2<<" "<<inter<<endl;
            cout<<"P "<<P.x<<" "<<P.y<<endl;
            cout<<"nP "<<nP.x<<" "<<nP.y<<endl;
            cout<<"nP2 "<<nP2.x<<" "<<nP2.y<<endl;*/
            nP3.x = nP.x * inter + (1-inter) * nP2.x;
            nP3.y = nP.y * inter + (1-inter) * nP2.y;
            nP3.z = nP.z * inter + (1-inter) * nP2.z;
            //cout<<"nP3 "<<nP3.x<<" "<<nP3.y<<endl;
        }
        else
        {
            nP3 = nP;
        }
        N = {nP3.x-P.x,nP3.y-P.y,nP3.z-P.z};
        double r3 = (r+r2)/2.0;
        N.x/=r3;
        N.y/=r3;
        N.z/=r3;
        double m3 = dist(N.x,N.y,N.z);
        //cout<<"m "<<m0<<" "<<m2<<" "<<m3<<endl;
        //cout<<"a "<<a<<" dt "<<dt<<" "<<dt0<<" "<<dt2<<" "<<dt3<<endl;
        /*N.x *= d;
        N.y *= d;
        N.z *= d;*/
    }
    return {nP3,N,{0,0,0}};
}

void courbeGradient(Point P, int dir, Point PrevN, vector<PointCourbe> *Pv, double dt, double off)
{
    int g = 0;
    double val = func(P.x,P.y,P.z);
    Point N = {0,0,0};
    float m = dist(PrevN.x,PrevN.y,PrevN.z);
    int st = 0;
    if(dir!=0)
    {
        int d = 1;
        cout<<endl<<endl<<endl;
        Point P2 = addP(P,PrevN,-dt);
        val = func(P.x,P.y,P.z);
        double val2 = func(P2.x,P2.y,P2.z);
        while(val!=-1)
        {
            //cout<<st<<" "<<P.x<<" "<<P.y<<endl<<endl;
            int d = 1;
            int iv = 0;
            while(iv<targets.size() && targets[iv]+0.002<val)
            {
                if(abs(targets[iv]-val)<0.003)
                {
                    glBegin(GL_POINTS);
                    glVertex3f(P.x,P.y,P.z);
                    glEnd();
                    //cout<<val<<" ";
                    if(dir==1)
                    {
                        PointCourbe Pc = Pv->at(Pv->size()-1);
                        //cout<<targets[iv]<<" "<<val<<endl;
                        if(targets[iv]!=Pc.val || dist2(P,Pc.P)>0.02001)
                        {
                            Pv->push_back({P,targets[iv],-1,-1,-1,-1});
                        }
                    }
                    else
                    {
                        PointCourbe Pc = Pv->at(0);
                        //cout<<targets[iv]<<" "<<val<<endl;
                        if(targets[iv]!=Pc.val || dist2(P,Pc.P)>0.02001)
                        {
                            Pv->insert(Pv->begin(),{P,targets[iv],-1,-1,-1,-1});
                        }
                    }
                }
                iv++;
            }
            Repere R = getVecteurGradient(P,PrevN,dt,off);
            Point nP3 = R.P;
            N = R.N;
            /*Repere R = getRepere(P);
            Point N = R.N;
            double ang1 = getAngle(PrevN,N);
            if(ang1>160)
            {
                //cout<<"a0 "<<ang1<<endl;
                glPointSize(20.0f);
                glColor3f(0,1,1);
                glBegin(GL_POINTS);
                glVertex3f(P.x,P.y,P.z);
                glEnd();
                d=-d;
            }
            if(d==-1)
                N = addP({0,0,0},N,d);
            Point nP3 = addP(P,N,1/(20*dt));
            Point N2 = addP({0,0,0},N,1/(20*dt));
            cout<<dist(N.x,N.y,N.z)<<" "<<dist(N2.x,N2.y,N2.z)<<" "<<dt<<" "<<1/(20*dt)<<endl;*/
            //cout<<"N00 "<<PrevN.x<<" "<<PrevN.y<<endl;
            PrevN=N;
            //cout<<"N01 "<<PrevN.x<<" "<<PrevN.y<<endl;
            glColor3f(color,color,color);
            glBegin(GL_LINES);
            glVertex3f(P.x,P.y,P.z);
            glVertex3f(nP3.x,nP3.y,nP3.z);
            glEnd();
            if(dir==1)
            {
                glColor3f(0,1,0);
            }
            else
            {
                glColor3f(0,0,1);
            }
            glPointSize(5.0f);
            glBegin(GL_POINTS);
            glVertex3f(P.x,P.y,P.z);
            glEnd();
            P2 = P;
            P = nP3;
            val2 = val;
            val = func(P.x,P.y,P.z);
            st++;
        }
        if(dir==1)
            Pv->push_back({P2,0,-1,-1,-1,-1});
        else
           Pv->insert(Pv->begin(),{P2,0,-1,-1,-1,-1});
    }
    else
    {
        cout<<"g"<<endl;
        if(val==0)
            return;
        glColor3f(0,1,0);
        glBegin(GL_POINTS);
        glVertex3f(P.x,P.y,P.z);
        glEnd();
        Repere R = getRepere(P);
        N = R.N;
        double m = dist(N.x,N.y,N.z);
        double r = dt * m + off;
        Point nP1 = {P.x+N.x*r,P.y+N.y*r,P.z+N.z*r};
        Point nP2 = {P.x-N.x*r,P.y-N.y*r,P.z-N.z*r};
        glBegin(GL_POINTS);
        glColor3f(1,1,0);
        glVertex3f(nP1.x,nP1.y,nP1.z);
        glVertex3f(nP2.x,nP2.y,nP2.z);
        glEnd();
        glBegin(GL_LINES);
        glColor3f(color,color,color);
        glVertex3f(P.x,P.y,P.z);
        glVertex3f(nP1.x,nP1.y,nP1.z);
        glVertex3f(P.x,P.y,P.z);
        glVertex3f(nP2.x,nP2.y,nP2.z);
        glEnd();
        Point N2 = {-N.x,-N.y,-N.z};
        courbeGradient(nP1,1,N,Pv,dt,off);
        courbeGradient(nP2,-1,N2,Pv,dt,off);
    }
    cout<<"H"<<endl; 
}

Point addP(Point p1, Point p2, double mult)
{
    return {p1.x+p2.x*mult,p1.y+p2.y*mult,p1.z+p2.z*mult};
}

float courbeBezier(Point P0, Point P1, Point P2, Point P3, bool dessin)
{
    Point PrevP = P0;
    float longueur = 0;
    
    glBegin(GL_LINE_STRIP);
    for(int i=0; i<=100; i++)
    {
        float t = i*1.0/100.0;
        Point P = {0,0,0};
        P = addP(P,P0,(1-t)*(1-t)*(1-t));
        P = addP(P,P1,3*t*(1-t)*(1-t));
        P = addP(P,P2,3*t*t*(1-t));
        P = addP(P,P3,t*t*t);
        longueur+=dist2(PrevP,P);
        PrevP = P;
        //glVertex3f(PrevP.x,PrevP.y,PrevP.z);
        if(dessin)
            glVertex3f(P.x,P.y,P.z);
    }
    glEnd();
    return longueur;
}

void stockerBezier(half_edge* he, vector<Point> *courbe)
{
    Point P0 = he->origine;
    Point P1 = he->controleO;
    Point P2 = he->controleI;
    Point P3 = he->incident;
    for(int i=0; i<=100; i++)
    {
        float t = i*1.0/100.0;
        Point P = {0,0,0};  
        P = addP(P,P0,(1-t)*(1-t)*(1-t));
        P = addP(P,P1,3*t*(1-t)*(1-t));
        P = addP(P,P2,3*t*t*(1-t));
        P = addP(P,P3,t*t*t);
        courbe->push_back(P);
    }
}


Point chercherBezier(half_edge* he, Point Pc, float start, float end, int step, float dist)
{
    Point P0 = he->origine;
    Point P1 = he->controleO;
    Point P2 = he->controleI;
    Point P3 = he->incident;
    float t1 = 0.75*start+0.25*end;
    Point Pr1 = {0,0,0};  
    Pr1 = addP(Pr1,P0,(1-t1)*(1-t1)*(1-t1));
    Pr1 = addP(Pr1,P1,3*t1*(1-t1)*(1-t1));
    Pr1 = addP(Pr1,P2,3*t1*t1*(1-t1));
    Pr1 = addP(Pr1,P3,t1*t1*t1);
    float d1 = dist2(Pr1,Pc);
    float t2 = 0.25*start+0.75*end;
    Point Pr2 = {0,0,0};  
    Pr2 = addP(Pr2,P0,(1-t2)*(1-t2)*(1-t2));
    Pr2 = addP(Pr2,P1,3*t2*(1-t2)*(1-t2));
    Pr2 = addP(Pr2,P2,3*t2*t2*(1-t2));
    Pr2 = addP(Pr2,P3,t2*t2*t2);
    float d2 = dist2(Pr2,Pc);
    float mid = 0.5*start + 0.5*end;
    float distZ = dist2(Pr1,Pr2);
    if(d1<d2)
    {
        if(d1<dist || distZ<0.0000001)
        {
            //cout<<d1<<" "<<dist<<" "<<distZ<<" "<<step<<endl;
            return Pr1;
        }
        else
        {
            return chercherBezier(he,Pc,start,0.5*start+0.5*end,step+1,dist);
        }
    }
    else
    {
        if(d2<dist || distZ<0.0000001)
        {
            //cout<<d1<<" "<<dist<<" "<<distZ<<" "<<step<<endl;
            return Pr2;
        }
        else
        {
            return chercherBezier(he,Pc,0.5*start+0.5*end,end,step+1,dist);
        }
    }
}

Point chercherBezier2(half_edge* he, Point Pc, double start, double end, double dist)
{
    string out = "";
    int step = 0;
    Point P0 = he->origine;
    Point P1 = he->controleO;
    Point P2 = he->controleI;
    Point P3 = he->incident;
    Point Pr = {0,0,0};
    double p = abs(end-start)/10.0f;
    double distMin = dist2(Pc,P0);
    do{
        double tmin = start;
        Point Pr3 =  addP({0,0,0},P0,(1-tmin)*(1-tmin)*(1-tmin));
        Pr3 = addP(Pr3,P1,3*tmin*(1-tmin)*(1-tmin));
        Pr3 = addP(Pr3,P2,3*tmin*tmin*(1-tmin));
        Pr3 = addP(Pr3,P3,tmin*tmin*tmin);
        Pr = Pr3;
        double distMin = dist2(Pr,Pc);
        p = abs(end-start)/10.0f;
        for(float t=start; t<=end; t+=p)
        {
            Pr3 = addP({0,0,0},P0,(1-t)*(1-t)*(1-t));
            Pr3 = addP(Pr3,P1,3*t*(1-t)*(1-t));
            Pr3 = addP(Pr3,P2,3*t*t*(1-t));
            Pr3 = addP(Pr3,P3,t*t*t);
            float d = dist2(Pr3,Pc);
            if(d<distMin)
            {
                Pr = Pr3;
                distMin = d;
                tmin = t;
            }
            step++;
        }
        start = tmin - p;
        end = tmin + p;
    }while(distMin>dist && p>0.00001);
    /*
    if(d1<=dist || d2<=dist)
        cout<<out<<endl;
        */
    /*
    if(d1>dist && d2>dist && (d1<dis1*0.05 || d2<dis2*0.05))
    {
        cout<<d1<<" "<<dis1<<" "<<d2<<" "<<dis2<<endl;
        //cout<<"Plop"<<endl;
        Point Pr = Pr1;
        Point Pr3 = {0,0,0};
        float distMin = d1;
        int i = 0;
        for(float t=start; t<=end; t+=abs(end-start)/30.0f)
        {
            //cout<<i<<endl;
            i++;
            Pr3 = addP({0,0,0},P0,(1-t)*(1-t)*(1-t));
            Pr3 = addP(Pr3,P1,3*t*(1-t)*(1-t));
            Pr3 = addP(Pr3,P2,3*t*t*(1-t));
            Pr3 = addP(Pr3,P3,t*t*t);
            float d = dist2(Pr3,Pc);
            if(d<distMin)
            {
                Pr = Pr3;
                distMin = d;
            }
        }
        return Pr3;
    }*/
    //cout<<step<<endl;
    return Pr;
}

Point distanceMinBezier(heFace f, Point P, float dist)
{
    half_edge *he = f.incidente;
    int id = he->id;
    Point ret = he->origine;
    float mindist = dist2(ret,P);
    do
    {
        Point Pr = chercherBezier(he,P,0,1,0,dist);
        float dist = dist2(P,Pr);
        if(dist<mindist)
        {
            mindist=dist;
            ret = Pr;
        }
        he = he->next;
    }while(he->id!=id);
    return ret;
}

Point distanceMinBezier2(heFace f, Point P, float dist)
{
    half_edge *he = f.incidente;
    int id = he->id;
    Point ret = he->origine;
    float mindist = dist2(ret,P);
    do
    {
        Point Pr = chercherBezier2(he,P,0,1,dist);
        float dist = dist2(P,Pr);
        if(dist<mindist)
        {
            mindist=dist;
            ret = Pr;
        }
        he = he->next;
    }while(he->id!=id);
    return ret;
}

RetourDensite calculDensité(float x, float y, bool in, Point P, float yMax, heFace f, bool cherche, int st, int last)
{
    RetourDensite rd = {0,0};
    if(y>=yMax)
    {
        if(in)
        {
            glColor3f(0,0,1);
            glBegin(GL_POINTS);
            while(st>=last)
            {
                rd.nb-=1;
                rd.somme-=func(x,y,0);
                st--;
                y-=0.01;
                if(dessinne)
                    glVertex3f(x,y,0);
            }
            glEnd();
        }
        return rd;
    }
    /*if(dessinne)
    {
        cout<<x<<"  \t"<<y<<"  \t"<<P.x<<"  \t"<<P.y<<"  \t"<<in<<endl;
        glPointSize(7.0f);
        glColor3f(0,1,1);
        glBegin(GL_POINTS);
        glVertex3f(P.x,P.y,0);
        glEnd();
    }*/
    if(P.y<=y || st%5==0)
    {
        if(dessinne)
            cout<<"Plop"<<y<<" "<<P.y<<" ch "<<cherche<<endl;
        P = distanceMinBezier(f, {x,y,0},0.05);
        //cout<<"Plop "<<dist<<endl;
        float dist = dist2({x,y,0},P);
        if(dessinne)
            cout<<dist<<" "<<P.y<<endl;
        if(y>=P.y && abs(x-P.x)<0.02 && cherche)
        {
            if(dessinne)
            {
                if(!in)
                    cout<<"IN"<<endl;
                else
                    cout<<"OUT"<<endl;
            }
            cherche = false;
            in = !in;
            last = st;
            /*if(dessinne)
            {
                glColor3f(0,1,0);
                glBegin(GL_POINTS);
                glVertex3f(x,y,0);
                glEnd();
                glColor3f(0,0,1);
                glBegin(GL_POINTS);
                glVertex3f(P.x,P.y,0);
                glEnd();
            }*/
        }
        if(y<P.y)
        {
            cherche = true;
        }
    }
    if(in)
    {
        rd.nb=1;
        rd.somme=func(x,y,0);
    }
    if(dessinne)
    {
        glPointSize(2.0f);
        if(in)
        {
            glColor3f(0,1,0);
        }
        else
        {
            glColor3f(0,0,1);
        }
        glBegin(GL_POINTS);
        glVertex3f(x,y,0);
        glEnd();
        glPointSize(2.5f);
        glColor3f(0,1,1);
        glBegin(GL_POINTS);
        glVertex3f(P.x,P.y,0);
        glEnd();
    }
    
    RetourDensite retour = {0,0};
    if(!in)
    {
        if(P.y<=y)
        {
            if(dessinne)
                cout<<"Skip 0.5"<<endl;
            retour = calculDensité(x,y+0.05,in,P,yMax,f,cherche,0,last);
        }
        else
        {
            if(dessinne)
                cout<<"Skip Long"<<endl;
            retour = calculDensité(x,P.y,in,P,yMax,f,cherche,0,last);
        }
    }
    else
    {
        if(dessinne)
            cout<<"In"<<endl;
        retour = calculDensité(x,y+0.01,in,P,yMax,f,cherche,st+1,last);
    }
    
    rd.nb+=retour.nb;
    rd.somme+=retour.somme;
    return rd;
}

float densite(heFace f)
{
    vector<vector<Point>> beziers;
    half_edge *he = f.incidente;
    float xMin = he->origine.x;
    float xMax = he->origine.x;
    float yMin = he->origine.y;
    float yMax = he->origine.y;
    //cout<<xMin<<" "<<xMax<<" "<<yMin<<" "<<yMax<<endl;
    int id = he->id;
    do
    {
        if(he->origine.x<xMin)
        {
            xMin=he->origine.x-0.02f;
        }
        if(he->origine.x>xMax)
        {
            xMax=he->origine.x+0.03f;
        }
        if(he->origine.y<yMin)
        {
            yMin=he->origine.y-0.02f;
        }
        if(he->origine.y>yMax)
        {
            yMax=he->origine.y+0.03f;
        }

        if(he->controleO.x<xMin)
        {
            xMin=he->controleO.x;
        }
        if(he->controleO.x>xMax)
        {
            xMax=he->controleO.x;
        }
        if(he->controleO.y<yMin)
        {
            yMin=he->controleO.y;
        }
        if(he->controleO.y>yMax)
        {
            yMax=he->controleO.y;
        }

        if(he->controleI.x<xMin)
        {
            xMin=he->controleI.x;
        }
        if(he->controleI.x>xMax)
        {
            xMax=he->controleI.x;
        }
        if(he->controleI.y<yMin)
        {
            yMin=he->controleI.y;
        }
        if(he->controleI.y>yMax)
        {
            yMax=he->controleI.y;
        }

        if(he->incident.x<xMin)
        {
            xMin=he->incident.x-0.02f;
        }
        if(he->incident.x>xMax)
        {
            xMax=he->incident.x+0.03f;
        }
        if(he->incident.y<yMin)
        {
            yMin=he->incident.y-0.02f;
        }
        if(he->incident.y>yMax)
        {
            yMax=he->incident.y+0.03f;
        }
        /*vector<Point> bezier;
        stockerBezier(he, &bezier);
        beziers.push_back(bezier);*/
        he = he->next;
    }while(he->id!=id);
    cout<<"Limite "<<xMin<<"\t"<<yMin<<"\t"<<xMax<<"\t"<<yMax<<endl;
    if(dessinne)
    {
        glColor3f(1,1,1);
        glBegin(GL_LINE_LOOP);
        glVertex3f(xMin,yMin,0);
        glVertex3f(xMin,yMax,0);
        glVertex3f(xMax,yMax,0);
        glVertex3f(xMax,yMin,0);
        glEnd();
    }
    Point P = distanceMinBezier(f,{xMin,yMin,0},0.05);
    RetourDensite total = {0,0};
    dessinne = false;
    int h=0;
    int nb = 0;
    for(float x = xMin; x<=xMax; x+=0.01)
    {
        //cout<<total.somme<<" "<<total.nb<<endl;
        if(x>(xMin+xMax)/2 && h==0)
        {
            dessinne=true;
            h=1;
        }
        else
        {
            dessinne=false;
        }
        RetourDensite rd = calculDensité(x,yMin,false,P,yMax,f,true,0,0);
        total.nb += rd.nb;
        total.somme += rd.somme;
        nb++;
        if(abs(x-P.x)<0.01 || nb == 10)
        {
            P = distanceMinBezier(f,{x,yMin,0},0.05);
            nb= 0;
        }
    }
    if(total.nb==0)
        return 0;
    return total.somme / (float)total.nb;
}

Point getVecteurBezier(half_edge* he, float t, float dif)
{
    Point P0 = he->origine;
    Point P1 = he->controleO;
    Point P2 = he->controleI;
    Point P3 = he->incident;
    Point P = {0,0,0};  
    P = addP(P,P0,(1-t)*(1-t)*(1-t));
    P = addP(P,P1,3*t*(1-t)*(1-t));
    P = addP(P,P2,3*t*t*(1-t));
    P = addP(P,P3,t*t*t);
    Point Pt = {0,0,0};
    t+=dif;
    Pt = addP(Pt,P0,(1-t)*(1-t)*(1-t));
    Pt = addP(Pt,P1,3*t*(1-t)*(1-t));
    Pt = addP(Pt,P2,3*t*t*(1-t));
    Pt = addP(Pt,P3,t*t*t);
    return {Pt.x-P.x,Pt.y-P.y,Pt.z-P.z};
}

float densite2(heFace f)
{
    cout<<"Start Densite"<<endl;
    half_edge* he1 = f.incidente;
    int id_ori = he1->id;
    half_edge* he2 = he1->previous;
    Point P = he1->origine;
    Point P2 = he1->controleO;
    Point P3 = he2->controleI;
    Point V1 = {P2.x-P.x,P2.y-P.y,P2.z-P.z};
    Point V2 = {P3.x-P.x,P3.y-P.y,P3.z-P.z};
    float n1 = dist(V1.x,V1.y,V1.z);
    float n2 = dist(V2.x,V2.y,V2.z);
    while(n1==0 || n2==0)
    {
        he2 = he1;
        he1 = he1->next;
        if(he1->id!=id_ori)
            return -1;
        P = he1->origine;
        P2 = he1->controleO;
        P3 = he2->controleI;
        V1 = {P2.x-P.x,P2.y-P.y,P2.z-P.z};
        V2 = {P3.x-P.x,P3.y-P.y,P3.z-P.z};
        n1 = dist(V1.x,V1.y,V1.z);
        n2 = dist(V2.x,V2.y,V2.z);
    }
    cout<<"Angle Min initial trouvé"<<endl;
    he1 = f.incidente;
    float l = 0;
    int nbC = 0;
    do
    {
        l+= courbeBezier(he1->origine,he1->controleO,he1->controleI,he1->incident,false);
        nbC++;
        he1 = he1->next;
    }while(he1->id!=id_ori);
    cout<<"Fin Calcul Longueur"<<endl;
    if(l<0.01*nbC)
    {
        glColor3f(0,0,1);
        he1 = f.incidente;
        do
        {
            l+= courbeBezier(he1->origine,he1->controleO,he1->controleI,he1->incident,true);
            he1 = he1->next;
        }while(he1->id!=id_ori);
        return -1;
    }
    float dp = V1.x*V2.x + V1.y*V2.y + V1.z*V1.z;
    float angleMin = abs((M_PI/2-M_PI/8)-acos(dp / (n1*n2)));
    while(angleMin<0)
    {
        angleMin+=M_PI*2;
    }
    while(angleMin>M_PI*2)
    {
        angleMin-=M_PI*2;
    }
    float limite = 0.005;
    half_edge* he = he1;
    do
    {
        cout<<angleMin<<endl;
        he2 = he1;
        he1 = he1->next;
        V1 = getVecteurBezier(he1,0,limite*2);
        V2 = getVecteurBezier(he2,1-limite*2,limite*2);
        n1 = dist(V1.x,V1.y,V1.z);
        n2 = dist(V2.x,V2.y,V2.z);
        if(n1!=0 && n2!=0)
        {
            dp = V1.x*V2.x + V1.y*V2.y + V1.z*V1.z;
            float angle = abs(M_PI/2-acos(dp / (n1*n2)));
            if(!isnan(angle))
            {
                while(angle<0)
                {
                    angle+=M_PI*2;
                }
                while(angle>M_PI*2)
                {
                    angle-=M_PI*2;
                }
                if(angle<angleMin)
                {
                    angleMin = angle;
                    he = he1;
                }
            }
            
        }
    }while(he1->id!=id_ori);
    cout<<"Angle Min initial"<<endl;
    he1 = he;
    he2 = he1->previous;
    cout<<angleMin<<endl;
    P = he1->origine;
    
    Point v0 = getVecteurBezier(he1,0,limite*2);
    Point v1 = getVecteurBezier(he2,1-limite*2,limite*2);
    cout<<"P  "<<P.x<<" "<<P.y<<" "<<P.z<<endl;
    cout<<"v0 "<<v0.x<<" "<<v0.y<<" "<<v0.z<<endl;
    cout<<"v1 "<<v1.x<<" "<<v1.y<<" "<<v1.z<<endl;
    
    glPointSize(10.0f);
    
    Point Po = {0,0,0};
    
    Po.x = P.x + v0.x - v1.x;
    Po.y = P.y + v0.y - v1.y;
    Po.z = P.z + v0.z - v1.z;
    cout<<"Po  "<<Po.x<<" "<<Po.y<<" "<<Po.z<<endl;
    Point Pc = distanceMinBezier2(f,Po,limite/2);
    //glColor3f(0,1,1);
    
    glBegin(GL_POINTS);

    //glVertex3f(Po.x,Po.y,Po.z);
    if(func(Po.x,Po.y,Po.z)==0)
    {
        cout<<"w"<<endl;
        Point V = {Pc.x-Po.x,Pc.y-Po.y,Pc.z-Po.z};
        float n = dist(V.x,V.y,V.z);
        V.x/=n;
        V.y/=n;
        V.z/=n;
        while(func(Po.x,Po.y,Po.z)==0)
        {
            cout<<"p"<<endl;
            Po.x+=V.x*0.01;
            Po.y+=V.y*0.01;
            Po.z+=V.z*0.01;
        }
        Po.x+=V.x*0.01;
        Po.y+=V.y*0.01;
        Po.z+=V.z*0.01;
        //glVertex3f(Po.x,Po.y,Po.z);
    }
    cout<<"Point init != 0"<<endl;
    float c = 1.0f;
    Pc = distanceMinBezier2(f,Po,limite/2);
    while(dist2(Pc,Po)<limite)
    {
        cout<<"c"<<endl;
        glColor3f(c,c,c);
        c*=0.8;
        Point V = {P.x-Po.x,P.y-Po.y,P.z-Po.z};
        Point V2 = {Pc.x-Po.x,Pc.y-Po.y,Pc.z-Po.z};
        float n = dist(V.x,V.y,V.z);
        V.x/=n;
        V.y/=n;
        V.z/=n;
        n = dist(V2.x,V2.y,V2.z);
        V2.x/=n;
        V2.y/=n;
        V2.z/=n;
        Po.x-=V.x*limite/5 + V2.x*limite/5;
        Po.y-=V.y*limite/5 + V2.y*limite/5;
        Po.z-=V.z*limite/5 + V2.z*limite/5;
        //glVertex3f(Po.x,Po.y,Po.z);
        Pc = distanceMinBezier2(f,Po,limite/2);
    }
    glEnd();
    cout<<"Point init trouvé"<<endl;
    
    marque.clear();
    pile.push_back({0,0,0});
    
    glColor3f(1,1,1);
    courbeBezier(he1->origine,he1->controleO,he1->controleI,he1->incident,true);
    courbeBezier(he2->origine,he2->controleO,he2->controleI,he2->incident,true);
    glColor3f(0,0,1);
    glPointSize(1.0f);
    glBegin(GL_POINTS);
    RetourDensite rd = {0,0};
    int xmin = 0;
    int xmax = 0;
    int ymin = 0;
    int ymax = 0;
    bool plopX = false;
    bool plopY = false;
    while(pile.size()!=0)
    {
        Point Pc = pile[0];
        int x = (int)Pc.x;
        int y = (int)Pc.y;
        if(x>xmax)
        {
            xmax = x;
            plopX = true;
        }
        if(x<xmin)
        {
            xmin = x;
            plopX = true;
        }
        if(y<ymin)
        {
            ymin = y;
            plopY = true;
        }
        if(y>ymax)
        {
            ymax = y;
            plopY = true;
        }
        if((x%100==0 && plopX) || (y%100==0 && plopY))
            cout<<xmin<<" "<<ymin<<" "<<xmax<<" "<<ymax<<" "<<x<<" "<<y<<endl;
        if(xmax - xmin > 1600 || ymax - ymin > 1600)
        {
            cout<<"Plop"<<endl;
            break;
        }
        plopX = false;
        plopY = false;
        string k = to_string(x)+":"+to_string(y);
        //cout<<k;
        if(marque.find(k) == marque.end())
        {
            Point P = {Po.x+x*limite/4,Po.y+y*limite/4,Po.z};
            //cout<<" new "<<P.x<<" "<<P.y<<" "<<P.z;
            marque[k]=true;
            glPointSize(0.5f);
            Point B = distanceMinBezier2(f,P,limite/2);
            float value = val(P.x,P.y,P.z);
            if(value!=0)
            {
                rd.nb += 1;
                rd.somme += value;
                if(dist2(B,P)>=limite)
                {
                    
                    //glColor3f(0,1,0);
                    pile.push_back({(float)x-1,(float)y,0});
                    pile.push_back({(float)x+1,(float)y,0});
                    pile.push_back({(float)x,(float)y-1,0});
                    pile.push_back({(float)x,(float)y+1,0});
                }
                else
                {
                    glPointSize(1.0f);
                    //glColor3f(0,0,1);
                    //cout<<" limite";
                }
            }
            else
            {
                glPointSize(1.0f);
                //glColor3f(1,1,0);
                //cout<<" dehors";
            }
            glVertex3f(P.x,P.y,P.z);
            
        }
        //cout<<" "<<pile.size();
        pile.erase(pile.begin());
        //cout<<" "<<pile.size()<<endl;
    }
    glEnd();
    glColor3f(0,1,0);
    glPointSize(5.0f);
    glBegin(GL_POINTS);
    glVertex3f(Po.x,Po.y,Po.z);
    glColor3f(1,1,0);
    glVertex3f(P.x,P.y,P.z);
    P2 = he2->incident;
    glColor3f(0,1,1);
    glVertex3f(P2.x,P2.y,P2.z);
    glEnd();
    if(rd.nb==0)
        return -1;
    return rd.somme / (float)rd.nb;
}

//------------------------------------------------------
void affichage(void)
{
    debug.open ("debug",ios::out | ios::trunc);
    Faces.clear();
    hes.clear();
    cout<<"1"<<endl;
    color = 0;
    surface.clear();
    listeFace.clear();
    targets.clear();
    courbes.clear();
    courbeFinale.clear();
    courbeGrad.clear();
    glPointSize(7.0f);
	glMatrixMode(GL_MODELVIEW);
  /* effacement de l'image avec la couleur de fond */
	glClear(GL_COLOR_BUFFER_BIT);
	glPushMatrix();
	glTranslatef(decX,decY,cameraDistance);
glRotatef(cameraAngleX,1.,0.,0.)	;
glRotatef(cameraAngleY,0.,1.,0.);
	affiche_repere();
    float n = 1000;

    for(int i=-10; i<=10; i++)
    {
        Point P = {-1.0f ,i*0.1f ,0};
        if(i==10)
        {
            pointsImportants.push_back(surface.size());
        }
        surface.push_back(P);
    }
    for(int i=-10; i<=10; i++)
    {
        Point P = {i*0.1f ,1.0f ,0};
        if(i==10)
        {
            pointsImportants.push_back(surface.size());
        }
        surface.push_back(P);
    }
    for(int i=10; i>=-10; i--)
    {
        Point P = {1.0f ,i*0.1f,0};
        if(i==-10)
        {
            pointsImportants.push_back(surface.size());
        }
        surface.push_back(P);
    }
    for(int i=10; i>=-10; i--)
    {
        Point P = {i*0.1f ,-1.0f ,0};
        if(i==-10)
        {
            pointsImportants.push_back(surface.size());
        }
        surface.push_back(P);
    }
    float min = 1;
    float max = 0;
    for(int i=-n; i<=n; i++)
    {
        for(int j=-n; j<=n; j++)
        {
            float var = func((float)(i)/(float)(n),(float)(j)/(float)(n),0);
            if(var<min)
            {
                min = var;
            }
            if(var>max)
            {
                max=var;
            }
        }
    }
    double GradMax = -1;
    double GradMin = -1;
    for(int i=-n; i<n; i++)
    {
        for(int j=-n; j<n; j++)
        {
            float xb = (float)(i)/(float)(n);
            float yb = (float)(j)/(float)(n);
            Point P = {xb,yb,0};
            Repere R = getRepere(P);
            double m = dist(R.N.x,R.N.y,R.N.z);
            if(GradMax==-1 || m>GradMax)
                GradMax = m;
            if(GradMin==-1 || m<GradMin)
                GradMin = m;
            float dn = 1.0f/(float)n;
            float var = func(xb,yb,0);
            glBegin(GL_QUADS);
            glColor3f((var-min)/(max-min),0,0);
            glVertex3f(xb,yb,0);

            var = func(xb+dn,yb,0);
            glColor3f((var-min)/(max-min),0,0);
            glVertex3f(xb+dn,yb,0);

            var = func(xb+dn,yb+dn,0);
            glColor3f((var-min)/(max-min),0,0);
            glVertex3f(xb+dn,yb+dn,0);

            var = func(xb,yb+dn,0);
            glColor3f((var-min)/(max-min),0,0);
            glVertex3f(xb,yb+dn,0);
            glEnd();
        }
    }
    double dt = 1;
    double off = 0;
    if(GradMax!=GradMin)
    {
        dt = (normeMin-normeMax)/(GradMin-GradMax);
        off = 0.01 - dt*GradMin;
    }
    else
    {
        dt = 1/(GradMax*10);
    }
    cout<<GradMin<<" "<<GradMax<<" "<<dt<<"x + "<<off<<endl;
    newMarchingSquare(0.5,0.1);
    /*glBegin(GL_LINES);
    glColor3f(1,1,0);
    glVertex3f(1,-1,0);
    glVertex3f(-1,1,0);
    glEnd();
    vector<PointCourbe> cG;
    Point P = {pX,pY,0};
    color = 1;
    courbeGradient(P,0,{0,0,0},&cG,dt,off);
    //Repere R = getRepere(P);
    //P = {0.001,-0.001,0};
    Repere R = getRepere(P);
    glBegin(GL_LINES);
    glColor3f(0,1,0);
    glVertex3f(P.x,P.y,P.z);
    glVertex3f(P.x+R.T.x*0.1,P.y+R.T.y*0.1,P.z+R.T.z*0.1);
    glColor3f(0,0,1);
    glVertex3f(P.x,P.y,P.z);
    glVertex3f(P.x+R.N.x*0.1,P.y+R.N.y*0.1,P.z+R.N.z*0.1);
    */
    glEnd();/*
    P = {0.001,0.001,0};
    R = getRepere(P);
    glBegin(GL_LINES);
    glColor3f(1,1,0);
    glVertex3f(P.x,P.y,P.z);
    glVertex3f(P.x+R.T.x*0.1,P.y+R.T.y*0.1,P.z+R.T.z*0.1);
    glColor3f(0,0,0);
    glVertex3f(P.x,P.y,P.z);
    glVertex3f(P.x+R.N.x*0.1,P.y+R.N.y*0.1,P.z+R.N.z*0.1);
    glEnd();
    P = {-0.001,-0.001,0};
    R = getRepere(P);
    glBegin(GL_LINES);
    glColor3f(0,1,1);
    glVertex3f(P.x,P.y,P.z);
    glVertex3f(P.x+R.T.x*0.1,P.y+R.T.y*0.1,P.z+R.T.z*0.1);
    glColor3f(1,1,1);
    glVertex3f(P.x,P.y,P.z);
    glVertex3f(P.x+R.N.x*0.1,P.y+R.N.y*0.1,P.z+R.N.z*0.1);
    glEnd();
    /*
    targets.push_back(0.1);
    targets.push_back(0.2);
    targets.push_back(0.3);
    targets.push_back(0.4);
    targets.push_back(0.5);
    targets.push_back(0.6);
    targets.push_back(0.7);
    targets.push_back(0.8);
    targets.push_back(0.9);
    //separation(0.05);
    cout<<"2"<<endl;
    for(auto it = targets.begin(); it!=targets.end(); it++)
    {
        float v = *it;
        Point C = {v*v,v*v,v*v};
        lancerMarchingSquare(v,C);
    }
    /*
    int nbC = 0;
    float minC = 1;
    float minI = -1;
    float maxC = 0;
    float maxI = -1;
    do{
        nbC=0;
        minC = 1;
        maxC = 0;
        int i = 0;
        for(auto it = targets.begin(); it!=targets.end(); it++)
        {
            float v = *it;
            vector<vector<Point>> c = courbes[v*1000];
            nbC+=c.size();
            if(c.size()!=0)
            {
                if(v<minC)
                {
                    minC = v;
                    minI = i;
                }
                if(v>maxC)
                {
                    maxC = v;
                    maxI = i;
                }
            }
            i++;
        }
        //cout<<endl<<"2.5 "<<nbC<<" "<<minI<<" "<<minC<<" "<<maxI<<" "<<maxC<<endl;
        if(minI==maxI)
        {
            minI--;
            maxI++;
        }
        if(nbC<6)
        {
            i=0;
            for(auto it = targets.begin(); it!=targets.end(); it++)
            {
                if(i>=minI && i<maxI)
                {
                    float v = *it;
                    it++;
                    float v2 = *it;
                    float nv = (v+v2)/2;
                    //cout<<v<<" "<<nv<<" "<<v2<<endl;
                    it = targets.insert(it,nv);
                    Point C = {nv*nv,nv*nv,nv*nv};
                    lancerMarchingSquare(nv,C);
                }
                i++;
            }
        }
        //cout<<endl;
        for(auto it = targets.begin(); it!=targets.end(); it++)
        {
            float v = *it;
            //cout<<v<<endl;
        }
    }while(nbC<6);
    cout<<"2 fin"<<endl;
    int s = 0;
    int nc = -1;
    float nt = -1;
    
    for(auto it = targets.begin(); it!=targets.end(); it++)
    {
        float v = *it;
        vector<vector<Point>> c = courbes[v*1000];
        int k=0;
        for(auto it2 = c.begin(); it2!=c.end(); it2++)
        {
            if(it2->size()>s)
            {
                s = it2->size();
                nc = k;
                nt = v;
            }
            k++;
        }
    }
    if(nc!=-1)
    {
        vector<vector<Point>> ct = courbes[nt*1000];
        vector<Point> c = ct[nc];
        int div = 10;
        if(s<20)
        {
            div = 4;
        }
        int ps = (s-2)/div;
        
        Point P = c[1];
        vector<PointCourbe> cG;
        //cout<<nt<<endl;
        cG.push_back({P,nt,-1,-1,-1,-1});
        courbeGradient(P,0,{0,0,0},&cG);
        courbeGrad.push_back(cG);
        for(int i=1; i<div ; i++)
        {
            Point P = c[ps*i];
            vector<PointCourbe> cG1;
            //cout<<nt<<endl;
            cG1.push_back({P,nt,-1,-1,-1,-1});
            courbeGradient(P,0,{0,0,0},&cG1);
            //cout << i << " " << cG1.size()<<endl;
            color+=0.05;
            courbeGrad.push_back(cG1);
            //cout<<endl<<endl;
        }
        
        P = c[s-2];
        vector<PointCourbe> cG2;
        //cout<<nt<<endl;
        cG2.push_back({P,nt,-1,-1,-1,-1});
        courbeGradient(P,0,{0,0,0},&cG2);
        courbeGrad.push_back(cG2);
        
    }
    /*
    //cout<<endl<<endl;;
    for(int i=0; i<courbeGrad.size(); i++)
    {
        vector<PointCourbe> cG = courbeGrad[i];
        for(int j=0; j<cG.size(); j++)
        {
            PointCourbe Pc = cG[j];
            Pc.graC = i;
            Pc.indGra = j;
            vector<vector<Point>> ct = courbes[Pc.val*1000];
            float min = -1;
            int nc = -1;
            Point Pmin = {0,0,0};
            for(int k=0; k<ct.size(); k++)
            {
                vector<Point> c = ct[k];
                for(int l=0; l<c.size(); l++)
                {
                    Point P = c[l];
                    float d = dist2(Pc.P,P);
                    if(min==-1 || d<min)
                    {
                        min = d;
                        nc = k;
                        Pmin = P;
                    }
                }
            }
            Pc.isoC = nc;
            vector<vector<PointCourbe>> ct2 = courbeFinale[Pc.val*1000];
            if(Pc.isoC==-1 && Pc.val==0)
            {
                Pc.isoC=0;
            }
            else
            {
                if(Pc.isoC!=-1)
                    Pc.P = Pmin;
            }
            //cout<<Pc.val<<" "<<Pc.isoC<<" "<<ct2.size()<<endl;
            while(ct2.size()<=Pc.isoC && Pc.isoC!=-1)
            {
                vector<PointCourbe> c;
                ct2.push_back(c);
            }
            if(Pc.val==0)
            {
                
                if(ct2[Pc.isoC].size()%2==0)
                {
                    Pc.isoInd = i;
                }
                else
                {
                    Pc.isoInd = 2*courbeGrad.size()-i-1;
                }
                auto it = ct2[Pc.isoC].begin();
                while(it!=ct2[Pc.isoC].end() && (*it).isoInd<Pc.isoInd)
                {
                    it++;
                }
                ct2[Pc.isoC].insert(it,Pc);
                
            }
            else
            {
                if(ct2[Pc.isoC].size()>0)
                {
                    PointCourbe Pc2 = ct2[Pc.isoC][ct2[Pc.isoC].size()-1];
                    if(Pc2.graC != Pc.graC || dist2(Pc.P,Pc2.P)>sqrt(0.02*0.02*2))
                    {
                        Pc.isoInd = ct2[Pc.isoC].size();
                        ct2[Pc.isoC].push_back(Pc);
                    }
                    else
                    {
                        Pc.isoInd=-1;
                    }
                }
                else
                {
                    Pc.isoInd = 0;
                    ct2[Pc.isoC].push_back(Pc);
                }
                
            }
            courbeFinale[Pc.val*1000] = ct2;
            cG[j] = Pc;
        }
        courbeGrad[i]=cG;
        //cout<<endl;
    }
    cout<<"3"<<endl;
    /*for(int i=0; i<targets.size(); i++)
    {
        vector<vector<PointCourbe>> cs = courbeFinale[targets[i]*1000];
        for(int j=0; j<cs.size(); j++)
        {
            vector<Point> oriC = courbes[targets[i]*1000][j];
            vector<PointCourbe> c = cs[j];
            PointCourbe Pc = c[0];
            Point P = oriC[1];
            int ind1 = Pc.graC;
            vector<PointCourbe> cG;
            Pc = {P,targets[i],j,0,ind1,0};
            c.insert(c.begin(),Pc);
            cG.push_back(Pc);
            courbeGradient(P,0,{0,0,0},&cG);
            cout<<targets[i]<<" "<<j<<" first "<<cG.size()<<endl;
            color+=0.05;
            PointCourbe pc = cG[0];
            for(int h=1; h<cG.size(); h++)
            {
                //cout<<"h "<<h<<endl;
                if(pc.val == cG[i].val && dist2(pc.P,cG[1].P)>sqrt(0.02*0.02*2))
                {
                    //cout<<"h1"<<endl;
                    pc = cG[i];
                    pc.isoInd = -1;
                    cG[i]=pc;
                }
                else
                {
                    //cout<<"h2"<<endl;
                    pc = cG[i];
                }
                //cout<<"hfin"<<endl;
            }
            auto it = courbeGrad.begin();
            for(int h=0; h<ind1; h++)
            {
                it++;
            }
            courbeGrad.insert(it,cG);
            Pc = c[c.size()-1];
            P = oriC[oriC.size()-2];
            int ind2 = Pc.graC;
            Pc = {P,targets[i],j,c.size(),ind2,0};
            c.push_back(Pc);
            vector<PointCourbe> cG2;
            cG2.push_back(Pc);
            courbeGradient(P,0,{0,0,0},&cG2);
            cout<<targets[i]<<" "<<j<<" second "<<cG2.size()<<endl;
            color+=0.05;
            PointCourbe pc2 = cG2[0];
            for(int h=1; h<cG2.size(); h++)
            {
                if(pc2.val == cG2[i].val && dist2(pc2.P,cG2[1].P)>sqrt(0.02*0.02*2))
                {
                    pc2 = cG2[i];
                    pc2.isoInd = -1;
                    cG2[i]=pc2;
                }
                else
                {
                    pc2 = cG2[i];
                }
            }
            it = courbeGrad.begin();
            for(int h=0; h<ind2; h++)
            {
                it++;
            }
            it++;
            courbeGrad.insert(it,cG2);
            for(int h=0; h<c.size(); h++)
            {
                Pc = c[h];
                Pc.isoInd = h;
                c[h]=Pc;
            }
            cs[j]=c;
            for(int k=0; k<targets.size(); k++)
            {
                vector<vector<PointCourbe>> cs2 = courbeFinale[targets[i]*1000];
                for(int l=0; l<cs.size(); l++)
                {
                    if(k!=i && l!=j)
                    {
                        vector<PointCourbe> c = cs2[j];
                        PointCourbe Pc = c[0];
                        if(Pc.graC>=ind1)
                        {
                            Pc.graC = Pc.graC+1;
                        }
                        if(Pc.graC>=ind2)
                        {
                            Pc.graC = Pc.graC+1;
                        }
                        Pc = c[c.size()-1];
                        if(Pc.graC>=ind1)
                        {
                            Pc.graC = Pc.graC+1;
                        }
                        if(Pc.graC>=ind2)
                        {
                            Pc.graC = Pc.graC+1;
                        }
                    }
                    
                }
            }
        }
        courbeFinale[targets[i]*1000] = cs;
    }*/
    /*
    cout<<"4"<<endl;
    for(int i=0; i<courbeGrad.size();i++)
    {
        vector<PointCourbe> cG = courbeGrad[i];
        int h = 0;
        for(auto it = cG.begin(); it!=cG.end(); it++)
        {
            PointCourbe Pc = *it;
            //cout<<Pc.val<<"\t"<<Pc.isoC<<"\t"<<Pc.isoInd<<"\t"<<Pc.graC<<"\t"<<Pc.indGra<<endl;
            if(Pc.isoInd==-1)
            {
                int ind = Pc.indGra;
                cG.erase(it);
                it--;
                h--;
            }
            else
            {
                if(Pc.indGra!=h)
                {
                    Pc.indGra=h;
                    cG[h]=Pc;
                    //cout<<Pc.indGra<<endl;
                }
            }
            h++;
        }
        courbeGrad[i]=cG;
    }
    cout<<"5"<<endl;
    /*for(int i=0; i<courbeGrad.size(); i++)
    {
        vector<PointCourbe> cG = courbeGrad[i];
        for(int j=0; j<cG.size(); j++)
        {
            PointCourbe Pc = cG[j];
            Pc.indGra = j;
            Pc.graC = i;
            if(Pc.val==0)
            {
                vector<vector<PointCourbe>> ct = courbeFinale[0];
                if(ct.size()==0)
                {
                    vector<PointCourbe> c;
                    ct.push_back(c);
                }
                if(ct[Pc.isoC].size()%2==0)
                {
                    Pc.isoInd = i;
                }
                else
                {
                    Pc.isoInd = 2*courbeGrad.size()-i-1;
                }
                auto it = ct[Pc.isoC].begin();
                while(it!=ct[Pc.isoC].end() && (*it).isoInd<Pc.isoInd)
                {
                    it++;
                }
                ct[Pc.isoC].insert(it,Pc);
                courbeFinale[0] = ct;
            }
            cG[j]=Pc;
        }
    }*/
    cout<<"6"<<endl;
    /*
    cout<<endl<<endl;
    //cout<<courbeFinale[0][0].size()<<endl<<endl;
    for(int i=0; i<courbeGrad.size();i++)
    {
        vector<PointCourbe> cG = courbeGrad[i];
        for(auto it = cG.begin(); it!=cG.end(); it++)
        {
            PointCourbe Pc = *it;
            cout<<Pc.val<<"\t"<<Pc.isoC<<"\t"<<Pc.isoInd<<"\t"<<Pc.graC<<"\t"<<Pc.indGra<<endl;
        }
        cout<<endl<<endl;
    }
    cout<<endl<<endl;
    */
    /*
    targets.insert(targets.begin(),0);
    for(int i=0; i<targets.size(); i++)
    {
        //cout<<endl<<targets[i]<<endl;
        vector<vector<PointCourbe>> cs = courbeFinale[targets[i]*1000];
        for(int j=0; j<cs.size(); j++)
        {
            //cout<<"Courbe "<<j<<endl;
            vector<PointCourbe> c = cs[j];
            glPointSize(8.0f);
            glBegin(GL_POINTS);
            for(int k=0; k<c.size(); k++)
            {
                PointCourbe Pc = c[k];
                glColor3f((float)k/(float)c.size(),(float)k/(float)c.size(),(float)k/(float)c.size());
                glVertex3f(Pc.P.x,Pc.P.y,Pc.P.z);
                //cout<<k<<" "<<Pc.val<<"\t"<<Pc.isoC<<"\t"<<Pc.isoInd<<"\t"<<Pc.graC<<"\t"<<Pc.indGra<<"\t"<<Pc.P.x<<"\t"<<Pc.P.y<<"\t"<<Pc.P.z<<endl;
            }
            glEnd();
        }
    }

    
    cout<<"7"<<endl;
    /*
    for(int i=0; i<courbeGrad.size()-1; i++)
    {
        vector<PointCourbe> cG = courbeGrad[i];
        for(int j=0; j<cG.size()-1; j++)
        {
            int nbNoP = 0;
            bool noP = false;
            PointCourbe Pc = cG[j];
            vector<PointCourbe> sf;
            sf.push_back(Pc);
            bool finIso = true;
            bool cont = true;
            float valFin = 0;
            int step = 0;
            //cout<<"\nFace\n"<<step<<"\t"<<Pc.val<<"\t"<<Pc.isoC<<"\t"<<Pc.isoInd<<"\t"<<Pc.graC<<"\t"<<Pc.indGra<<endl;
            PointCourbe Pc2 = Pc;
            do
            {
                noP = false;
                if(step==0 || (step>=2 && cont))
                {
                    //cout<<"Gra";
                    vector<PointCourbe> grad = courbeGrad[Pc2.graC];
                    if(step<2)
                    {
                        //cout<<"0";
                        if(Pc2.indGra<grad.size()-1)
                        {
                            //cout<<" sens";
                            Pc2 = grad[Pc2.indGra+1];
                        }
                        else
                        {
                            //cout<<" opposé";
                            PointCourbe Pc3 = Pc2;
                            vector<PointCourbe> iso = courbeFinale[Pc2.val*1000][Pc3.isoC];
                            while(Pc3.isoInd==iso.size()-1)
                            {

                                Pc3 = grad[Pc3.indGra-1];
                                iso = courbeFinale[Pc2.val*1000][Pc3.isoC];
                                //cout<<"Pc3\t"<<Pc3.val<<"\t"<<Pc3.isoC<<"\t"<<Pc3.isoInd<<"\t"<<Pc3.graC<<"\t"<<Pc3.indGra<<endl;
                            }
                            Pc2 = iso[Pc3.isoInd+1];
                        }
                    }
                    else
                    {
                        //cout<<"1";
                        if(Pc2.indGra>0)
                        {
                            //cout<<" sens";
                            Pc2 = grad[Pc2.indGra-1];
                        }
                        else
                        {
                            //cout<<" opposé";
                            /*vector<PointCourbe> iso = courbeFinale[Pc2.val*1000][Pc2.isoC];
                            PointCourbe Pc3 = iso[Pc2.isoInd-1];
                            grad = courbeGrad[Pc3.graC];
                            Pc2 = grad[Pc3.indGra-1];*/
    /*
                            Pc2 = grad[grad.size()-1];
                        }
                    }
                    //cout<<endl;
                }
                else
                {
                    //cout<<"iso";
                    vector<PointCourbe> iso = courbeFinale[Pc2.val*1000][Pc2.isoC];
                    if(step<2)
                    {
                        //cout<<"0";
                        if(Pc2.val==0)
                        {
                            if(Pc2.isoInd>0)
                            {
                                //cout<<" sens";
                                Pc2 = iso[Pc2.isoInd-1];
                            }
                            else
                            {
                                //cout<<" opposé";
                                noP = true;
                                nbNoP++;
                                step-=2;
                            }
                        }
                        else
                        {
                            if(Pc2.isoInd<iso.size()-1)
                            {
                                //cout<<" sens";
                                Pc2 = iso[Pc2.isoInd+1];
                            }
                            else
                            {
                                //cout<<" opposé";
                                noP = true;
                                nbNoP++;
                                step-=2;
                            }
                        }
                        
                    }
                    else
                    {
                        //cout<<"1";
                        if(Pc2.isoInd>0)
                        {
                            //cout<<" sens";
                            Pc2 = iso[Pc2.isoInd-1];
                        }
                        else
                        {
                            //cout<<" opposé";
                            Pc2 = iso[iso.size()-1];
                        }
                    }
                    //cout<<endl;
                }
                if(!finIso)
                {
                    if(Pc2.val==valFin)
                    {
                        cont=false;
                    }
                    if(Pc2.graC==Pc.graC)
                        break;
                }
                if(step==1)
                {
                    bool verif =false;
                    vector<PointCourbe> cgv = courbeGrad[Pc2.graC];
                    for(int l=0; l<cgv.size(); l++)
                    {
                        PointCourbe Pc3 = cgv[l];
                        if(abs(Pc3.val-Pc.val)<0.001)
                        {
                            verif=true;
                        }
                    }
                    if(!verif)
                    {
                        finIso=false;
                        valFin = cG[j-1].val;
                    }
                }
                step++;
                if(!noP)
                {
                    sf.push_back(Pc2);
                    //cout<<step<<"\t"<<Pc2.val<<"\t"<<Pc2.isoC<<"\t"<<Pc2.isoInd<<"\t"<<Pc2.graC<<"\t"<<Pc2.indGra<<endl;
                }
            }while((Pc2.isoC!=Pc.isoC || (Pc.val==0 && abs(Pc.isoInd-Pc2.isoInd)>1)) || (Pc2.val != Pc.val && finIso));
            if(Pc.isoInd==Pc2.isoInd)
            {
                sf.pop_back();
            }
            Face f;
            f.sommets = sf;
            listeFace.push_back(f);
            j+=nbNoP;
        }
    }
    cout<<"8"<<endl;
    Face f0;
    Face f1;
    vector<PointCourbe> c0 = courbeGrad[0];
    for(int i=c0.size()-1; i>=0; i--)
    {
        PointCourbe Pc = c0[i];
        f0.sommets.push_back(Pc);
    }
    PointCourbe Pc0 = c0[0];
    PointCourbe Pc1 = c0[c0.size()-1];
    float min0 = dist2(Pc0.P,surface[0]);
    int i0 = 0;
    float min02 = dist2(Pc0.P,surface[1]);
    int i02 = 1;
    if(min02<min0)
    {
        i0=1;
        i02=0;
        float t = min0;
        min0 = min02;
        min02 = t;
    }
    float min1 = dist2(Pc1.P,surface[0]);
    int i1 = 0;
    float min12 = dist2(Pc1.P,surface[1]);
    int i12 = 1;
    if(min12<min1)
    {
        i1=1;
        i12=0;
        float t = min1;
        min1 = min12;
        min12 = t;
    }
    for(int i=0; i<surface.size();i++)
    {
        Point P = surface[i];
        float d0 = dist2(P,Pc0.P);
        if(d0<min0)
        {
            min02 = min0;
            i02 = i0;
            i0 = i;
            min0 = d0;
        }
        else
        {
            if(d0<min02)
            {
                i02 = i;
                min02 = d0;
            }
        }
        float d1 = dist2(P,Pc1.P);
        if(d1<min1)
        {
            min12 = min1;
            i12 = i1;
            i1 = i;
            min1 = d1;
        }
        else
        {
            if(d1<min12)
            {
                i12 = i;
                min12 = d1;
            }
        }
    }
    cout<<"face0 Debut : "<<i0<<" "<<i02<<endl;
    cout<<"face0 Fin : "<<i1<<" "<<i12<<endl;
    int dt = 1;
    int nb = 0;
    for(int i=i0; i!=i1; i++)
    {
        if(i==surface.size())
        {
            i=0;
        }
        nb++;
    }
    int nb2 = 0;
    for(int i=i0; i!=i1; i--)
    {
        if(i==-1)
        {
            i=surface.size()-1;
        }  
        nb2++;
    }
    if(nb2<nb)
    {
        dt=-1;
    }
    int no = 0;
    for(int i=i0; i!=i1; i+=dt)
    {
        cout<<i<<endl;
        if(i==-1)
        {
            i=surface.size()-1;
        }
        if(i==surface.size())
        {
            i=0;
        }
        for(int j=0; j<pointsImportants.size(); j++)
        {
            if(i==pointsImportants[j])
            {
                 f0.sommets.push_back({surface[i],0,1,no,-1,-1});
                no++;
            }
        }
       
    }
    /*Point P0 = surface[i0];
    Point P1 = surface[i02];
    float dist0 = dist2(P0,Pc1.P);
    float dist1 = dist2(P1,Pc1.P);
    if(dist0<dist1)
    {
        int dt = 1;
        cout<<"i0 plus proche"<<endl;
        Point P10 = surface[i1];
        Point P11 = surface[i12];
        float dist10 = dist2(P10,Pc0.P);
        float dist11 = dist2(P11,Pc0.P);
        int target =i1;
        if(dist11<dist10)
            target = i12;
        int nb = 0;
        for(int i=i0; i!=target; i++)
        {
            if(i==surface.size())
            {
                i=0;
            }
            nb++;
        }
        int nb2 = 0;
        for(int i=i0; i!=target; i--)
        {
            if(i==-1)
            {
                i=surface.size()-1;
            }  
            nb2++;
        }
        if(nb2<nb)
        {
            dt=-1;
        }
        
        cout<<"target "<<target<<endl;
        int no = 0;
        for(int i=i0; i!=target; i+=dt)
        {
            cout<<i<<endl;
            if(i==-1)
            {
                i=surface.size()-1;
            }
            if(i==surface.size())
            {
                i=0;
            }
            for(int j=0; j<pointsImportants.size(); j++)
            {
                if(i==pointsImportants[j])
                {
                     f0.sommets.push_back({surface[i],0,1,no,-1,-1});
                    no++;
                }
            }
           
        }
    }
    else
    {
        int dt = 1;
        cout<<"i02 plus proche"<<endl;
        Point P10 = surface[i1];
        Point P11 = surface[i12];
        float dist10 = dist2(P10,Pc0.P);
        float dist11 = dist2(P11,Pc0.P);
        int target =i1;
        if(dist11<dist10)
            target = i12;
        int nb = 0;
        for(int i=i02; i!=target; i++)
        {
            if(i==surface.size())
            {
                i=0;
            }
            nb++;
        }
        int nb2 = 0;
        for(int i=i02; i!=target; i--)
        {
            if(i==-1)
            {
                i=surface.size()-1;
            }
            nb2++;
        }
        if(nb2<nb)
        {
            dt=-1;
        }
        
        cout<<"target "<<target<<endl;
        int no = 0;
        for(int i=i02; i!=target; i+=dt)
        {
            cout<<i<<endl;
            if(i==-1)
            {
                i=surface.size()-1;
            }
            if(i==surface.size())
            {
                i=0;
            }
            for(int j=0; j<pointsImportants.size(); j++)
            {
                if(i==pointsImportants[j])
                {
                     f0.sommets.push_back({surface[i],0,1,no,-1,-1});
                    no++;
                }
            }
        }
    }
    */
    /*
    listeFace.insert(listeFace.begin(),f0);
    c0 = courbeGrad[courbeGrad.size()-1];
    for(int i=0; i<c0.size(); i++)
    {
        PointCourbe Pc = c0[i];
        f1.sommets.push_back(Pc);
    }
    
    Pc0 = c0[0];
    Pc1 = c0[c0.size()-1];
    min0 = dist2(Pc0.P,surface[0]);
    i0 = 0;
    min02 = dist2(Pc0.P,surface[1]);
    i02 = 1;
    if(min02<min0)
    {
        i0=1;
        i02=0;
        float t = min0;
        min0 = min02;
        min02 = t;
    }
    min1 = dist2(Pc1.P,surface[0]);
    i1 = 0;
    min12 = dist2(Pc1.P,surface[1]);
    i12 = 1;
    if(min12<min1)
    {
        i1=1;
        i12=0;
        float t = min1;
        min1 = min12;
        min12 = t;
    }
    for(int i=0; i<surface.size();i++)
    {
        Point P = surface[i];
        float d0 = dist2(P,Pc0.P);
        if(d0<min0)
        {
            min02 = min0;
            i02 = i0;
            i0 = i;
            min0 = d0;
        }
        else
        {
            if(d0<min02)
            {
                i02 = i;
                min02 = d0;
            }
        }
        float d1 = dist2(P,Pc1.P);
        if(d1<min1)
        {
            min12 = min1;
            i12 = i1;
            i1 = i;
            min1 = d1;
        }
        else
        {
            if(d1<min12)
            {
                i12 = i;
                min12 = d1;
            }
        }
    }
    cout<<"face0 Debut : "<<i1<<" "<<i12<<endl;
    cout<<"face0 Fin : "<<i0<<" "<<i02<<endl;
    dt = 1;
    nb = 0;
    for(int i=i1; i!=i0; i++)
    {
        if(i==surface.size())
        {
            i=0;
        }
        nb++;
    }
    nb2 = 0;
    for(int i=i1; i!=i0; i--)
    {
        if(i==-1)
        {
            i=surface.size()-1;
        }
        nb2++;
    }
    if(nb2<nb)
    {
        dt=-1;
    }
    no = 0;
    for(int i=i1; i!=i0; i+=dt)
    {
        cout<<i<<endl;
        if(i==-1)
        {
            i=surface.size()-1;
        }
        if(i==surface.size())
        {
            i=0;
        }
        for(int j=0; j<pointsImportants.size(); j++)
        {
            if(i==pointsImportants[j])
            {
                 f1.sommets.push_back({surface[i],0,1,no,-1,-1});
                no++;
            }
        }
    }
    /*P0 = surface[i1];
    P1 = surface[i12];
    dist0 = dist2(P0,Pc0.P);
    dist1 = dist2(P1,Pc0.P);
    if(dist0<dist1)
    {
        cout<<"i1 plus proche"<<endl;
        int dt = 1;
        Point P10 = surface[i0];
        Point P11 = surface[i02];
        float dist10 = dist2(P10,Pc1.P);
        float dist11 = dist2(P11,Pc1.P);
        int target =i0;
        if(dist11<dist10)
            target = i02;
        int nb = 0;
        for(int i=i1; i!=target; i++)
        {
            if(i==surface.size())
            {
                i=0;
            }
            nb++;
        }
        int nb2 = 0;
        for(int i=i1; i!=target; i--)
        {
            if(i==-1)
            {
                i=surface.size()-1;
            }
            nb2++;
        }
        if(nb2<nb)
        {
            dt=-1;
        }
        cout<<"target "<<target<<endl;
        int no = 0;
        for(int i=i1; i!=target; i+=dt)
        {
            cout<<i<<endl;
            if(i==-1)
            {
                i=surface.size()-1;
            }
            if(i==surface.size())
            {
                i=0;
            }
            for(int j=0; j<pointsImportants.size(); j++)
            {
                if(i==pointsImportants[j])
                {
                     f1.sommets.push_back({surface[i],0,1,no,-1,-1});
                    no++;
                }
            }
        }
    }
    else
    {
        cout<<"i12 plus proche"<<endl;
        int dt = -1;
        Point P10 = surface[i0];
        Point P11 = surface[i02];
        float dist10 = dist2(P10,Pc1.P);
        float dist11 = dist2(P11,Pc1.P);
        int target =i0;
        if(dist11<dist10)
            target = i02;
        int nb = 0;
        for(int i=i12; i!=target; i++)
        {
            if(i==surface.size())
            {
                i=0;
            }
            nb++;
        }
        int nb2 = 0;
        for(int i=i12; i!=target; i--)
        {
            if(i==-1)
            {
                i=surface.size()-1;
            }
            nb2++;
        }
        if(nb2<nb)
        {
            dt=-1;
        }
        cout<<"target "<<target<<endl;
        int no = 0;
        for(int i=i12; i!=target; i+=dt)
        {
            cout<<i<<endl;
            if(i==-1)
            {
                i=surface.size()-1;
            }
            if(i==surface.size())
            {
                i=0;
            }
            for(int j=0; j<pointsImportants.size(); j++)
            {
                if(i==pointsImportants[j])
                {
                     f1.sommets.push_back({surface[i],0,1,no,-1,-1});
                    no++;
                }
            }
        }
    }
    */
    //listeFace.push_back(f1);
    //cout<<nbF<<" / "<<listeFace.size()<<endl;
    
    int h = 60;
    /*
    for(int i=0; i<nbF;i++)
    {
        Face f = listeFace[i];
        Point c = hsv2rgb(h*i,1.0,1.0);
        glBegin(GL_POLYGON);
        glColor3f(c.x,c.y,c.z);
        for(int j=0; j<f.sommets.size(); j++)
        {
            PointCourbe Pc = f.sommets[j];
            if(i==nbF-1)
            {
                cout<<Pc.val<<"\t"<<Pc.isoC<<"\t"<<Pc.isoInd<<"\t"<<Pc.graC<<"\t"<<Pc.indGra<<"\t"<<Pc.P.x<<"\t"<<Pc.P.y<<"\t"<<Pc.P.z<<endl;
            }
            glVertex3f(Pc.P.x,Pc.P.y,Pc.P.z);
        }
        glEnd();
        cout<<f.sommets.size()<<endl;
        for(int j=0; j<f.sommets.size()-1;j++)
        {
            PointCourbe Pc = f.sommets[j];
            PointCourbe Pc2 = f.sommets[j+1];
            Repere R = getRepere(Pc.P);
            Repere R2 = getRepere(Pc2.P);
            int dir = 1;
            int dir2 = 1;
            float dp = dist2(Pc.P,Pc2.P);
            if(Pc.val == Pc2.val && Pc.isoC==Pc2.isoC)
            {
                float d1 = dist2(addP(Pc.P,R.T,dp*0.5),Pc2.P);
                float d2 = dist2(addP(Pc.P,R.T,-dp*0.5),Pc2.P);
                if(d2<d1)
                    dir=-1;
                d1 = dist2(addP(Pc2.P,R2.T,dp*0.5),Pc.P);
                d2 = dist2(addP(Pc2.P,R2.T,-dp*0.5),Pc.P);
                if(d2<d1)
                    dir2=-1;
                //courbeBezier(Pc.P,addP(Pc.P,R.T,dir*dp*0.5),addP(Pc2.P,R2.T,dir2*dp*0.5),Pc2.P,true);
            }
            else
            {
                float d1 = dist2(addP(Pc.P,R.N,dp*0.5),Pc2.P);
                float d2 = dist2(addP(Pc.P,R.N,-dp*0.5),Pc2.P);
                if(d2<d1)
                    dir=-1;
                d1 = dist2(addP(Pc2.P,R2.N,dp*0.5),Pc.P);
                d2 = dist2(addP(Pc2.P,R2.N,-dp*0.5),Pc.P);
                if(d2<d1)
                    dir2=-1;
                //courbeBezier(Pc.P,addP(Pc.P,R.N,dir*dp*0.5),addP(Pc2.P,R2.N,dir2*dp*0.5),Pc2.P,true);
            }
        }
        PointCourbe Pc = f.sommets[0];
        PointCourbe Pc2 = f.sommets[f.sommets.size()-1];
        Repere R = getRepere(Pc.P);
        Repere R2 = getRepere(Pc2.P);
        int dir = 1;
        int dir2 = 1;
        float dp = dist2(Pc.P,Pc2.P);
        if(Pc.val == Pc2.val && Pc.isoC==Pc2.isoC)
        {
            float d1 = dist2(addP(Pc.P,R.T,dp*0.5),Pc2.P);
            float d2 = dist2(addP(Pc.P,R.T,-dp*0.5),Pc2.P);
            if(d2<d1)
                dir=-1;
            d1 = dist2(addP(Pc2.P,R2.T,dp*0.5),Pc.P);
            d2 = dist2(addP(Pc2.P,R2.T,-dp*0.5),Pc.P);
            if(d2<d1)
                dir2=-1;
            //courbeBezier(Pc.P,addP(Pc.P,R.T,dir*dp*0.5),addP(Pc2.P,R2.T,dir2*dp*0.5),Pc2.P,true);
        }
        else
        {
            float d1 = dist2(addP(Pc.P,R.N,dp*0.5),Pc2.P);
            float d2 = dist2(addP(Pc.P,R.N,-dp*0.5),Pc2.P);
            if(d2<d1)
                dir=-1;
            d1 = dist2(addP(Pc2.P,R2.N,dp*0.5),Pc.P);
            d2 = dist2(addP(Pc2.P,R2.N,-dp*0.5),Pc.P);
            if(d2<d1)
                dir2=-1;
            //courbeBezier(Pc.P,addP(Pc.P,R.N,dir*dp*0.5),addP(Pc2.P,R2.N,dir2*dp*0.5),Pc2.P,true);
        }
        if(i==nbF-1)
            cout<<endl;
    }
    */
    /*
    int tid=0;
    cout<<"Plop0"<<endl;
    float o = -M_PI/2;
    
    for(int i=0; i<listeFace.size();i++)
    {
        float somme = 0;
        float divid = 0;
        Face f = listeFace[i];
        PointCourbe Pc0 = f.sommets[0];
        PointCourbe Pc = f.sommets[f.sommets.size()-1];
        int ik = tid;
        string str1 = to_string((int)(Pc0.val*1000))+to_string(Pc0.isoC)+to_string(Pc0.isoInd)+to_string(Pc0.graC)+to_string(Pc0.indGra);
        string str0 = to_string((int)(Pc.val*1000))+to_string(Pc.isoC)+to_string(Pc.isoInd)+to_string(Pc.graC)+to_string(Pc.indGra);
        string str = str0+":"+str1;
        string str_ori = str;
        //half_edge he0 = {ik,Pc0.P,{0,0,0},{0,0,0},Pc.P,i,&nullHe,&nullHe,&nullHe};
        hes[str] = new half_edge;
        hes[str]->id = ik;
        hes[str]->origine = Pc.P;
        hes[str]->densityO = val(Pc.P.x,Pc.P.y,Pc.P.z);
        hes[str]->incident = Pc0.P;
        hes[str]->densityI = val(Pc0.P.x,Pc0.P.y,Pc0.P.z);
        hes[str]->face = i;
        hes[str]->next = &nullHe;
        hes[str]->previous = &nullHe;
        hes[str]->opposite = &nullHe;
        heFace hef = {i,0,hes[str]};
        Repere R = getRepere(Pc0.P);
        if(Pc0.val==0)
        {
            R = getRepere0(R);
        }
        if(hes[str]->densityI == 0)
        {
            float t = 0.01;
            float dt = -0.01;
            float val = 0;
            float dir = -1;
            Point N = {0,0,0};
            N.x = R.T.x * cos(o) - R.T.y * sin(o);
            N.y = R.T.x * sin(o) + R.T.y * cos(o);
            N.z = R.T.z;
            if(dist(Pc0.P.x+t*N.x,Pc0.P.y+t*N.y,Pc0.P.z+t*N.z)<dist(Pc0.P.x,Pc0.P.y,Pc0.P.z))
            {
                dt*=-1;;
            }
            glBegin(GL_POINTS);
            glPointSize(1.0f);
            float c = 0;
            while(val==0 && abs(t)<1)
            {
                c+=0.1;
                glColor3f(c,c,c);
                val=func(Pc0.P.x+t*N.x,Pc0.P.y+t*N.y,Pc0.P.z+t*N.z);
                glVertex3f(Pc0.P.x+t*N.x,Pc0.P.y+t*N.y,Pc0.P.z+t*N.z);
                t+=dt;
            }
            glPointSize(2.0f);
            glVertex3f(Pc0.P.x+t*N.x,Pc0.P.y+t*N.y,Pc0.P.z+t*N.z);
            glEnd();
            cout<<t<<endl;
            hes[str]->densityI = val;
        }
        Repere R2 = getRepere(Pc.P);
        if(Pc.val==0)
        {
            R2 = getRepere0(R2);
        }
        if(hes[str]->densityO == 0)
        {
            float t = 0.01;
            float dt = -0.01;
            float val = 0;
            float dir = -1;
            Point N = {0,0,0};
            N.x = R2.T.x * cos(o) - R2.T.y * sin(o);
            N.y = R2.T.x * sin(o) + R2.T.y * cos(o);
            N.z = R2.T.z;
            if(dist(Pc.P.x+t*N.x,Pc.P.y+t*N.y,Pc.P.z+t*N.z)<dist(Pc.P.x,Pc.P.y,Pc.P.z))
            {
                dt*=-1;
            }
            glBegin(GL_POINTS);
            glPointSize(1.0f);
            float c = 0;
            while(val==0 && abs(t)<1)
            {
                c+=0.1;
                glColor3f(c,c,c);
                val=func(Pc.P.x+t*N.x,Pc.P.y+t*N.y,Pc.P.z+t*N.z);
                glVertex3f(Pc.P.x+t*N.x,Pc.P.y+t*N.y,Pc.P.z+t*N.z);
                t+=dt;
            }
            glPointSize(2.0f);
            glVertex3f(Pc.P.x+t*N.x,Pc.P.y+t*N.y,Pc.P.z+t*N.z);
            glEnd();
            cout<<t<<endl;
            hes[str]->densityO = val;
        }
        glColor3f(0,1,0);
        //glColor3f((float)i/(float)listeFace.size(),1,0);
        int dir = 1;
        int dir2 = 1;
        float dp = dist2(Pc0.P,Pc.P);
        if(Pc0.val == Pc.val && (Pc0.isoC==Pc.isoC || Pc.val==0))
        {
            float d1 = dist2(addP(Pc0.P,R.T,dp*0.5),Pc.P);
            float d2 = dist2(addP(Pc0.P,R.T,-dp*0.5),Pc.P);
            if(d2<d1)
                dir=-1;
            d1 = dist2(addP(Pc.P,R2.T,dp*0.5),Pc0.P);
            d2 = dist2(addP(Pc.P,R2.T,-dp*0.5),Pc0.P);
            if(d2<d1)
                dir2=-1;
            float mult1 = 0.4;
            float mult2 = 0.4;
            if(Pc0.val==0)
            {
                mult1 = 0.05;
            } 
            if(Pc.val==0)
            {
                mult2 = 0.05;
            }
            
            glColor3f(0,1,0);
            //glColor3f((float)i/(float)listeFace.size(),1,0);
            hes[str]->controleI = addP(Pc0.P,R.T,dir*dp*mult1);
            hes[str]->controleO = addP(Pc.P,R2.T,dir2*dp*mult2);
            float l = courbeBezier(Pc0.P,hes[str]->controleI,hes[str]->controleO,Pc.P,true);
            float avg = (hes[str]->densityO + hes[str]->densityI)/2;
            somme+=avg*l;
            divid+=l;
            
        }
        else
        {
            float d1 = dist2(addP(Pc0.P,R.N,dp*0.5),Pc.P);
            float d2 = dist2(addP(Pc0.P,R.N,-dp*0.5),Pc.P);
            if(d2<d1)
                dir=-1;
            d1 = dist2(addP(Pc.P,R2.N,dp*0.5),Pc0.P);
            d2 = dist2(addP(Pc.P,R2.N,-dp*0.5),Pc0.P);
            if(d2<d1)
                dir2=-1;
            float mult1 = 0.4;
            float mult2 = 0.4;
            if(Pc0.val==0)
            {
                mult1 = 0.1;
            }
            if(Pc.val==0)
            {
                mult2 = 0.1;
            }
            
            glColor3f(0,1,0);
            //glColor3f((float)i/(float)listeFace.size(),1,0);
            hes[str]->controleI = addP(Pc0.P,R.N,dir*dp*mult1);
            hes[str]->controleO = addP(Pc.P,R2.N,dir2*dp*mult2);
            float l = courbeBezier(Pc0.P,hes[str]->controleI,hes[str]->controleO,Pc.P,true);
            float avg = (hes[str]->densityO + hes[str]->densityI)/2;
            somme+=avg*l;
            divid+=l;
            
        }
        hes[str]->densityCI = val(hes[str]->controleI.x,hes[str]->controleI.y,hes[str]->controleI.z);
        hes[str]->densityCO = val(hes[str]->controleO.x,hes[str]->controleO.y,hes[str]->controleO.z);
        tid++;
        //cout<<i<<endl<<"\t"<<hes[str]->id<<endl;
        for(int j=1; j<f.sommets.size(); j++)
        {

            PointCourbe Pc1 = f.sommets[j-1];
            string str1 = to_string((int)(Pc1.val*1000))+to_string(Pc1.isoC)+to_string(Pc1.isoInd)+to_string(Pc1.graC)+to_string(Pc1.indGra);
            PointCourbe Pc2 = f.sommets[j];
            string str2 = to_string((int)(Pc2.val*1000))+to_string(Pc2.isoC)+to_string(Pc2.isoInd)+to_string(Pc2.graC)+to_string(Pc2.indGra);
            int ik2 = tid;
            hes[str1+":"+str2] = new half_edge;
            hes[str1+":"+str2]->id = ik2;
            hes[str1+":"+str2]->origine = Pc1.P;
            hes[str1+":"+str2]->densityO = val(Pc1.P.x,Pc1.P.y,Pc1.P.z);
            hes[str1+":"+str2]->incident = Pc2.P;
            hes[str1+":"+str2]->densityI = val(Pc2.P.x,Pc2.P.y,Pc2.P.z);
            hes[str1+":"+str2]->face = i;
            hes[str1+":"+str2]->next = &nullHe;
            hes[str1+":"+str2]->previous = hes[str];
            hes[str1+":"+str2]->opposite = &nullHe;
            Repere R = getRepere(Pc1.P);
            if(Pc1.val==0)
            {
                R = getRepere0(R);
            }
            if(hes[str1+":"+str2]->densityI == 0)
            {
                float t = 0.01;
                float dt = -0.01;
                float val = 0;
                float dir = -1;
                Point N = {0,0,0};
                N.x = R.T.x * cos(o) - R.T.y * sin(o);
                N.y = R.T.x * sin(o) + R.T.y * cos(o);
                N.z = R.T.z;
                if(dist(Pc1.P.x+t*N.x,Pc1.P.y+t*N.y,Pc1.P.z+t*N.z)<dist(Pc1.P.x,Pc1.P.y,Pc1.P.z))
                {
                    dt*=-1;
                }
                float c = 0,
                glBegin(GL_POINTS);
                glPointSize(1.0f);
                while(val==0 && abs(t)<1)
                {
                    c+=0.1;
                    glColor3f(c,c,c);
                    val=func(Pc1.P.x+t*N.x,Pc1.P.y+t*N.y,Pc1.P.z+t*N.z);
                    glVertex3f(Pc1.P.x+t*N.x,Pc1.P.y+t*N.y,Pc1.P.z+t*N.z);
                    t+=dt;
                }
                glPointSize(2.0f);
                glVertex3f(Pc1.P.x+t*N.x,Pc1.P.y+t*N.y,Pc1.P.z+t*N.z);
                glEnd();
                cout<<t<<endl;
                hes[str1+":"+str2]->densityI = val;
            }
            Repere R2 = getRepere(Pc2.P);
            if(Pc2.val==0)
            {
                R2 = getRepere0(R2);
            }
            if(hes[str1+":"+str2]->densityO == 0)
            {
                float t = 0.01;
                float dt = -0.01;
                float val = 0;
                float dir = -1;
                Point N = {0,0,0};
                N.x = R2.T.x * cos(o) - R2.T.y * sin(o);
                N.y = R2.T.x * sin(o) + R2.T.y * cos(o);
                N.z = R2.T.z;
                if(dist(Pc2.P.x+t*N.x,Pc2.P.y+t*N.y,Pc2.P.z+t*N.z)<dist(Pc2.P.x,Pc2.P.y,Pc2.P.z))
                {
                    dt *=-1;
                }
                float c = 0;
                glBegin(GL_POINTS);
                glPointSize(1.0f);
                while(val==0 && abs(t)<1)
                {
                    c+=0.1;
                    glColor3f(c,c,c);
                    val=func(Pc2.P.x+t*N.x,Pc2.P.y+t*N.y,Pc2.P.z+t*N.z);
                    glVertex3f(Pc2.P.x+t*N.x,Pc2.P.y+t*N.y,Pc2.P.z+t*N.z);
                    t+=dt;
                }
                glPointSize(2.0f);
                glVertex3f(Pc2.P.x+t*N.x,Pc2.P.y+t*N.y,Pc2.P.z+t*N.z);
                glEnd();
                cout<<t<<endl;
                hes[str1+":"+str2]->densityO = val;
            }
            glColor3f(0,1,0);
            //glColor3f((float)i/(float)listeFace.size(),1,0);
            int dir = 1;
            int dir2 = 1;
            float dp = dist2(Pc1.P,Pc2.P);
            if(Pc1.val == Pc2.val && (Pc2.isoC==Pc1.isoC || Pc1.val==0))
            {
                float d1 = dist2(addP(Pc1.P,R.T,dp*0.5),Pc2.P);
                float d2 = dist2(addP(Pc1.P,R.T,-dp*0.5),Pc2.P);
                if(d2<d1)
                    dir=-1;
                d1 = dist2(addP(Pc2.P,R2.T,dp*0.5),Pc1.P);
                d2 = dist2(addP(Pc2.P,R2.T,-dp*0.5),Pc1.P);
                if(d2<d1)
                    dir2=-1;
                float mult1 = 0.4;
                float mult2 = 0.4;
                if(Pc1.val==0)
                {
                    mult1 = 0.05;
                } 
                if(Pc2.val==0)
                {
                    mult2 = 0.05;
                }
                
                glColor3f(0,1,0);
                //glColor3f((float)i/(float)listeFace.size(),1,0);
                hes[str1+":"+str2]->controleO = addP(Pc1.P,R.T,dir*dp*mult1);
                hes[str1+":"+str2]->controleI = addP(Pc2.P,R2.T,dir2*dp*mult2);
                float l = courbeBezier(Pc1.P,hes[str1+":"+str2]->controleI,hes[str1+":"+str2]->controleO,Pc2.P,true);
                float avg = (hes[str1+":"+str2]->densityO + hes[str1+":"+str2]->densityI)/2;
                somme+=avg*l;
                divid+=l;
                
            }
            else
            {
                float d1 = dist2(addP(Pc1.P,R.N,dp*0.5),Pc2.P);
                float d2 = dist2(addP(Pc1.P,R.N,-dp*0.5),Pc2.P);
                if(d2<d1)
                    dir=-1;
                d1 = dist2(addP(Pc2.P,R2.N,dp*0.5),Pc1.P);
                d2 = dist2(addP(Pc2.P,R2.N,-dp*0.5),Pc1.P);
                if(d2<d1)
                    dir2=-1;
                float mult1 = 0.4;
                float mult2 = 0.4;
                if(Pc1.val==0)
                {
                    mult1 = 0.1;
                } 
                if(Pc2.val==0)
                {
                    mult2 = 0.1;
                }
                
                glColor3f(0,1,0);
                //glColor3f((float)i/(float)listeFace.size(),1,0);
                hes[str1+":"+str2]->controleO = addP(Pc1.P,R.N,dir*dp*mult1);
                hes[str1+":"+str2]->controleI = addP(Pc2.P,R2.N,dir2*dp*mult2);
                float l = courbeBezier(Pc1.P,hes[str1+":"+str2]->controleI,hes[str1+":"+str2]->controleO,Pc2.P,true);
                float avg = (hes[str1+":"+str2]->densityO + hes[str1+":"+str2]->densityI)/2;
                somme+=avg*l;
                divid+=l;
                
            }
            hes[str1+":"+str2]->densityCI = val(hes[str1+":"+str2]->controleI.x,hes[str1+":"+str2]->controleI.y,hes[str1+":"+str2]->controleI.z);
            hes[str1+":"+str2]->densityCO = val(hes[str1+":"+str2]->controleO.x,hes[str1+":"+str2]->controleO.y,hes[str1+":"+str2]->controleO.z);
            if(hes.find(str2+":"+str1)!= hes.end())
            {
                half_edge* opp = hes[str2+":"+str1];
                opp->opposite = hes[str1+":"+str2];
                hes[str1+":"+str2]->opposite = opp;
            }
            
            hes[str]->next = hes[str1+":"+str2];
            //cout<<"\t"<<hes[str1+":"+str2]->id<<"\t"<<hes[str1+":"+str2]->previous->id<<"\t"<<hes[str1+":"+str2]->opposite->id<<"\t"<<hes[str1+":"+str2]->next->id<<endl;
            //cout<<"\t"<<hes[str]->id<<"\t"<<hes[str]->previous->id<<"\t"<<hes[str]->opposite->id<<"\t"<<hes[str]->next->id<<endl;

            tid++;
            str = str1+":"+str2;
            
        }
        hes[str]->next = hef.incidente;
        hef.incidente->previous = hes[str];
        //cout<<hef.incidente->id<<"\t"<<hef.incidente->previous->id<<"\t"<<hef.incidente->opposite->id<<"\t"<<hef.incidente->next->id<<endl;
        
        if(hes.find(str1+":"+str0)!= hes.end())
        {
            half_edge* opp = hes[str1+":"+str0];
            opp->opposite = hef.incidente;
            hef.incidente->opposite = opp;
        }
        if(divid!=0)
        {
            hef.density = somme/divid;
        }
        Faces.push_back(hef);
    }
    cout<<"Plop1"<<endl;
    
        
    
    glPointSize(1.0f);
    /*if(nbF!=-1)
    {
        /*heFace f = Faces[nbF];
        half_edge *he = f.incidente;
        int id = he->id;
        do{
            Point P0 = he->origine;
            Point P1 = he->controleO;
            Point P2 = he->controleI;
            Point P3 = he->incident;
            glBegin(GL_LINES);
            glColor3f(1,1,1);
            glVertex3f(P0.x,P0.y,P0.z);
            glVertex3f(P1.x,P1.y,P1.z);
            glColor3f(0,0,1);
            glVertex3f(P2.x,P2.y,P2.z);
            glVertex3f(P3.x,P3.y,P3.z);
            glEnd();
            he = he->next;
        }while(he->id!=id);*/
        /*
        heFace f = Faces[nbF];
        float t = f.density;
        f.density = densite2(f);
        cout<<t<<" > "<<f.density<<endl;
        /*half_edge *he = f.incidente;
        int id = he->id;
        glColor3f(1,1,1);
        do{
            courbeBezier(he->origine,he->controleO,he->controleI,he->incident,true);
            he = he->next;
        }while(he->id!=id);*/
    //}
    /*
    if(nbF!=-1)
    {
        for(int i=0; i<Faces.size(); i++)
        {
            heFace f = Faces[i];
            float t = f.density;
            float c = (float)i/(float)Faces.size();
            if(i%2==0)
                glColor3f(0,c,0);
            else
                glColor3f(0,0,c);
            float g = densite2(f);
            if(g!=-1)
                f.density = g;
            cout<<"Face "<<i<<" : "<<t<<" > "<<f.density<<endl;
        }
    }
    
    
    cout<<"Plop2"<<endl;
    /*
    for(int i=0; i<Faces.size(); i++)
    {
        heFace f = Faces[i];
        float val = f.density;
        int h = (int)(val*360);
        Point c = hsv2rgb(h,1.,1.);
        glColor3f(c.x,c.y,c.z);
        half_edge* he = f.incidente;
        int id_ori = he->id;
        glBegin(GL_POLYGON);
        do
        {
            //glVertex3f(he->origine.x,he->origine.y,he->origine.z);
            glVertex3f(he->incident.x,he->incident.y,he->incident.z);
            he = he->next;
        }while(he->id != id_ori);
        glEnd();
    }
    glColor3f(1,0,1);
    for(int i=0; i<Faces.size(); i++)
    {
        heFace f = Faces[i];
        half_edge* he = f.incidente;
        int id_ori = he->id;
        cout<<i<<endl<<"\t"<<he->id<<endl;
        do
        {
            //glBegin(GL_LINES);
            //glVertex3f(he->origine.x,he->origine.y,he->origine.z);
            //glVertex3f(he->incident.x,he->incident.y,he->incident.z);
            //glEnd();
            he = he->next;
            cout<<"\t"<<he->id<<endl;
        }while(he->id != id_ori); 
    }
    */
    /*
    ofstream myfile;
    myfile.open ("output.json",ios::out | ios::trunc);
    myfile << "{" <<endl<<"\t\"Faces\":["<<endl;
    for(int i=0;i<Faces.size();i++)
    {
        //cout<<i<<endl;
        heFace f = Faces[i];
        myfile<<"\t\t{"<<endl<<"\t\t\t\"Density\":"<<f.density<<","<<endl<<"\t\t\t\"Edges\":["<<endl;
        half_edge* he = f.incidente;
        int id_ori = he->id;
        do
        {
            //cout<<"Plop"<<endl;
            //cout<<he->id<<endl;
            myfile<<"\t\t\t\t{"<<endl<<"\t\t\t\t\t\"Points\":["<<endl;
            Point P = he->origine;
            myfile<<"\t\t\t\t\t\t{"<<endl<<"\t\t\t\t\t\t\t\"x\":"<<P.x<<","<<endl<<"\t\t\t\t\t\t\t\"y\":"<<P.y<<","<<endl<<"\t\t\t\t\t\t\t\"z\":"<<P.z<<","<<endl<<"\t\t\t\t\t\t\t\"Density\":"<<he->densityO<<endl<<"\t\t\t\t\t\t},"<<endl;
            P = he->controleO;
            myfile<<"\t\t\t\t\t\t{"<<endl<<"\t\t\t\t\t\t\t\"x\":"<<P.x<<","<<endl<<"\t\t\t\t\t\t\t\"y\":"<<P.y<<","<<endl<<"\t\t\t\t\t\t\t\"z\":"<<P.z<<","<<endl<<"\t\t\t\t\t\t\t\"Density\":"<<he->densityCO<<endl<<"\t\t\t\t\t\t},"<<endl;
            P = he->controleI;
            myfile<<"\t\t\t\t\t\t{"<<endl<<"\t\t\t\t\t\t\t\"x\":"<<P.x<<","<<endl<<"\t\t\t\t\t\t\t\"y\":"<<P.y<<","<<endl<<"\t\t\t\t\t\t\t\"z\":"<<P.z<<","<<endl<<"\t\t\t\t\t\t\t\"Density\":"<<he->densityCI<<endl<<"\t\t\t\t\t\t},"<<endl;
            P = he->incident;
            myfile<<"\t\t\t\t\t\t{"<<endl<<"\t\t\t\t\t\t\t\"x\":"<<P.x<<","<<endl<<"\t\t\t\t\t\t\t\"y\":"<<P.y<<","<<endl<<"\t\t\t\t\t\t\t\"z\":"<<P.z<<","<<endl<<"\t\t\t\t\t\t\t\"Density\":"<<he->densityI<<endl<<"\t\t\t\t\t\t}"<<endl;
            myfile<<"\t\t\t\t\t],"<<endl;
            //ADJACENCE
            half_edge* opp = he->opposite;
            myfile<<"\t\t\t\t\t\"Adjacency\":"<<endl;
            if(opp->id!=-1)
            {
                myfile<<"\t\t\t\t\t{"<<endl<<"\t\t\t\t\t\t\"Face\":"<<opp->face<<","<<endl;
                int num = 0;
                heFace f2 = Faces[opp->face];
                half_edge* he2 = f2.incidente;
                int id_ori2 = he2->id;
                while(opp->id != he2->id)
                {
                    num++;
                    he2 = he2->next;
                    if(he2->id == id_ori2)
                    {
                        break;
                    }
                }
                if(he2->id == id_ori2 && num!=0)
                {
                    num = -1;
                }
                myfile<<"\t\t\t\t\t\t\"Edge\":"<<num<<endl;
            }
            else
            {
                myfile<<"\t\t\t\t\t{"<<endl<<"\t\t\t\t\t\t\"Face\":-1,"<<endl;
                myfile<<"\t\t\t\t\t\t\"Edge\":-1"<<endl;
            }
            myfile<<"\t\t\t\t\t}"<<endl;
            myfile<<"\t\t\t\t}";
            he = he->next;
            if(he->id!=id_ori)
                myfile<<","<<endl;
            else
                myfile<<endl;
        }while(he->id != id_ori); 
        myfile<<"\t\t\t]"<<endl;
        myfile<<"\t\t}";
        if(i!=Faces.size()-1)
            myfile<<","<<endl;
        else
            myfile<<endl;
    }
    myfile<<"\t]"<<endl<<"}";
    myfile.close();
    cout<<"Plop3"<<endl;
    /*
    for(int i=0; i<listeFace.size();i++)
    {
        Face f = listeFace[i];
        for(int j=0; j<f.sommets.size()-1;j++)
        {
            PointCourbe Pc = f.sommets[j];
            PointCourbe Pc2 = f.sommets[j+1];
            Repere R = getRepere(Pc.P);
            Repere R2 = getRepere(Pc2.P);
            int dir = 1;
            int dir2 = 1;
            float dp = dist2(Pc.P,Pc2.P);
            if(Pc.val == Pc2.val && Pc.isoC==Pc2.isoC)
            {
                float d1 = dist2(addP(Pc.P,R.T,dp*0.5),Pc2.P);
                float d2 = dist2(addP(Pc.P,R.T,-dp*0.5),Pc2.P);
                if(d2<d1)
                    dir=-1;
                d1 = dist2(addP(Pc2.P,R2.T,dp*0.5),Pc.P);
                d2 = dist2(addP(Pc2.P,R2.T,-dp*0.5),Pc.P);
                if(d2<d1)
                    dir2=-1;
                courbeBezier(Pc.P,addP(Pc.P,R.T,dir*dp*0.5),addP(Pc2.P,R2.T,dir2*dp*0.5),Pc2.P,true);
            }
            else
            {
                float d1 = dist2(addP(Pc.P,R.N,dp*0.5),Pc2.P);
                float d2 = dist2(addP(Pc.P,R.N,-dp*0.5),Pc2.P);
                if(d2<d1)
                    dir=-1;
                d1 = dist2(addP(Pc2.P,R2.N,dp*0.5),Pc.P);
                d2 = dist2(addP(Pc2.P,R2.N,-dp*0.5),Pc.P);
                if(d2<d1)
                    dir2=-1;
                courbeBezier(Pc.P,addP(Pc.P,R.N,dir*dp*0.5),addP(Pc2.P,R2.N,dir2*dp*0.5),Pc2.P,true);
            }
        }
        PointCourbe Pc = f.sommets[0];
        PointCourbe Pc2 = f.sommets[f.sommets.size()-1];
        Repere R = getRepere(Pc.P);
        Repere R2 = getRepere(Pc2.P);
        int dir = 1;
        int dir2 = 1;
        float dp = dist2(Pc.P,Pc2.P);
        if(Pc.val == Pc2.val && Pc.isoC==Pc2.isoC)
        {
            float d1 = dist2(addP(Pc.P,R.T,dp*0.5),Pc2.P);
            float d2 = dist2(addP(Pc.P,R.T,-dp*0.5),Pc2.P);
            if(d2<d1)
                dir=-1;
            d1 = dist2(addP(Pc2.P,R2.T,dp*0.5),Pc.P);
            d2 = dist2(addP(Pc2.P,R2.T,-dp*0.5),Pc.P);
            if(d2<d1)
                dir2=-1;
            courbeBezier(Pc.P,addP(Pc.P,R.T,dir*dp*0.5),addP(Pc2.P,R2.T,dir2*dp*0.5),Pc2.P,true);
        }
        else
        {
            float d1 = dist2(addP(Pc.P,R.N,dp*0.5),Pc2.P);
            float d2 = dist2(addP(Pc.P,R.N,-dp*0.5),Pc2.P);
            if(d2<d1)
                dir=-1;
            d1 = dist2(addP(Pc2.P,R2.N,dp*0.5),Pc.P);
            d2 = dist2(addP(Pc2.P,R2.N,-dp*0.5),Pc.P);
            if(d2<d1)
                dir2=-1;
            courbeBezier(Pc.P,addP(Pc.P,R.N,dir*dp*0.5),addP(Pc2.P,R2.N,dir2*dp*0.5),Pc2.P,true);
        }
    }*/
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
        nbF++;
        if(nbF>=listeFace.size()) nbF = listeFace.size()-1;
      glutPostRedisplay();
      break;
    case '-': //* affichage du carre plein 
        nbF--;
	if (nbF < 1 ) nbF=0;
      glutPostRedisplay();
      break;
  case 'l': //* affichage du carre plein 
        nbE++;
        cout<<nbE<<endl;
      glutPostRedisplay();
      break;
    case 'm': //* affichage du carre plein 
        nbE--;
        cout<<nbE<<endl;
      glutPostRedisplay();
      break;
  case 'x':
    pX+=0.1;
    glutPostRedisplay();
    break;
case 'X':
    pX-=0.1;
    glutPostRedisplay();
    break;
case 'y':
    pY+=0.1;
    glutPostRedisplay();
    break;
case 'Y':
    pY-=0.1;
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

    
    
