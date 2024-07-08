#include <stdlib.h>
#include <cmath>
#include <GL/glut.h>
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <ctime>
#include <cstdlib> 
#include <iostream> 
#include <chrono>

using namespace std;
using namespace chrono;

void affichage(void);

void clavier(unsigned char touche,int x,int y);
void affiche_repere(void);
void initMap();
void mouse(int, int, int, int);
void mouseMotion(int, int);
//void reshape(int,int);

int nbF = 0;

struct Point{
    float x;
    float y;
    float z;
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
    Point controle;
    int face;
    half_edge* next;
    half_edge* previous;
    half_edge* opposite;
};

struct heFace{
    int id;
    half_edge* incidente;
};

half_edge nullHe = {-1,{0,0,0},{0,0,0},-1,nullptr,nullptr,nullptr};

map<string,half_edge> hes;
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

int nb = 0;
//vector<Point> P = {};
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

//-----------------------------------------------------

Point hsv2rgb(int h,float s, float v)
{
    h = h%360;
    float c = s*v;
    float x = c * ( 1 -abs((h/60)%2 - 1));
    float m = v-c;
    Point C = {0,0,0};
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

float func(float x, float y, float z)
{
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
    if(x>=-1 && y>=-1 && x<=1 && y<=1)
        return (v1+v2)/2;
    return 0;
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
    float dt = 0.001;
    float target = func(P.x,P.y,P.z);
    float p1 = target - func(P.x-dt,P.y-dt,P.z);
    float p2 = target - func(P.x+dt,P.y-dt,P.z);
    float p3 = target - func(P.x+dt,P.y+dt,P.z);
    float p4 = target - func(P.x-dt,P.y+dt,P.z);
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
    /*if(abs(p1-target)<0.0001 || abs(p2-target)<0.0001 || abs(p3-target)<0.0001 || abs(p4-target)<0.0001)
    {
        return {P,{1,1,1},T}; 
    }*/
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
            if(p1==target || p2==target || p3==target || p4==target)
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
        glColor3f(r,g,b);
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


void courbeGradient(Point P, int dir, Point PrevN, vector<PointCourbe> *Pv)
{
    int g = 0;
    float val = func(P.x,P.y,P.z);
    if(val==0)
        return;
    Point N = {0,0,0};
    float dt = 0.005;
    if(dir!=0)
    {
        float val = func(P.x,P.y,P.z);
        while(val!=0)
        {
            
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
                        //cout<<(targets[iv]!=Pc.val)<<" "<<dist2(P,Pc.P)<<endl;
                        if(targets[iv]!=Pc.val || dist2(P,Pc.P)>0.02001)
                        {
                            Pv->push_back({P,targets[iv],-1,-1,-1,-1});
                        }
                    }
                    else
                    {
                        PointCourbe Pc = Pv->at(0);
                        //cout<<(targets[iv]!=Pc.val)<<" "<<dist2(P,Pc.P)<<endl;
                        if(targets[iv]!=Pc.val || dist2(P,Pc.P)>0.02001)
                        {
                            Pv->insert(Pv->begin(),{P,targets[iv],-1,-1,-1,-1});
                        }
                    }
                }
                iv++;
            }
            Repere R = getRepere(P);
            N = R.N;
            Point nP3 = {0,0,0};
            if(abs(dist(N.x,N.y,N.z)-1)>0.0001)
            {
                N = PrevN;
                nP3 = {P.x+N.x*dt*dir,P.y+N.y*dt*dir,P.z+N.z*dt*dir};
            }
            else
            {
                Point nP = {P.x+N.x*dt*dir,P.y+N.y*dt*dir,P.z+N.z*dt*dir};
                Repere R2 = getRepere(nP);
                N = R2.N;
                Point nP2 = {P.x+N.x*dt*dir,P.y+N.y*dt*dir,P.z+N.z*dt*dir};
                float val1 = val-func(nP.x,nP.y,nP.z);
                float val2 = val-func(nP2.x,nP2.y,nP2.z);
                float inter = interpolation(val2,val1);
                nP3.x = nP.x * inter + (1-inter) * nP2.x;
                nP3.y = nP.y * inter + (1-inter) * nP2.y;
                nP3.z = nP.z * inter + (1-inter) * nP2.z;
                N = {nP3.x-P.x,nP3.y-P.y,nP3.z-P.z};
                N.x *= dir;
                N.y *= dir;
                N.z *= dir;
                float m = dist(N.x,N.y,N.z);
                N.x/=m;
                N.y/=m;    
                N.z/=m;
            }
            PrevN=N;
            glBegin(GL_LINES);
            glVertex3f(P.x,P.y,P.z);
            glVertex3f(nP3.x,nP3.y,nP3.z);
            glEnd();
            P = nP3;
            val = func(P.x,P.y,P.z);
        }
        if(dir==1)
            Pv->push_back({P,0,-1,-1,-1,-1});
        else
           Pv->insert(Pv->begin(),{P,0,-1,-1,-1,-1}); 
    }
    else
    {
        glColor3f(0,1,0);
        glBegin(GL_POINTS);
        glVertex3f(P.x,P.y,P.z);
        glEnd();
        Repere R = getRepere(P);
        N = R.N;
        Point nP1 = {P.x+N.x*dt,P.y+N.y*dt,P.z+N.z*dt};
        Point nP2 = {P.x-N.x*dt,P.y-N.y*dt,P.z-N.z*dt};
        glBegin(GL_LINES);
        glColor3f(1,0,1);
        glVertex3f(P.x,P.y,P.z);
        glVertex3f(nP1.x,nP1.y,nP1.z);
        glVertex3f(P.x,P.y,P.z);
        glVertex3f(nP2.x,nP2.y,nP2.z);
        glEnd();
        courbeGradient(nP1,1,N,Pv);
        courbeGradient(nP2,-1,N,Pv);
    }    
}


//------------------------------------------------------
void affichage(void)
{
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
	glTranslatef(0,0,cameraDistance);
glRotatef(cameraAngleX,1.,0.,0.)	;
glRotatef(cameraAngleY,0.,1.,0.);
	affiche_repere();
    float n = 100;

    for(int i=-n; i<n; i++)
    {
        for(int j=-n; j<n; j++)

        {
            glBegin(GL_QUADS);
            glColor3f(abs(func((float)(i)/n,(float)(j)/n,0)),0,0);
            glVertex3f((float)(i)/n,(float)(j)/n,0);

            glColor3f(abs(func((float)(i+1)/n,(float)(j)/n,0)),0,0);
            glVertex3f((float)(i+1)/n,(float)(j)/n,0);

            glColor3f(abs(func((float)(i+1)/n,(float)(j+1)/n,0)),0,0);
            glVertex3f((float)(i+1)/n,(float)(j+1)/n,0);

            glColor3f(abs(func((float)(i)/n,(float)(j+1)/n,0)),0,0);
            glVertex3f((float)(i)/n,(float)(j+1)/n,0);
            glEnd();
        }
    }
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
    for(auto it = targets.begin(); it!=targets.end(); it++)
    {
        float v = *it;
        Point C = {v*v,v*v,v*v};
        lancerMarchingSquare(v,C);
    }

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
        int ps = s/div;
        for(int i=0; i<=div ; i++)
        {
            Point P = c[ps*i];
            vector<PointCourbe> cG;
            //cout<<nt<<endl;
            cG.push_back({P,nt,-1,-1,-1});
            courbeGradient(P,0,{0,0,0},&cG);
            courbeGrad.push_back(cG);
            //cout<<endl<<endl;
        }
    }
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
                    Pc.isoInd = courbeGrad.size()+i;
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
                    if(Pc2.graC != Pc.graC || dist2(Pc.P,Pc2.P)>0.02)
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
    }
    for(int i=0; i<courbeGrad.size();i++)
    {
        vector<PointCourbe> cG = courbeGrad[i];
        for(auto it = cG.begin(); it!=cG.end(); it++)
        {
            PointCourbe Pc = *it;
            cout<<Pc.val<<"\t"<<Pc.isoC<<"\t"<<Pc.isoInd<<"\t"<<Pc.graC<<"\t"<<Pc.indGra<<endl;
        }
    }
    /*
    targets.insert(targets.begin(),0);
    for(int i=0; i<targets.size(); i++)
    {
        cout<<endl<<targets[i]<<endl;
        vector<vector<PointCourbe>> cs = courbeFinale[targets[i]*1000];
        for(int j=0; j<cs.size(); j++)
        {
            cout<<"Courbe "<<j<<endl;
            vector<PointCourbe> c = cs[j];
            for(int k=0; k<c.size(); k++)
            {
                PointCourbe Pc = c[k];
                cout<<Pc.val<<"\t"<<Pc.isoC<<"\t"<<Pc.isoInd<<"\t"<<Pc.graC<<"\t"<<Pc.indGra<<endl;
            }
        }
    }
    cout<<endl<<endl;*/
    for(int i=0; i<courbeGrad.size()-1; i++)
    {
        vector<PointCourbe> cG = courbeGrad[i];
        for(int j=0; j<cG.size()-1; j++)
        {
            bool noP = false;
            PointCourbe Pc = cG[j];
            vector<PointCourbe> sf;
            sf.push_back(Pc);
            bool finIso = true;
            bool cont = true;
            float valFin = 0;
            int step = 0;
            //cout<<"Face\n"<<step<<"\t"<<Pc.val<<"\t"<<Pc.isoC<<"\t"<<Pc.isoInd<<"\t"<<Pc.graC<<"\t"<<Pc.indGra<<endl;
            PointCourbe Pc2 = Pc;
            do
            {
                if(step==0 || (step>=2 && cont))
                {
                    vector<PointCourbe> grad = courbeGrad[Pc2.graC];
                    if(step<2)
                    {
                        if(Pc2.indGra<grad.size()-1)
                        {
                            Pc2 = grad[Pc2.indGra+1];
                        }
                        else
                        {
                            Pc2 = grad[0];
                        }
                    }
                    else
                    {
                        if(Pc2.indGra>0)
                        {
                            Pc2 = grad[Pc2.indGra-1];
                        }
                        else
                        {
                            Pc2 = grad[grad.size()-1];
                        }
                    }
                }
                else
                {
                    vector<PointCourbe> iso = courbeFinale[Pc2.val*1000][Pc2.isoC];
                    if(step<2)
                    {
                        if(Pc2.isoInd<iso.size()-1)
                        {
                            Pc2 = iso[Pc2.isoInd+1];
                        }
                        else
                        {
                            noP = true;
                            step-=2;
                        }
                    }
                    else
                    {
                        if(Pc2.isoInd>0)
                        {
                            Pc2 = iso[Pc2.isoInd-1];
                        }
                        else
                        {
                            Pc2 = iso[iso.size()-1];
                        }
                    }
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
                    sf.push_back(Pc2);
                //cout<<step<<"\t"<<Pc2.val<<"\t"<<Pc2.isoC<<"\t"<<Pc2.isoInd<<"\t"<<Pc2.graC<<"\t"<<Pc2.indGra<<endl;
            }while((Pc2.isoC!=Pc.isoC || Pc2.val != Pc.val && finIso));
            if(Pc.isoInd==Pc2.isoInd)
            {
                sf.pop_back();
            }
            Face f;
            f.sommets = sf;
            listeFace.push_back(f);
        }
    }
    cout<<nbF<<" / "<<listeFace.size()<<endl;
    int h = 2000 / listeFace.size();
    int tid=0;
    for(int i=0; i<nbF;i++)
    {
        Face f = listeFace[i];
        PointCourbe Pc0 = f.sommets[0];
        half_edge he = {tid,Pc0.P,{0,0,0},i,&nullHe,&nullHe,&nullHe};
        heFace hef = {i,&he};
        tid++;
        for(int j=1; j<f.sommets.size(); j++)
        {
            PointCourbe Pc2 = f.sommets[j];
            half_edge he2 = {tid,Pc2.P,{0,0,0},i,&nullHe,&he,&nullHe};
            he.next =&he2;
            tid++;
            he = he2;
        }
        he.next = hef.incidente;
        hef.incidente->previous = &he;
        Faces.push_back(hef);
        /*Point c = hsv2rgb(h*i,1,1);
        glBegin(GL_POLYGON);
        glColor3f(c.x,c.y,c.z);
        for(int j=0; j<f.sommets.size(); j++)
        {
            PointCourbe Pc = f.sommets[j];
            if(i==nbF-1)
            {
                cout<<Pc.val<<"\t"<<Pc.isoC<<"\t"<<Pc.isoInd<<"\t"<<Pc.graC<<"\t"<<Pc.indGra<<endl;
            }
            glVertex3f(Pc.P.x,Pc.P.y,Pc.P.z);
        }
        glEnd();
        if(i==nbF-1)
            cout<<endl;
        */
    }
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
        if(nbF>listeFace.size()) nbF = listeFace.size()-1;
      glutPostRedisplay();
      break;
    case '-': //* affichage du carre plein 
        nbF--;
	if (nbF < 1 ) nbF=0;
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

    
    
