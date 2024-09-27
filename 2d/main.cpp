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
void initMap();
void mouse(int, int, int, int);
void mouseMotion(int, int);


double normeMin = .001;
double normeMax = .005;

float color=0;

int nbF = 0;
int nbE = 0;

struct Point{
    double x;
    double y;
    double z;
};

struct PointCourbe{
    Point P = {0,0,0};
    float val = -1;
    int isoC = -1;
    int isoInd = -1;
    int graC = -1;
    int indGra = -1;
};


struct Repere{
    Point P;
    Point N;
    Point T;
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
double distProj(Point P0, Point P1, Point Pt);

map<string,bool> marque;
vector<Point> pile;

half_edge nullHe = {-1,{0,0,0},0,{0,0,0},0,{0,0,0},0,{0,0,0},0,-1,nullptr,nullptr,nullptr};

map<string,half_edge*> hes;
vector<heFace> Faces;

vector<float> targets;
map<int,vector<vector<PointCourbe*>>> newCourbeFinale; //Iso-lignes
map<string,int> msTable;
vector<vector<PointCourbe*>> NewCourbeGrad; //Courbes de gradient


float msMoy = 0;
int nbMs = 0;
float cbMoy = 0;
int nbCb = 0;

float moy = 0;
int ng=0;
float moy2=0;
int n2=0;



bool dessinne = true;

int nb = 0;

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

double pX = 0;
double pY = 0;
float dp = 0.0001;

//------------------------------------------------------
int main(int argc,char **argv)
{
    initMap();
    srand(time(0));
  /* initialisation de glut et creation
     de la fenetre */
  glutInit(&argc,argv);
  glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
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
	 gluPerspective(60.0,(GLfloat)200.0/(GLfloat)200.0,0.0001,100.0);
	glMatrixMode(GL_MODELVIEW);
 gluLookAt(0.,0.,2, 0.,0.,0., 0.,1.,0.);

/* Entree dans la boucle principale glut */
  glutMainLoop();
  return 0;
}


//-----------------------------------------------------
//Calcul la valeur RGB d'une couleur HSV passée en paramètre
/**
 * s (saturation) et v (valeur) sont compris entre 0 et 1
 * h (teinte) prends des valeurs tel que 0 correspond à 0° et 1 correspond à 360°
 * Si la valeur de h est inférieur à 0 ou supérieur à 1 alors elle sera modifié pour y être (1.0001 -> 0.0001 ; 2 -> 1...)
*/
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

//Calcul le produit scalaire entre deux vecteurs
float produitScalaire(Point V1, Point V2)
{
    return V1.x*V2.x+V1.y*V2.y+V1.z*V2.z;
}

//Donne la pondération d'un point pour faire une moyenne pondérée entre deux points pour trouver le point d'équilibre entre deux points dont les valeurs sont entrées en paramètre
/*
Exemple : Un point de valeur -2 et un point de valeur 3, on cherche le point de valeur 0 entre les deux. Il est plus près de celui de valeur -2.
Le point de valeur 0 est égal à 3/5 de celui de valeur -2 et 2/5 de celui de valeur 3. Cette fonction renverrait 3/5.
*/
float interpolation(float a, float b)
{
    float aa = abs(a);
    float ab = abs(b);
    if(aa+ab!=0)
        return ab/(aa+ab);
    return 0.5;
}

//Calcul le point interpolé linéairement entre deux points
/**
 * P0 est le premier point (pour t=0)
 * P1 est le deuxième point (pour t=1)
 * t est la valeur qui donne l'interpolation
*/
Point interPoint(Point P0, Point P1, double t)
{
    Point P = {0,0,0};
    P.x = t * P0.x + (1-t)*P1.x;
    P.y = t * P0.y + (1-t)*P1.y;
    P.z = t * P0.z + (1-t)*P1.z;
    return P;
}

//Donne la distance entre 2 points
float dist2(Point P1, Point P2)
{
    float x = P2.x - P1.x;
    float y = P2.y - P1.y;
    float z = P2.z - P1.z;
    return sqrt(x*x+y*y+z*z);
}

//Renvoie la distance d'un point par rapport au point {0,0,0}
//Donne aussi la norme d'un vecteur
float dist(float x, float y, float z)
{
    return sqrt(x*x+y*y+z*z);
}

//Renvoie la valeur de la fonction pour un point donné
double val(double x, double y, double z)
{  
    double v = 0.5*cos(x+2*y)+0.5*cos(3*x+2*y);
    v/=2;
    v+=0.5;
    return v;
}

//Applique un masque sur la valeur de la fonction pour un point donné
//Le masque transforme la valeur de la fonction en -1 si le point est en dehors de la surface 
double func(float x, float y, float z)
{
    if(x>=-1 && y>=-1 && x<=1 && y<=1)
        return val(x,y,z);
    else
        return -1;
}

//Remplit la map qui associe à chaque configuration de Marching Square un numéro représentant le traitement du carré
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

//Calcul un repère orthogonal local sur un point donné
//Ce repère est donné par le gradient de la fonction au point donné ainsi que la tangente à l'iso-ligne qui passerait par ce point (le gradient et la tangente à l'iso-ligne sont orthogonaux pour un point donnée)
Repere getRepere(Point P)
{
    double dtO = 0.00000001;
    double dt = dtO;
    double target = val(P.x,P.y,P.z);
    double p1 = val(P.x+dt,P.y,P.z);
    double p2 = val(P.x,P.y+dt,P.z);
    double p3 = val(P.x-dt,P.y,P.z);
    double p4 = val(P.x,P.y-dt,P.z);
    bool p = true;
    //On cherche une distance entre le point central et les points qui servent à calculer le gradient tel que les valeurs sont différentes afin d'avoir un gradient qui ne soit pas nul dans la plupart des cas
    //Même si le gradient local au point est nul, on utilise la tendance proche
    while(abs(p1-p3)<dtO/10 && abs(p2-p4)<dtO/10)
    {
        if(p)
            dt/=10;
        else
            dt*=1.1;
        p1 = val(P.x+dt,P.y,P.z);
        p2 = val(P.x,P.y+dt,P.z);
        p3 = val(P.x-dt,P.y,P.z);
        p4 = val(P.x,P.y-dt,P.z);
        if(dt<dtO/1000)
        {
            p=false;
            dt = dtO;
        }
        if(dt>normeMin && (abs(target-p1)>dtO/10 || abs(target-p3)>dtO/10 || abs(target-p2)>dtO/10 || abs(target-p4)>dtO/10))
        {
            break;
        }
    }
    double dFx = 0;
    double dFy = 0;
    double dX = 1;
    double dY = 1;
    //De plus, si le point dont on cherche le gradient est un extremum local dont les valeurs autour changent à la même fréquence, alors les points avant et après (en x ou y) peuvent avoir les mêmes valeurs
    //Dans ce cas, on ne prend qu'une des deux valeurs de cette direction
    if(abs(p1-p3)>dtO/10)
    {
        dFx = p1-p3;
        dX = 2*dt;
    }
    else
    {
        if(abs(target-p1)>dtO/10)
        {
            dFx = p1 - target;
            dX = dt;
        }
        else
        {
            dFx = target - p3;
            dX = dt;
        }
    }
    if(abs(p2-p4)>dtO/10)
    {
        dFy = p2-p4;
        dY = 2*dt;
    }
    else
    {
        if(abs(target-p2)>dtO/10)
        {
            dFy = p2 - target;
            dY = dt;
        }
        else
        {
            dFy = target - p4;
            dY = dt;
        }
    }
    float o = -M_PI/2;
    Point N = {0,0,0};
    Point T = {0,0,0};
    N.x = dFx/dX;
    N.y = dFy/dY;
    //La tangente à l'iso-ligne qui passerait par ce point est orthogonale au gradient
    T.x = N.x * cos(o) - N.y * sin(o);
    T.y = N.x * sin(o) + N.y * cos(o);
    T.z = N.z;
    //On normalise la tangente puisqu'on a pas d'information sur la courbure de l'iso-ligne
    float m2 = dist(T.x,T.y,T.z);
    T.x/=m2;
    T.y/=m2;
    T.z/=m2;
    return {P,N,T};
}

//Lance Marching Square pour obtenir une liste des points des courbes de gradient étant sur cette iso-ligne
//Cette fonction ne récupère donc que les points d'intersections entre iso-lignes et courbes de gradient
/**
 * target est la valeur de l'iso-ligne approximée
 * pas est la moitié des bords des carrés de Marching Square (c'est la distance minimal entre le point central et les bords du carré)
 * iso est la liste des points des courbes de gradient faisant partie de cette iso-ligne
 * P est le point de départ, c'est à dire le point central du premier carré
 * nb indique le numéro de l'iso-ligne pour une valeur
 * dir indique la direction de déplacement le long de l'iso-ligne
 * prev indique d'où vient le carré précédent afin de ne pas faire demi-tour (ou de ne pas dédoubler continuellement le nombre de carré)
*/
void findIso(double target, double pas, vector<PointCourbe*> *iso, Point P, int nb, int dir, int prev)
{
    double r=0.25, g=0.25, b=0;
    double t = func(P.x,P.y,P.z);
    double dmin = 0.0000000001;
    vector<PointCourbe*> add;
    bool fin = false;
    bool adding = false;
    //Cherche les points des courbes de gradient qui sont dans le carré
    //Si un de ces points a déjà été ajouté dans une iso-ligne (normalement celle qui est entrain d'être créer), alors on s'arrête (parce que ça veux dire que l'algorithme a fait une boucle)
    for(int i = 0; i<NewCourbeGrad.size(); i++)
    {
        vector<PointCourbe*> c = NewCourbeGrad[i];
        for(int j = 0; j<c.size(); j++)
        {
            Point Pc = c[j]->P;
            double val = c[j]->val;
            int isoC = c[j]->isoC;
            double x = abs(P.x-Pc.x);
            double y = abs(P.y-Pc.y);
            if(x<=pas && y<=pas && val==target && !fin)
            {
                if(isoC!=-1)
                    fin = true;
                else
                    add.push_back(c[j]);
            }
        }
    }
    if(add.size()!=0)
    {
        adding =true;
        r+= 0.25;
        g+= 0.25;
    }
    if(fin)
        r += 0.25;
    
    if(!fin)
    {
        double l0 = func(P.x-pas,P.y-pas,P.z);
        double l1 = func(P.x+pas,P.y-pas,P.z);
        double l2 = func(P.x+pas,P.y+pas,P.z);
        double l3 = func(P.x-pas,P.y+pas,P.z);
        double v0 = target - val(P.x-pas,P.y-pas,P.z);
        double v1 = target - val(P.x+pas,P.y-pas,P.z);
        double v2 = target - val(P.x+pas,P.y+pas,P.z);
        double v3 = target - val(P.x-pas,P.y+pas,P.z);
        //On traite les cas où un ou plusieurs coins du carré se trouve sur l'iso-ligne
        if(v0==0)
        {
            if((v1>0 && v3<0) || (v1<0 && v3>0))
            {
                if(abs(v1)<abs(v3))
                    v0 = copysign(dmin,v1);
                else
                    v0 = copysign(dmin,v3);
            }
            else
            {
                if(abs(v1)>abs(v3))
                    v0 = copysign(dmin,-v1);
                else
                    v0 = copysign(dmin,-v3);
            }
        }
        if(v1==0)
        {
            if((v0>0 && v2<0) || (v0<0 && v2>0))
            {
                if(abs(v0)<abs(v2))
                    v1 = copysign(dmin,v0);
                else
                    v1 = copysign(dmin,v2);
            }
            else
            {
                v1 = copysign(dmin,-v0);
            }
        }
        if(v2==0)
        {
            if((v1>0 && v3<0) || (v1<0 && v3>0))
            {
                if(abs(v1)<abs(v3))
                    v2 = copysign(dmin,v1);
                else
                    v2 = copysign(dmin,v3);
            }
            else
            {
                v2 = copysign(dmin,-v1);
            }
        }
        if(v3==0)
        {
            if((v0>0 && v2<0) || (v0<0 && v2>0))
            {
                if(abs(v0)<abs(v2))
                    v3 = copysign(dmin,v0);
                else
                    v3 = copysign(dmin,v2);
            }
            else
            {
                v3 = copysign(dmin,-v0);
            }
        }
        Point P0 = {P.x-pas,P.y-pas,P.z};
        Point P1 = {P.x+pas,P.y-pas,P.z};
        Point P2 = {P.x+pas,P.y+pas,P.z};
        Point P3 = {P.x-pas,P.y+pas,P.z};
        Point P01 = {P.x,P.y-2*pas,P.z};
        Point P12= {P.x+2*pas,P.y,P.z};
        Point P23 = {P.x,P.y+2*pas,P.z};
        Point P30 = {P.x-2*pas,P.y,P.z};
        //On identifie la configuration des coins du carré
        string s = to_string(v0>0)+to_string(v1>0)+to_string(v2>0)+to_string(v3>0);
        string s2 = to_string(l0==-1)+to_string(l1==-1)+to_string(l2==-1)+to_string(l3==-1);
        int chx = msTable[s];
        Point Pt0 = P;
        Point Pt1 = P;
        Point Pt2 = P;
        Point Pt3 = P;
        Point Pf1 = P;
        Point Pf2 = P;
        switch(chx)
        {
            case 1:

                if(prev==1)
                {
                    Pt0 = P2;
                    Pt1 = P3;
                    Pt2 = P0;
                    Pt3 = P3;
                }
                else
                {
                    Pt0 = P0;
                    Pt1 = P3;
                    Pt2 = P2;
                    Pt3 = P3;
                }
                
                break;
            case 2:
                if(prev==4)
                {
                    Pt0 = P1;
                    Pt1 = P2;
                    Pt2 = P3;
                    Pt3 = P2;
                }
                else
                {
                    Pt0 = P3;
                    Pt1 = P2;
                    Pt2 = P1;
                    Pt3 = P2;
                }
                
                break;
            case 3:
                if(prev==4)
                {
                    Pt0 = P1;
                    Pt1 = P2;
                    Pt2 = P0;
                    Pt3 = P3;
                }
                else
                {
                    Pt0 = P0;
                    Pt1 = P3;
                    Pt2 = P1;
                    Pt3 = P2;
                }
                break;
            case 4:
                if(prev==4)
                {
                    Pt0 = P1;
                    Pt1 = P2;
                    Pt2 = P0;
                    Pt3 = P1;
                }
                else
                {
                    Pt0 = P0;
                    Pt1 = P1;
                    Pt2 = P1;
                    Pt3 = P2;
                }
                break;
            case 6:
                if(prev==3)
                {
                    Pt0 = P1;
                    Pt1 = P0;
                    Pt2 = P2;
                    Pt3 = P3;
                }
                else
                {
                    Pt0 = P2;
                    Pt1 = P3;
                    Pt2 = P1;
                    Pt3 = P0;
                }
                break;
            case 7:
                if(prev==3)
                {
                    Pt0 = P1;
                    Pt1 = P0;
                    Pt2 = P0;
                    Pt3 = P3;
                }
                else
                {
                    Pt0 = P0;
                    Pt1 = P3;
                    Pt2 = P1;
                    Pt3 = P0;
                }
                break;
        }
        //On identifie les points de l'iso-ligne sur les bords du carré (les bords concernés dépendent de la configuration)
        for(int i=0; i<1000; i++)
        {
            double t = i*0.001;
            Pf1 = interPoint(Pt0,Pt1,t);
            if(abs(val(Pf1.x,Pf1.y,Pf1.z)-target)<0.00001)
                break;
        }
        for(int i=0; i<1000; i++)
        {
            double t = i*0.001;
            Pf2 = interPoint(Pt2,Pt3,t);
            if(abs(val(Pf2.x,Pf2.y,Pf2.z)-target)<0.00001)
                break;
        }
        //On ordonne les points des courbes de gradient selon la distance de leur projection entre les deux points de l'iso-ligne de l'étape précédente
        vector<PointCourbe*> adding;
        for(int i=0; i<add.size(); i++)
        {
            PointCourbe* Pct = add[i];
            Point Pt = Pct->P;
            double distP = distProj(Pf1,Pf2,Pt);
            bool added = false;
            for(auto it = adding.begin(); it!=adding.end() && !added; it++)
            {
                Point Pa = (*it)->P;
                double distP2 = distProj(Pf1,Pf2,Pa);
                if(distP<distP2)
                {
                    added = true;
                    adding.insert(it,Pct);
                }
            }
            if(!added)
                adding.push_back(Pct);
        }
        //On ajoute les points des courbes de gradient dans l'iso-ligne (dans l'ordre estimé à l'étape précédente)
        for(int i=0; i<adding.size(); i++)
        {
            PointCourbe* Pca = adding[i];
            Pca->isoC = nb;
            if(dir == -1)
            {
                iso->insert(iso->begin(),Pca);
            }
            else
            {
                iso->push_back(Pca);
            }
        }
        //On vérifie si le carré touche le bord de la surface
        int chx2 = msTable[s2];
        if(chx2!=0)
        {
            g += 0.5;
            //On vérifie si l'iso-ligne sort de la surface
            if(func(Pf2.x,Pf2.y,Pf2.z)==-1 || func(Pf1.x,Pf1.y,Pf1.z)==-1)
            {
                //On identifie le point de l'iso-ligne le plus proche du bord (ou proche en dehors de la surface)
                fin = true;
                Point Pf = Pf2;
                Point Po = Pf1;
                if(func(Pf2.x,Pf2.y,Pf2.z)!=-1)
                {
                    Po = Pf2;
                    Pf = Pf1;
                }
                Point V = addP(Pf,Po,-1);
                int st=0;
                int lim = (pas*4)/0.001;
                Repere R = getRepere(Po);
                double ao = abs(getAngle(V,R.T));
                if(ao>90)
                    ao=abs(ao-180);
                ao/=180;
                while(func(Po.x,Po.y,Po.z)!=-1 && dist(V.x,V.y,V.z)>0.000001 && dist2(Po,P)<pas*1.5)
                {
                    double aMin = -110;
                    double ang = (0.01*aMin)*M_PI*ao;
                    Point Vmi = V;
                    Vmi.x = V.x * cos(ang) - V.y * sin(ang);
                    Vmi.y = V.x * sin(ang) + V.y * cos(ang);
                    Vmi.z = V.z;
                    Point Pmi = addP(Po,Vmi,min(0.001,0.9*dist(V.x,V.y,V.z)));
                    double difMin = abs(val(Pmi.x,Pmi.y,Pmi.z)-target);
                    for(int i=-110; i<=110; i++)
                    {
                        double a = (i*0.01)*M_PI*ao;
                        Point Vm = V;
                        Vm.x = V.x * cos(a) - V.y * sin(a);
                        Vm.y = V.x * sin(a) + V.y * cos(a);
                        Vm.z = V.z;
                        Point Pm = addP(Po,Vm,min(0.001,0.9*dist(V.x,V.y,V.z)));
                        double diff = abs(val(Pm.x,Pm.y,Pm.z)-target);
                        if(diff<0.00001)
                        {
                            aMin = i;
                            difMin = diff;
                            Pmi = Pm;
                            break;
                        }
                        else
                        {
                            if(diff<difMin)
                            {
                                aMin = i;
                                difMin = diff;
                                Pmi = Pm;
                            }
                        }
                    }
                    Po = Pmi;
                    V = addP(Pf,Po,-1);
                    st++;
                    
                }
                if(dir!=-1)
                {
                    iso->push_back(new PointCourbe);
                    iso->back()->P = Po;
                    iso->back()->val = target;
                    iso->back()->isoC = nb;
                }
                else
                {
                   iso->insert(iso->begin(),new PointCourbe);
                   iso->front()->P = Po;
                   iso->front()->val = target;
                   iso->front()->isoC = nb;
                }
            }
            
        }
        if(chx2==0 && l0==-1)
            fin = true;
        //On regarde si l'algorithme devrait se terminer
        if(!fin)
        {
            //On continu l'algorithme
            if(dir==0)
            {
                //Si c'est le premier carré, alors on lance le carré dans les deux directions de l'iso-ligne
                switch(chx)
                {
                    case 1:
                        findIso(target,pas,iso,P23,nb,-1,3);
                        findIso(target,pas,iso,P30,nb,1,4);
                        break;
                    case 2:
                        findIso(target,pas,iso,P12,nb,-1,2);
                        findIso(target,pas,iso,P23,nb,1,3);
                        break;
                    case 3:
                        findIso(target,pas,iso,P12,nb,-1,2);
                        findIso(target,pas,iso,P30,nb,1,4);
                        break;
                    case 4:
                        findIso(target,pas,iso,P01,nb,-1,1);
                        findIso(target,pas,iso,P12,nb,1,2);
                        break;
                    case 6:
                        findIso(target,pas,iso,P01,nb,-1,1);
                        findIso(target,pas,iso,P23,nb,1,3);
                        break;
                    case 7:
                        findIso(target,pas,iso,P01,nb,-1,1);
                        findIso(target,pas,iso,P30,nb,1,4);
                        break;
                }
            }
            else
            {
                //Sinon, on ne continu que dans la direction qui ne fait pas revenir en arrière
                switch(chx)
                {
                    case 1:
                        if(prev!=1)
                            findIso(target,pas,iso,P23,nb,dir,3);
                        else
                            findIso(target,pas,iso,P30,nb,dir,4);
                        break;
                    case 2:
                        if(prev!=4)
                            findIso(target,pas,iso,P12,nb,dir,2);
                        else
                            findIso(target,pas,iso,P23,nb,dir,3);
                        break;
                    case 3:
                        if(prev!=4)
                            findIso(target,pas,iso,P12,nb,dir,2);
                        else
                            findIso(target,pas,iso,P30,nb,dir,4);
                        break;
                    case 4:
                        if(prev!=3)
                            findIso(target,pas,iso,P01,nb,dir,1);
                        else
                            findIso(target,pas,iso,P12,nb,dir,2);
                        break;
                    case 6:
                        if(prev!=3)
                            findIso(target,pas,iso,P01,nb,dir,1);
                        else
                            findIso(target,pas,iso,P23,nb,dir,3);
                        break;
                    case 7:
                        if(prev!=3)
                            findIso(target,pas,iso,P01,nb,dir,1);
                        else
                            findIso(target,pas,iso,P30,nb,dir,4);
                        break;
                }
            }
        }
        else
        {
            //Si l'algorithme devait finir lors de sa première itération (typiquement parce qu'il commence au bord de la surface), alors on avance le carré dans la direction de l'iso-ligne qui ne sort pas de la surface
            if(dir==0)
            {
                b+=0.75;
                double temp=0;
                Point Pt1 = P3;
                Point Pt2 = P3;
                Point Po01 = P3;
                Point Po02 = P3;
                Point Po1 = P0;
                Point Po2 = P2;
                bool b1 = func(Pt1.x,Pt1.y,Pt1.z)==-1;
                bool b2 = func(Pt2.x,Pt2.y,Pt2.z)==-1;
                bool b1b = abs(val(Pt1.x,Pt1.y,Pt1.z)-target)<0.001;
                bool b2b = abs(val(Pt2.x,Pt2.y,Pt2.z)-target)<0.001;
                switch(chx)
                {
                    case 1:
                        temp=0;
                        Po01 = P3;
                        Po02 = P3;
                        Po1 = P0;
                        Po2 = P2;
                        Pt1 = Po01;
                        Pt2 = Po02;
                        b1 = func(Pt1.x,Pt1.y,Pt1.z)==-1;
                        b2 = func(Pt2.x,Pt2.y,Pt2.z)==-1;
                        b1b = abs(val(Pt1.x,Pt1.y,Pt1.z)-target)<0.001;
                        b2b = abs(val(Pt2.x,Pt2.y,Pt2.z)-target)<0.001;
                        while(!(b1 && b1b) && !(b2 && b2b) && temp<=1)
                        {
                            Pt1 = interPoint(Po1,Po01,temp);
                            Pt2 = interPoint(Po2,Po02,temp);
                            b1 = func(Pt1.x,Pt1.y,Pt1.z)==-1;
                            b2 = func(Pt2.x,Pt2.y,Pt2.z)==-1;
                            b1b = abs(val(Pt1.x,Pt1.y,Pt1.z)-target)<0.001;
                            b2b = abs(val(Pt2.x,Pt2.y,Pt2.z)-target)<0.001;
                            temp+=pas*0.01;
                        }
                        if(b1 && b1b)
                        {
                            findIso(target,pas,iso,P30,nb,1,4);
                        }
                        else
                        {
                            if(b2 && b2b)
                            {
                                findIso(target,pas,iso,P23,nb,-1,3);
                            }
                        }
                        break;
                    case 2:
                        temp=0;
                        Po01 = P2;
                        Po02 = P2;
                        Po1 = P3;
                        Po2 = P1;
                        Pt1 = Po01;
                        Pt2 = Po02;
                        b1 = func(Pt1.x,Pt1.y,Pt1.z)==-1;
                        b2 = func(Pt2.x,Pt2.y,Pt2.z)==-1;
                        b1b = abs(val(Pt1.x,Pt1.y,Pt1.z)-target)<0.001;
                        b2b = abs(val(Pt2.x,Pt2.y,Pt2.z)-target)<0.001;
                        while(!(b1 && b1b) && !(b2 && b2b) && temp<=1)
                        {
                            Pt1 = interPoint(Po1,Po01,temp);
                            Pt2 = interPoint(Po2,Po02,temp);
                            b1 = func(Pt1.x,Pt1.y,Pt1.z)==-1;
                            b2 = func(Pt2.x,Pt2.y,Pt2.z)==-1;
                            b1b = abs(val(Pt1.x,Pt1.y,Pt1.z)-target)<0.001;
                            b2b = abs(val(Pt2.x,Pt2.y,Pt2.z)-target)<0.001;
                            temp+=pas*0.01;
                        }
                        if(b1 && b1b)
                        {
                            findIso(target,pas,iso,P23,nb,1,3);
                        }
                        else
                        {
                            if(b2 && b2b)
                            {
                                findIso(target,pas,iso,P12,nb,-1,2);
                            }
                        }
                        break;
                    case 3:
                        temp=0;
                        Po01 = P3;
                        Po02 = P2;
                        Po1 = P0;
                        Po2 = P1;
                        Pt1 = Po01;
                        Pt2 = Po02;
                        b1 = func(Pt1.x,Pt1.y,Pt1.z)==-1;
                        b2 = func(Pt2.x,Pt2.y,Pt2.z)==-1;
                        b1b = abs(val(Pt1.x,Pt1.y,Pt1.z)-target)<0.001;
                        b2b = abs(val(Pt2.x,Pt2.y,Pt2.z)-target)<0.001;
                        while(!(b1 && b1b) && !(b2 && b2b) && temp<=1)
                        {
                            Pt1 = interPoint(Po1,Po01,temp);
                            Pt2 = interPoint(Po2,Po02,temp);
                            b1 = func(Pt1.x,Pt1.y,Pt1.z)==-1;
                            b2 = func(Pt2.x,Pt2.y,Pt2.z)==-1;
                            b1b = abs(val(Pt1.x,Pt1.y,Pt1.z)-target)<0.001;
                            b2b = abs(val(Pt2.x,Pt2.y,Pt2.z)-target)<0.001;
                            temp+=pas*0.01;
                        }
                        if(b1 && b1b)
                        {
                            findIso(target,pas,iso,P30,nb,1,4);
                        }
                        else
                        {
                            if(b2 && b2b)
                            {
                                findIso(target,pas,iso,P12,nb,-1,2);
                            }
                        }
                        break;
                    case 4:
                        temp=0;
                        Po01 = P1;
                        Po02 = P1;
                        Po1 = P2;
                        Po2 = P0;
                        Pt1 = Po01;
                        Pt2 = Po02;
                        b1 = func(Pt1.x,Pt1.y,Pt1.z)==-1;
                        b2 = func(Pt2.x,Pt2.y,Pt2.z)==-1;
                        b1b = abs(val(Pt1.x,Pt1.y,Pt1.z)-target)<0.001;
                        b2b = abs(val(Pt2.x,Pt2.y,Pt2.z)-target)<0.001;
                        while(!(b1 && b1b) && !(b2 && b2b) && temp<=1)
                        {
                            Pt1 = interPoint(Po1,Po01,temp);
                            Pt2 = interPoint(Po2,Po02,temp);
                            b1 = func(Pt1.x,Pt1.y,Pt1.z)==-1;
                            b2 = func(Pt2.x,Pt2.y,Pt2.z)==-1;
                            b1b = abs(val(Pt1.x,Pt1.y,Pt1.z)-target)<0.001;
                            b2b = abs(val(Pt2.x,Pt2.y,Pt2.z)-target)<0.001;
                            temp+=pas*0.01;
                        }
                        if(b1 && b1b)
                        {
                            findIso(target,pas,iso,P12,nb,1,2);
                        }
                        else
                        {
                            if(b2 && b2b)
                            {
                                findIso(target,pas,iso,P01,nb,-1,1);
                            }
                        }
                        break;
                    case 6:
                        temp=0;
                        Po01 = P3;
                        Po02 = P0;
                        Po1 = P2;
                        Po2 = P1;
                        Pt1 = Po01;
                        Pt2 = Po02;
                        b1 = func(Pt1.x,Pt1.y,Pt1.z)==-1;
                        b2 = func(Pt2.x,Pt2.y,Pt2.z)==-1;
                        b1b = abs(val(Pt1.x,Pt1.y,Pt1.z)-target)<0.001;
                        b2b = abs(val(Pt2.x,Pt2.y,Pt2.z)-target)<0.001;
                        while(!(b1 && b1b) && !(b2 && b2b) && temp<=1)
                        {
                            Pt1 = interPoint(Po1,Po01,temp);
                            Pt2 = interPoint(Po2,Po02,temp);
                            b1 = func(Pt1.x,Pt1.y,Pt1.z)==-1;
                            b2 = func(Pt2.x,Pt2.y,Pt2.z)==-1;
                            b1b = abs(val(Pt1.x,Pt1.y,Pt1.z)-target)<0.001;
                            b2b = abs(val(Pt2.x,Pt2.y,Pt2.z)-target)<0.001;
                            temp+=pas*0.01;
                        }
                        if(b1 && b1b)
                        {
                            findIso(target,pas,iso,P23,nb,1,3);
                        }
                        else
                        {
                            if(b2 && b2b)
                            {
                                (target,pas,iso,P01,nb,-1,1);
                            }
                        }
                        break;
                    case 7:
                        temp=0;
                        Po01 = P3;
                        Po02 = P0;
                        Po1 = P0;
                        Po2 = P1;
                        Pt1 = Po01;
                        Pt2 = Po02;
                        b1 = func(Pt1.x,Pt1.y,Pt1.z)==-1;
                        b2 = func(Pt2.x,Pt2.y,Pt2.z)==-1;
                        b1b = abs(val(Pt1.x,Pt1.y,Pt1.z)-target)<0.001;
                        b2b = abs(val(Pt2.x,Pt2.y,Pt2.z)-target)<0.001;
                        while(!(b1 && b1b) && !(b2 && b2b) && temp<=1)
                        {
                            Pt1 = interPoint(Po1,Po01,temp);
                            Pt2 = interPoint(Po2,Po02,temp);
                            b1 = func(Pt1.x,Pt1.y,Pt1.z)==-1;
                            b2 = func(Pt2.x,Pt2.y,Pt2.z)==-1;
                            b1b = abs(val(Pt1.x,Pt1.y,Pt1.z)-target)<0.001;
                            b2b = abs(val(Pt2.x,Pt2.y,Pt2.z)-target)<0.001;
                            temp+=pas*0.01;
                        }
                        if(b1 && b1b)
                        {
                            findIso(target,pas,iso,P30,nb,1,4);
                        }
                        else
                        {
                            if(b2 && b2b)
                            {
                                findIso(target,pas,iso,P01,nb,-1,1);
                            }
                        }
                        break;
                }
            }
        }
    }
    //On dessine le carré, la couleur change en fonction des cas rencontrés pour ce carré (voir les variables r, g et b)
    Point P0 = {P.x-pas,P.y-pas,P.z};
    Point P1 = {P.x+pas,P.y-pas,P.z};
    Point P2 = {P.x+pas,P.y+pas,P.z};
    Point P3 = {P.x-pas,P.y+pas,P.z};
    glColor3f(r,g,b);
    glBegin(GL_LINE_LOOP);
    glVertex3f(P0.x,P0.y,P0.z);
    glVertex3f(P1.x,P1.y,P1.z);
    glVertex3f(P2.x,P2.y,P2.z);
    glVertex3f(P3.x,P3.y,P3.z);
    glEnd();
}

//Lance Marching Square pour obtenir une liste de points proche d'une iso-ligne
//Cette fonction est utilisé pour trouver des points de départ pour les courbes de gradient autour des extremums locaux
/**
 * target est la valeur de l'iso-ligne approximée
 * pas est la moitié des bords des carrés de Marching Square (c'est la distance minimal entre le point central et les bords du carré)
 * iso est la liste des points centraux des carrés de Marching Square
 * P est le point de départ, c'est à dire le point central du premier carré
 * dir indique la direction de déplacement le long de l'iso-ligne
 * prev indique d'où vient le carré précédent afin de ne pas faire demi-tour (ou de ne pas dédoubler continuellement le nombre de carré)
*/
void createIso(double target, double pas, vector<Point> *iso, Point P, int dir, int prev)
{
    double r=0.25, g=0.25, b=0;
    double t = func(P.x,P.y,P.z);
    bool ajout = true;
    double min = 0.0000000001;
    bool fin = false;
    for(int i=0; i<iso->size() && ajout; i++)
    {
        Point Pi = iso->at(i);
        if(dist2(P,Pi)<pas)
        {
            ajout = false;
        }
    }
    if(ajout)
        iso->push_back(P);
    else
        fin = true;
    
    if(!fin)
    {
        double l0 = func(P.x-pas,P.y-pas,P.z);
        double l1 = func(P.x+pas,P.y-pas,P.z);
        double l2 = func(P.x+pas,P.y+pas,P.z);
        double l3 = func(P.x-pas,P.y+pas,P.z);
        double v0 = target - val(P.x-pas,P.y-pas,P.z);
        double v1 = target - val(P.x+pas,P.y-pas,P.z);
        double v2 = target - val(P.x+pas,P.y+pas,P.z);
        double v3 = target - val(P.x-pas,P.y+pas,P.z);
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
        Point P0 = {P.x-pas,P.y-pas,P.z};
        Point P1 = {P.x+pas,P.y-pas,P.z};
        Point P2 = {P.x+pas,P.y+pas,P.z};
        Point P3 = {P.x-pas,P.y+pas,P.z};
        Point P01 = {P.x,P.y-2*pas,P.z};
        Point P12= {P.x+2*pas,P.y,P.z};
        Point P23 = {P.x,P.y+2*pas,P.z};
        Point P30 = {P.x-2*pas,P.y,P.z};
        string s = to_string(v0>0)+to_string(v1>0)+to_string(v2>0)+to_string(v3>0);
        string s2 = to_string(l0==-1)+to_string(l1==-1)+to_string(l2==-1)+to_string(l3==-1);
        int chx = msTable[s];
        int chx2 = msTable[s2];
        if(chx2!=0)
        {
            g += 0.5;
            Point Pf0 = P;
            Point Pf1 = P;
            Point Pt0 = P;
            Point Pt1 = P;
            Point Pt2 = P;
            Point Pt3 = P;
            double temp = 0;
            switch(chx2)
            {
                case 1:
                    Pt0 = P0;
                    Pt1 = P3;
                    Pt2 = P2;
                    Pt3 = P3;
                    break;
                case 2:
                    Pt0 = P1;
                    Pt1 = P2;
                    Pt2 = P3;
                    Pt3 = P2;
                    break;
                case 3:
                    Pt0 = P0;
                    Pt1 = P3;
                    Pt2 = P1;
                    Pt3 = P2;
                    break;
                case 4:
                    Pt0 = P0;
                    Pt1 = P1;
                    Pt2 = P2;
                    Pt3 = P1;
                    break;
                case 6:
                    Pt0 = P0;
                    Pt1 = P1;
                    Pt2 = P2;
                    Pt3 = P3;
                    break;
                case 7:
                    Pt0 = P3;
                    Pt1 = P0;
                    Pt2 = P1;
                    Pt3 = P0;
                    break;
            }
            Point Pt = Pt1;
            bool p = func(Pt.x,Pt.y,Pt.z)==-1;
            while(temp!=1 && (func(Pt.x,Pt.y,Pt.z)==-1)==p)
            {
                Pt.x = temp*Pt0.x + (1-temp)*Pt1.x;
                Pt.y = temp*Pt0.y + (1-temp)*Pt1.y;
                Pt.z = temp*Pt0.z + (1-temp)*Pt1.z;
                temp+=pas*0.1;
            }
            Pf0 = Pt;
            temp = 0;
            Pt = Pt3;
            p = func(Pt.x,Pt.y,Pt.z)==-1;
            while(temp!=1 && (func(Pt.x,Pt.y,Pt.z)==-1)==p)
            {
                Pt.x = temp*Pt2.x + (1-temp)*Pt3.x;
                Pt.y = temp*Pt2.y + (1-temp)*Pt3.y;
                Pt.z = temp*Pt2.z + (1-temp)*Pt3.z;
                temp+=pas*0.1;
            }
            Pf1 = Pt;
            temp = 0;
            while(temp!=1 && abs(val(Pt.x,Pt.y,Pt.z)-target)<0.001)
            {
                Pt.x = temp*Pf0.x + (1-temp)*Pf1.x;
                Pt.y = temp*Pf0.y + (1-temp)*Pf1.y;
                Pt.z = temp*Pf0.z + (1-temp)*Pf1.z;
                temp+=pas*0.001;
            }
            if(temp!=1)
            {
                fin = true;
                if(dir!=-1)
                {
                    iso->push_back(Pt);
                }
                else
                {
                   iso->insert(iso->begin(), Pt);
                }
            }
        }
        if(chx2==0 && l0==-1)
            fin = true;
        if(!fin)
        {
            if(dir==0)
            {
                switch(chx)
                {
                    case 1:
                        createIso(target,pas,iso,P23,-1,3);
                        createIso(target,pas,iso,P30,1,4);
                        break;
                    case 2:
                        createIso(target,pas,iso,P12,-1,2);
                        createIso(target,pas,iso,P23,1,3);
                        break;
                    case 3:
                        createIso(target,pas,iso,P12,-1,2);
                        createIso(target,pas,iso,P30,1,4);
                        break;
                    case 4:
                        createIso(target,pas,iso,P01,-1,1);
                        createIso(target,pas,iso,P12,1,2);
                        break;
                    case 6:
                        createIso(target,pas,iso,P01,-1,1);
                        createIso(target,pas,iso,P23,1,3);
                        break;
                    case 7:
                        createIso(target,pas,iso,P01,-1,1);
                        createIso(target,pas,iso,P30,1,4);
                        break;
                }
            }
            else
            {
                switch(chx)
                {
                    case 1:
                        if(prev!=1)
                            createIso(target,pas,iso,P23,dir,3);
                        else
                            createIso(target,pas,iso,P30,dir,4);
                        break;
                    case 2:
                        if(prev!=4)
                            createIso(target,pas,iso,P12,dir,2);
                        else
                            createIso(target,pas,iso,P23,dir,3);
                        break;
                    case 3:
                        if(prev!=4)
                            createIso(target,pas,iso,P12,dir,2);
                        else
                            createIso(target,pas,iso,P30,dir,4);
                        break;
                    case 4:
                        if(prev!=3)
                            createIso(target,pas,iso,P01,dir,1);
                        else
                            createIso(target,pas,iso,P12,dir,2);
                        break;
                    case 6:
                        if(prev!=3)
                            createIso(target,pas,iso,P01,dir,1);
                        else
                            createIso(target,pas,iso,P23,dir,3);
                        break;
                    case 7:
                        if(prev!=3)
                            createIso(target,pas,iso,P01,dir,1);
                        else
                            createIso(target,pas,iso,P30,dir,4);
                        break;
                }
            }
        }
        else
        {
            if(dir==0)
            {
                b+=0.75;
                double temp=0;
                Point Pt1 = P3;
                Point Pt2 = P3;
                Point Po01 = P3;
                Point Po02 = P3;
                Point Po1 = P0;
                Point Po2 = P2;
                bool b1 = func(Pt1.x,Pt1.y,Pt1.z)==-1;
                bool b2 = func(Pt2.x,Pt2.y,Pt2.z)==-1;
                bool b1b = abs(val(Pt1.x,Pt1.y,Pt1.z)-target)<0.001;
                bool b2b = abs(val(Pt2.x,Pt2.y,Pt2.z)-target)<0.001;
                switch(chx)
                {
                    case 1:
                        temp=0;
                        Po01 = P3;
                        Po02 = P3;
                        Po1 = P0;
                        Po2 = P2;
                        Pt1 = Po01;
                        Pt2 = Po02;
                        b1 = func(Pt1.x,Pt1.y,Pt1.z)==-1;
                        b2 = func(Pt2.x,Pt2.y,Pt2.z)==-1;
                        b1b = abs(val(Pt1.x,Pt1.y,Pt1.z)-target)<0.001;
                        b2b = abs(val(Pt2.x,Pt2.y,Pt2.z)-target)<0.001;
                        while(!(b1 && b1b) && !(b2 && b2b) && temp<=1)
                        {
                            Pt1 = interPoint(Po1,Po01,temp);
                            Pt2 = interPoint(Po2,Po02,temp);
                            b1 = func(Pt1.x,Pt1.y,Pt1.z)==-1;
                            b2 = func(Pt2.x,Pt2.y,Pt2.z)==-1;
                            b1b = abs(val(Pt1.x,Pt1.y,Pt1.z)-target)<0.001;
                            b2b = abs(val(Pt2.x,Pt2.y,Pt2.z)-target)<0.001;
                            temp+=pas*0.01;
                        }
                        if(b1 && b1b)
                        {
                            createIso(target,pas,iso,P30,1,4);
                        }
                        else
                        {
                            if(b2 && b2b)
                            {
                                createIso(target,pas,iso,P23,-1,3);
                            }
                        }
                        break;
                    case 2:
                        temp=0;
                        Po01 = P2;
                        Po02 = P2;
                        Po1 = P3;
                        Po2 = P1;
                        Pt1 = Po01;
                        Pt2 = Po02;
                        b1 = func(Pt1.x,Pt1.y,Pt1.z)==-1;
                        b2 = func(Pt2.x,Pt2.y,Pt2.z)==-1;
                        b1b = abs(val(Pt1.x,Pt1.y,Pt1.z)-target)<0.001;
                        b2b = abs(val(Pt2.x,Pt2.y,Pt2.z)-target)<0.001;
                        while(!(b1 && b1b) && !(b2 && b2b) && temp<=1)
                        {
                            Pt1 = interPoint(Po1,Po01,temp);
                            Pt2 = interPoint(Po2,Po02,temp);
                            b1 = func(Pt1.x,Pt1.y,Pt1.z)==-1;
                            b2 = func(Pt2.x,Pt2.y,Pt2.z)==-1;
                            b1b = abs(val(Pt1.x,Pt1.y,Pt1.z)-target)<0.001;
                            b2b = abs(val(Pt2.x,Pt2.y,Pt2.z)-target)<0.001;
                            temp+=pas*0.01;
                        }
                        if(b1 && b1b)
                        {
                            createIso(target,pas,iso,P23,1,3);
                        }
                        else
                        {
                            if(b2 && b2b)
                            {
                                createIso(target,pas,iso,P12,-1,2);
                            }
                        }
                        break;
                    case 3:
                        temp=0;
                        Po01 = P3;
                        Po02 = P2;
                        Po1 = P0;
                        Po2 = P1;
                        Pt1 = Po01;
                        Pt2 = Po02;
                        b1 = func(Pt1.x,Pt1.y,Pt1.z)==-1;
                        b2 = func(Pt2.x,Pt2.y,Pt2.z)==-1;
                        b1b = abs(val(Pt1.x,Pt1.y,Pt1.z)-target)<0.001;
                        b2b = abs(val(Pt2.x,Pt2.y,Pt2.z)-target)<0.001;
                        while(!(b1 && b1b) && !(b2 && b2b) && temp<=1)
                        {
                            Pt1 = interPoint(Po1,Po01,temp);
                            Pt2 = interPoint(Po2,Po02,temp);
                            b1 = func(Pt1.x,Pt1.y,Pt1.z)==-1;
                            b2 = func(Pt2.x,Pt2.y,Pt2.z)==-1;
                            b1b = abs(val(Pt1.x,Pt1.y,Pt1.z)-target)<0.001;
                            b2b = abs(val(Pt2.x,Pt2.y,Pt2.z)-target)<0.001;
                            temp+=pas*0.01;
                        }
                        if(b1 && b1b)
                        {
                            createIso(target,pas,iso,P30,1,4);
                        }
                        else
                        {
                            if(b2 && b2b)
                            {
                                createIso(target,pas,iso,P12,-1,2);
                            }
                        }
                        break;
                    case 4:
                        temp=0;
                        Po01 = P1;
                        Po02 = P1;
                        Po1 = P2;
                        Po2 = P0;
                        Pt1 = Po01;
                        Pt2 = Po02;
                        b1 = func(Pt1.x,Pt1.y,Pt1.z)==-1;
                        b2 = func(Pt2.x,Pt2.y,Pt2.z)==-1;
                        b1b = abs(val(Pt1.x,Pt1.y,Pt1.z)-target)<0.001;
                        b2b = abs(val(Pt2.x,Pt2.y,Pt2.z)-target)<0.001;
                        while(!(b1 && b1b) && !(b2 && b2b) && temp<=1)
                        {
                            Pt1 = interPoint(Po1,Po01,temp);
                            Pt2 = interPoint(Po2,Po02,temp);
                            b1 = func(Pt1.x,Pt1.y,Pt1.z)==-1;
                            b2 = func(Pt2.x,Pt2.y,Pt2.z)==-1;
                            b1b = abs(val(Pt1.x,Pt1.y,Pt1.z)-target)<0.001;
                            b2b = abs(val(Pt2.x,Pt2.y,Pt2.z)-target)<0.001;
                            temp+=pas*0.01;
                        }
                        if(b1 && b1b)
                        {
                            createIso(target,pas,iso,P12,1,2);
                        }
                        else
                        {
                            if(b2 && b2b)
                            {
                                createIso(target,pas,iso,P01,-1,1);
                            }
                        }
                        break;
                    case 6:
                        temp=0;
                        Po01 = P3;
                        Po02 = P0;
                        Po1 = P2;
                        Po2 = P1;
                        Pt1 = Po01;
                        Pt2 = Po02;
                        b1 = func(Pt1.x,Pt1.y,Pt1.z)==-1;
                        b2 = func(Pt2.x,Pt2.y,Pt2.z)==-1;
                        b1b = abs(val(Pt1.x,Pt1.y,Pt1.z)-target)<0.001;
                        b2b = abs(val(Pt2.x,Pt2.y,Pt2.z)-target)<0.001;
                        while(!(b1 && b1b) && !(b2 && b2b) && temp<=1)
                        {
                            Pt1 = interPoint(Po1,Po01,temp);
                            Pt2 = interPoint(Po2,Po02,temp);
                            b1 = func(Pt1.x,Pt1.y,Pt1.z)==-1;
                            b2 = func(Pt2.x,Pt2.y,Pt2.z)==-1;
                            b1b = abs(val(Pt1.x,Pt1.y,Pt1.z)-target)<0.001;
                            b2b = abs(val(Pt2.x,Pt2.y,Pt2.z)-target)<0.001;
                            temp+=pas*0.01;
                        }
                        if(b1 && b1b)
                        {
                            createIso(target,pas,iso,P23,1,3);
                        }
                        else
                        {
                            if(b2 && b2b)
                            {
                                (target,pas,iso,P01,-1,1);
                            }
                        }
                        break;
                    case 7:
                        temp=0;
                        Po01 = P3;
                        Po02 = P0;
                        Po1 = P0;
                        Po2 = P1;
                        Pt1 = Po01;
                        Pt2 = Po02;
                        b1 = func(Pt1.x,Pt1.y,Pt1.z)==-1;
                        b2 = func(Pt2.x,Pt2.y,Pt2.z)==-1;
                        b1b = abs(val(Pt1.x,Pt1.y,Pt1.z)-target)<0.001;
                        b2b = abs(val(Pt2.x,Pt2.y,Pt2.z)-target)<0.001;
                        while(!(b1 && b1b) && !(b2 && b2b) && temp<=1)
                        {
                            Pt1 = interPoint(Po1,Po01,temp);
                            Pt2 = interPoint(Po2,Po02,temp);
                            b1 = func(Pt1.x,Pt1.y,Pt1.z)==-1;
                            b2 = func(Pt2.x,Pt2.y,Pt2.z)==-1;
                            b1b = abs(val(Pt1.x,Pt1.y,Pt1.z)-target)<0.001;
                            b2b = abs(val(Pt2.x,Pt2.y,Pt2.z)-target)<0.001;
                            temp+=pas*0.01;
                        }
                        if(b1 && b1b)
                        {
                            createIso(target,pas,iso,P30,1,4);
                        }
                        else
                        {
                            if(b2 && b2b)
                            {
                                createIso(target,pas,iso,P01,-1,1);
                            }
                        }
                        break;
                }
            }
        }
    }
    Point P0 = {P.x-pas,P.y-pas,P.z};
    Point P1 = {P.x+pas,P.y-pas,P.z};
    Point P2 = {P.x+pas,P.y+pas,P.z};
    Point P3 = {P.x-pas,P.y+pas,P.z};
}

//Donne l'angle (en degrés) entre deux vecteurs
double getAngle(Point V0, Point V1)
{
    double d = dist(V0.x,V0.y,V0.z)*dist(V1.x,V1.y,V1.z);
    if(d==0)
        return 0;
    double c = (V0.x*V1.x+V0.y*V1.y+V0.z*V1.z);
    double o = acos(clamp(c/d,-1.0,1.0));
    return o*180.0/M_PI;
}

//Calcul le gradient avec une norme recadrée en un point donné
Repere getVecteurGradient(Point P, Point PrevN, double dt, double off)
{
    double val = func(P.x,P.y,P.z);
    int d = 1;
    Repere R = getRepere(P);
    Point N = R.N;
    double ang1 = abs(getAngle(PrevN,N));
    double ang2 = 0;
    bool plop = false;
    if(ang1>160)
    {
        if(ang1>260)
        {
            plop=true;
        }
        d=-1;
    }
    else
    {
        if(ang1>60 && ang1<110)
        {
            plop=true;
            double p1 = func(P.x+0.0001,P.y,P.z);
            double p2 = func(P.x,P.y+0.0001,P.z);
        }
    }
    int d2 = d;
    Point nP3 = {0,0,0};
    if(isnan(R.N.x))
    {
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
        Point nP = {P.x+N.x*r/m*d,P.y+N.y*r/m*d,P.z+N.z*r/m*d};
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
            
            double ang2 = abs(getAngle(R2.N,N));
            if(ang2>160)
            {
                d2=-d;
            }
            else
            {
                if(ang2>8)
                {
                    plop=true;
                    double p1 = func(P.x+0.0001,P.y,P.z);
                    double p2 = func(P.x,P.y+0.0001,P.z);
                    double p3 = func(nP.x+0.0001,nP.y,nP.z);
                    double p4 = func(nP.x,nP.y+0.0001,nP.z);
                    double v = func(nP.x,nP.y,nP.z);
                }
            }
            N = R2.N;
            m2 = dist(N.x,N.y,N.z);
            r2 = m2 * dt + off;
            Point nP2 = {P.x+N.x*r2/m2*d2,P.y+N.y*r2/m2*d2,P.z+N.z*r2/m2*d2};
            double mp2 = dist2(nP2,P);
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
            float inter = interpolation((double)nb2,(double)nb0);
            nP3.x = nP.x * inter + (1-inter) * nP2.x;
            nP3.y = nP.y * inter + (1-inter) * nP2.y;
            nP3.z = nP.z * inter + (1-inter) * nP2.z;
        }
        else
        {
            nP3 = nP;
        }
        N = {nP3.x-P.x,nP3.y-P.y,nP3.z-P.z};
        double c1 = 1;
        double c2 = 1;
        if(d==-1)
            c1 = 0;
        if(d2==-1)
            c2 = 0;
        glColor3f(0.5,c1,c2);
        glBegin(GL_LINES);
        glVertex3f(P.x,P.y,P.z);
        glVertex3f(nP3.x,nP3.y,nP3.z);
        glEnd();
        double l = dist(N.x,N.y,N.z);
        double m3 = (m+m2)/2.0;
        double r3 = m3/l;
        N.x*=r3;
        N.y*=r3;
        N.z*=r3;
    }
    return {nP3,N,{0,0,0}};
}

//Lance une courbe de gradient dans les deux sens en partant du point de départ donné
//La courbe s'arrête quand elle rencontre un extremum local (sauf si cet extremum fait partie d'une ligne, dans ce cas, elle repart de l'autre côté)
//La courbe stocke aussi les points qui ont une valeur de fonction pour laquelle on construit les iso-lignes
/**
 * P le point de départ de la courbe
 * dir la direction (mettre 0 pour lancer dans les deux directions)
 * PrevN le gradient au point précédent de la courbe (si dir est 0, cette valeur n'a pas d'intérêt en entrée)
 * vPc est la liste des points de la courbe qui sont soit sur une iso-ligne soit sur un bord de la surface
 * dt et off sont les valeurs permettant le recadrage des normes des gradients entre normeMin et normeMax tel que nouvelleNorme = norme * dt + off et normeMin<=nouvelleNorme<=normeMax
*/
void courbeGradient(Point P, int dir, Point PrevN, vector<PointCourbe*> *vPc, double dt, double off)
{
    int g = 0;
    double valP = func(P.x,P.y,P.z);
    Point N = {0,0,0};
    float m = dist(PrevN.x,PrevN.y,PrevN.z);
    int st = 0;
    if(dir!=0)
    {
        double m = dist(PrevN.x,PrevN.y,PrevN.z);
        double r = dt * m + off;
        Point P2 = addP(P,PrevN,-r/m);
        valP = func(P.x,P.y,P.z);
        double val2 = func(P2.x,P2.y,P2.z);
        int d=1;
        if(val2>valP)
        {
            d=-1;
        }
        bool change = false;
        while(valP>=0)
        {
            glPointSize(3.0f);
            int iv = 0;
            for(int iv=0; iv<targets.size()-1; iv++)
            {
                double v = targets[iv];
                if((val2<=v && valP>=v) || (val2>=v && valP<=v))
                {
                    Point Pv = P2;
                    double inter = 0;
                    double iMin = 0;
                    double vMin = abs(v-func(Pv.x,Pv.y,Pv.z));
                    while(inter<=1)
                    {
                        Pv.x = P.x * inter + P2.x * (1-inter);
                        Pv.y = P.y * inter + P2.y * (1-inter);
                        Pv.z = P.z * inter + P2.z * (1-inter);
                        if(abs(v-func(Pv.x,Pv.y,Pv.z)) < vMin)
                        {
                            iMin = inter;
                            vMin = abs(v-func(Pv.x,Pv.y,Pv.z));
                        }
                        inter+=0.0001;
                    }
                    Pv.x = P.x * iMin + P2.x * (1-iMin);
                    Pv.y = P.y * iMin + P2.y * (1-iMin);
                    Pv.z = P.z * iMin + P2.z * (1-iMin);
                    Point c = hsv2rgb(v,1,1);
                    glColor3f(c.x,c.y,c.z);
                    glBegin(GL_POINTS);
                    glVertex3f(Pv.x,Pv.y,Pv.z);
                    glEnd();
                    if(dir==1)
                    {
                        vPc->push_back(new PointCourbe);
                        vPc->back()->P = Pv;
                        vPc->back()->val = v;
                    }
                    else
                    {
                       vPc->insert(vPc->begin(),new PointCourbe);
                       vPc->front()->P = Pv;
                       vPc->front()->val = v;

                    }
                }
                else
                {
                    double v2 = targets[iv+1];
                    if(valP<=v2 && valP>=v && val2<=v2 && val2>=v)
                    {
                        double d0 = abs(valP-v);
                        double d1 = abs(val2-v);
                        double d2 = abs(valP-v2);
                        double d3 = abs(val2-v2);
                        if(d0<0.05 && d1<0.05)
                        {
                            Point P0 = P;
                            Point P1 = P2;
                            while(dist2(P0,P1)>normeMin*0.0001)
                            {
                                int iMin = 0;
                                double min = abs(func(P0.x,P0.y,P0.z)-v);
                                for(int i=1; i<=10; i++)
                                {
                                    double inter = (double)i/10.0;
                                    Point Pv = {0,0,0};
                                    Pv.x = P0.x * inter + P1.x * (1-inter);
                                    Pv.y = P0.y * inter + P1.y * (1-inter);
                                    Pv.z = P0.z * inter + P1.z * (1-inter);
                                    double vP = func(Pv.x,Pv.y,Pv.z);
                                    if(abs(vP-v)<min)
                                    {
                                        iMin = i;
                                        min = abs(vP-v);
                                    }
                                }
                                if(min<0.0000001)
                                {
                                    Point c = hsv2rgb(v,1,1);
                                    glColor3f(c.x,c.y,c.z);
                                    double inter = (double)iMin/10.0;
                                    Point Pv = {0,0,0};
                                    Pv.x = P0.x * inter + P1.x * (1-inter);
                                    Pv.y = P0.y * inter + P1.y * (1-inter);
                                    Pv.z = P0.z * inter + P1.z * (1-inter);
                                    glBegin(GL_POINTS);
                                    glVertex3f(Pv.x,Pv.y,Pv.z);
                                    glEnd();
                                    if(dir==1)
                                    {
                                        vPc->push_back(new PointCourbe);
                                        vPc->back()->P = Pv;
                                        vPc->back()->val = v;
                                    }
                                    else
                                    {
                                       vPc->insert(vPc->begin(),new PointCourbe);
                                       vPc->front()->P = Pv;
                                       vPc->front()->val = v;

                                    }
                                    break;
                                }
                                else
                                {
                                    if(iMin==0)
                                        iMin = 1;
                                    if(iMin==10)
                                        iMin = 9;
                                    double inter = (double)(iMin-1)/10.0;
                                    Point Pv = {0,0,0};
                                    Pv.x = P0.x * inter + P1.x * (1-inter);
                                    Pv.y = P0.y * inter + P1.y * (1-inter);
                                    Pv.z = P0.z * inter + P1.z * (1-inter);
                                    double inter2 = (double)(iMin-1)/10.0;
                                    Point Pv2 = {0,0,0};
                                    Pv2.x = P0.x * inter2 + P1.x * (1-inter2);
                                    Pv2.y = P0.y * inter2 + P1.y * (1-inter2);
                                    Pv2.z = P0.z * inter2 + P1.z * (1-inter2);
                                    P0 = Pv;
                                    P1 = Pv2;
                                }
                            }
                        }
                        if(d2<0.05 && d3<0.05)
                        {
                            Point P0 = P;
                            Point P1 = P2;
                            while(dist2(P0,P1)>normeMin*0.0001)
                            {
                                int iMin = 0;
                                double min = abs(func(P0.x,P0.y,P0.z)-v2);
                                for(int i=1; i<=10; i++)
                                {
                                    double inter = (double)i/10.0;
                                    Point Pv = {0,0,0};
                                    Pv.x = P1.x * inter + P0.x * (1-inter);
                                    Pv.y = P1.y * inter + P0.y * (1-inter);
                                    Pv.z = P1.z * inter + P0.z * (1-inter);
                                    double vP = func(Pv.x,Pv.y,Pv.z);
                                    if(abs(vP-v2)<min)
                                    {
                                        iMin = i;
                                        min = abs(vP-v2);
                                    }
                                }
                                if(min<0.0000001)
                                {
                                    Point c = hsv2rgb(v,1,1);
                                    glColor3f(c.x,c.y,c.z);
                                    double inter = (double)iMin/10.0;
                                    Point Pv = {0,0,0};
                                    Pv.x = P1.x * inter + P0.x * (1-inter);
                                    Pv.y = P1.y * inter + P0.y * (1-inter);
                                    Pv.z = P1.z * inter + P0.z * (1-inter);
                                    glBegin(GL_POINTS);
                                    glVertex3f(Pv.x,Pv.y,Pv.z);
                                    glEnd();
                                    if(dir==1)
                                    {
                                        vPc->push_back(new PointCourbe);
                                        vPc->back()->P = Pv;
                                        vPc->back()->val = v2;
                                    }
                                    else
                                    {
                                       vPc->insert(vPc->begin(),new PointCourbe);
                                       vPc->front()->P = Pv;
                                       vPc->front()->val = v2;

                                    }
                                    break;
                                }
                                else
                                {
                                    if(iMin==0)
                                        iMin = 1;
                                    if(iMin==10)
                                        iMin = 9;
                                    double inter = (double)(iMin-1)/10.0;
                                    Point Pv = {0,0,0};
                                    Pv.x = P1.x * inter + P0.x * (1-inter);
                                    Pv.y = P1.y * inter + P0.y * (1-inter);
                                    Pv.z = P1.z * inter + P0.z * (1-inter);
                                    double inter2 = (double)(iMin+1)/10.0;
                                    Point Pv2 = {0,0,0};
                                    Pv2.x = P1.x * inter2 + P0.x * (1-inter2);
                                    Pv2.y = P1.y * inter2 + P0.y * (1-inter2);
                                    Pv2.z = P1.z * inter2 + P0.z * (1-inter2);
                                    P0 = Pv;
                                    P1 = Pv2;
                                }
                            }
                        }
                    }

                }
            }
            Repere R = getVecteurGradient(P,PrevN,dt,off);
            Point nP3 = R.P;
            N = R.N;
            PrevN=N;
            P2 = P;
            P = nP3;
            val2 = valP;
            valP = func(P.x,P.y,P.z);
            if(val2>valP)
            {
                if(d==1)
                    change = true; 
            }
            else
            {
                if(d==-1)
                    change = true;   
            }
            if(change)
            {
                P2 = P;
                float mult = 0.1;
                double p1 = val(P.x+mult*normeMin,P.y,P.z);
                double p2 = val(P.x,P.y+mult*normeMin,P.z);
                double p3 = val(P.x-mult*normeMin,P.y,P.z);
                double p4 = val(P.x,P.y-mult*normeMin,P.z);
                bool fin = false;
                if(d==1)
                {
                    if(valP<=p1 && valP<=p2 && valP<=p3 && valP<=p4)
                    {
                        fin = true;
                    }
                }
                else
                {
                    if(valP>=p1 && valP>=p2 && valP>=p3 && valP>=p4)
                    {
                        fin = true;
                    }
                }
                int nc = 0;
                while(!fin){
                    int dirW = 0;
                    if(d==1)
                    {
                        if(p1>p2 && p1>p3 && p1>p4)
                        {
                            dirW = 1;
                        }
                        if(p2>p1 && p2>p3 && p2>p4)
                        {
                            dirW = 2;
                        }
                        if(p3>p1 && p3>p2 && p3>p4)
                        {
                            dirW = 3;
                        }
                        if(p4>p1 && p4>p2 && p4>p3)
                        {
                            dirW = 4;
                        }
                    }
                    else
                    {
                        if(p1<p2 && p1<p3 && p1<p4)
                        {
                            dirW = 1;
                        }
                        if(p2<p1 && p2<p3 && p2<p4)
                        {
                            dirW = 2;
                        }
                        if(p3<p1 && p3<p2 && p3<p4)
                        {
                            dirW = 3;
                        }
                        if(p4<p1 && p4<p2 && p4<p3)
                        {
                            dirW = 4;
                        }
                    }
                    if((dirW==1 && nc==3) || (dirW==3 && nc==1) || (dirW==2 && nc==4) || (dirW==4 && nc==2))
                        mult*=0.1;
                    switch(dirW)
                    {
                        case 1:
                            P.x+=normeMin*mult;
                            break;
                        case 2:
                            P.y+=normeMin*mult;
                            break;
                        case 3:
                            P.x-=normeMin*mult;
                            break;
                        case 4:
                            P.y-=normeMin*mult;
                            break;
                    }
                    nc=dirW;
                    valP = func(P.x,P.y,P.z);
                    p1 = val(P.x+mult*normeMin,P.y,P.z);
                    p2 = val(P.x,P.y+mult*normeMin,P.z);
                    p3 = val(P.x-mult*normeMin,P.y,P.z);
                    p4 = val(P.x,P.y-mult*normeMin,P.z);
                    if(dirW==0)
                        fin = true;
                }
                bool plop = false;
                double mindif = 1;
                for(int i=0; i<10000; i++)
                {
                    double a = i*0.0002*M_PI;
                    Point Pv = P;
                    Pv.x += normeMax * cos(a);
                    Pv.y += normeMax * sin(a);
                    double dv = abs(val(Pv.x,Pv.y,Pv.z)-val(P.x,P.y,P.z));
                    if(dv<0.000001)
                    {
                        if(dv<mindif)
                            mindif = dv;
                        plop = true;
                    }
                }
                if(!plop)
                {
                    valP=-2;
                    glPointSize(5.0f);
                    glColor3f(0,0,0);
                }
                else
                {
                    glPointSize(15.0f);
                    glColor3f(0,1,1);
                    d*=-1;
                    valP = func(P.x,P.y,P.z);
                    val2 = valP;
                    Point PrevN2 = addP(P,P2,-1);
                    if(dist(PrevN2.x,PrevN2.y,PrevN2.z)>0)
                    {
                        PrevN = PrevN2;
                    }
                    P2 = P;
                    P = addP(P,PrevN,1);
                    change =false;
                }
                glBegin(GL_POINTS);
                glVertex3f(P.x,P.y,P.z);
                glEnd();
                
            }
            st++;
        }
        if(valP==-1)
        {
            if(dir==1)
            {
                vPc->push_back(new PointCourbe);
                vPc->back()->P = P2;
            }
            else
            {
               vPc->insert(vPc->begin(),new PointCourbe);
               vPc->front()->P = P2;
            }
        }

    }
    else
    {
        if(valP==-1)
            return;
        Repere R = getRepere(P);
        N = R.N;
        double m = dist(N.x,N.y,N.z);
        double r = dt * m + off;
        Point nP1 = {P.x+N.x*r/m,P.y+N.y*r/m,P.z+N.z*r/m};
        Point nP2 = {P.x-N.x*r/m,P.y-N.y*r/m,P.z-N.z*r/m};
        glBegin(GL_POINTS);
        glColor3f(color,color,color);
        glVertex3f(P.x,P.y,P.z);
        glVertex3f(nP1.x,nP1.y,nP1.z);
        glVertex3f(nP2.x,nP2.y,nP2.z);
        glEnd();
        Point N2 = {-N.x,-N.y,-N.z};
        courbeGradient(nP1,1,N,vPc,dt,off);
        courbeGradient(nP2,-1,N2,vPc,dt,off);
    } 
}

//Fais la somme entre un point et un deuxième multiplié par un scalaire (les points en paramètres et le résultat peuvent être des points ou des vecteurs)
//Exemples : Point 2 + (Point 1 * -1) donne le vecteur Point 1 -> Point 2 ; Point 0 + (Vecteur * 1) donne le point à la pointe du Vecteur en partant du Point 0
Point addP(Point p1, Point p2, double mult)
{
    return {p1.x+p2.x*mult,p1.y+p2.y*mult,p1.z+p2.z*mult};
}

//Calcul la longueur d'une courbe de Bézier définie par 4 points de contrôle donnés et peut la déssiner
/**
 * P0, P1, P2, P3 sont les points de contrôle dans l'ordre
 * dessin indique si la courbe doit être déssinée (true pour dessiner)
*/ 
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
        if(dessin)
            glVertex3f(P.x,P.y,P.z);
    }
    glEnd();
    return longueur;
}

//Cherche le point sur une courbe de Bézier le plus proche du point donné
/**
 * start et end sont les valeurs de t qui donnent la délimitations de la zone de recherche sur la courbe de Bézier
 * dist est la distance recherché entre le point donné et le point le plus proche souhaité
 * On arrête la recherche quand la différence entre start et end devient suffisament faible
*/
Point chercherBezier(half_edge* he, Point Pc, double dist)
{
    int step = 0;
    Point P0 = he->origine;
    Point P1 = he->controleO;
    Point P2 = he->controleI;
    Point P3 = he->incident;
    Point Pr = {0,0,0};
    double start = 0;
    double end = 1;
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
        for(double t=start; t<=end; t+=p)
        {
            Pr3 = addP({0,0,0},P0,(1-t)*(1-t)*(1-t));
            Pr3 = addP(Pr3,P1,3*t*(1-t)*(1-t));
            Pr3 = addP(Pr3,P2,3*t*t*(1-t));
            Pr3 = addP(Pr3,P3,t*t*t);
            double d = dist2(Pr3,Pc);
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
        if(start<0)
            start = 0;
        if(start>1)
            start = 1;
        if(end<0)
            end = 0;
        if(end>1)
            end = 1;
    }while(distMin>dist && abs(end-start)>0.000001);
    return Pr;
}

//Cherche le point sur un des bord de la cellule le plus proche du point donné
/***
 * dist est la distance recherché entre le point donné et le point le plus proche souhaité
*/
Point distanceMinBezier(heFace f, Point P, float dist)
{
    half_edge *he = f.incidente;
    int id = he->id;
    Point ret = he->origine;
    float mindist = dist2(ret,P);
    do
    {
        Point Pr = chercherBezier(he,P,dist);
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

//Donne le vecteur entre deux points (le premier est C(t) et le deuxième est C(t+dif)) de la courbe de Bézier 
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

//Applique une rotation d'un angle (en radians) à un vecteur (stocké dans la structure Point)
Point rotation(Point V, double angle)
{
    Point R = {0,0,0};
    R.x = cos(angle)*V.x - sin(angle)*V.y;
    R.y = sin(angle)*V.x + cos(angle)*V.y;
    R.z = V.z;
    return R;
}

//Calcul la densité moyenne dans une cellule
double densite(heFace f)
{
    half_edge* he1 = f.incidente;
    int id_ori = he1->id;
    half_edge* he2 = he1->previous;
    double minX = he1->origine.x;
    double maxX = he1->origine.x;
    double minY = he1->origine.y;
    double maxY = he1->origine.y;
    Point Somme = {0,0,0};
    int nbE = 0;
    float limite = normeMax/2;
    half_edge* he = he1;
    int pos = 0;
    int neg = 0;
    //On trouve le rectangle englobant la cellule et on trouve le sens des angles des bords de la cellule afin d'identifier ceux qui font que la cellule est concave (si elle l'est)
    do
    {
        he2 = he1;
        he1 = he1->next;
        Point P0 = he1->origine;
        Point P1 = he1->controleO;
        Point P2 = he1->controleI;
        Point P3 = he1->incident;
        if(P0.x<minX)
            minX = P0.x-limite*2;
        if(P0.x>maxX)
            maxX = P0.x+limite*2;
        if(P1.x<minX)
            minX = P1.x-limite*2;
        if(P1.x>maxX)
            maxX = P1.x+limite*2;
        if(P2.x<minX)
            minX = P2.x-limite*2;
        if(P2.x>maxX)
            maxX = P2.x+limite*2;
        if(P3.x<minX)
            minX = P3.x-limite*2;
        if(P3.x>maxX)
            maxX = P3.x+limite*2;
        if(P0.y<minY)
            minY = P0.y-limite*2;
        if(P0.y>maxY)
            maxY = P0.y+limite*2;
        if(P1.y<minY)
            minY = P1.y-limite*2;
        if(P1.y>maxY)
            maxY = P1.y+limite*2;
        if(P2.y<minY)
            minY = P2.y-limite*2;
        if(P2.y>maxY)
            maxY = P2.y+limite*2;
        if(P3.y<minY)
            minY = P3.y-limite*2;
        if(P3.y>maxY)
            maxY = P3.y+limite*2;
        Somme = addP(Somme,P0,1);
        nbE++;
        Point V1 = getVecteurBezier(he1,0,limite*2);
        Point V2 = getVecteurBezier(he2,1,-limite*2);
        double n1 = dist(V1.x,V1.y,V1.z);
        double n2 = dist(V2.x,V2.y,V2.z);
        if(n1!=0 && n2!=0)
        {
            double ang = getAngle(V1,V2);
            if(ang!=0)
            {
                Point Vr = rotation(V1,(ang/180.0)*M_PI);
                double ang2 = getAngle(Vr,V2);
                if(ang2<ang)
                {
                    pos++;
                }
                else
                {
                    neg++;
                }
            }
        }
    }while(he1->id!=id_ori);
    Somme = addP({0,0,0},Somme,1/(double)nbE);
    int sens = -1;
    if(pos>neg)
    {
        sens = 1;
    }
    he1 = f.incidente;
    he2 = he1->previous;
    double longu = -1;
    //On choisi un coin de la cellule tel qu'il n'est pas concave par rapport au reste de la cellue, qu'il ait un angle élevé et que ses côtés adjacents soient de longueur suffisante
    //On calcule aussi une approximation du barycentre de la face (à partir des coins)
    do{
        Point P = he1->origine;
        Point P2 = he1->controleO;
        Point P3 = he2->controleI;
        Point V0 = addP(he1->controleO,he1->origine,-1);
        Point V1 = addP(he2->controleI,he1->origine,-1);
        double a1 = getAngle(V0,V1);
        Point V01 = rotation(V0,M_PI*a1/180.0*sens);
        double a01 = getAngle(V01,V1);
        if(a01<a1 && (longu==-1 || (dist2(P,P2)+dist2(P,P3))*a1*a1>longu))
        {
            he = he1;
            longu = (dist2(P,P2)+dist2(P,P3))*a1*a1;
        }
        he2 = he1;
        he1 = he1->next;
    }while(he1->id!=id_ori);
    he1 = he;
    he2 = he1->previous;
    Point P = he1->origine;
    Point v0 = getVecteurBezier(he1,0,limite*2);
    Point v1 = getVecteurBezier(he2,1,-limite*2);
    double ang = getAngle(v0,v1);
    Point v2 = rotation(v0,(ang/360.0)*sens*M_PI);
    Point Po = {0,0,0};
    //On prend notre point de départ depuis le coin choisi précédemment, sur la bissectrice de l'angle du coin
    Po = addP(P,v2,1);
    Point Pc = distanceMinBezier(f,Po,limite/2);
    //On déplace notre point le long de la bissectrice jusqu'à l'éloigner suffisamment des bords de la cellue
    while(dist2(Po,Pc)<limite)
    {
        Po = addP(Po,v2,1);
        Pc = distanceMinBezier(f,Po,limite/2);
    }
    //Si le point de départ est hors de la surface, on le rapproche du barycentre (approximation) de la face
    if(func(Po.x,Po.y,Po.z)==-1)
    {
        Point V = addP(Somme, Po, -1.0);
        float n = dist(V.x,V.y,V.z);
        V.x/=n;
        V.y/=n;
        V.z/=n;
        V = addP({0,0,0},V,0.5*limite);
        while(func(Po.x,Po.y,Po.z)==-1)
        {
            Po = addP(Po,V,1);
        }
        Po = addP(Po,V,1);
        Pc = distanceMinBezier(f,Po,limite/2);
    }
    //Si le point se retrouve trop proche d'un bord, on l'éloigne du bord et de son coin d'origine
    while(dist2(Pc,Po)<limite)
    {
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
        Pc = distanceMinBezier(f,Po,limite/2);
    }

    //On propage les points via une pile
    marque.clear();
    pile.push_back({0,0,0});
    RetourDensite rd = {0,0};
    while(pile.size()!=0)
    {
        Point Pc = pile[0];
        int x = (int)Pc.x;
        int y = (int)Pc.y;
        Point P = {Po.x+x*limite/4,Po.y+y*limite/4,Po.z};
        //On vérifie si on se trouve toujours dans le rectangle englobant (on en sort seulement en cas d'erreur, généralement un mauvais placement ou déplacements du point de départ)
        bool nope = P.x<minX || P.x>maxX || P.y<minY || P.y>maxY;
        string k = to_string(x)+":"+to_string(y);
        //Si le point n'a pas déjà été traité et qu'il n'est pas sorti du rectangle englobant, alors on le traite
        if(marque.find(k) == marque.end() && !nope)
        {
            //On marque le point
            marque[k]=true;
            //On ajoute sa valeur dans le calcul de la densité moyenne
            float value = val(P.x,P.y,P.z);
            rd.nb += 1;
            rd.somme += value;
            //On ne propage notre point que s'il est suffisament éloigné du bord
            Point B = distanceMinBezier(f,P,limite/2);
            if(dist2(B,P)>=limite)
            {
                pile.push_back({(float)x-1,(float)y,0});
                pile.push_back({(float)x+1,(float)y,0});
                pile.push_back({(float)x,(float)y-1,0});
                pile.push_back({(float)x,(float)y+1,0});
            }
            
        }
        pile.erase(pile.begin());
        if(nope)
        {
            debug<<"\t\tErreur : sortie de la cellue"<<endl;
            return -1.0;
        }
    }
    if(rd.nb==0)
        return -1.0;
    return rd.somme / (float)rd.nb;
}

//Cherche un extremum local
/**
 * Po est le point de départ de la recherche
 * ordre indique si on cherche un maximum ou un minimum local (1 pour un maximum, autres pour un minimum)
 * extrem est la liste où sont stockés les extremums
 * lim est le nombre maximum d'itération de la recherche avec le pas maximum (utilisé pour éviter les chevauchements des zones de recherches)
*/
void findExtrem(Point Po, int ordre, vector<Point> *extrem, int lim)
{
    int st = 0;
    Point P = Po;
    double valP = func(P.x,P.y,P.z);
    float mult = 0.1;
    double p1 = val(P.x+mult*normeMin,P.y,P.z);
    double p2 = val(P.x,P.y+mult*normeMin,P.z);
    double p3 = val(P.x-mult*normeMin,P.y,P.z);
    double p4 = val(P.x,P.y-mult*normeMin,P.z);
    bool fin = false;
    if(ordre==1)
    {
        if(valP<=p1 && valP<=p2 && valP<=p3 && valP<=p4)
        {
            fin = true;
        }
    }
    else
    {
        if(valP>=p1 && valP>=p2 && valP>=p3 && valP>=p4)
        {
            fin = true;
        }
    }
    int nc = 0;
    while(!fin){
        int dirW = 0;
        if(ordre==1)
        {
            if(p1>p2 && p1>p3 && p1>p4)
            {
                dirW = 1;
            }
            if(p2>p1 && p2>p3 && p2>p4)
            {
                dirW = 2;
            }
            if(p3>p1 && p3>p2 && p3>p4)
            {
                dirW = 3;
            }
            if(p4>p1 && p4>p2 && p4>p3)
            {
                dirW = 4;
            }
        }
        else
        {
            if(p1<p2 && p1<p3 && p1<p4)
            {
                dirW = 1;
            }
            if(p2<p1 && p2<p3 && p2<p4)
            {
                dirW = 2;
            }
            if(p3<p1 && p3<p2 && p3<p4)
            {
                dirW = 3;
            }
            if(p4<p1 && p4<p2 && p4<p3)
            {
                dirW = 4;
            }
        }
        double f1 = func(P.x+mult*normeMin,P.y,P.z);
        double f2 = func(P.x,P.y+mult*normeMin,P.z);
        double f3 = func(P.x-mult*normeMin,P.y,P.z);
        double f4 = func(P.x,P.y-mult*normeMin,P.z);
        while((dirW==1 && f1==-1) || (dirW==3 && f3==-1) || (dirW==2 && f2==-1) || (dirW==4 && f4==-1)){
            if((dirW==1 && f1==-1) || (dirW==3 && f3==-1) || (dirW==2 && f2==-1) || (dirW==4 && f4==-1))
                mult*=0.1;
            f1 = func(P.x+mult*normeMin,P.y,P.z);
            f2 = func(P.x,P.y+mult*normeMin,P.z);
            f3 = func(P.x-mult*normeMin,P.y,P.z);
            f4 = func(P.x,P.y-mult*normeMin,P.z);
            if(dirW==0 || mult<=0.0000001)
            {
                fin = true;
                break;
            }
        }
        if((dirW==1 && nc==3) || (dirW==3 && nc==1) || (dirW==2 && nc==4) || (dirW==4 && nc==2))
            mult*=0.1;
        switch(dirW)
        {
            case 1:
                P.x+=normeMin*mult;
                break;
            case 2:
                P.y+=normeMin*mult;
                break;
            case 3:
                P.x-=normeMin*mult;
                break;
            case 4:
                P.y-=normeMin*mult;
                break;
        }
        nc=dirW;
        valP = func(P.x,P.y,P.z);
        p1 = val(P.x+mult*normeMin,P.y,P.z);
        p2 = val(P.x,P.y+mult*normeMin,P.z);
        p3 = val(P.x-mult*normeMin,P.y,P.z);
        p4 = val(P.x,P.y-mult*normeMin,P.z);
        st++;
        if(mult<0.1)
            st=0;
        if(dirW==0 || mult<0.00000000001 || st>=lim)
        {
            fin = true;
        }
    }
    if(st<lim)
    {
        bool fin = false;
        
        for(int i=0; i<100; i++)
        {
            double a = i*M_PI/50;
            Point Pg = {P.x+0.001*cos(a),P.y+0.001*sin(a),P.z};
            double f = func(Pg.x,Pg.y,Pg.z);
            if(f!=-1)
            {
                if(ordre==1)
                {
                    if(f>=valP)
                        fin = true;
                }
                else
                {
                    if(f<=valP)
                        fin = true;
                }
            }
        }
        if(!fin)
            extrem->push_back(P);
    }
}

//Vérifie si deux PointCourbes contiennent les mêmes données
bool equal(PointCourbe Pc, PointCourbe Pc2)
{
    bool b = Pc.isoC==Pc2.isoC;
    b = b && Pc.isoInd==Pc2.isoInd;
    b = b && Pc.graC==Pc2.graC;
    b = b && Pc.indGra==Pc2.indGra;
    b = b && (abs(Pc.val-Pc2.val)<0.00001);
    return b;
}

//Calcul la distance entre la projection d'un point sur la droite passant par deux points et l'un des deux points
double distProj(Point P0, Point P1, Point Pt)
{
    Point V = addP(P1,P0,-1);
    Point Vt = addP(Pt,P0,-1);
    double ang = abs(getAngle(V,Vt));
    double rad = ang * (M_PI/180);
    double distP = cos(rad)*dist2(P0,Pt);
    return distP;
}

//------------------------------------------------------
void affichage(void)
{
    debug.open ("debug",ios::out | ios::trunc);
    Faces.clear();
    hes.clear();
    color = 0;
    surface.clear();
    targets.clear();
    newCourbeFinale.clear();
    NewCourbeGrad.clear();
    glPointSize(7.0f);
	glMatrixMode(GL_MODELVIEW);
    /* effacement de l'image avec la couleur de fond */
	glClear(GL_COLOR_BUFFER_BIT);
	glPushMatrix();
    	glTranslatef(decX,decY,cameraDistance);
        glRotatef(cameraAngleX,1.,0.,0.);
        glRotatef(cameraAngleY,0.,1.,0.);
        cout<<"Commencement (détails dans le fichier debug)"<<endl;
        debug<<"Commencement"<<endl;
        //Placer la légende des couleurs représentant les différentes valeurs de la fonction (pour les points des courbes de gradient) de rouge (valeur 0) vers le vert puis bleu et finalement retour au rouge à la valeur 1
        for(int i=0; i<=10; i++)
        {
            double v = (double)i/10.0;
            Point C = hsv2rgb(v,1,1);
            glColor3f(C.x,C.y,C.z);
            glPointSize(5.0f);
            glBegin(GL_POINTS);
            glVertex3f(1.1f,v,0);
            glEnd();
        }
        float n = 1000;

        debug<<"Bords de la surface"<<endl;
        //Créer les bords de la surface
        for(int i=-10; i<10; i++)
        {
            Point P = {-1.0f ,i*0.1f ,0};
            surface.push_back(P);
        }
        for(int i=-10; i<10; i++)
        {
            Point P = {i*0.1f ,1.0f ,0};
            surface.push_back(P);
        }
        for(int i=10; i>-10; i--)
        {
            Point P = {1.0f ,i*0.1f,0};
            surface.push_back(P);
        }
        for(int i=10; i>-10; i--)
        {
            Point P = {i*0.1f ,-1.0f ,0};
            surface.push_back(P);
        }
        debug<<"\t"<<surface.size()<<" points"<<endl;

        debug<<"Valeurs maximale et minimale de la fonction"<<endl;
        //On trouve les valeurs minimales et maximales de la fonction (après discrétisation)
        float Vmin = 1;
        float Vmax = 0;
        for(int i=-n; i<=n; i++)
        {
            for(int j=-n; j<=n; j++)
            {
                float var = func((float)(i)/(float)(n),(float)(j)/(float)(n),0);
                if(var<Vmin)
                {
                    Vmin = var;
                }
                if(var>Vmax)
                {
                    Vmax=var;
                }
            }
        }
        debug<<"\tMaximum = "<<Vmax<<" et Minimum = "<<Vmin<<endl;

        //On dessine la surface avec des couleurs pour représenter la fonction : noir pour 0 et blanc pour 1
        //On identifie aussi la norme de gradient maximum et minimum
        //On trouve aussi les extremums locaux
        vector<Point> extremum0;
        vector<Point> extremum1;
        double GradMax = -1;
        double GradMin = -1;
        int lim = ((int)floor((2.0/2001.0)/(normeMin*0.06)));
        cout<<"Recherche des extremums"<<endl;
        debug<<"Recherche des extremums"<<endl;
        double augm = (1.0/(double)(2*n))*100;
        for(int i=-n; i<=n; i++)
        {
            double prec = ((double)(i+n)/(double)(2*n))*100;
            if(floor(prec)-(((int)(round(prec)))%10)>prec-augm)
                debug<<"\t"<<floor(prec)-(((int)(round(prec)))%10)<<"%"<<endl;
            for(int j=-n; j<=n; j++)
            {
                float xb = (float)(i)/(float)(n);
                float yb = (float)(j)/(float)(n);
                Point P = {xb,yb,0};
                //Trouver les extremums
                findExtrem(P,-1,&extremum0,lim);
                findExtrem(P,1,&extremum1,lim);
                //Calculer les normes de gradients
                Repere R = getRepere(P);
                double m = dist(R.N.x,R.N.y,R.N.z);
                if(GradMax==-1 || m>GradMax)
                    GradMax = m;
                if(GradMin==-1 || m<GradMin)
                    GradMin = m;
                //Affichage
                if(i!=n && j!=n)
                {
                    float dn = 1.0f/(float)n;
                    float var = func(xb,yb,0);
                    glBegin(GL_QUADS);
                    glColor3f((var-Vmin)/(Vmax-Vmin),(var-Vmin)/(Vmax-Vmin),(var-Vmin)/(Vmax-Vmin));
                    glVertex3f(xb,yb,0);

                    var = func(xb+dn,yb,0);
                    glColor3f((var-Vmin)/(Vmax-Vmin),(var-Vmin)/(Vmax-Vmin),(var-Vmin)/(Vmax-Vmin));
                    glVertex3f(xb+dn,yb,0);

                    var = func(xb+dn,yb+dn,0);
                    glColor3f((var-Vmin)/(Vmax-Vmin),(var-Vmin)/(Vmax-Vmin),(var-Vmin)/(Vmax-Vmin));
                    glVertex3f(xb+dn,yb+dn,0);

                    var = func(xb,yb+dn,0);
                    glColor3f((var-Vmin)/(Vmax-Vmin),(var-Vmin)/(Vmax-Vmin),(var-Vmin)/(Vmax-Vmin));
                    glVertex3f(xb,yb+dn,0);
                    glEnd();
                }
            }
        }

        debug<<"Recadrage de la norme des gradient"<<endl;
        //On calcule la fonction de recadrage la norme du gradient
        double dt = 1;
        double off = 0;
        if(GradMax!=GradMin)
        {
            dt = (normeMin-normeMax)/(GradMin-GradMax);
            off = normeMin - dt*GradMin;
        }
        else
        {
            dt = 1/(GradMax*10);
        }
        debug<<"Fonction : nouvelleNorme = "<<dt<<" * norme + "<<off<<endl;


        //On supprime les extremums en double et ceux qui ne sont pas un point extrême (comparé à une ligne de valeur extrême)
        //extremum0 contient les minimums locaux et extremum1 les maximums locaux
        vector<Point> e0;
        vector<Point> e1;
        bool plopU = true;
        cout<<"Traitement des extremums"<<endl;
        debug<<"Traitement des extremums (suppression des doublons et autres valeurs non conformes)"<<endl;
        debug<<"Origalement :\n\tN° minimums locaux = "<<extremum0.size()<<endl<<"\tN° maximums locaux = "<<extremum1.size()<<endl;
        for(int i=0; i<extremum0.size(); i++)
        {
            Point P0 = extremum0[i];
            bool plop = true;
            for(int j=0; j<i; j++)
            {
                Point P1 = extremum0[j];
                if(dist2(P0,P1)<0.0001)
                    plop = false;
                else
                {
                    if(abs(val(P0.x,P0.y,P0.z)-val(P1.x,P1.y,P1.z))<0.001)
                    {
                        bool plop2 = true;
                        for(int k=1; k<10 && plop2; k++)
                        {
                            double t = 0.1*k;
                            Point P2 = interPoint(P0,P1,t);
                            if(abs(val(P0.x,P0.y,P0.z)-val(P2.x,P2.y,P2.z))>0.001)
                                plop2 = false;
                        }
                        if(plop2)
                        {
                            plop = false;
                            if(j==extremum0.size()-2 && i==extremum0.size()-1)
                                plopU = false;
                        }
                    }
                }
                
            }
            if(plop && plopU)
                e0.push_back(P0);
        }
        extremum0.clear();
        for(int i=0; i<e0.size(); i++)
        {
            Point P0 = e0[i];
            extremum0.push_back(e0[i]);
        }
        e0.clear();
        plopU = true;
        for(int i=0; i<extremum1.size(); i++)
        {
            Point P0 = extremum1[i];
            bool plop = true;
            for(int j=0; j<i; j++)
            {
                Point P1 = extremum1[j];
                if(dist2(P0,P1)<0.0001)
                    plop = false;
                else
                {
                    if(abs(val(P0.x,P0.y,P0.z)-val(P1.x,P1.y,P1.z))<0.001)
                    {
                        bool plop2 = true;
                        for(int k=1; k<10 && plop2; k++)
                        {
                            double t = 0.1*k;
                            Point P2 = interPoint(P0,P1,t);
                            if(abs(val(P0.x,P0.y,P0.z)-val(P2.x,P2.y,P2.z))>0.001)
                                plop2 = false;
                        }
                        if(plop2)
                        {
                            plop = false;
                            if(j==extremum0.size()-2 && i==extremum0.size()-1)
                                plopU = false;
                        }
                    }
                }
                
            }
            if(plop && plopU)
                e1.push_back(P0);
        }
        extremum1.clear();
        for(int i=0; i<e1.size(); i++)
        {
            Point P0 = e1[i];
            extremum1.push_back(e1[i]);
        }
        e1.clear();
        debug<<"Finalement :\n\tN° minimums locaux = "<<extremum0.size()<<endl<<"\tN° maximums locaux = "<<extremum1.size()<<endl;

        //On ordonne les extremums par ordre croissant de différence de valeur avec les maximums et minimums globaux (0 et 1)
        //Exemple : 0.2 est à 0.2 de 0 et 0.85 est à 0.15 de 1. On mets donc 0.85 avant 0.2

        debug<<"Traitement des extremums (traitement des valeurs des extremums et recherche d'iso-lignes proches)"<<endl;
        vector<int> ordre;
        map<int,vector<Point>> points;
        map<string,double> mapVals;
        glPointSize(3.0f);
        glBegin(GL_POINTS);
        glColor3f(0,1,0);
        for(int i=0; i<extremum0.size(); i++)
        {
            Point Em = extremum0[i];
            double v0 = val(Em.x,Em.y,Em.z);
            double v1 = 1-v0;
            double v = v1;
            if(v0<v1)
                v = v0;
            int v2 = (int)(v*1000);
            bool ok = false;
            for(auto it = ordre.begin(); it!=ordre.end() && !ok; it++)
            {
                int x = *it;
                if(x==v2)
                {
                    ok = true;
                    break;
                }
                if(x>v2)
                {
                    ok = true;
                    ordre.insert(it,v2);
                }
            }
            if(ordre.size()==0)
            {
                ordre.push_back(v2);
                ok = true;
            }
            if(!ok)
                ordre.push_back(v2);
            vector<Point> c = points[v2];
            c.push_back(Em);
            points[v2] = c;
            glVertex3f(Em.x,Em.y,Em.z);

        }
        glColor3f(0,0,1);
        for(int i=0; i<extremum1.size(); i++)
        {
            Point Em = extremum1[i];
            double v0 = val(Em.x,Em.y,Em.z);
            double v1 = 1-v0;
            double v = v1;
            if(v0<v1)
                v = v0;
            int v2 = (int)(v*1000);
            bool plop = false;
            for(auto it = ordre.begin(); it!=ordre.end(); it++)
            {
                int x = *it;
                if(x==v2)
                {
                    plop = true;
                    break;
                }
                if(x>v2)
                {
                    plop = true;
                    ordre.insert(it,v2);
                    break;
                }
            }
            if(ordre.size()==0)
            {
                ordre.push_back(v2);
                plop = true;
            }
            if(!plop)
                ordre.push_back(v2);
            vector<Point> c = points[v2];
            c.push_back(Em);
            points[v2] = c;
            glVertex3f(Em.x,Em.y,Em.z);
        }
        glEnd();


        //On choisie la valeur d'un point autour de chaque extremum pour l'ajouter à la liste des valeurs d'iso-lignes (sauf si une valeur proche à déjà été ajoutée)
        vector<double> vals;
        for(int i=0; i<ordre.size(); i++)
        {
            double rayon = normeMax*10;
            int x = ordre[i];
            vector<Point> c = points[x];
            for(int j=0; j<c.size(); j++)
            {
                double max = 0;
                double val = 0;
                Point P = c[j];
                double vP = func(P.x,P.y,P.z);
                glColor3f(1,0,0);
                glBegin(GL_LINE_LOOP);
                for(int k=0; k<100; k++)
                {
                    double a = k*M_PI/50;
                    Point Pg = {P.x+rayon*cos(a),P.y+rayon*sin(a),P.z};
                    double va = func(Pg.x,Pg.y,Pg.z);
                    double v = abs(va-vP);
                    if(v>max && va!=-1)
                    {
                        max = v;
                        val = va;
                    }
                    glVertex3f(Pg.x,Pg.y,Pg.z);
                }
                mapVals[to_string(x)+":"+to_string(j)] = val;
                glEnd();
                bool deja = false;
                for(int k=0; k<vals.size() && !deja; k++)
                {
                    double v = vals[k];
                    if(abs(val-v)<0.000001)
                    {
                        deja = true;
                    }
                }
                if(!deja)
                {
                    vals.push_back(val);
                }
            }
        }
        
        debug<<"Remplissage des valeurs pour les iso-lignes"<<endl;
        //On ajoute les valeurs précédentes à la liste des valeurs d'iso-lignes dans l'ordre croissant
        for(int i=0; i<=10; i++)
        {
            double v = i*0.1;
            targets.push_back(v);
            debug<<"\t"<<v<<" (valeur de base)"<<endl;
            for(int k=0; k<vals.size(); k++)
            {
                double va = vals[k];
                if(va>v && va<(v+0.1))
                {
                    targets.push_back(va);
                    debug<<"\t"<<va<<" (valeur pour les extremums)"<<endl;
                }
            }
        }

        cout<<"Lancement des courbes de gradient"<<endl;
        debug<<"Lancement des courbes de gradient autour des extremums"<<endl;
        //On lance des courbes de gradients le long de l'iso-lignes des valeurs autour des extremums (l'iso-ligne est approximé par une liste des centres des carrrés de Marching Square, et elle n'est pas gardée)
        double pasIso = 0.01;
        for(int i=0; i<ordre.size(); i++)
        {
            int x = ordre[i];
            vector<Point> c = points[x];
            for(int j=0; j<c.size(); j++)
            {
                double val = mapVals[to_string(x)+":"+to_string(j)];
                debug<<"\tCourbes de gradient pour l'iso-ligne n°"<<j<<" pour la valeur "<<val<<endl;
                Point P = c[j];
                double vP = func(P.x,P.y,P.z);
                glColor3f(1,0,0);
                int nb = 100;
                Point start = P;
                for(int k=0; k<nb; k++)
                {
                    double a = k*M_PI/((double)nb/2.0);
                    double rayon = normeMax*10;
                    Point Pg = {P.x+rayon*cos(a),P.y+rayon*sin(a),P.z};
                    double va = func(Pg.x,Pg.y,Pg.z);
                    if(abs(va-val)<0.0001)
                        start = Pg;
                }
                vector<Point> isoT;
                createIso(val,pasIso,&isoT,start, 0, 0);
                int nbP = isoT.size();
                int nbD = nbP / 3.0;
                double pas = (double)(nbP)/((double)nbD);
                if(pas<1)
                    pas=1;
                int nbG = 0;
                for(double k=0; k<nbP; k+=pas)
                {
                    int ind = (int)round(k);
                    Point Pg = isoT[ind];
                    vector<PointCourbe*> cG;
                    bool ok = true;
                    double distMin = 10;
                    for(int l=0; l<NewCourbeGrad.size();l++)
                    {
                        vector<PointCourbe*> c = NewCourbeGrad[l];
                        for(int m=0; m<c.size(); m++)
                        {
                            PointCourbe* Pc = c[m];
                            if(dist2(Pc->P, Pg)<distMin)
                                distMin = dist2(Pc->P, Pg);
                            if(dist2(Pc->P, Pg)<pasIso*2)
                            {
                                ok = false;
                            }
                        }
                    }
                    if(ok && func(Pg.x,Pg.y,Pg.z)>=0)
                    {
                        courbeGradient(Pg,0,{0,0,0},&cG,dt,off);
                        PointCourbe* Pc0 = cG[0];
                        PointCourbe* Pc1 = cG[cG.size()-1];
                        if(cG.size()>2)
                        {
                            NewCourbeGrad.push_back(cG);
                            nbG++;
                        }
                        
                    }
                }
                debug<<"\t\t"<<nbG<<" courbes lancés"<<endl;
            }
        }
        
        debug<<"Lancement des courbes de gradient le long des bords"<<endl;
        //On lance des courbes de gradients depuis les bords de la surface, s'il n'y a aucune courbes de gradient proches
        double dmin = dist2(surface[0],surface[1])/2;
        for(int i=0; i<surface.size();i++)
        {
            debug<<"\t Depuis le point "<<i<<" du bord : ";
            vector<PointCourbe*> cG;
            Point P = surface[i];
            bool ok = true;
            for(int j=0; j<NewCourbeGrad.size();j++)
            {
                vector<PointCourbe*> c = NewCourbeGrad[j];
                Point P0 = c[0]->P;
                Point P1 = c[c.size()-1]->P;
                double d0 = dist2(P,P0);
                double d1 = dist2(P,P1);
                if(d0<dmin || d1<dmin)
                    ok=false;
            }
            if(ok)
            {
                debug<<"Courbe lancée"<<endl;
                courbeGradient(P,0,{0,0,0},&cG,dt,off);
                if(cG.size()>2)
                {
                    NewCourbeGrad.push_back(cG);
                }
                glColor3f(0,1,0);
            }
            else
            {
                debug<<"Pas de courbe"<<endl;
                glColor3f(0,0,1);
            }
            glBegin(GL_POINTS);
            glVertex3f(P.x,P.y,P.z);
            glEnd();
        }
        
        cout<<"Construction des iso-lignes"<<endl;
        debug<<"Construction des iso-lignes"<<endl;
        //Pour chaque point (sauf ceux sur les bords) de chaque courbe de gradient, s'il n'est pas déjà sur une iso-ligne identifiée, on calcule l'iso-ligne dont il fait partie
        for(int i = 0; i<NewCourbeGrad.size(); i++)
        {
            vector<PointCourbe*> c = NewCourbeGrad[i];
            debug<<"\tLe long de la courbe de gradient n°"<<i<<endl;
            for(int j = 0; j<c.size(); j++)
            {
                c[j]->graC = i;
                c[j]->indGra = j;
                if(c[j]->isoC==-1)
                {
                    double val = c[j]->val;
                    if(val!=-1)
                    {
                        Point P = c[j]->P;
                        vector<PointCourbe*> iso;
                        int d = newCourbeFinale[(int)(round(val*1000))].size();
                        findIso(val, pasIso, &iso, P, d, 0,0);
                        glPointSize(5.0f);
                        glColor3f(1,0,1);
                        glBegin(GL_POINTS);
                        for(int k=0;k<iso.size(); k++)
                        {
                            Point Po = iso[k]->P;
                            glVertex3f(Po.x,Po.y,Po.z);
                        }
                        glEnd();
                        debug<<"\t\tAjout de l'iso-ligne n°"<<d<<" pour la valeur "<<val<<". Nombre de points : "<<iso.size()<<endl;
                        vector<vector<PointCourbe*>> c = newCourbeFinale[(int)(round(val*1000))];
                        c.push_back(iso);
                        newCourbeFinale[(int)(round(val*1000))] = c;
                    }
                }
            }
        }

        //Les courbes de gradient et les iso-lignes ne contiennent que les points d'intersection entre elles et avec le bord
        
        //On numérote les points et courbes afin de retrouver facilement les voisins de chaque points
        for(int i=0; i<targets.size(); i++)
        {
            int v = (int)(round(targets[i]*1000));
            vector<vector<PointCourbe*>> c = newCourbeFinale[v];
            for(int j=0; j<c.size(); j++)
            {
                vector<PointCourbe*> iso = c[j];
                for(int k=0; k<iso.size(); k++)
                {
                    PointCourbe* Pc = iso[k];
                    Pc->isoInd = k;
                }
            }
        }
        
        debug<<"Détection des iso-lignes proches"<<endl;
        //On identifie les iso-lignes proches (si leurs valeurs sont proches aussi) et on en marque une (par couple proche) pour suppression
        for(int i=1; i<targets.size(); i++)
        {
            int v = (int)(round(targets[i]*1000));
            double v0 = targets[i-1];
            debug<<"\tPour les iso-lignes de valeurs "<<targets[i]<<" et "<<v0<<endl;
            vector<vector<PointCourbe*>> c = newCourbeFinale[v];
            for(int j=0; j<c.size(); j++)
            {
                vector<PointCourbe*> iso = c[j];
                double distAvg = 0;
                int nbI = 0;
                vector<int> numIso;
                for(int k=0; k<iso.size(); k++)
                {
                    PointCourbe Pc = *(iso[k]);
                    if(Pc.graC!=-1 && Pc.indGra!=0)
                    {
                        vector<PointCourbe*> grad = NewCourbeGrad[Pc.graC];
                        PointCourbe Pc2 = *(grad[Pc.indGra-1]);
                        if(abs(Pc2.val-v0)<0.00001)
                        {
                            nbI++;
                            distAvg+=dist2(Pc.P,Pc2.P);
                            numIso.push_back(Pc2.isoC);
                        }
                        else
                        {
                            if(Pc.indGra+1<grad.size())
                            {
                                Pc2 = *(grad[Pc.indGra+1]);
                                if(abs(Pc2.val-v0)<0.00001)
                                {
                                    nbI++;
                                    distAvg+=dist2(Pc.P,Pc2.P);
                                    numIso.push_back(Pc2.isoC);
                                }
                            }
                            
                        }
                    }
                }
                if(nbI!=0)
                {
                    distAvg=(distAvg/(double)nbI);
                    if(distAvg<normeMax*5)
                    {
                        int nb = iso.size();
                        map<int,int> nbIso;
                        int Nmax = 0;
                        int val = -1;
                        for(int k=0; k<numIso.size(); k++)
                        {
                            int n = numIso[k];
                            if(nbIso.find(n) != nbIso.end())
                            {
                                nbIso[n]=1+ nbIso[n];
                            }
                            else
                            {
                                nbIso[n]=1;
                            }
                            if(nbIso[n]>Nmax)
                            {
                                Nmax =  nbIso[n];
                                val = n;
                            }
                        }
                        vector<vector<PointCourbe*>> c0 = newCourbeFinale[(int)(round(v0*1000))];
                        if(val!=-1)
                        {
                            vector<PointCourbe*> iso0 = c0[val];
                            int nb0 = iso0.size();
                            if(nb0>nb)
                            {
                                debug<<"\t\tSuppression de l'iso-ligne n°"<<j<<" de valeur "<<targets[i]<<" à cause de l'iso-ligne n°"<<val<<" de valeur "<<v0<<endl;
                                for(int k=0; k<iso.size(); k++)
                                {
                                    PointCourbe* Pc = iso[k];
                                    Pc->isoC=-1;
                                }
                            }
                            else
                            {
                                debug<<"\t\tSuppression de l'iso-ligne n°"<<val<<" de valeur "<<v0<<" à cause de l'iso-ligne n°"<<j<<" de valeur "<<targets[i]<<endl;
                                for(int k=0; k<iso0.size(); k++)
                                {
                                    PointCourbe* Pc = iso0[k];
                                    Pc->isoC=-1;
                                }
                            }
                            
                        }
                    }
                }
            }
        }

        //On supprime les iso-lignes qui ont été marquées à l'étape précédente
        for(int i=0; i<targets.size(); i++)
        {
            int v = (int)(round(targets[i]*1000));
            vector<vector<PointCourbe*>> c = newCourbeFinale[v];
            for(auto it = c.begin(); it!=c.end(); it++)
            {
                vector<PointCourbe*> iso = *(it);
                if(iso[0]->isoC==-1)
                {
                    auto it2 = it;
                    it2--;
                    c.erase(it);
                    it=it2;
                }
            }
            newCourbeFinale[v] = c;
        }

        debug<<"Suppression des points des iso-lignes supprimées dans les courbes de gradient"<<endl;
        //On supprime les points des courbes de gradient qui sont sur les iso-lignes supprimées
        for(int i=0; i<NewCourbeGrad.size(); i++)
        {
            vector<PointCourbe*> grad = NewCourbeGrad[i];
            for(auto it = grad.begin(); it!=grad.end(); it++)
            {
                PointCourbe* Pc = *(it);
                if(Pc->val!=-1 && Pc->isoC==-1)
                {
                    auto it2 = it;
                    it2--;
                    grad.erase(it);
                    it=it2;
                }
            }
            NewCourbeGrad[i] = grad;
        }

         cout<<"Construction des bords"<<endl;
         debug<<"Construction des bords"<<endl;
        //On fait la liste des points (d'iso-lignes et de courbes de gradient) sur les bords de la surface
        vector<PointCourbe*> tour0;
        for(int i=0; i<targets.size(); i++)
        {
            int v = (int)(round(targets[i]*1000));
            vector<vector<PointCourbe*>> c = newCourbeFinale[v];
            for(int j = 0; j<c.size(); j++)
            {
                vector<PointCourbe*> iso = c[j];
                for(int k=0; k<iso.size(); k++)
                {
                    PointCourbe* Pc = iso[k];
                    if(Pc->graC==-1)
                        tour0.push_back(Pc);
                }
            }
        }
        for(int i=0; i<NewCourbeGrad.size(); i++)
        {
            vector<PointCourbe*> grad = NewCourbeGrad[i];
            for(int j=0; j<grad.size(); j++)
            {
                PointCourbe* Pc = grad[j];
                if(Pc->val==-1)
                {
                    tour0.push_back(Pc);
                }
            }
        }

        debug<<"Identification des points des bords présentant un angle important avec ses voisins"<<endl;
        //On identifie les points du bords de la surface qu'il faut garder pour maintenir la forme
        for(int i=0; i<surface.size();i++)
        {
            int inP = i-1;
            if(inP<0)
                inP=surface.size()-1;
            int inN = i+1;
            if(inN>=surface.size())
                inN=0;
            Point Pp = surface[inP];
            Point P = surface[i];
            Point Pn = surface[inN];
            Point Vp = addP(P,Pp,-1);
            Point Vn = addP(Pn,P,-1);
            double ang = abs(getAngle(Vp,Vn));
            if(ang>30)
                pointsImportants.push_back(i);
        }
        debug<<"\t"<<pointsImportants.size()<<" points trouvés"<<endl;

        debug<<"Construction de la liste ordonnée des points des courbes de gradient, iso-lignes et ceux de l'étape précédente le long du bord"<<endl;
        //On ordonne les points (d'iso-lignes, de courbes de gradients et ceux identifiés à l'étape précédente) le long du bord 
        vector<PointCourbe*> tour;
        for(int i=0; i<surface.size(); i++)
        {
            int indN = i+1;
            if(indN>=surface.size())
                indN=0;
            Point P = surface[i];
            Point Pn = surface[indN];
            PointCourbe* Important;
            bool impo = false;
            for(int j=0; j<pointsImportants.size() && !impo; j++)
            {
                if(indN==pointsImportants[j])
                {
                    Important = new PointCourbe;
                    Important->P = Pn;
                    Important->val=-1;
                    Important->isoC=-2;
                    Important->graC=-2;
                    impo = true;
                }
            }
            double distSur = dist2(P,Pn);
            double x = P.x-Pn.x;
            double y = P.y-Pn.y;
            int xy = 0;
            int dir = 0;
            if(abs(x)>abs(y))
            {
                xy=1;
                if(x<0)
                    dir = 1;
                else
                    dir = -1;
            }
            else
            {
                xy=-1;
                if(y<0)
                    dir = 1;
                else
                    dir = -1;
            }
            vector<int> added;
            for(int j=0; j<10000; j++)
            {
                double t = j*0.0001;
                Point Pi = interPoint(Pn,P,t);
                for(int k=0; k<tour0.size(); k++)
                {
                    Point Pt = tour0[k]->P;
                    double d = dist2(Pi,Pt);
                    if(d<normeMax)
                    {
                        double distP = distProj(P,Pn,Pt);
                        if(distP>=0)
                        {
                            bool already = false;
                            for(int l=0; l<added.size() && !already; l++)
                            {

                                if(added[l]==k)
                                {
                                    already=true;
                                }
                            }
                            if(!already)
                            {
                                if(abs(distP-distSur)<normeMin && impo)
                                {
                                    impo=false;
                                }
                                bool add = false;
                                for(auto it = added.begin(); it!=added.end() && !add; it++)
                                {
                                    int l = *it;
                                    Point Pa = tour0[l]->P;
                                    double distP2 = distProj(P,Pn,Pa);
                                    if(distP<distP2)
                                    {
                                        add=true;
                                        added.insert(it,k);
                                    }

                                }
                                if(!add)
                                    added.push_back(k);
                            }
                        }
                    }
                }
            }
            for(int j=0; j<added.size(); j++)
            {
                int k = added[j];
                PointCourbe* Pct = tour0[k];
                if(Pct->isoC!=-2 && Pct->graC!=-2)
                {
                    Point Pt = Pct->P;
                    double distP = distProj(P,Pn,Pt);
                    if(distP>distSur && impo)
                    {
                        Important->isoInd = tour.size();
                        Important->indGra = tour.size();
                        tour.push_back(Important);
                        impo = false;
                    }
                    if(Pct->isoC==-1)
                    {
                        Pct->isoC=-2;
                        Pct->isoInd=tour.size();
                    }
                    else
                    {
                        Pct->graC=-2;
                        Pct->indGra=tour.size();
                    }
                    tour.push_back(Pct);
                }
            }
            if(impo)
            {
                Important->isoInd = tour.size();
                Important->indGra = tour.size();
                tour.push_back(Important);
            }
            
        }
        debug<<"\t"<<tour.size()<<" points le long du bord"<<endl;
        
         cout<<"Fin de la phase de préparation"<<endl;
         debug<<"Fin de la phase de préparation"<<endl;
    glPopMatrix();
    /* on force l'affichage du resultat */
    glFlush();
    //On efface tout pour avoir un rendu final propre
    //L'affichage des étapes précédentes est gardé afin de montrer l'évolution
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f );
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glPushMatrix();
    glTranslatef(0,0,0);
    glRotatef(cameraAngleX,1.,0.,0.)    ;
    glRotatef(cameraAngleY,0.,1.,0.);
    map<int,vector<vector<PointCourbe*>>> newCourbeFinale2;
    vector<vector<PointCourbe*>> NewCourbeGrad2;

    //On redéssine la surface avec les couleurs pour la fonction
    for(int i=-n; i<=n; i++)
    {
        for(int j=-n; j<=n; j++)
        {

            float xb = (float)(i)/(float)(n);
            float yb = (float)(j)/(float)(n);
            if(i!=n && j!=n)
            {
                float dn = 1.0f/(float)n;
                float var = func(xb,yb,0);
                glBegin(GL_QUADS);
                glColor3f((var-Vmin)/(Vmax-Vmin),(var-Vmin)/(Vmax-Vmin),(var-Vmin)/(Vmax-Vmin));
                glVertex3f(xb,yb,0);

                var = func(xb+dn,yb,0);
                glColor3f((var-Vmin)/(Vmax-Vmin),(var-Vmin)/(Vmax-Vmin),(var-Vmin)/(Vmax-Vmin));
                glVertex3f(xb+dn,yb,0);

                var = func(xb+dn,yb+dn,0);
                glColor3f((var-Vmin)/(Vmax-Vmin),(var-Vmin)/(Vmax-Vmin),(var-Vmin)/(Vmax-Vmin));
                glVertex3f(xb+dn,yb+dn,0);

                var = func(xb,yb+dn,0);
                glColor3f((var-Vmin)/(Vmax-Vmin),(var-Vmin)/(Vmax-Vmin),(var-Vmin)/(Vmax-Vmin));
                glVertex3f(xb,yb+dn,0);
                glEnd();
            }
        }
    }

    cout<<"Vérification du sens de parcours des iso-lignes et des bords"<<endl;
    debug<<"Vérification du sens de parcours des iso-lignes"<<endl;

    //On fait une première étape de tri des points des iso-lignes
    for(int i=0; i<targets.size(); i++)
    {
        int v = (int)(round(targets[i]*1000));
        vector<vector<PointCourbe*>> c = newCourbeFinale[v];
        for(int j=0; j<c.size(); j++)
        { 
            
            vector<PointCourbe*> iso = c[j];
            vector<PointCourbe*> newIso;
            bool stop = true;
            int x = 0;
            do{
                stop = true;
                newIso.push_back(iso[0]);
                for(int k=0; k<iso.size()-1; k++)
                {
                    PointCourbe* Pc0 = iso[k];
                    PointCourbe* Pc1 = iso[k+1];
                    
                    bool plop = false;
                    while(Pc1->isoC==-1 && k<iso.size()-1)
                    {
                        
                        Pc1 = iso[k+1];
                        k++;
                        plop = true;
                    }
                    if(plop)
                        k--;
                    
                    if(Pc0->isoC!=-1 && Pc1->isoC!=-1)
                    {
                        double d = dist2(Pc0->P,Pc1->P);
                        for(int l = k+2; l<iso.size() && d>normeMax; l++)
                        {
                            PointCourbe* Pc = iso[l];
                            double dc = distProj(Pc0->P,Pc1->P,Pc->P);
                            Point Pi= interPoint(Pc0->P,Pc1->P,dc/d);
                            double da = dist2(Pi,Pc->P);
                            if(dc>0 && dc<d && da<normeMax)
                            {
                                if(Pc->isoC!=-1)
                                {
                                    Pc->isoC=-1;
                                    
                                    newIso.push_back(Pc);
                                    stop = false;
                                }
                                else
                                {
                                    
                                }                              
                            }
                        }
                        
                        newIso.push_back(Pc1);
                    }
                }
                
                
                
                for(int k=0; k<newIso.size(); k++)
                {
                    newIso[k]->isoC = j;
                    newIso[k]->isoInd=k;
                    
                    
                }
                
                
                iso = newIso;
                newIso.clear();
                x++;
            }while(!stop);
            c[j]=iso;
        }
        newCourbeFinale[v]=c;
    }
    
    //On identifie le sens des iso-lignes et on les inverse si nécessaire
    //On fait en sorte que toutes les iso-lignes soient dans un même sens par rapport au gradient local
    for(int i=0; i<targets.size(); i++)
    {
        int v = (int)(round(targets[i]*1000));
        vector<vector<PointCourbe*>> c = newCourbeFinale[v];
        for(int j=0; j<c.size(); j++)
        {
            vector<PointCourbe*> iso = c[j];
            if(iso.size()>2)
            {
                
                PointCourbe* Pc = iso[1];
                Point Po = Pc->P;
                int ind = 1;
                double distMult = -1;
                for(int k=1; k<iso.size()-1; k++)
                {
                    Point P0 = iso[k-1]->P;
                    Point P1 = iso[k]->P;
                    Point P2 = iso[k+1]->P;
                    double d0 = dist2(P0,P1);
                    double d2 = dist2(P0,P2);
                    double da = (d0+d2)/2;
                    if(distMult==-1)
                    {
                        distMult = min(d0,d2);
                        ind = k;
                    }
                    else
                    {
                        double dd = min(d0,d2);
                        if(dd>distMult)
                        {
                            distMult = dd;
                            ind = k;
                        }
                    }
                }
                Po = iso[ind]->P;
                bool fin = false;
                int sens = 0;
                double distMinT = -1;
                int indMinT = -1;
                double distMinT2 = -1;
                int indMinT2 = -1;
                while(!fin)
                {
                    Repere R = getRepere(Po);
                    Point To = R.T;
                    Point Pi = addP(Po,To,distMult*0.01);
                    double diffV = abs(val(Pi.x,Pi.y,Pi.z)-val(Po.x,Po.y,Po.z));
                    if(diffV>0.000001)
                    {
                        for(int k=1; k<250; k++)
                        {
                            double a1 = k*0.001*M_PI;
                            double a2 = -k*0.001*M_PI;
                            Point T1 = {0,0,0};
                            T1.x = To.x * cos(a1) - To.y * sin(a1);
                            T1.y = To.x * sin(a1) + To.y * cos(a1);
                            T1.z = To.z;
                            Point T2 = {0,0,0};
                            T2.x = To.x * cos(a2) - To.y * sin(a2);
                            T2.y = To.x * sin(a2) + To.y * cos(a2);
                            T2.z = To.z;
                            Point Pr1 = addP(Po,T1,distMult*0.01);
                            Point Pr2 = addP(Po,T2,distMult*0.01);
                            double v1 = abs(val(Pr1.x,Pr1.y,Pr1.z)-val(Po.x,Po.y,Po.z));
                            double v2 = abs(val(Pr2.x,Pr2.y,Pr2.z)-val(Po.x,Po.y,Po.z));
                            if(v1<diffV)
                            {
                                diffV = v1;
                                Pi = Pr1;
                            }
                            if(v2<diffV)
                            {
                                diffV = v2;
                                Pi = Pr2;
                            }
                            if(diffV<0.000001)
                                break;
                        }
                    }
                    double distMin = -1;
                    double distMin2 = -1;
                    int indMin = -1;
                    int indMin2 = -1;
                    Po = Pi;
                    for(int k=0; k<iso.size(); k++)
                    {
                        Point P=iso[k]->P;
                        double d = dist2(P,Po);
                        if(distMin==-1 || d<distMin)
                        {
                            distMin2=distMin;
                            indMin2=indMin;
                            distMin=d;
                            indMin=k;
                        }
                        else
                        {
                            if(distMin2==-1 || d<distMin2)
                            {
                                distMin2=d;
                                indMin2=k;
                            }
                        }
                    }
                    if(distMin>distMinT && distMin2>distMinT2 && ((indMin==0 && indMin2==1)|| (indMin==iso.size()-1 && indMin2==iso.size()-2)))
                        fin = true;
                    else
                    {
                        if(indMin == ind)
                        {
                            if(distMin2<distMult*0.05)
                            {
                                if(indMin2<ind)
                                    sens=-1;
                                else
                                    sens=1;
                                fin = true;
                            }
                        }
                        else
                        {
                            if(distMin<distMult*0.05)
                            {
                                if(indMin<ind)
                                    sens=-1;
                                else
                                    sens=1;
                                fin = true;
                            }
                        }
                    }
                    distMinT = distMin;
                    indMinT  = indMin;
                    distMinT2= distMin2;
                    indMinT2 = indMin2;
                    
                }
                
                if(sens==-1)
                {
                    
                    vector<PointCourbe*> iso2;
                    for(int k=iso.size()-1; k>=0; k--)
                    {
                        iso2.push_back(iso[k]);
                    }
                    
                    for(int k=0; k<iso2.size(); k++)
                    {
                        PointCourbe* Pc = iso2[k];
                        
                        Pc->isoInd = k;
                    }
                    
                    c[j]=iso2;
                }
            }
        }
        newCourbeFinale[v] = c;
    }

    //On renumérote les points
    for(int i=0; i<targets.size(); i++)
    {
        int v = (int)(round(targets[i]*1000));
        vector<vector<PointCourbe*>> c = newCourbeFinale[v];
        for(int j=0; j<c.size(); j++)
        {
            vector<PointCourbe*> iso = c[j];
            for(int k=0; k<iso.size(); k++)
            {
                PointCourbe* Pc = iso[k];
                Pc->isoC = j;
                Pc->isoInd = k;
            }
        }
    }

    for(int i=0; i<NewCourbeGrad.size(); i++)
    {
        vector<PointCourbe*> grad = NewCourbeGrad[i];
        for(int j=0; j<grad.size(); j++)
        {
            PointCourbe* Pc = grad[j];
            Pc->graC=i;
            Pc->indGra=j;
        }
    }

    debug<<"Vérification du sens de parcours des points sur les bords"<<endl;
    //On identifie le sens des points du bords (par rapport à l'angle d'arrivé des courbes/iso-lignes sur le bord) et on l'inverse si nécéssaire
    double dMax = -1;
    int indMax = -1;
    for(int i=0; i<tour.size(); i++)
    {
        PointCourbe* Pct = tour[i];
        int indN = i+1;
        if(indN>=tour.size())
            indN=0;
        PointCourbe* Pcn = tour[indN];
        double d = dist2(Pct->P,Pcn->P);
        if(d>dMax)
        {
            dMax = d;
            indMax = i;
        }
    }
    PointCourbe* PctM = tour[indMax];
    int indN = indMax+1;
    if(indN>=tour.size())
        indN=0;
    PointCourbe* PcnM = tour[indN];
    PointCourbe* PcaM;
    if(PctM->isoC!=-2)
    {
        int v = (int)(round(PctM->val*1000));
        vector<PointCourbe*> iso = newCourbeFinale[v][PctM->isoC];
        if(PctM->isoInd==0)
        {
            PcaM = iso[1];
        }
        else
        {
            PcaM = iso[iso.size()-2];
        }
    }
    else
    {
        vector<PointCourbe*> grad = NewCourbeGrad[PctM->graC];
        if(PctM->indGra==0)
        {
            PcaM = grad[1];
        }
        else
        {
            PcaM = grad[grad.size()-2];
        }
    }
    Point PtM = PctM->P;
    Point PnM = PcnM->P;
    Point PaM = PcaM->P;
    Point VnM = addP(PnM,PtM,-1);
    Point VaM = addP(PaM,PtM,-1);
    double angVM = abs(getAngle(VnM,VaM));
    double oM = angVM * M_PI/180;
    Point Vrot = {0,0,0};
    Vrot.x = VaM.x * cos(oM) - VaM.y * sin(oM);
    Vrot.y = VaM.x * sin(oM) + VaM.y * cos(oM);
    Vrot.z = VaM.z;
    double newAngVM = abs(getAngle(VnM,Vrot));
    bool tourInverse = false;
    if(newAngVM<angVM)
    {
        tourInverse = true;
    }
    if(tourInverse)
    {
        vector<PointCourbe*> tourI;
        for(int i=0 ; i<tour.size(); i++)
        {
            tourI.insert(tourI.begin(),tour[i]);
        }
        tour.clear();
        tour = tourI;
        for(int i=0 ; i<tour.size(); i++)
        {
            PointCourbe* Pc = tour[i];
            if(Pc->isoC==-2)
            {
                Pc->isoInd = i;
                if(Pc->graC==-2)
                    Pc->indGra = i;
            }
            else
            {
                Pc->indGra = i;
            }
        }
    }


     cout<<"Gestion des courbes de gradients se rejoignant"<<endl;
     debug<<"Gestion des courbes de gradients se rejoignant"<<endl;
    //On a des courbes de gradients très proches les unes des autres sur certaines iso-lignes (plus les iso-lignes sont proches d'un extremum plus les courbes de gradient se rejoignent)
    //Cela cause de très petite cellules
    //C'est pour cette raison qu'on cherche à en supprimer afin de créer de plus grosses cellules
    //Cependant, lors de ces croisements et dû à la nature des approximations des courbes, celles-ci peuvent se croiser ou simplement être trop proche pour identifié l'ordre des intersections
    //Pour cela, on identifie tous les points qui font partie de ces groupements
    //Ensuite, pour chaque iso-ligne (et le bord), on identifie l'ordre des points qui ne font pas partie de ces groupements
    //Cela nous donne l'ordre des courbes de gradients par groupe. En utilisant ces ordres et les courbes de gradients des points groupés, on peut trouver leur ordre dans le groupe
    vector<vector<PointCourbe*>> allProche;
    vector<vector<PointCourbe*>> ordreGra;
    int sta = -1;
    bool isLink = false;
    vector<PointCourbe*> proche;
    bool fusion = false;
    vector<PointCourbe*> ordreP;
    vector<PointCourbe*> ordreP2;
    Point Prev;
    debug<<"Identification des groupes de points le long des bords"<<endl;
    //On identifie les groupes le long du bord
    for(int k=0; k<tour.size(); k++)
    {
        PointCourbe* Pc0 = tour[k];
        if(Pc0->isoC==-2)
        {
            if(proche.empty())
            {
                proche.push_back(Pc0);
                sta = k;
                fusion = true;
            }
            else
            {
                double distMax = dist2(Pc0->P,proche[0]->P);
                double distMin = dist2(Pc0->P,proche[0]->P);
                for(int l=1; l<proche.size(); l++)
                {
                    double d = dist2(Pc0->P,proche[l]->P);
                    if(d>distMax)
                        distMax = d;
                    if(d<distMin)
                        distMin = d;
                }
                if(distMax<=pasIso)
                {
                    proche.push_back(Pc0);
                }
                else
                {
                    if(proche.size()==1)
                    {
                        proche[0]=Pc0;
                        sta = k;
                        fusion = false;
                    }
                    else
                    {
                        
                        if(fusion && distMin<pasIso*2)
                        {
                            for(int m=0; m<proche.size(); m++)
                            {
                                allProche.back().push_back(proche[m]);
                            }

                        }
                        else
                        {
                            allProche.push_back(proche);
                        }
                        proche.clear();
                        proche.push_back(Pc0);
                        sta = k;
                        fusion = true;
                        
                    }
                }
            }
        }
    }
    if(proche.size()>1)
        allProche.push_back(proche);

    debug<<"Identification des points le long des bords dont l'ordre est sûr"<<endl;
    //On fait l'ordre des points non groupés le long du bord
    for(int k=0; k<tour.size(); k++)
    {
        PointCourbe* Pc0 = tour[k];
        if(Pc0->isoC==-2 && Pc0->graC!=-2)
        {
            bool ok = true;
            int ind = k;
            PointCourbe* Pc1 = tour[k];
            do{
                ind--;
                if(ind<0)
                    ind = tour.size()-1;
                Pc1 = tour[ind];
            }while(Pc1->graC==-2);
            if(dist2(Pc0->P,Pc1->P)<pasIso)
            {
                ok = false;
            }
            ind = k;
            PointCourbe* Pc2 = tour[k];
            do{
                ind++;
                if(ind>=tour.size())
                    ind = 0;
                Pc2 = tour[ind];
            }while(Pc2->graC==-2);
            
            if(dist2(Pc0->P,Pc2->P)<pasIso)
            {
                ok = false;
            }
            if(ok)
            {
                if(Pc0->indGra==0)
                {
                    if(ordreP.size()>0)
                    {
                        auto it = ordreP.end();
                        do{
                            it--;
                        }while((*it)->indGra!=0 && (*it)->graC!=Pc0->graC && it!=ordreP.begin());
                        if((*it)->graC==Pc0->graC || (*it)->indGra==0)
                        {
                            it++;
                            ordreP.insert(it,Pc0);
                        }
                        else
                        {
                            ordreP.insert(it,Pc0);
                        }
                    }
                    else
                    {
                        ordreP.push_back(Pc0);
                    }
                }
                else
                {
                    if(ordreP.size()>0)
                    {
                        auto it = ordreP.end();
                        do{
                            it--;
                        }while((*it)->indGra!=0 && it!=ordreP.begin());
                        if((*it)->indGra!=0)
                            ordreP.insert(it,Pc0);
                        else
                        {
                            it++;
                            ordreP.insert(it,Pc0);
                        }
                    }
                    else
                    {
                        ordreP.push_back(Pc0);
                    }
                    
                }
            }
            else
            {
                if(ordreP.size()>1)
                    ordreGra.push_back(ordreP);
                ordreP.clear();
            }
        }
    }
    if(ordreP.size()>1)
        ordreGra.push_back(ordreP);
    ordreP.clear();
    
    proche.clear();

    debug<<"Identification des groupes de points sur les iso-lignes"<<endl;
    //On identifie les groupes le long des iso-lignes
    for(int i=0; i<targets.size(); i++)
    {
        int v = (int)(round(targets[i]*1000));
        vector<vector<PointCourbe*>> c = newCourbeFinale[v];
        for(int j=0; j<c.size(); j++)
        { 
            
            vector<PointCourbe*> iso = c[j];
            fusion = false;
            int sta = -1;
            vector<vector<PointCourbe*>> pr;
            for(int k=0; k<iso.size(); k++)
            {
                PointCourbe* Pc0 = iso[k];
                if(Pc0->graC!=-2)
                {
                    int already = -1;
                    for(int l=0; l<pr.size() && already==-1; l++)
                    {
                        for(int m=0; m<pr[l].size() && already==-1; m++)
                        {
                            PointCourbe* Pcc = pr[l][m];
                            if(equal(*Pc0,*Pcc))
                                already = l;
                        }
                    }
                    for(int l=0; l<iso.size(); l++)
                    {
                        if(k!=l)
                        {
                            PointCourbe* Pc1 = iso[l];
                            if(dist2(Pc0->P,Pc1->P)<pasIso && Pc1->graC!=-2)
                            {
                                bool ok = true;
                                for(int l1=0; l1<pr.size() && ok; l1++)
                                {
                                    for(int m=0; m<pr[l1].size() && ok; m++)
                                    {
                                        PointCourbe* Pcc = pr[l1][m];
                                        if(equal(*Pc1,*Pcc))
                                            ok=false;;
                                    }
                                }
                                if(ok)
                                {
                                    if(already==-1)
                                    {
                                        already = pr.size();
                                        vector<PointCourbe*> pro;
                                        pro.push_back(Pc0);
                                        pro.push_back(Pc1);
                                        pr.push_back(pro);
                                    }
                                    else
                                    {
                                        pr[already].push_back(Pc1);
                                    }
                                }
                            }
                        }
                    }
                }

            }
            for(int k=0; k<pr.size(); k++)
            {
                allProche.push_back(pr[k]);
            }
        }
    }
    
    debug<<"Identification des pointssur les iso-lignes dont l'ordre est sûr"<<endl;
    //On fait l'ordre des points non groupés le long des iso-lignes
    for(int i=0; i<targets.size(); i++)
    {
        int v = (int)(round(targets[i]*1000));
        vector<vector<PointCourbe*>> c = newCourbeFinale[v];
        for(int j=0; j<c.size(); j++)
        {
            
            vector<PointCourbe*> iso = c[j];
            int sta = -1;
            for(int k=0; k<iso.size(); k++)
            {
                PointCourbe* Pc0 = iso[k];
                bool ok = true;
                for(int l=0; l<iso.size() && ok; l++)
                {
                    if(l!=k)
                    {
                        PointCourbe* Pc1 = iso[l];
                        ok = dist2(Pc0->P,Pc1->P)>pasIso;
                    }
                }
                if(ok)
                {
                    if(Pc0->graC!=-2)
                    {
                        ordreP.push_back(Pc0);
                    }
                }
                else
                {
                    if(ordreP.size()>1)
                        ordreGra.push_back(ordreP);
                    ordreP.clear();
                }
            }
            if(ordreP.size()>1)
                ordreGra.push_back(ordreP);
            ordreP.clear();

            
        }
    }

    //On supprime les doublons dans les ordres (à cause des ordres le long du bord)
    //Les courbes de gradient peuvent commencer et finir le long d'un bord, donc apparaître deux fois dans les ordres
    //La façon dont elles sont traitées lors de la formation des ordres fait que les doublons se suivent
    for(int i=0; i<ordreGra.size(); i++)
    {
        vector<PointCourbe*> o = ordreGra[i];
        int prev = -1;
        for(auto it = o.begin(); it!=o.end(); it++)
        {
            PointCourbe* Pc = *it;
            
            if(Pc->graC == prev)
            {
                auto it2 = it;
                it2--;
                o.erase(it);
                it=it2;
            }
            else
            {
                prev = Pc->graC;
            }
        }
        ordreGra[i] = o;
    }


    //On identifie si les ordres le long du bord ou des iso-lignes formant une boucle peuvent être fusionner (si un ordre commence un début de la liste de points et un autre finit à la fin de la liste)
    int prem = -1;
    int dern = 0;
    for(int i=0; i<tour.size(); i++)
    {
        PointCourbe* Pc = tour[i];
        if(Pc->isoC==-2 && Pc->graC!=-2)
        {
            if(prem==-1)
                prem=i;
            dern = i;
        }
    }

    vector<vector<PointCourbe*>> ordreGra2;
    vector<int> fused;
    for(int i=0; i<ordreGra.size(); i++)
    {
        vector<PointCourbe*> o = ordreGra[i];
        PointCourbe* Pc0 = o[0];
        for(int j=i+1; j<ordreGra.size(); j++)
        {
            vector<PointCourbe*> o1 = ordreGra[j];
            PointCourbe* Pc1 = o1[o1.size()-1];
            if(Pc1->isoC==Pc0->isoC && (Pc0->isoC==-2 || abs(Pc0->val-Pc1->val)<0.0001))
            {
                if(Pc0->isoC==-2)
                {
                    bool ind0 = false;
                    bool indF = false;
                    for(int k=0; k<o.size() && !ind0; k++)
                    {
                        PointCourbe* Pc = o[k];
                        if(Pc->isoInd==prem)
                            ind0 = true;
                    }
                    for(int k=0; k<o1.size() && !indF && ind0; k++)
                    {
                        PointCourbe* Pc = o1[k];
                        if(Pc->isoInd==dern)
                            indF = true;
                    }
                    if(ind0 && indF)
                    {
                        
                        fused.push_back(i);
                        fused.push_back(j);
                        if(tour[dern]->indGra!=0)
                        {
                            int numP = 0;
                            for(int k=0; k<o.size(); k++)
                            {
                                PointCourbe* Pc = o[k];
                                if(Pc->indGra!=0)
                                {
                                    auto it = o1.end();
                                    do{
                                        it--;
                                    }while((*it)->indGra!=0 && it!=o1.begin());
                                    if((*it)->indGra!=0)
                                        o1.insert(it,Pc);
                                    else
                                    {
                                        it++;
                                        for(int l=0; l<numP; l++)
                                            it++;
                                        o1.insert(it,Pc);
                                    }
                                    numP++;
                                }
                                else
                                {
                                    numP=0;
                                    bool doublon = false;
                                    for(int l=0; l<o1.size() && !doublon; l++)
                                    {
                                        if(o1[l]->graC==Pc->graC)
                                            doublon = true;
                                    }
                                    if(!doublon)
                                        o1.push_back(Pc);
                                }
                            }
                            ordreGra2.push_back(o1);
                        }
                        else
                        {
                            o1.insert(o1.end(),o.begin(),o.end());
                            ordreGra2.push_back(o1);
                        }
                    }
                }
                else
                {
                    int v = (int)(round(Pc0->val*1000));
                    vector<PointCourbe*> iso = newCourbeFinale[v][Pc0->isoC];
                    if(Pc0->isoInd==0 && Pc1->isoInd==iso.size()-1)
                    {
                        
                        o1.insert(o1.end(),o.begin(),o.end());
                        ordreGra2.push_back(o1);
                        fused.push_back(i);
                        fused.push_back(j);
                    }
                }
            }
        }
        bool fuse = false;
        for(int j=0; j<fused.size(); j++)
        {
            if(i==fused[j])
                fuse=true;
        }
        if(!fuse)
            ordreGra2.push_back(o);
    }
    ordreGra = ordreGra2;

    debug<<"Identification de "<<allProche.size()<<" groupes de points et "<<ordreGra.size()<<" suite de points dont l'ordre est sûr"<<endl;

    //On ordonne les ordres par taille décroissante (afin de commencer le traitement des groupes avec la plus grande quantité d'information)
    vector<vector<PointCourbe*>> newOrdre;
    for(int i=0; i<ordreGra.size(); i++)
    {
        vector<PointCourbe*> o = ordreGra[i];
        auto it = newOrdre.begin();
        bool placed = false;
        while(it!=newOrdre.end() && !placed)
        {
            if(o.size()>it->size())
            {
                placed = true;
                newOrdre.insert(it,o);
            }
            it++;
        }
        if(!placed)
            newOrdre.push_back(o);
    } 

    debug<<"Traitement des groupes pour réordonner les points et supprimer des courbes de gradient"<<endl;
    //On traite les groupes. On ordonne les points de chaque groupe et on identifie les extrémités de ces groupes (premier et dernier points du groupe)
    //Ces extrémités forment les bords de la cellule fusionnée
    int minIsoC = -3;
    for(int i=0; i<allProche.size(); i++)
    {
        vector<PointCourbe*> pr = allProche[i];
        PointCourbe* Pc = pr[0];
        debug<<"\tGroupe n°"<<i<<" de taille "<<pr.size();

        int iP = Pc->isoInd;
        int iN = Pc->isoInd;
        for(int j=1; j<pr.size(); j++)
        {
            PointCourbe* Pcp = pr[j];
            if(Pcp->isoInd<iP)
                iP=Pcp->isoInd;
            if(Pcp->isoInd>iN)
                iN=Pcp->isoInd;
        }
        PointCourbe* PcP;
        PointCourbe* PcN;
        if(Pc->isoC!=-2)
        {
            debug<<" sur l'iso-ligne n°"<<Pc->isoC<<" de valeur "<<Pc->val<<endl;
            int v = (int)(round(Pc->val*1000));
            
            vector<PointCourbe*> iso = newCourbeFinale[v][Pc->isoC];
            
            PcP = iso[iP];
            if(iP==0)
            {
                if(PcP->graC!=-2)
                {
                    PcP = iso.back();
                }
            }
            else
            {
                PcP = iso[iP-1];
            }
            PcN = iso[iN];
            if(iN==iso.size()-1)
            {
                if(PcN->graC!=-2)
                {
                    PcN = iso.front();
                }
            }
            else
            {
                PcN = iso[iN+1];
            }
        }
        else
        {
            debug<<" le long du bord"<<endl;
            if(iP==0)
                PcP=tour.back();
            else
            
            if(iN==tour.size()-1)
                PcN = tour.front();
            else
                    PcP = tour[iP-1];
            PcN = tour[iN+1];
        }
        
        
        
        vector<vector<int>> ordres;
        bool fin = false;
        for(int j=0; j<newOrdre.size() && !fin; j++)
        {
            
            vector<int> ordre;
            bool stop = false;
            bool last = false;
            vector<PointCourbe*> o = newOrdre[j];
            for(int k=0; k<o.size() && !stop; k++)
            {
                
                if(o[k]->graC==PcP->graC)
                {
                    if(ordre.size()>0)
                        stop=true;
                    else
                        ordre.push_back(PcP->graC);
                }

                for(int l=0; l<pr.size() && !stop; l++)
                {
                    if(o[k]->graC==pr[l]->graC)
                    {
                        if(last)
                            stop=true;
                        else
                            ordre.push_back(pr[l]->graC);
                    }
                }

                if(o[k]->graC==PcN->graC)
                {
                    last=true;
                    ordre.push_back(PcN->graC);
                }

            }
            if(ordre.size()==pr.size()+2)
                fin = true;
            if(ordre.size()>1)
            {
                
                ordres.push_back(ordre);
                ordre.clear();
            }
        }
        
        vector<int> ordreFinal;
        int numIt=1;
        fin = false;
        
        int nb0 = 0;
        while(!fin)
        {
            numIt=0;
            
            for(int j=0; j<ordres.size(); j++)
            {
                
                vector<int> o = ordres[j];
                if(ordreFinal.size()==0)
                {
                    for(int k=0; k<o.size(); k++)
                    {
                        ordreFinal.push_back(o[k]);
                        numIt++;
                    }
                }
                else
                {
                    for(int k=0; k<o.size(); k++)
                    {
                        int num = o[k];
                        bool exist = false;
                        for(int l=0; l<ordreFinal.size() && !exist; l++)
                        {
                            if(ordreFinal[l]==num)
                                exist=true;
                        }
                        if(!exist)
                        {
                            numIt++;
                            if(num==PcP->graC)
                            {
                                ordreFinal.insert(ordreFinal.begin(),num);
                            }
                            else
                            {
                                if(num==PcN->graC)
                                {
                                    ordreFinal.push_back(num);
                                }
                                else
                                {
                                    int numP = -1;
                                    int numN = -1;
                                    if(k>0)
                                    {
                                        numP = o[k-1];
                                    }
                                    if(k<o.size()-1)
                                    {
                                        numN = o[k+1];
                                    }
                                    for(auto it=ordreFinal.begin(); it!=ordreFinal.end() && !exist; it++)
                                    {
                                        int no = *it;
                                        if(no==numN)
                                        {
                                            exist=true;
                                            ordreFinal.insert(it,num);
                                        }
                                        if(no==numP)
                                        {
                                            it++;
                                            exist=true;
                                            ordreFinal.insert(it,num);
                                        }
                                    }
                                }
                            }
                            
                        }
                    }
                }
                
            }
            
            if(ordreFinal.size()>0)
            {
                if(numIt==0)
                    nb0++;
                else
                    nb0=0;
                if((ordreFinal.back()==PcN->graC || nb0==2) && numIt==0)
                    fin = true;
            }
            
        }
        vector<PointCourbe*> newPr;
        int st = pr[0]->isoInd;
        for(int k=1; k<pr.size(); k++)
        {
            if(pr[k]->isoInd<st)
                st=pr[k]->isoInd;
        }
        for(int j=0; j<ordreFinal.size(); j++)
        {
            for(int k=0; k<pr.size(); k++)
            {
                if(pr[k]->graC==ordreFinal[j])
                    newPr.push_back(pr[k]);
            }
        }
        
        for(int j=0; j<newPr.size(); j++)
        {
            PointCourbe* nPc = newPr[j];
            if(j!=0 && j!=newPr.size()-1)
            {
                nPc->isoC=(-4-nPc->isoC);
                if(nPc->isoC<minIsoC)
                    minIsoC = nPc->isoC;
            }
            int nIn = st+j;
            
            nPc->isoInd=nIn;
        }
        for(int k=0; k<pr.size(); k++)
        {
            PointCourbe* Pc = pr[k];
            if(Pc->isoC>-4)
            {
                if(Pc->isoInd != st && Pc->isoInd != st+newPr.size()-1)
                    Pc->isoC = -3;
            }
        }
        
        allProche[i]=newPr;
    }

    //On réordonne les points des iso-lignes en fonction de l'ordre trouvé à l'étape précédente
    for(int i=0; i<targets.size(); i++)
    {
        int v = (int)(round(targets[i]*1000));
        vector<vector<PointCourbe*>> c = newCourbeFinale[v];
        for(int j=0; j<c.size(); j++)
        {
            vector<PointCourbe*> iso = c[j];
            vector<PointCourbe*> isoN;
            for(int k=0; k<iso.size(); k++)
            {
                for(int l=0; l<iso.size(); l++)
                {
                    PointCourbe* Pc = iso[l];
                    if(Pc->isoInd==k)
                    {
                        isoN.push_back(Pc);
                    }
                }
            }
            if(isoN.size()==iso.size())
            {
                c[j]=isoN;
            }
        }
        newCourbeFinale[v] = c;
    }

    //On continue le traitement de groupes afin de couper les cellules fusionnés si celles-ci sont maintenant trop grandes
    //Pour cela, on réajoute un maximum de points des mêmes courbes de gradient de façon à éviter des morceaux isolés de courbes de gradient
    /*
    Exemple : La même courbe de gradient a des points qui forment un groupe sur 4 iso-lignes succéssives.
    Si on réajoute les points sans faire attention aux courbes de gradient, 
    alors on pourrait ajouter un morceau de cette courbe de gradient seulement entre les 2e et 3e iso-lignes ainsi qu'avant la 1e et après la 4e
    créant trois morceaux séparés de la même courbe de gradient.
    Cela aurait plusieurs effets, casser l'alignement des cellules et compliquer la formation des cellules si on ne considère pas les morceaux comme des courbes de gradient différentes
    */

    //On ordonne les groupes de points proches par ordre croissant de taille
    /*
    Parce qu'en ajoutant des courbes de gradient en commençant par les plus grand groupes
    On court le risque qu'aucune de celles ajoutées ne passent dans les petites groupes
    Et que lors du traitement des petites groupes on rajoute une autre courbe de gradient qui passerait dans un grand groupe
    Faisant une séparation de plus que prévue dans la cellule fusionné à partir du grand groupe
    */
    vector<vector<PointCourbe*>> allP;
    for(int i=0; i<allProche.size(); i++)
    {
        vector<PointCourbe*> p = allProche[i];
        auto it = allP.begin();
        bool placed = false;
        while(it!=allP.end() && !placed)
        {
            if(p.size()<it->size())
            {
                placed = true;
                allP.insert(it,p);
            }
            it++;
        }
        if(!placed)
            allP.push_back(p);
    }

    debug<<"Choix des courbes de gradient à garder"<<endl;
    //On identifie les courbes de gradient à réajouter en fonction de la longueur d'un bord de la cellule fusionnée le long d'une autre iso-lignes (ou le bord)
    double longLim = 0.4;
    for(int i=0; i<allP.size(); i++)
    {
        vector<PointCourbe*> p = allP[i];
        vector<int> graOk;
        for(int j=1; j<p.size()-1; j++)
        {
            PointCourbe* Pc = p[j];
            if(Pc->isoC>=-2 || Pc->isoC<minIsoC)
            {
                
                graOk.push_back(Pc->graC);
            }
        }
        PointCourbe* Pc0 = p[0];
        
        PointCourbe* Pc1 = p.back();
        vector<PointCourbe*> gra0 = NewCourbeGrad[Pc0->graC];
        vector<PointCourbe*> gra1 = NewCourbeGrad[Pc1->graC];
        vector<int> inGra0;
        vector<int> inGra1;
        for(int j=0; j<gra0.size(); j++)
        {
            if(j!=Pc0->indGra)
            {
                PointCourbe* Pg0 =gra0[j];
                for(int k=0; k<gra1.size(); k++)
                {
                    if(k!=Pc1->indGra)
                    {
                        PointCourbe* Pg1 =gra1[k];
                        if(abs(Pg0->val-Pg1->val)<0.0001 && Pg0->isoC==Pg1->isoC && Pg0->isoC>=-2)
                        {
                            inGra0.push_back(j);
                            inGra1.push_back(k);
                        }
                    }
                }
            }
        }
        if(inGra0.size()>0)
        {
            
            int iG0 = inGra0[0];
            int iG1 = inGra1[0];
            int dM = abs(Pc0->indGra-iG0)+abs(Pc1->indGra-iG1);
            for(int j=1; j<inGra0.size(); j++)
            {
                int i0 = inGra0[j];
                int i1 = inGra1[j];
                int dif = abs(Pc0->indGra-i0)+abs(Pc1->indGra-i1);
                if(dif<dM)
                {
                    iG0=i0;
                    iG1=i1;
                    dM=dif;
                }
            }
            
            PointCourbe* Pg0 = gra0[iG0];
            int num0 = Pg0->isoC;
            if(num0<=-4)
            {
                if(num0>=minIsoC)
                    num0+=4;
                else
                    num0+=(minIsoC+1);
                num0*=-1;
            }
            PointCourbe* Pg1 = gra1[iG1];
            
            vector<PointCourbe*> c;
            double l = 0;
            map<int,double> longu;
            if(num0!=-2)
            {
                int v = (int)(round(Pg0->val*1000));
                c = newCourbeFinale[v][num0];
            }
            else
                c = tour;
            int ind = Pg0->isoInd;
            int nb0 = 0;
            int nb1 = 0;
            if((c[0]->graC!=-2 && c.back()->graC!=-2) || num0==-2)
            {
                for(int j=ind; j!=Pg1->isoInd; j++)
                {
                    if(j==c.size()-1)
                        j=-1;
                    nb0++;
                }
                for(int j=ind; j!=Pg1->isoInd; j--)
                {
                    if(j==0)
                        j=c.size();
                    nb1++;
                }
            }
            bool inv = nb1<nb0;
            vector<int> indOk;
            longu[ind]=0;
            while(ind!=Pg1->isoInd)
            {
                PointCourbe* Pi0 = c[ind];
                PointCourbe* Pi1;
                do{
                    if(inv)
                    {
                        ind--;
                        if(ind==-1)
                            ind=c.size()-1;
                        Pi1 = c[ind];
                    }
                    else
                    {
                        ind++;
                        if(ind==c.size())
                            ind=0;
                        Pi1 = c[ind];
                    }
                }while(Pi1->graC==-2);
                
                for(int j=0; j<graOk.size(); j++)
                {
                    if(Pi1->graC==graOk[j])
                        indOk.push_back(Pi1->isoInd);
                }
                double d = dist2(Pi0->P,Pi1->P);
                l+=d;
                longu[Pi1->isoInd] = l;
            }
            int nbR = l/longLim;
            if(abs(nbR*longLim-l)<0.00001)
                nbR--;
            
            vector<int> newCoupe;
            double lC = l/(double)(nbR+1);
            indOk.push_back(Pg1->isoInd);
            
            while(nbR>=indOk.size())
            {
                double prevC = 0;
                double maxC = -1;
                int maxInd = -1;
                for(int j=0; j<indOk.size(); j++)
                {
                    double coupe = longu[indOk[j]];
                    double longCoupe = coupe-prevC;
                    
                    if(longCoupe>maxC)
                    {
                        maxC = longCoupe;
                        maxInd = j;
                    }
                    prevC = coupe;
                }
                
                int nbC = maxC/lC;
                if(abs(nbC*lC-maxC)<0.00001)
                    nbC--;
                double lNC = maxC/(nbC+1);
                maxC = lNC;
                int indStart = Pg0->isoInd;
                int indFin = indOk[maxInd];
                if(maxInd!=0)
                {
                    indStart = indOk[indMax-1];
                    maxC+=longu[maxInd];
                }
                
                int nbCE = 0;
                while(nbCE<nbC && indStart!=indFin)
                {
                    if(inv)
                    {
                        indStart--;
                        if(indStart==-1)
                            indStart=c.size()-1;
                    }
                    else
                    {
                        indStart++;
                        if(indStart==c.size())
                            indStart=0;
                    }
                    double lP = longu[indStart];
                    if(lP>=maxC)
                    {
                        
                        newCoupe.push_back(indStart);
                        for(auto it=indOk.begin(); it!=indOk.end(); it++)
                        {
                            if(((*it)>indStart && !inv) || ((*it)<indStart && inv))
                            {
                                indOk.insert(it,indStart);
                                break;
                            }
                        }
                        nbCE++;
                        maxC+=lNC;
                    }
                }
                
            }
            for(int j=0; j<newCoupe.size(); j++)
            {
                int indC = newCoupe[j];
                PointCourbe* Pc = c[indC];
                int iG = Pc->graC;
                debug<<"\tCourbe de gradient n°"<<iG<<" gardée pour un groupe sur l'iso-ligne n°"<<Pc0->isoC<<" de valeur "<<Pc0->val<<endl;
                vector<PointCourbe*> mGra = NewCourbeGrad[iG];
                for(int k=0; k<mGra.size(); k++)
                {
                    PointCourbe* PgM = mGra[k];
                    if(PgM->isoC<=-4)
                    {
                        int numG = (PgM->isoC+4)*-1; 
                        PgM->isoC= numG;
                    }
                }
            }
        }
    }

    cout<<"Phase de suppression de courbes de gradient et des points affectés"<<endl;
    debug<<"Phase de suppression de courbes de gradient et des points affectés"<<endl;
    //On supprime les points (des groupes) qui n'ont pas été réajoutés
    //On supprime aussi les iso-lignes ou courbes de gradient devenues trop courtes (en nombre de points)
    //On marque les points de ces iso-lignes et courbes de gradient pour suppression
    //On recommence tant qu'il y a des suppressions prévues
    int nx =0;
    bool supT = false;
    do{
        supT = false;
        for(int i=0; i<targets.size(); i++)
        {
            int v = (int)(round(targets[i]*1000));
            vector<vector<PointCourbe*>> c = newCourbeFinale[v];
            for(int j=0; j<c.size(); j++)
            {
                vector<PointCourbe*> iso = c[j];
                
                for(auto it=iso.begin(); it!=iso.end(); it++)
                {
                    PointCourbe* Pc = *it;
                    if(Pc->isoC<-2)
                    {
                        
                        
                        auto it2 = it;
                        it2--;
                        iso.erase(it);
                        it=it2;
                    }
                }
                
                for(int k=0; k<iso.size(); k++)
                {
                    iso[k]->isoInd=k;
                }
                c[j]=iso;
            }
            newCourbeFinale[v]=c;
        }
        for(int i=0; i<NewCourbeGrad.size(); i++)
        {
            vector<PointCourbe*> gra = NewCourbeGrad[i];
            
            for(auto it=gra.begin(); it!=gra.end(); it++)
            {
                PointCourbe* Pc = *it;
                if(Pc->isoC<-2)
                {
                    
                    
                    auto it2 = it;
                    it2--;
                    gra.erase(it);
                    it=it2;
                }
            }
            
            for(int j=0; j<gra.size(); j++)
            {
                gra[j]->indGra=j;
            }
            NewCourbeGrad[i] = gra;
        }
        bool finSup = false;
        vector<PointCourbe*> iso = newCourbeFinale[77][0];
        while(!finSup)
        {
            bool sup1 = false;
            for(auto it=tour.begin(); it!=tour.end() && !sup1; it++)
            {
                PointCourbe* Pc = *it;
                if(Pc->isoC<-2)
                {
                    
                    
                    sup1 = true;
                    tour.erase(it);
                    it = tour.begin();
                }
            }
            if(!sup1)
                finSup=true;
        }
        
        for(int i=0; i<tour.size(); i++)
        {
            PointCourbe* Pc = tour[i];
            if(Pc->isoC==-2)
                Pc->isoInd=i;
            if(Pc->graC==-2)
                Pc->indGra=i;
        }
        
        for(int i=0; i<targets.size(); i++)
        {
            int v = (int)(round(targets[i]*1000));
            vector<vector<PointCourbe*>> c = newCourbeFinale[v];
            bool fin = false;
            while(!fin)
            {
                bool sup1 = false;
                for(auto it=c.begin(); it!=c.end() && !sup1; it++)
                {
                    if((*it).size()<=2)
                    {
                        sup1=true;
                        supT = true;
                        vector<PointCourbe*> iso = *it;
                        for(int k=0; k<iso.size(); k++)
                        {
                            
                            iso[k]->isoC=-3;
                        }
                        c.erase(it);
                        it=c.begin();
                    }
                }
                if(!sup1)
                    fin =true;
            }
            for(int j=0; j<c.size(); j++)
            {
                vector<PointCourbe*> iso = c[j];
                for(int k=0; k<iso.size(); k++)
                {
                    iso[k]->isoC=j;
                }
            }
            newCourbeFinale[v]=c;
        }
        
        finSup = false;
        while(!finSup)
        {
            bool sup1 = false;
            for(auto it=NewCourbeGrad.begin(); it!=NewCourbeGrad.end() && !sup1; it++)
            {
                if((*it).size()<=1)
                {
                    sup1=true;
                    supT=true;
                    vector<PointCourbe*> gra = *it;
                    for(int k=0; k<gra.size(); k++)
                    {
                        
                        gra[k]->isoC=-3;
                    }
                    NewCourbeGrad.erase(it);
                    it=NewCourbeGrad.begin();
                }
            }
            if(!sup1)
                finSup=true;
        }
        for(int i=0; i<NewCourbeGrad.size(); i++)
        {
            vector<PointCourbe*> gra = NewCourbeGrad[i];
            for(int j=0; j<gra.size(); j++)
            {
                gra[j]->graC=i;
            }
        }
        glEnd();
        nx++;
    }while(supT);

    
    vector<vector<PointCourbe*>> allFaces;

    //On forment les listes de sommets des cellules

    cout<<"Construction des cellules"<<endl;
    debug<<"Construction des cellules"<<endl;
    //On forment celle le long du bord qui ne seront pas créer par les étapes d'après
    bool skip = false;
    for(int i=0; i<tour.size(); i++)
    {
        
        PointCourbe* Po = tour[i];
        bool start = !skip;
        skip = false;
        if(Po->graC==-2)
        {
            if(Po->isoC==-2)
            {
                start=false;
            }
            else
            {
                if(Po->isoInd!=0)
                {
                    vector<PointCourbe*> isoT = newCourbeFinale[(int)(round(Po->val*1000))][Po->isoC];
                    int indT = Po->isoInd-1;
                    PointCourbe* Pi = isoT[indT];
                    while(Pi->indGra==0 && Pi->isoInd!=0)
                    {
                        indT--;
                        Pi = isoT[indT];
                    }
                    if(Pi->isoInd!=0)
                    {
                        vector<PointCourbe*> graT = NewCourbeGrad[Pi->graC];
                        PointCourbe* Pg = graT[Pi->indGra-1];
                        if(Pg->isoC!=-2)
                        {
                            start = false;
                        }
                        else
                            skip = true;
                    }
                }
                else
                {
                    start = false;
                }
            }
        }
        else
        {
            if(Po->indGra!=0)
            {
                start =false;
            }
        }
        if(start)
        {
            vector<PointCourbe*> face;
            PointCourbe* Pa = Po;
            int dir = 5;
            int dGra = 1;
            
            do{
                bool modif = false;
                face.push_back(Pa);
                bool continu = false;
                if(Pa->graC!=-2 && (dir==2 || dir==4 || dir==5))
                {
                    vector<PointCourbe*> gra = NewCourbeGrad[Pa->graC];
                    if((dir==2 && Pa->indGra==0) || (dir==4 && Pa->indGra==gra.size()-1))
                    {
                        continu = true;
                    }
                    else
                    {
                        if(dir==2 || (dir==5 && Pa->indGra!=0))
                        {
                            modif=true;
                            dir=3;
                            Pa = gra[Pa->indGra-1];
                        }
                        else
                        {
                            modif=true;
                            dir=1;
                            Pa = gra[Pa->indGra+1];
                        }
                    }
                }
                if(Pa->isoC==-2 && Pa->graC==-2 && dir==5 && !modif)
                {
                    modif=true;
                    int indT = Pa->indGra+1;
                    if(indT==tour.size())
                        indT=0;
                    Pa = tour[indT];
                }
                if((Pa->graC==-2 || Pa->isoC==-2) && dir!=5 && !modif)
                {
                    dir=5;
                    int indT = -1;
                    if(Pa->graC==-2)
                        indT = Pa->indGra+1;
                    else
                        indT = Pa->isoInd+1;
                    if(indT==tour.size())
                        indT=0;
                    Pa = tour[indT];
                    modif=true;
                }
                if(continu || ((dir==1 || dir==3 || dir==5) && Pa->isoC!=-2 && !modif))
                {
                    vector<PointCourbe*> iso = newCourbeFinale[(int)(round(Pa->val*1000))][Pa->isoC];
                    if(continu)
                    {
                        modif=true;
                        if(dir==2)
                        {
                            int indI = Pa->isoInd-1;
                            if(indI<0)
                                indI=iso.size()-1;
                            Pa = iso[indI];
                        }
                        else
                        {
                            int indI = Pa->isoInd+1;
                            if(indI>=iso.size())
                                indI=0;
                            Pa = iso[indI];
                        }
                    }
                    else
                    {
                        if(dir==1 || (dir==5 && Pa->isoInd!=0))
                        {
                            modif=true;
                            dir=2;
                            int indI = Pa->isoInd-1;
                            if(indI<0)
                                indI=iso.size()-1;
                            Pa = iso[indI];
                        }
                        else
                        {
                            modif=true;
                            dir=4;
                            int indI = Pa->isoInd+1;
                            if(indI>=iso.size())
                                indI=0;
                            Pa = iso[indI];
                        }
                    }
                }
                
            }while(!equal(*Pa,*Po));
            allFaces.push_back(face);
            
        }
    }
    debug<<"\tAprès la première phase : "<<allFaces.size()<<" faces crées"<<endl;

    //On forment les cellules qui sont le long des iso-lignes du côtés où la fonction est plus hautes
    //Pour une iso-ligne de valeur 0.1, on forme les cellules ayant cette iso-ligne et une de valeur 0.2 comme bordure (s'il n'y a pas d'iso-ligne intermédiaire entre les deux)
    for(int i=0; i<targets.size(); i++)
    {
        int v = (int)(round(targets[i]*1000));
        vector<vector<PointCourbe*>> c = newCourbeFinale[v];
        for(int k=0; k<c.size(); k++)
        {
            vector<PointCourbe*> isoOri = c[k];
            for(int j=0; j<isoOri.size(); j++)
            {
                vector<PointCourbe*> face;
                PointCourbe* Po = isoOri[j];
                if(j<isoOri.size()-1 || Po->graC!=-2)
                {
                    PointCourbe* Pa = Po;
                    int dir = 3;
                    if(Pa->graC==-2)
                        dir = 5;
                    int dGra = 1;
                    bool Isostart = true;
                    
                    do{
                        bool modif = false;
                        face.push_back(Pa);
                        bool continu = false;
                        if(Pa->graC!=-2 && (dir==2 || dir==4 || dir==5))
                        {
                            vector<PointCourbe*> gra = NewCourbeGrad[Pa->graC];
                            if((dir==2 && Pa->indGra==0) || (dir==4 && Pa->indGra==gra.size()-1))
                            {
                                continu = true;
                                if(Isostart)
                                {
                                    j++;
                                }
                            }
                            else
                            {
                                Isostart = false;
                                if(dir==2 || (dir==5 && Pa->indGra!=0))
                                {
                                    modif=true;
                                    dir=3;
                                    Pa = gra[Pa->indGra-1];
                                }
                                else
                                {
                                    modif=true;
                                    dir=1;
                                    Pa = gra[Pa->indGra+1];
                                }
                            }
                        }
                        if(Pa->isoC==-2 && Pa->graC==-2 && dir==5 && !modif)
                        {
                            modif=true;
                            int indT = Pa->indGra+1;
                            if(indT==tour.size())
                                indT=0;
                            Pa = tour[indT];
                        }
                        if((Pa->graC==-2 || Pa->isoC==-2) && dir!=5 && !modif)
                        {
                            Isostart = false;
                            dir=5;
                            int indT = -1;
                            if(Pa->graC==-2)
                                indT = Pa->indGra+1;
                            else
                                indT = Pa->isoInd+1;
                            if(indT==tour.size())
                                indT=0;
                            Pa = tour[indT];
                            modif=true;
                        }
                        if(continu || ((dir==1 || dir==3 || dir==5) && Pa->isoC!=-2 && !modif))
                        {
                            vector<PointCourbe*> iso = newCourbeFinale[(int)(round(Pa->val*1000))][Pa->isoC];
                            if(continu)
                            {
                                modif=true;
                                if(dir==2)
                                {
                                    int indI = Pa->isoInd-1;
                                    if(indI<0)
                                        indI=iso.size()-1;
                                    Pa = iso[indI];
                                }
                                else
                                {
                                    int indI = Pa->isoInd+1;
                                    if(indI>=iso.size())
                                        indI=0;
                                    Pa = iso[indI];
                                }
                            }
                            else
                            {
                                if(dir==1 || (dir==5 && Pa->isoInd!=0))
                                {
                                    modif=true;
                                    dir=2;
                                    int indI = Pa->isoInd-1;
                                    if(indI<0)
                                        indI=iso.size()-1;
                                    Pa = iso[indI];
                                }
                                else
                                {
                                    modif=true;
                                    dir=4;
                                    int indI = Pa->isoInd+1;
                                    if(indI>=iso.size())
                                        indI=0;
                                    Pa = iso[indI];
                                }
                            }
                        }
                        
                    }while(!equal(*Pa,*Po));
                    allFaces.push_back(face);
                    
                }
            }
        }
    }

    debug<<"\tAprès la deuxième phase : "<<allFaces.size()<<" faces crées"<<endl;

    //On forment les cellules le long des iso-lignes n'ayant pas d'iso-ligne de valeur inférieure et ne touchant pas le bord
    //Ces iso-lignes sont celles qui sont autour des minimums locaux et qui ne touchent pas le bord
    for(int i=0; i<targets.size() && !prem; i++)
    {
        int v = (int)(round(targets[i]*1000));
        vector<vector<PointCourbe*>> c = newCourbeFinale[v];
        for(int k=0; k<c.size(); k++)
        {
            vector<PointCourbe*> isoOri = c[k];
            bool ok = true;
            for(int j=0; j<isoOri.size(); j++)
            {
                if(isoOri[j]->graC==-2 || isoOri[j]->indGra!=0)
                    ok = false;
            }
            if(ok)
            {
                for(int j=1; j<isoOri.size(); j++)
                {
                    vector<PointCourbe*> face;
                    PointCourbe* Po = isoOri[j];
                    PointCourbe* Pa = Po;
                    int dir = 1;
                    if(Pa->graC==-2)
                        dir = 5;
                    int dGra = 1;
                    bool Isostart = true;
                    
                    do{
                        bool modif = false;
                        face.push_back(Pa);
                        bool continu = false;
                        if(Pa->graC!=-2 && (dir==2 || dir==4 || dir==5))
                        {
                            vector<PointCourbe*> gra = NewCourbeGrad[Pa->graC];
                            if((dir==2 && Pa->indGra==0) || (dir==4 && Pa->indGra==gra.size()-1))
                            {
                                continu = true;
                                if(Isostart)
                                {
                                    
                                    j++;
                                }
                            }
                            else
                            {
                                Isostart = false;
                                if(dir==2 || (dir==5 && Pa->indGra!=0))
                                {
                                    modif=true;
                                    dir=3;
                                    Pa = gra[Pa->indGra-1];
                                }
                                else
                                {
                                    modif=true;
                                    dir=1;
                                    Pa = gra[Pa->indGra+1];
                                }
                            }
                        }
                        if(Pa->isoC==-2 && Pa->graC==-2 && dir==5 && !modif)
                        {
                            modif=true;
                            int indT = Pa->indGra+1;
                            if(indT==tour.size())
                                indT=0;
                            Pa = tour[indT];
                        }
                        if((Pa->graC==-2 || Pa->isoC==-2) && dir!=5 && !modif)
                        {
                            Isostart = false;
                            dir=5;
                            int indT = -1;
                            if(Pa->graC==-2)
                                indT = Pa->indGra+1;
                            else
                                indT = Pa->isoInd+1;
                            if(indT==tour.size())
                                indT=0;
                            Pa = tour[indT];
                            modif=true;
                        }
                        if(continu || ((dir==1 || dir==3 || dir==5) && Pa->isoC!=-2 && !modif))
                        {
                            vector<PointCourbe*> iso = newCourbeFinale[(int)(round(Pa->val*1000))][Pa->isoC];
                            if(continu)
                            {
                                modif=true;
                                if(dir==2)
                                {
                                    int indI = Pa->isoInd-1;
                                    if(indI<0)
                                        indI=iso.size()-1;
                                    Pa = iso[indI];
                                }
                                else
                                {
                                    int indI = Pa->isoInd+1;
                                    if(indI>=iso.size())
                                        indI=0;
                                    Pa = iso[indI];
                                }
                            }
                            else
                            {
                                if(dir==1 || (dir==5 && Pa->isoInd!=0))
                                {
                                    modif=true;
                                    dir=2;
                                    int indI = Pa->isoInd-1;
                                    if(indI<0)
                                        indI=iso.size()-1;
                                    Pa = iso[indI];
                                }
                                else
                                {
                                    modif=true;
                                    dir=4;
                                    int indI = Pa->isoInd+1;
                                    if(indI>=iso.size())
                                        indI=0;
                                    Pa = iso[indI];
                                }
                            }
                        }
                        
                    }while(!equal(*Pa,*Po));
                    allFaces.push_back(face);
                    
                }
            }
            
        }
    }

    debug<<"\tAprès la troisième phase : "<<allFaces.size()<<" faces crées"<<endl;

    /*
    cout<<"Affichage des courbes de Bézier approximant les iso-lignes (bleu) et courbes de gradient (vert)"<<endl;
    debug<<"Affichage des courbes de Bézier approximant les iso-lignes (bleu) et courbes de gradient (vert)"<<endl;
    //On déssine les iso-lignes approximées par des courbes de Béziers
    for(int i=0; i<targets.size(); i++)
    {
        int v = (int)(round(targets[i]*1000));
        
        vector<vector<PointCourbe*>> c = newCourbeFinale[v];
        for(int j=0; j<c.size(); j++)
        {
            vector<PointCourbe*> iso = c[j];
            for(int k=0; k<iso.size()-1; k++)
            {
                Point P0 = iso[k]->P;
                Point P1 = iso[k+1]->P;
                
                double distP = dist2(P0,P1);
                Repere R0 = getRepere(P0);
                Repere R1 = getRepere(P1);
                Point Pc0 = addP(P0,R0.T,distP*0.4);
                Point Pc1 = addP(P1,R1.T,-distP*0.4);
                glColor3f(0,0,1);
                courbeBezier(P0,Pc0,Pc1,P1,true);
                
            }
            if(iso[0]->graC!=-2)
            {
                Point P0 = iso[iso.size()-1]->P;
                Point P1 = iso[0]->P;
                double distP = dist2(P0,P1);
                Repere R0 = getRepere(P0);
                Repere R1 = getRepere(P1);
                Point Pc0 = addP(P0,R0.T,distP*0.4);
                Point Pc1 = addP(P1,R1.T,-distP*0.4);
                glColor3f(0,0,1);
                courbeBezier(P0,Pc0,Pc1,P1,true);
            }
        }
    }

    //On déssine les courbes de gradients approximées par des courbes de Béziers
    for(int i=0; i<NewCourbeGrad.size(); i++)
    {
        bool sensCorrect = true;
        vector<PointCourbe*> c = NewCourbeGrad[i];
        for(int j=0; j<c.size()-1; j++)
        {
            PointCourbe* Pcc0 = c[j];
            PointCourbe* Pcc1 = c[j+1];
            Point P0 = Pcc0->P;
            Point P1 = Pcc1->P;
            double distP = dist2(P0,P1);
            Repere R0 = getRepere(P0);
            Repere R1 = getRepere(P1);
            double dt2 = (0.3*distP-0.5*distP)/(normeMin-normeMax);
            double off2 = 0.3*distP - dt*normeMin;
            double m0 = dist(R0.N.x,R0.N.y,R0.N.z);
            double r0 =  m0* dt + off;
            double r02 = r0*dt2 + off2;
            double m1 = dist(R1.N.x,R1.N.y,R1.N.z);
            double r1 =  m1* dt + off;
            double r12 = r1*dt2 + off2;
            Point Pc0 = addP(P0,R0.N,r02);
            Point Pc1 = addP(P1,R1.N,-r12);
            int v=-2;
            glColor3f(0,1,0);
            courbeBezier(P0,Pc0,Pc1,P1,true);
        }
    }
    */

    cout<<"Remplissage de la structure half-edge"<<endl;
    debug<<"Remplissage de la structure half-edge"<<endl;
    //On remplit la strucutre half-edge avec les données des cellules
    int numHe = 0;
    for(int i=0; i<allFaces.size(); i++)
    {
        vector<PointCourbe*> face = allFaces[i];
        if(i%50==0 || i==allFaces.size()-1)
        {
            debug<<"\t"<<i<<"/"<<allFaces.size()-1<<endl;
        }
        string strP = "";
        string strSt = "";
        for(int j=0; j<face.size(); j++)
        {
            PointCourbe* Pcs = face[j];
            int iN = j+1;
            if(iN>=face.size())
                iN=0;
            PointCourbe* Pcn = face[iN];
            string str1 = to_string((int)(Pcs->val*1000))+to_string(Pcs->isoC)+to_string(Pcs->isoInd)+to_string(Pcs->graC)+to_string(Pcs->indGra);
            string str0 = to_string((int)(Pcn->val*1000))+to_string(Pcn->isoC)+to_string(Pcn->isoInd)+to_string(Pcn->graC)+to_string(Pcn->indGra);
            string str = str0+":"+str1;
            Point Po = Pcs->P;
            Point Pi = Pcn->P;
            hes[str]=new half_edge;
            hes[str]->id=numHe;
            hes[str]->incident = Pi;
            hes[str]->origine = Po;
            hes[str]->face = i;
            hes[str]->densityI = val(Pi.x,Pi.y,Pi.z);
            hes[str]->densityO = val(Po.x,Po.y,Po.z);
            hes[str]->next = &nullHe;
            hes[str]->previous = &nullHe;
            hes[str]->opposite = &nullHe;
            double distP = dist2(Po,Pi);
            if(abs(Pcs->val-Pcn->val)<0.0001 && Pcs->isoC==Pcn->isoC && Pcs->isoC!=-2)
            {
                vector<PointCourbe*> iso = newCourbeFinale[(int)(round(Pcs->val*1000))][Pcs->isoC];
                int sens = -1;
                if(Pcn->isoInd>Pcs->isoInd)
                {
                    sens = 1;
                }
                if((Pcs->isoInd==0 && Pcn->isoInd==iso.size()-1) || (Pcn->isoInd==0 && Pcs->isoInd==iso.size()-1))
                    sens*=-1;
                Repere Ro = getRepere(Po);
                Repere Ri = getRepere(Pi);
                Point Pco = addP(Po,Ro.T,0.4*distP*sens);
                Point Pci = addP(Pi,Ri.T,0.4*distP*sens*-1);
                hes[str]->controleI = Pci;
                hes[str]->controleO = Pco;
                hes[str]->densityCI = val(Pci.x,Pci.y,Pci.z);
                hes[str]->densityCO = val(Pco.x,Pco.y,Pco.z);
            }
            else
            {
                if(Pcs->graC==Pcn->graC && Pcs->graC!=-2)
                {
                    int sens = -1;
                    if(hes[str]->densityI>hes[str]->densityO)
                    {
                        sens = 1;
                    }
                    Repere Ro = getRepere(Po);
                    Repere Ri = getRepere(Pi);
                    double dt2 = (0.3*distP-0.5*distP)/(normeMin-normeMax);
                    double off2 = 0.3*distP - dt*normeMin;
                    double m0 = dist(Ro.N.x,Ro.N.y,Ro.N.z);
                    double r0 =  m0* dt + off;
                    double r02 = r0*dt2 + off2;
                    double m1 = dist(Ri.N.x,Ri.N.y,Ri.N.z);
                    double r1 =  m1* dt + off;
                    double r12 = r1*dt2 + off2;
                    Point Pco = addP(Po,Ro.N,r02*sens);
                    Point Pci = addP(Pi,Ri.N,-r12*sens);
                    hes[str]->controleI = Pci;
                    hes[str]->controleO = Pco;
                    hes[str]->densityCI = val(Pci.x,Pci.y,Pci.z);
                    hes[str]->densityCO = val(Pco.x,Pco.y,Pco.z);
                }
                else
                {
                    Point Vo = {0,0,0};
                    Point Vi = {0,0,0};
                    int trouve = 0;
                    for(int k=0; k<surface.size() && trouve<2; k++)
                    {
                        Point P0 = surface[k];
                        int iN = k+1;
                        if(iN>=surface.size())
                            iN=0;
                        Point P1 = surface[iN];
                        double distS = dist2(P0,P1);
                        double distO = distProj(P0,P1,Po);
                        double dAdj=0;
                        for(int l=0; l<pointsImportants.size() && dAdj<0.000001; l++)
                        {
                            if(iN=pointsImportants[l])
                                dAdj=normeMax;
                        }
                        if(distO>=0 && distO<distS+dAdj)
                        {
                            Point PProj= interPoint(P1,P0,distO/distS);
                            double dSp = dist2(P0,PProj);
                            double distProj = dist2(PProj,Po);
                            if(distProj<normeMax)
                            {
                                if(dist(Vo.x,Vo.y,Vo.z)<=0.0001)
                                    trouve++;
                                Vo = addP(P1,P0,-1);
                            }
                        }
                        double distI = distProj(P0,P1,Pi);
                        if(distI>=0 && distI<distS+dAdj)
                        {
                            Point PProj= interPoint(P1,P0,distI/distS);
                            double dSp = dist2(P0,PProj);
                            double distProj = dist2(PProj,Pi);
                            if(distProj<normeMax)
                            {
                                if(dist(Vi.x,Vi.y,Vi.z)<=0.0001)
                                    trouve++;
                                Vi = addP(P1,P0,-1);
                            }
                        }
                    }
                    int sens = 1;
                    if(tourInverse)
                    {
                        sens = -1;
                    }
                    double normO = dist(Vo.x,Vo.y,Vo.z);
                    Vo = addP({0,0,0},Vo,sens*1.0/normO);
                    double normI = dist(Vi.x,Vi.y,Vi.z);
                    Vi = addP({0,0,0},Vi,sens*1.0/normI);
                    int indO = -1;
                    int indI = -1;
                    if(Pcs->isoC==-2)
                    {
                        indO=Pcs->isoInd;
                    }
                    else
                    {
                        indO=Pcs->indGra;
                    }
                    if(Pcn->isoC==-2)
                    {
                        indI=Pcn->isoInd;
                    }
                    else
                    {
                        indI=Pcn->indGra;
                    }
                    sens = 1;
                    int nb1 = abs(indO-indI);
                    int nb2 = tour.size()-nb1;
                    nb2 = abs(nb2);
                    if(indO>indI)
                        sens*=-1;
                    if(nb2<nb1)
                        sens*=-1;
                    Point Pco = addP(Po,Vo,0.4*distP*sens);
                    Point Pci = addP(Pi,Vi,0.4*distP*sens*-1);
                    hes[str]->controleI = Pci;
                    hes[str]->controleO = Pco;
                    hes[str]->densityCI = val(Pci.x,Pci.y,Pci.z);
                    hes[str]->densityCO = val(Pco.x,Pco.y,Pco.z);
                }
            }
            if(j!=0)
            {
                hes[strP]->next = hes[str];
                hes[str]->previous = hes[strP];
            }
            else
            {
                strSt = str;
            }
            if(j==face.size()-1)
            {
                hes[str]->next = hes[strSt];
                hes[strSt]->previous = hes[str];
            }
            string strOpp = str1+":"+str0;
            if(hes.find(strOpp)!=hes.end())
            {
                hes[str]->opposite=hes[strOpp];
                hes[strOpp]->opposite=hes[str];
            }
            strP = str;
            numHe++;
        }
        if(face.size()>0)
        {
            heFace hef = {i,0,hes[strSt]};
            Faces.push_back(hef);
        }
    }

    cout<<"Traitement sur les cellules (Fusion des points et des cellules)"<<endl;
    debug<<"Traitement sur les cellules (Fusion des points et des cellules)"<<endl;

    debug<<"\tFusion des groupes de points"<<endl;
    //On fusionne les groupes de points proches afin de limiter la taille minimale des bords des cellules finales
    glColor3f(1,0,0);
    int nbX = 0;
    for(int i=0; i<Faces.size(); i++)
    {
        if(i%50==0 || i==Faces.size()-1)
        {
            debug<<"\t\t"<<i<<"/"<<Faces.size()-1<<endl;
        }
        heFace hef = Faces[i];
        
        half_edge* heo = hef.incidente;
        half_edge* he = heo;
        bool change = false;
        do{
            double longE = courbeBezier(he->origine,he->controleO,he->controleI,he->incident,false);
            if(longE<normeMax)
            {
               
                
                nbX++;
                change = true;
                
                vector<half_edge*> edges;
                half_edge* heop = he;
                edges.insert(edges.begin(),heop);
                double l = longE;
                bool dessin = false;
                while(heop->previous->opposite->id!=-1)
                {
                    
                    glColor3f(1,0,0);
                    l = courbeBezier(heop->origine,heop->controleO,heop->controleI,heop->incident,dessin);
                    
                    if(l>normeMax)
                    {
                        break;
                    }
                    
                    half_edge* hep = heop->previous;
                    Point Pp = hep->controleI;
                    Point P = heop->origine;
                    Point Pn = heop->controleO;
                    Point Vp = addP(P,Pp,-1);
                    Point Vn = addP(Pn,P,-1);
                    double ang = abs(getAngle(Vp,Vn));
                    if((ang>80 && ang<100) || heop->opposite->id==-1)
                    {
                        heop = heop->previous->opposite->previous;
                    }
                    else
                    {
                        heop = heop->previous;
                    }
                    if(heop->id == he->id)
                    {
                        break;
                    }
                    else
                        edges.insert(edges.begin(),heop);
                }
                glColor3f(1,0,0);
                courbeBezier(heop->origine,heop->controleO,heop->controleI,heop->incident,dessin);
                
                half_edge* heon = he;
                double l2 = longE;
                while(heon->next->opposite->id!=-1)
                {
                    glColor3f(1,1,0);
                    l2 = courbeBezier(heon->origine,heon->controleO,heon->controleI,heon->incident,dessin);
                    
                    if(l2>normeMax)
                    {
                        break;
                    }
                    half_edge* hen = heon->next;
                    Point Pp = heon->controleI;
                    Point P = hen->origine;
                    Point Pn = hen->controleO;
                    Point Vp = addP(P,Pp,-1);
                    Point Vn = addP(Pn,P,-1);
                    double ang = abs(getAngle(Vp,Vn));
                    if((ang>80 && ang<100) || heop->opposite->id==-1)
                    {
                        heon = heon->next->opposite->next;
                    }
                    else
                    {
                        heon = heon->next;
                    }
                    if(heon->id == he->id)
                    {
                        break;
                    }
                    else
                        edges.push_back(heon);
                    
                }
                glColor3f(1,1,0);
                courbeBezier(heon->origine,heon->controleO,heon->controleI,heon->incident,dessin);
                
                Point somme = {0,0,0};
                int nb = 0;
                if(l>normeMax)
                {
                    somme = addP(somme, heop->incident,1);
                    nb++;
                }
                else
                {
                    somme = addP(somme, heop->incident,1);
                    somme = addP(somme, heop->origine,1);
                    nb+=2;
                }
                if(l2>normeMax)
                {
                    somme = addP(somme, heon->origine,1);
                    nb++;
                }
                else
                {
                    somme = addP(somme, heon->incident,1);
                    somme = addP(somme, heon->origine,1);
                    nb+=2;
                }
                for(int j=1; j<edges.size()-1; j++)
                {
                    somme = addP(somme, edges[j]->incident,1);
                    somme = addP(somme, edges[j]->origine,1);
                    nb+=2;
                }
                somme = addP({0,0,0},somme,1/(double)nb);
                
                if(l>normeMax)
                {
                    heop->incident = somme;
                    heop->opposite->origine = somme;
                }
                else
                {
                    heop->next->previous = heop->previous;
                    heop->previous->next = heop->next;
                    heop->previous->incident=somme;
                    heop->next->origine=somme;
                    heFace hef = Faces[heop->face];
                    hef.incidente = heop->next;
                    Faces[heop->face]=hef;
                    half_edge* opp = heop->opposite;
                    if(opp->id!=-1)
                    {
                        opp->next->previous = opp->previous;
                        opp->previous->next = opp->next;
                        opp->previous->incident=somme;
                        opp->next->origine=somme;
                        hef = Faces[opp->face];
                        hef.incidente = opp->next;
                        Faces[opp->face]=hef;
                    }
                }
                if(l2>normeMax)
                {
                    heon->origine = somme;
                    heon->opposite->incident = somme;
                }
                else
                {
                    heon->next->previous = heon->previous;
                    heon->previous->next = heon->next;
                    heon->previous->incident=somme;
                    heon->next->origine=somme;
                    heFace hef = Faces[heon->face];
                    hef.incidente = heon->next;
                    Faces[heon->face]=hef;
                    half_edge* opp = heon->opposite;
                    if(opp->id!=-1)
                    {
                        opp->next->previous = opp->previous;
                        opp->previous->next = opp->next;
                        opp->previous->incident=somme;
                        opp->next->origine=somme;
                        hef = Faces[opp->face];
                        hef.incidente = opp->next;
                        Faces[opp->face]=hef;
                    }
                }
                
                for(int j=1; j<edges.size()-1; j++)
                {
                    half_edge* hea = edges[j];
                    
                    hea->next->previous = hea->previous;
                    hea->previous->next = hea->next;
                    hea->previous->incident=somme;
                    hea->next->origine=somme;
                    heFace hef = Faces[hea->face];
                    hef.incidente = hea->next;
                    Faces[hea->face]=hef;
                    half_edge* opp = hea->opposite;
                    if(opp->id!=-1)
                    {
                        
                        opp->next->previous = opp->previous;
                        opp->previous->next = opp->next;
                        opp->previous->incident=somme;
                        opp->next->origine=somme;
                        hef = Faces[opp->face];
                        hef.incidente = opp->next;
                        Faces[opp->face]=hef;
                    }
                }
                
                break;
            }
            he = he->next;
        }while(he->id!=heo->id);
    }

    debug<<"\tFusion des cellules de petites tailles"<<endl;
    //On identifie les cellules de petites tailles et on les fusionnent avec une cellule voisine
    vector<int> faceSup;
    glColor3f(1,0,0);
    for(int i=0; i<Faces.size(); i++)
    {
        if(i%50==0 || i==Faces.size()-1)
        {
            debug<<"\t\t"<<i<<"/"<<Faces.size()-1<<endl;
        }
        heFace hef = Faces[i];
        half_edge* heO = hef.incidente;
        half_edge* he = heO;
        int nbC = 0;
        double longT = 0;
        do{
            double longB = courbeBezier(he->origine,he->controleO,he->controleI,he->incident,false);
            longT+=longB;
            he = he->next;
            nbC++;
        }while(he->id!=heO->id);
        if(longT<nbC*normeMax*5)
        {
            
            heO = hef.incidente;
            he = heO;
            do{
                
                
                he = he->next;
            }while(he->id!=heO->id && he->opposite->id==-1);
            
            half_edge* opp = he->opposite;
            if(opp->id!=-1)
            {
                faceSup.push_back(hef.id);
                
                heFace hef = Faces[opp->face];
                hef.incidente = opp->next;
                Faces[opp->face]=hef;
                opp->next->previous = he->previous;
                he->previous->next = opp->next;
                opp->previous->next = he->next;
                he->next->previous = opp->previous;
            }
        }
        if(nbC<3)
        {
            
            faceSup.push_back(hef.id);
        }
    }

    //On supprime les données de cellules fusionnées
    for(auto it=Faces.begin(); it!=Faces.end() && !supT; it++)
    {
        for(int i=0; i<faceSup.size(); i++)
        {
            if(it->id == faceSup[i])
            {
                Faces.erase(it);
                it = Faces.begin();
                break;
            }
        }
    }

    //On renumérote les faces et les faces incidentes des demi-arrêtes
    for(int i=0; i<Faces.size(); i++)
    {
        heFace hef = Faces[i];
        hef.id = i;
        half_edge* heo = hef.incidente;
        half_edge* he = heo;
        do{
            he->face = i;
            Point milieu = {0,0,0};
            milieu = addP(milieu,he->origine,0.5);
            milieu = addP(milieu,he->previous->incident,0.5);
            he->origine = milieu;
            he->previous->incident = milieu;
            he = he->next;
        }while(he->id!=heo->id);
        Faces[i]=hef;
    }

    cout<<"Affichage des cellules (rouge) et approximation de la densité moyenne"<<endl;
    debug<<"Affichage des cellules (rouge) et approximation de la densité moyenne"<<endl;
    //On déssine les cellules finales et on calcule une approximation de leur densité
    glColor3f(1,0,0);
    for(int i=0; i<Faces.size(); i++)
    {
        if(i%50==0 || i==Faces.size()-1)
        {
            debug<<"\t"<<i<<"/"<<Faces.size()-1<<endl;
        }
        heFace hef = Faces[i];
        half_edge* heO = hef.incidente;
        half_edge* he = heO;
        double longT = 0;
        double dens = 0;
        do{
            double longB = courbeBezier(he->origine,he->controleO,he->controleI,he->incident,true);
            longT+=longB;
            double de = (he->densityO+he->densityCO+he->densityCI+he->densityI)/4.0;
            dens += de*longB;
            he = he->next;
        }while(he->id!=heO->id);
        if(longT>0)
            hef.density = dens/longT;
        else
            hef.density = 0;
        Faces[i]=hef;
    }

    cout<<"Deuxième approximation de la densité moyenne des cellules"<<endl;
    debug<<"Deuxième approximation de la densité moyenne des cellules"<<endl;
    //On calcule une meilleure approximation de la densité des cellules par propagation de points
    for(int i=0; i<Faces.size(); i++)
    {
        heFace f = Faces[i];
        float t = f.density;
        f.density = densite(f);
        debug<<"\tFace n°"<<i<<" : première approximation = "<<t<<" ; deuxième approximation = "<<f.density<<endl;
        Faces[i]=f;
    }
    cout<<"Ecriture du fichier JSON"<<endl;
    debug<<"Ecriture du fichier JSON"<<endl;
    //On écrit les données de cellules dans un fichier JSON
    ofstream myfile;
    myfile.open ("output.json",ios::out | ios::trunc);
    myfile << "{" <<endl<<"\t\"Faces\":["<<endl;
    for(int i=0;i<Faces.size();i++)
    {
        if(i%50==0 || i==Faces.size()-1)
        {
            debug<<"\t\t"<<i<<"/"<<Faces.size()-1<<endl;
        }
        heFace f = Faces[i];
        myfile<<"\t\t{"<<endl<<"\t\t\t\"Density\":"<<f.density<<","<<endl<<"\t\t\t\"Edges\":["<<endl;
        half_edge* he = f.incidente;
        int id_ori = he->id;
        do
        {
            
            
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
    debug.close();
    cout<<"Fin"<<endl;
    glPopMatrix();
  /* on force l'affichage du resultat */
  //glFlush();
    glutSwapBuffers();
}

//------------------------------------------------------


//------------------------------------------------------
void clavier(unsigned char touche,int x,int y)
{

  switch (touche)
    {
    case '+': //* affichage du carre plein 
        nbF++;
        //if(nbF>=listeFace.size()) nbF = listeFace.size()-1;
      glutPostRedisplay();
      break;
    case '-': //* affichage du carre plein 
        nbF--;
	//if (nbF < 1 ) nbF=0;
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

    
    
