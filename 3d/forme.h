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

#include "config.h"

#ifndef FORME_H
#define FORME_H

using namespace std;

struct Point{
    double x;
    double y;
    double z;
    Point mult(double v)
    {
    	Point P;
    	P.x = x*v;
    	P.y = y*v;
    	P.z = z*v;
    	return P;
    }
    Point add(Point P)
    {
    	Point R = {P.x,P.y,P.z};
    	R.x += x;
    	R.y += y;
    	R.z += z;
    	return R;
    }
};

struct Triangle{
	Point P0;
	Point P1;
	Point P2;
};

struct Zone{
	vector<Triangle> faces; 
};

class Forme
{
	private :
		int type;
		bool points[8];
		vector<Forme> sousForme;
	public :
		Forme(int t, int p0, int p1, int p2, int p3);
		Forme();
		vector<Zone> dessin(Point p[8], double val[8]);
		string toString();
};

#endif