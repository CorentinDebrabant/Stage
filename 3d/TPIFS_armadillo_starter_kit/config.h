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

#ifndef CONFIG_H
#define CONFIG_H

using namespace std;

class Config{
    public :
        vector<bool> p;
        int max;
        Config();
        int val();
        Config add(int n);
        Config next();
        string toString();
        void inv();
        Config inv(int n);
};

#endif