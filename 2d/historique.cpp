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

Point findGradient(Point P)
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
    moy+=elapsed.count();
    ng++;
    return R1;
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
    
    if(ind<((x_size/pas)*(y_size/pas))/4)
    {
        int ix = round(x*1000);
        int iy = round(y*1000);
        int ip = round(pas*1000);
        string k = to_string(ix)+to_string(iy)+to_string(ip);
        
        if(Passage.find(k) == Passage.end())
        {
            Passage[k]=true;
            float p1 = target-func(x,y,0);
            float p2 = target-func(x+pas,y,0);
            float p3 = target-func(x+pas,y+pas,0);
            float p4 = target-func(x,y+pas,0);
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
                
                
                
                if(p1<0)
                {
                    if(p2<0)
                    {
                        if(p3<0)
                        {
                            if(p4<0)
                            {
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
                                
                            }
                        }
                        else
                        {
                            if(p4<0)
                            {
                                
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
                                
                            }
                            else
                            {
                                
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
                if(ch==0)
                {
                    if(abs(v0)<0.05 && abs(v1)<0.05 && abs(v2)<0.05 && abs(v3)<0.05)
                    {
                        bool p01 = false;
                        bool p03 = false;
                        bool p12 = false;
                        bool p23 = false;
                        for(int i=-9; i<10; i++)
                        {
                            double p = (double)i/10.0;
                            Point P01 = {P.x-Pi.z,P.y+p*Pi.z,P.z};
                            Point P12 = {P.x+p*Pi.z,P.y+Pi.z,P.z};
                            Point P23 = {P.x+Pi.z,P.y+p*Pi.z,P.z};
                            Point P03 = {P.x+p*Pi.z,P.y-Pi.z,P.z};
                            double v01 = abs(target-func(P01.x,P01.y,P01.z));
                            double v12 = abs(target-func(P12.x,P12.y,P12.z));
                            double v23 = abs(target-func(P23.x,P23.y,P23.z));
                            double v03 = abs(target-func(P03.x,P03.y,P03.z));
                            if(v01<abs(v0) && v01<abs(v1) && v01<0.01)
                                p01 = true;
                            if(v12<abs(v1) && v12<abs(v2) && v12<0.01)
                                p12 = true;
                            if(v23<abs(v2) && v23<abs(v3) && v23<0.01)
                                p23 = true;
                            if(v03<abs(v0) && v03<abs(v3) && v03<0.01)
                                p03 = true;
                        }
                        if(p01 && p12)
                        {
                            v1*=-1;
                        }
                        if(p01 && p03)
                        {
                            v0*=-1;
                        }
                        if(p01 && p23)
                        {
                            v0*=-1;
                            v3*=-1;
                        }
                        if(p12 && p23)
                        {
                            v2*=-1;
                        }
                        if(p23 && p03)
                        {
                            v3*=-1;
                        }
                        if(p12 && p03)
                        {
                            v0*=-1;
                            v1*=-1;
                        }
                    }
                }
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
                if(abs(Pi.z-pas)<min)
                {
                    pile.push_back({Pi.x-1,Pi.y,Pi.z});
                    pile.push_back({Pi.x+1,Pi.y,Pi.z});
                    pile.push_back({Pi.x,Pi.y-1,Pi.z});
                    pile.push_back({Pi.x,Pi.y+1,Pi.z});
                }
            }
        }
        pile.erase(pile.begin());
    }
    map<int,int> transfo;
    vector<int> plop;
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
            
            return Pr2;
        }
        else
        {
            return chercherBezier(he,Pc,0.5*start+0.5*end,end,step+1,dist);
        }
    }
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
    if(P.y<=y || st%5==0)
    {
        P = distanceMinBezier(f, {x,y,0},0.05);
        
        float dist = dist2({x,y,0},P);
        if(y>=P.y && abs(x-P.x)<0.02 && cherche)
        {
            cherche = false;
            in = !in;
            last = st;
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
            retour = calculDensité(x,y+0.05,in,P,yMax,f,cherche,0,last);
        }
        else
        {
            retour = calculDensité(x,P.y,in,P,yMax,f,cherche,0,last);
        }
    }
    else
    {
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
        
        
        
        he = he->next;
    }while(he->id!=id);
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

float densite2(heFace f)
{
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
    he1 = f.incidente;
    float l = 0;
    int nbC = 0;
    do
    {
        l+= courbeBezier(he1->origine,he1->controleO,he1->controleI,he1->incident,false);
        nbC++;
        he1 = he1->next;
    }while(he1->id!=id_ori);
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
    he1 = he;
    he2 = he1->previous;
    P = he1->origine;
    
    Point v0 = getVecteurBezier(he1,0,limite*2);
    Point v1 = getVecteurBezier(he2,1-limite*2,limite*2);
    
    
    
    
    glPointSize(10.0f);
    
    Point Po = {0,0,0};
    
    Po.x = P.x + v0.x - v1.x;
    Po.y = P.y + v0.y - v1.y;
    Po.z = P.z + v0.z - v1.z;
    
    Point Pc = distanceMinBezier2(f,Po,limite/2);
    
    
    glBegin(GL_POINTS);

    
    if(func(Po.x,Po.y,Po.z)==0)
    {
        
        Point V = {Pc.x-Po.x,Pc.y-Po.y,Pc.z-Po.z};
        float n = dist(V.x,V.y,V.z);
        V.x/=n;
        V.y/=n;
        V.z/=n;
        while(func(Po.x,Po.y,Po.z)==0)
        {
            
            Po.x+=V.x*0.01;
            Po.y+=V.y*0.01;
            Po.z+=V.z*0.01;
        }
        Po.x+=V.x*0.01;
        Po.y+=V.y*0.01;
        Po.z+=V.z*0.01;
        
    }
    
    float c = 1.0f;
    Pc = distanceMinBezier2(f,Po,limite/2);
    while(dist2(Pc,Po)<limite)
    {
        
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
        
        Pc = distanceMinBezier2(f,Po,limite/2);
    }
    glEnd();
    
    
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
            
        if(xmax - xmin > 1600 || ymax - ymin > 1600)
        {
            
            break;
        }
        plopX = false;
        plopY = false;
        string k = to_string(x)+":"+to_string(y);
        
        if(marque.find(k) == marque.end())
        {
            Point P = {Po.x+x*limite/4,Po.y+y*limite/4,Po.z};
            
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
                    
                    
                    pile.push_back({(float)x-1,(float)y,0});
                    pile.push_back({(float)x+1,(float)y,0});
                    pile.push_back({(float)x,(float)y-1,0});
                    pile.push_back({(float)x,(float)y+1,0});
                }
                else
                {
                    glPointSize(1.0f);
                    
                    
                }
            }
            else
            {
                glPointSize(1.0f);
                
                
            }
            glVertex3f(P.x,P.y,P.z);
            
        }
        
        pile.erase(pile.begin());
        
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