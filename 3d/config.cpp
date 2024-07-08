#include "config.h"

Config::Config()
{
    max = 8;
    p.push_back(false);
}

int Config::val()
{
    int v = 0;
    int p2 = 1;
    for(int i=0; i<p.size(); i++)
    {
        v+= p[i]*p2;
        p2*=2;
    }
    return v;
}

Config Config::add(int n)
{
    Config r;
    r.p = p;
    r.max = max;
    if(r.p.size()==0)
    	r.p.push_back(false);
    //temporaire
    for(int i=0; i<n; i++)
    {
        r = r.next();
    }
    return r;
}

Config Config::next()
{
    Config r;
    r.p = p;
    r.max = max;
    if(r.p.size()==0)
    {
        r.p.push_back(false);
    }
    else
    {
        bool aug = false;
        if(r.p[0])
        {
            r.p[0]=false;
            aug=true;
        }
        else
        {
            r.p[0]=true;
        }
        for(int i=1; i<p.size() && aug; i++)
        {
            if(r.p[i])
            {
                r.p[i]=false;
            }
            else
            {
                r.p[i]=true;
                aug=false;
            }
        }
        if(aug && r.p.size()<max)
        {
            r.p.push_back(true);
        }
    }
    return r;
}

string Config::toString()
{
    string s = "";
    for(int i=0; i<p.size(); i++)
    {
        s=to_string(p[i])+s;
    }
    if(p.size()<max)
    {
        for(int i=p.size(); i<max; i++)
        {
            s="0"+s;
        }
    }
    return s;
}

void Config::inv()
{
    while(p.size()<max)
    {
        p.push_back(false);
    }
    for(int i=0; i<p.size(); i++)
    {
        p[i]=!p[i];
    }
}

Config Config::inv(int n)
{
	Config r;
	r.p=p;
	r.max=max;
	if(n>max)
		return r;
	else
	{
		if(n>=r.p.size())
		{
			for(int i=r.p.size(); i<=n; i++)
			{
				r.p.push_back(false);
			}
		}
		r.p[n]=!r.p[n];
		return r;
	}
}