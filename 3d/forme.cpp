#include "forme.h"


using namespace std;

double Forme::interpolation(double a, double b)
{
    double aa = abs(a);
    double ab = abs(b);
    if(aa+ab!=0)
        return ab/(aa+ab);
    return 0.5;
}

Forme::Forme(int t, int p0, int p1, int p2, int p3)
{
	type = t;
	for(int i=0; i<8; i++)
	{
		if(i==p0 || i==p1 || i==p2 || i==p3)
			points[i]=true;
		else
			points[i]=false;
	}

	switch(type)
	{
		case 3:
			sousForme.push_back(Forme(1,p0,-1,-1,-1));
			sousForme.push_back(Forme(1,p1,-1,-1,-1));
			break;
		case 6:
			sousForme.push_back(Forme(1,p0,-1,-1,-1));
			sousForme.push_back(Forme(4,p1,p2,p3,-1));
			break;
		case 7:
			sousForme.push_back(Forme(1,p0,-1,-1,-1));
			sousForme.push_back(Forme(1,p1,-1,-1,-1));
			sousForme.push_back(Forme(1,p2,-1,-1,-1));
			sousForme.push_back(Forme(1,p3,-1,-1,-1));
			break;
		case 10:
			sousForme.push_back(Forme(1,p0,-1,-1,-1));
			sousForme.push_back(Forme(1,p1,-1,-1,-1));
			break;
		case 11:
			sousForme.push_back(Forme(1,p0,-1,-1,-1));
			sousForme.push_back(Forme(2,p1,p2,-1,-1));
			break;
		case 12:
			sousForme.push_back(Forme(1,p0,-1,-1,-1));
			sousForme.push_back(Forme(1,p1,-1,-1,-1));
			sousForme.push_back(Forme(1,p2,-1,-1,-1));
			break;
		case 13:
			sousForme.push_back(Forme(2,p0,p1,-1,-1));
			sousForme.push_back(Forme(2,p2,p3,-1,-1));
			break;
	}
}

Forme::Forme()
{
	type = -1;
	for(int i=0; i<8; i++)
	{
		points[i]=false;
	}
}

string Forme::toString()
{
	string s = "Forme de type "+to_string(type)+" avec les points : (";
	for(int i=0; i<8; i++)
	{
		if(points[i])
		{
			s+=to_string(i)+", ";
		}
	}
	s+=") et "+to_string(sousForme.size())+" sous formes";
	return s;
}

vector<Zone> Forme::dessin(Point p[8], double val[8])
{
	vector<Zone> vz;
	if(sousForme.size()!=0)
	{
		for(int i=0; i<sousForme.size(); i++)
		{
			vector<Zone> vz0 = sousForme[i].dessin(p,val);
			vz.insert(vz.end(),vz0.begin(),vz0.end());
		}
	}
	else
	{
		Zone z;
		switch(type)
		{
			case 1:
				{
					Config c;
					c.max = 3;
					for(int i=0; i<8; i++)
					{
						if(points[i])
							c = c.add(i);
					}
					Point P0 = p[c.val()];
					double v0 = val[c.val()];
					Point P1 = p[c.inv(0).val()];
					double v1 = val[c.inv(0).val()];
					Point P2 = p[c.inv(1).val()];
					double v2 = val[c.inv(1).val()];
					Point P3 = p[c.inv(2).val()];
					double v3 = val[c.inv(2).val()];
					double inter01 = interpolation(v0,v1);
					double inter02 = interpolation(v0,v2);
					double inter03 = interpolation(v0,v3);
					inter01 = 0.5;
					inter02 = 0.5;
					inter03 = 0.5;
					Point P01 = P0.mult(inter01).add(P1.mult(1-inter01));
					Point P02 = P0.mult(inter02).add(P2.mult(1-inter02));
					Point P03 = P0.mult(inter03).add(P3.mult(1-inter03));
					Triangle t = {P01,P02,P03};
					Zone z;
					z.faces.push_back(t);
					vz.push_back(z);
				}
				break;
			case 2:
				{
					Config c;
					c.max=3;
					bool b = false;
					for(int i=0; i<8; i++)
					{
						if(points[i] && !b)
						{
							c = c.add(i);
							b = true;
						}
					}
					double v0, v1, v01, v02, v11, v12;
					v0 = val[c.val()];
					Point P0 = p[c.val()];
					Config c01 = c.inv(0);
					Config c02 = c.inv(1);
					Config c03 = c.inv(2);
					int i1 = c01.val();
					int i2 = c02.val();
					int i3 = c03.val();
					Config c1;
					Point P01 = {0,0,0};
					Point P02 = {0,0,0};
					Point P1 = {0,0,0};
					if(points[i1])
					{
						P1 = p[i1];
						P01 = p[i2];
						P02 = p[i3];
						c1 = c01;
						v1 = val[i1];
						v01 = val[i2];
						v02 = val[i3];
					}
					else
					{
						if(points[i2])
						{
							P1 = p[i2];
							P01 = p[i1];
							P02 = p[i3];
							c1 = c02;
							v1 = val[i2];
							v01 = val[i1];
							v02 = val[i3];
						}
						else
						{
							P1 = p[i3];
							P01 = p[i1];
							P02 = p[i2];
							c1 = c03;
							v1 = val[i3];
							v01 = val[i1];
							v02 = val[i2];
						}
					}
					c01 = c1.inv(0);
					c02 = c1.inv(1);
					c03 = c1.inv(2);
					i1 = c01.val();
					i2 = c02.val();
					i3 = c03.val();
					Point P11 = {0,0,0};
					Point P12 = {0,0,0};
					if(i1==c.val())
					{
						P11 = p[i2];
						P12 = p[i3];
						v11 = val[i2];
						v12 = val[i3];
					}
					else
					{
						if(i2==c.val())
						{
							P11 = p[i1];
							P12 = p[i3];
							v11 = val[i1];
							v12 = val[i3];
						}
						else
						{
							P11 = p[i1];
							P12 = p[i2];
							v11 = val[i1];
							v12 = val[i2];
						}
					}
					//P01 = P0.mult(inter01).add(P1.mult(1-inter01));
					double inter01 = interpolation(v0,v01);
					double inter02 = interpolation(v0,v02);
					double inter11 = interpolation(v1,v11);
					double inter12 = interpolation(v1,v12);
					P01 = (P01.add(P0)).mult(0.5);
					P02 = (P02.add(P0)).mult(0.5);
					P11 = (P11.add(P1)).mult(0.5);
					P12 = (P12.add(P1)).mult(0.5);
					Triangle t0 = {P01,P02,P11};
					Triangle t1 = {P11,P12,P02};
					Zone z;
					z.faces.push_back(t0);
					z.faces.push_back(t1);
					vz.push_back(z);
				}
				break;
			case 4:
				{
					vector<Config> cs;
					for(int i=0; i<8; i++)
					{
						if(points[i])
						{
							Config c;
							c.max = 3;
							c = c.add(i);
							cs.push_back(c);
						}
					}
					int i0 = cs[0].val();
					int i1 = cs[1].val();
					int i2 = cs[2].val();
					int c0 = -1;
					int c1 = -1;
					int c2 = -1;
					for(int i=0; i<3; i++)
					{
						int k=0;
						Point P = p[cs[i].val()];
						for(int j=0; j<3; j++)
						{
							int v = cs[i].inv(j).val();
							if(v!=i0 && v!=i1 && v!=i2)
							{
								k++;
							}
						}
						if(k==1)
							c0 = i;
						else
						{
							if(c1==-1)
								c1=i;
							else
								c2=i;
						}
					}
					int axe = -1;
					Point P0 = p[cs[c0].val()];
					Point P0axe;
					for(int j=0; j<3; j++)
					{
						int v = cs[c0].inv(j).val();
						if(v!=i0 && v!=i1 && v!=i2)
						{
							axe = j;
							P0axe = (p[v].add(P0)).mult(0.5);
						}
					}
					Point P1 = p[cs[c1].val()];
					Point P1axe;
					Point P1autre;
					for(int j=0; j<3; j++)
					{
						int v = cs[c1].inv(j).val();
						if(v!=i0 && v!=i1 && v!=i2)
						{
							if(j==axe)
							{
								P1axe = (p[v].add(P1)).mult(0.5);
							}
							else
							{
								P1autre = (p[v].add(P1)).mult(0.5);
							}
							
						}
					}
					Point P2 = p[cs[c2].val()];
					Point P2axe;
					Point P2autre;
					for(int j=0; j<3; j++)
					{
						int v = cs[c2].inv(j).val();
						if(v!=i0 && v!=i1 && v!=i2)
						{
							if(j==axe)
							{
								P2axe = (p[v].add(P2)).mult(0.5);
							}
							else
							{
								P2autre = (p[v].add(P2)).mult(0.5);
							}
							
						}
					}
					Triangle t0 = {P0axe,P1axe,P2axe};
					Triangle t1 = {P1autre,P1axe,P2axe};
					Triangle t2 = {P1autre,P2axe,P2autre};
					Zone z;
					z.faces.push_back(t0);
					z.faces.push_back(t1);
					z.faces.push_back(t2);
					vz.push_back(z);
				}
				break;
			case 5:
				{
					vector<Config> cs;
					for(int i=0; i<8; i++)
					{
						if(points[i])
						{
							Config c;
							c.max = 3;
							c = c.add(i);
							cs.push_back(c);
						}
					}
					int i0 = cs[0].val();
					int i1 = cs[1].val();
					int i2 = cs[2].val();
					int axe = -1;
					for(int i=0; i<3; i++)
					{
						int v = cs[3].inv(i).val();
						if(v!=i0 && v!=i1 && v!=i2)
							axe = i;
					}
					Point P0 = p[cs[0].val()];
					Point P01 = p[cs[0].inv(axe).val()];
					Point P1 = p[cs[1].val()];
					Point P11 = p[cs[1].inv(axe).val()];
					Point P2 = p[cs[2].val()];
					Point P21 = p[cs[2].inv(axe).val()];
					Point P3 = p[cs[3].val()];
					Point P31 = p[cs[3].inv(axe).val()];
					P01 = (P01.add(P0)).mult(0.5);
					P11 = (P11.add(P1)).mult(0.5);
					P21 = (P21.add(P2)).mult(0.5);
					P31 = (P31.add(P3)).mult(0.5);
					Triangle t0 = {P01,P11,P21};
					Triangle t1 = {P11,P21,P31};
					Zone z;
					z.faces.push_back(t0);
					z.faces.push_back(t1);
					vz.push_back(z);
				}
				break;
			case 8:
				{
					vector<Config> cs;
					for(int i=0; i<8; i++)
					{
						if(points[i])
						{
							Config c;
							c.max = 3;
							c = c.add(i);
							cs.push_back(c);
						}
					}
					int i0 = cs[0].val();
					int i1 = cs[1].val();
					int i2 = cs[2].val();
					int i3 = cs[3].val();
					map<int,int> corr;
					map<int,int> corr2;
					vector<int> ps;
					
					for(int i=0; i<4; i++)
					{
						Config c = cs[i];
						int k=0;
						for(int j=0; j<3; j++)
						{
							int v = c.inv(j).val();
							if(v!=i0 && v!=i1 && v!=i2 && v!=i3)
							{
								if(corr.find(v)!=corr.end())
								{
									corr2[v]=i;
								}
								else
								{
									corr[v]=i;
									ps.push_back(v);
								}
							}
						}
					}
					int a0, b0, c0;
					int a1, b1, c1;
					int a2, b2, c2;
					c0 = ps[0];
					b0 = corr[c0];
					a0 = corr2[c0];
					b0 = cs[b0].val();
					a0 = cs[a0].val();
					c1 = ps[1];
					b1 = corr[c1];
					a1 = corr2[c1];
					b1 = cs[b1].val();
					a1 = cs[a1].val();
					c2 = ps[2];
					b2 = corr[c2];
					a2 = corr2[c2];
					b2 = cs[b2].val();
					a2 = cs[a2].val();
					int p0, p1, p2, pc01, pc12, pc20;
					p0 = a0;
					p1 = b0;
					pc01 = c0;
					if(a1==p1)
					{
						p2 = b1;
						pc12 = c1;
						pc20 = c2;
					}
					else
					{
						if(a1==p0)
						{
							p2 = b1;
							pc20 = c1;
							pc12 = c2;
						}
						else
						{
							p2 = a1;
							if(b1==p0)
							{
								pc20 = c1;
								pc12 = c2;
							}
							else
							{
								pc20 = c2;
								pc12 = c1;
							}
						}
					}
					Point P0, P1, P2, P010, P011, P120, P121, P200, P201;
					P0 = p[p0];
					P1 = p[p1];
					P2 = p[p2];
					P010 = (P0.add(p[pc01])).mult(0.5);
					P011 = (P1.add(p[pc01])).mult(0.5);
					P120 = (P1.add(p[pc12])).mult(0.5);
					P121 = (P2.add(p[pc12])).mult(0.5);
					P200 = (P0.add(p[pc20])).mult(0.5);
					P201 = (P2.add(p[pc20])).mult(0.5);
					Triangle t0 = {P010,P011,P200};
					Triangle t1 = {P011,P200,P120};
					Triangle t2 = {P200,P120,P201};
					Triangle t3 = {P120,P201,P121};
					Zone z;
					z.faces.push_back(t0);
					z.faces.push_back(t1);
					z.faces.push_back(t2);
					z.faces.push_back(t3);
					vz.push_back(z);
				}
				break;
			case 9:
				{
					vector<Config> cs;
					for(int i=0; i<8; i++)
					{
						if(points[i])
						{
							Config c;
							c.max = 3;
							c = c.add(i);
							cs.push_back(c);
						}
					}
					int i0 = cs[0].val();
					int i1 = cs[1].val();
					int i2 = cs[2].val();
					int i3 = cs[3].val();
					map<int,int> corr;
					map<int,int> corr2;
					vector<int> ps;
					for(int i=0; i<4; i++)
					{
						Config c = cs[i];
						int k=0;
						for(int j=0; j<3; j++)
						{
							int v = c.inv(j).val();
							if(v!=i0 && v!=i1 && v!=i2 && v!=i3)
							{
								if(corr.find(i)!=corr.end())
								{
									corr2[i]=v;
								}
								else
								{
									corr[i]=v;
								}
							}
						}
					}
					bool a0 = false;
					bool b0 = false;
					Point A0, A1, B01, B02, B11, B12;
					Point PA0, PA1, PB0, PB1;
					int a1,a2,b01,b02,b11,b12;
					int oa0, oa1, ob0, ob1;
					for(int i=0; i<4; i++)
					{
						if(corr2.find(i)!=corr2.end())
						{
							int x0 = corr[i];
							int x1 = corr2[i];
							if(b0)
							{
								B11 = p[x0];
								b11 = x0;
								B12 = p[x1];
								b12 = x1;
								PB1 = p[cs[i].val()];
								ob1 = cs[i].val();
							}
							else
							{
								b0=true;
								B01 = p[x0];
								b01 = x0;
								B02 = p[x1];
								b02 = x1;
								PB0 = p[cs[i].val()];
								ob0 = cs[i].val();
							}
							
						}
						else
						{
							int x0 = corr[i];
							if(a0)
							{
								A1 = p[x0];
								a2 = x0;
								PA1 = p[cs[i].val()];
								oa1 = cs[i].val();
							}
							else
							{
								a0=true;
								A0 = p[x0];
								a1 = x0;
								PA0 = p[cs[i].val()];
								oa0 = cs[i].val();
							}
							
						}
					}
					A0 = A0.mult(0.5).add(PA0.mult(0.5));
					A1 = A1.mult(0.5).add(PA1.mult(0.5));
					B01 = B01.mult(0.5).add(PB0.mult(0.5));
					B02 = B02.mult(0.5).add(PB0.mult(0.5));
					B11 = B11.mult(0.5).add(PB1.mult(0.5));
					B12 = B12.mult(0.5).add(PB1.mult(0.5));
					Triangle t0 = {B01,B02,A0};
					Triangle t1 = {B11,B12,A0};
					Point P1;
					Point P2;
					if(a1==b01)
					{
						P1 = B02;
					}
					if(a1==b02)
					{
						P1 = B01;
					}
					if(a1==b11)
					{
						P1 = B12;
					}
					if(a1==b12)
					{
						P1 = B11;
					}
					if(a2==b01)
					{
						P2 = B01;
					}
					if(a2==b02)
					{
						P2 = B02;
					}
					if(a2==b11)
					{
						P2 = B11;
					}
					if(a2==b12)
					{
						P2 = B12;
					}
					Triangle t2 = {A0,P1,P2};
					Triangle t3 = {A1,P1,P2};
					Zone z;
					z.faces.push_back(t0);
					z.faces.push_back(t1);
					z.faces.push_back(t2);
					z.faces.push_back(t3);
					vz.push_back(z);
				}
				break;
		}
	}
	return vz;
}