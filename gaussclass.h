#include<iostream>
#include<fstream>
#include<vector>
#include"factorial.h"
#include<math.h>
using namespace std;
double pi=3.1415926535;
double e=2.718281828459;
double cut=1E-10;

class gaussjifen
{
public:
	double A,a,b,c;
	double x1,y1,z1;
	int xm,ym,zm;
	void out()
	{
		cout<<A<<','<<a<<','<<b<<','<<c<<','<<x1<<','<<y1<<','<<z1<<','<<xm<<','<<ym<<','<<zm<<endl;
	}
	gaussjifen()
	{
		A=0;
		a=0,b=0,c=0;
		x1=0,y1=0,z1=0;
		xm=0,ym=0,zm=0;
	}
	gaussjifen(double ad,double a1,double b1,double c1,double x,double y,double z,int xn,int yn,int zn)
	{
		A=ad;
		a=a1,b=b1,c=c1;
		x1=x,y1=y,z1=z;
		xm=xn,ym=yn,zm=zn;
	}
	double jifen()	//对单高斯函数求积分
	{
		if (abs(A)<cut)
		{
			return 0;
		}
		double const si=4;
		double Rx=0;
		for(int l=0;l<=xm/2;l++)
		{
			Rx=Rx+factor1(xm)*pow(pi,0.5)/factor1(l)/factor1(xm-2*l)/pow(si,l)*pow(x1,xm-2*l)*pow(a,-0.5-l);
		}
		double Ry=0;
		for(int l=0;l<=ym/2;l++)
		{
			Ry=Ry+factor1(ym)*pow(pi,0.5)/factor1(l)/factor1(ym-2*l)/pow(si,l)*pow(y1,ym-2*l)*pow(b,-0.5-l);
		}
		double Rz=0;
		for(int l=0;l<=zm/2;l++)
		{
			Rz=Rz+factor1(zm)*pow(pi,0.5)/factor1(l)/factor1(zm-2*l)/pow(si,l)*pow(z1,zm-2*l)*pow(c,-0.5-l);
		}
		return A*Rx*Ry*Rz;
	}
	double jifenx()	//对单高斯函数的x部分求积分
	{
		double const si=4;
		double Rx=0;
		for(int l=0;l<=xm/2;l++)
		{
			Rx=Rx+factor1(xm)*pow(pi,0.5)/factor1(l)/factor1(xm-2*l)/pow(si,l)*pow(x1,xm-2*l)*pow(a,-0.5-l);
		}
		return A*Rx;
	}
	double jifeny()	//对单高斯函数y部分求积分
	{
		double const si=4;
		double Ry=0;
		for(int l=0;l<=ym/2;l++)
		{
			Ry=Ry+factor1(ym)*pow(pi,0.5)/factor1(l)/factor1(ym-2*l)/pow(si,l)*pow(y1,ym-2*l)*pow(b,-0.5-l);
		}
		return A*Ry;
	}
	double jifenz()	//对单高斯函数z部分求积分
	{
		double const si=4;
		double Rz=0;
		for(int l=0;l<=zm/2;l++)
		{
			Rz=Rz+factor1(zm)*pow(pi,0.5)/factor1(l)/factor1(zm-2*l)/pow(si,l)*pow(z1,zm-2*l)*pow(c,-0.5-l);
		}
		return A*Rz;
	}
	double jifent(double a1c,double b1c,double a2c,double b2c,double a3c,double b3c,int n1,int n2,int n3)	//对高斯函数乘一自变量的线性函数的幂次方做积分
	{
		double Rx=0;
		for(int jx=0;jx<=xm;jx++)
		{
			for(int ix=0;ix<=n1;ix++)
			{
				gaussjifen jifenx(1,a,0,0,0,0,0,ix+jx,0,0);
				Rx=Rx+factor1(xm)/factor1(jx)/factor1(xm-jx)*factor1(n1)/factor1(ix)/factor1(n1-ix)*pow(a1c*x1+b1c,n1-ix)*pow(x1,xm-jx)*pow(a1c,ix)*jifenx.jifenx();
			}
		}
		double Ry=0;
		for(int jy=0;jy<=ym;jy++)
		{
			for(int iy=0;iy<=n2;iy++)
			{
				gaussjifen jifeny(1,0,b,0,0,0,0,0,iy+jy,0);
				Ry=Ry+factor1(ym)/factor1(jy)/factor1(ym-jy)*factor1(n2)/factor1(iy)/factor1(n2-iy)*pow(a2c*y1+b2c,n2-iy)*pow(y1,ym-jy)*pow(a2c,iy)*jifeny.jifeny();
			}
		}
		double Rz=0;
		for(int jz=0;jz<=zm;jz++)
		{
			for(int iz=0;iz<=n3;iz++)
			{	gaussjifen jifenz(1,0,0,c,0,0,0,0,0,iz+jz);
				Rz=Rz+factor1(zm)/factor1(jz)/factor1(zm-jz)*factor1(n3)/factor1(iz)/factor1(n3-iz)*pow(a3c*z1+b3c,n3-iz)*pow(z1,zm-jz)*pow(a3c,iz)*jifenz.jifenz();
			}
		}
		return Rx*Ry*Rz*A;
	}
	gaussjifen operator*(gaussjifen &chengji)	//对两高斯函数求乘积
	{
		double a1,b1,c1;
		double B,a2,b2,c2;
		double x2,y2,z2;
		double C,a3,b3,c3;
		double x3,y3,z3;
		int xm3,ym3,zm3;
		a1=a,a2=chengji.a;
		b1=b,b2=chengji.b;
		c1=c,c2=chengji.c;
		B=chengji.A,x2=chengji.x1,y2=chengji.y1,z2=chengji.z1;

		C=A*B*exp(-(a1*a2*(x1-x2)*(x1-x2)/(a1+a2)+b1*b2*(y1-y2)*(y1-y2)/(b1+b2)+c1*c2*(z1-z2)*(z1-z2)/(c1+c2)));
		a3=a1+a2;
		b3=b1+b2;
		c3=c1+c2;
		x3=(a1*x1+a2*x2)/(a1+a2);
		y3=(b1*y1+b2*y2)/(b1+b2);
		z3=(c1*z1+c2*z2)/(c1+c2);
		xm3=xm+chengji.xm;
		ym3=ym+chengji.ym;
		zm3=zm+chengji.zm;
		
		gaussjifen C3(C,a3,b3,c3,x3,y3,z3,xm3,ym3,zm3);
		return C3;
	}
	gaussjifen operator*(double ch)	//对高斯函数和常数求乘积
	{
		double C;
		C=A*ch;
		gaussjifen C3(C,a,b,c,x1,y1,z1,xm,ym,zm);
		return C3;
	}
	double operator/(gaussjifen &chengji)	//对高斯函数乘另一高斯函数的二阶导数(dx^2+dy^2+dz^2)求积分
	{
		if (abs(A)<cut)
		{
			return 0;
		}
		if (abs(chengji.A)<cut)
		{
			return 0;
		}
		double a1,b1,c1;
		double B,a2,b2,c2;
		double x2,y2,z2;
		a1=a,a2=chengji.a;
		b1=b,b2=chengji.b;
		c1=c,c2=chengji.c;
		B=chengji.A,x2=chengji.x1,y2=chengji.y1,z2=chengji.z1;
		double Qx[5]={4*a2*a2,-8*a2*a2*x2,4*a2*a2*x2*x2-4*a2*chengji.xm-2*a2,4*a2*x2*chengji.xm,chengji.xm*(chengji.xm-1)};
		double Qy[5]={4*b2*b2,-8*b2*b2*y2,4*b2*b2*y2*y2-4*b2*chengji.ym-2*b2,4*b2*y2*chengji.ym,chengji.ym*(chengji.ym-1)};
		double Qz[5]={4*c2*c2,-8*c2*c2*z2,4*c2*c2*z2*z2-4*c2*chengji.zm-2*c2,4*c2*z2*chengji.zm,chengji.zm*(chengji.zm-1)};
		int ic=0;
		int jc=0;
		int kc=0;
		double sum=0;
		for(int i=2+chengji.xm;i>=0;i--)
		{
			gaussjifen g1(A,a1,b1,c1,x1,y1,z1,xm,ym,zm);
			gaussjifen g2(B*Qx[ic],a2,b2,c2,x2,y2,z2,i,chengji.ym,chengji.zm);
			sum=sum+(g1*g2).jifen();
			ic++;
		}
		for(int j=2+chengji.ym;j>=0;j--)
		{
			gaussjifen g1(A,a1,b1,c1,x1,y1,z1,xm,ym,zm);
			gaussjifen g2(B*Qy[jc],a2,b2,c2,x2,y2,z2,chengji.xm,j,chengji.zm);
			sum=sum+(g1*g2).jifen();
			jc++;
		}
		for(int k=2+chengji.zm;k>=0;k--)
		{
			gaussjifen g1(A,a1,b1,c1,x1,y1,z1,xm,ym,zm);
			gaussjifen g2(B*Qz[kc],a2,b2,c2,x2,y2,z2,chengji.xm,chengji.ym,k);
			sum=sum+(g1*g2).jifen();
			kc++;
		}

		return sum;
	}
};

class gaussgroup
{
public:
	int n;
	vector<gaussjifen> G;

	gaussgroup(){};

	gaussgroup(int ng)
	{
		n=ng;
		G.reserve(n);
	}

	void out()
	{
		for (int i=0;i<n;i++)
		{
			G[i].out();
		}
	}

	gaussgroup (double ad,double a1,double b1,double c1,double x,double y,double z,int xn,int yn,int zn)
	{
		n=(xn+1)*(yn+1)*(zn+1);
		G.reserve(n);
		for (int i=0;i<xn+1;i++)
			for (int j=0;j<yn+1;j++)
				for (int k=0;k<zn+1;k++)
				{
					G.push_back(gaussjifen((factor1(xn)/factor1(i)/factor1(xn-i)*pow(-x,xn-i))*(factor1(yn)/factor1(j)/factor1(yn-j)*pow(-y,yn-j))*(factor1(zn)/factor1(k)/factor1(zn-k)*pow(-z,zn-k))*ad,a1,b1,c1,x,y,z,i,j,k));
				}
	}

	double jifen()	//对单高斯函数求积分
	{
		double sum=0;
		for (int i=0;i<n;i++)
		{
			sum=sum+G[i].jifen();
		}
		return sum;
	}

	gaussgroup operator*(double ch)	//对高斯函数和常数求乘积
	{
		gaussgroup G2(n);
		for (int i=0;i<n;i++)
		{
			G2.G.push_back(G[i]*ch);
		}
		return G2;
	}

	gaussgroup operator*(gaussgroup &chengji)	//对两高斯函数求乘积
	{
		gaussgroup G2(n*chengji.n);
		for (int i=0;i<n;i++)
		{
			for (int j=0;j<chengji.n;j++)
			{
				G2.G.push_back(G[i]*chengji.G[j]);
			}
		}
		return G2;
	}

	double operator/(gaussgroup &chengji)	//对高斯函数乘另一高斯函数的二阶导数(dx^2+dy^2+dz^2)求积分
	{
		double sum=0;
		for (int i=0;i<n;i++)
		{
			for (int j=0;j<chengji.n;j++)
			{
				sum=sum+G[i]/(chengji.G)[j];
			}
		}
		return sum;
	}
};

	double djifen(gaussjifen &j,gaussjifen &g1,gaussjifen &g2)	//对两高斯函数与无变量系数的高斯函数乘积求积分
	{
		if (abs(g1.A)<cut)
		{
			return 0;
		}
		if (abs(g2.A)<cut)
		{
			return 0;
		}
		double A,a1,b1,c1;
		double x1,y1,z1;
		int xm1,ym1,zm1;
		double B,a2,b2,c2;
		double x2,y2,z2;
		int xm2,ym2,zm2;
		double C,a3,b3,c3;
		double x3,y3,z3;
		a1=g1.a,a2=g2.a;
		b1=g1.b,b2=g2.b;
		c1=g1.c,c2=g2.c;
		A=g1.A,x1=g1.x1,y1=g1.y1,z1=g1.z1;
		B=g2.A,x2=g2.x1,y2=g2.y1,z2=g2.z1;
		xm1=g1.xm,ym1=g1.ym,zm1=g1.zm;
		xm2=g2.xm,ym2=g2.ym,zm2=g2.zm;

		C=j.A,a3=j.a,b3=j.b,c3=j.c,x3=j.x1,y3=j.y1,z3=j.z1;
		double const si=4;
		double R=0;
		for(int i=0;i<=xm1/2;i++)
		{
			for(int j=0;j<=ym1/2;j++)
			{
				for(int k=0;k<=zm1/2;k++)
					{
						gaussjifen gA(1,a3*a1/(a3+a1),b3*b1/(b3+b1),c3*c1/(c3+c1),x1,y1,z1,0,0,0);
						gaussjifen gB(1,a2,b2,c2,x2,y2,z2,xm2,ym2,zm2);
						R=R+(pow(pi,0.5)*factor1(xm1)/factor1(i)/factor1(xm1-2*i)/pow(si,i)*pow(a3+a1,-0.5-i))*(pow(pi,0.5)*factor1(ym1)/factor1(j)/factor1(ym1-2*j)/pow(si,j)*pow(b3+b1,-0.5-j))*(pow(pi,0.5)*factor1(zm1)/factor1(k)/factor1(zm1-2*k)/pow(si,k)*pow(c3+c1,-0.5-k))*(gA*gB).jifent(a3/(a3+a1),a1/(a3+a1)*x1,b3/(b3+b1),b1/(b3+b1)*y1,c3/(c3+c1),c1/(c3+c1)*z1,xm1-2*i,ym1-2*j,zm1-2*k);
					}
			}
		}
		return A*B*C*R;
	}

	double dgjifen(gaussgroup &jg,gaussgroup &g1,gaussgroup &g2)	//对两高斯函数与无变量系数的高斯函数乘积求积分
	{
		double sum=0;
		for (int i=0;i<jg.n;i++)
			for (int j=0;j<g1.n;j++)
				for (int k=0;k<g2.n;k++)
				{
					sum=sum+djifen(jg.G[i],g1.G[j],g2.G[k]);
				}
		return sum;
	}

	double dsjifen(gaussjifen &j,gaussjifen &g1,gaussjifen &g2)	//对两高斯函数与无变量系数的高斯函数乘积求积分（轨道为球对称高斯函数）
	{
		double A,a1,b1,c1;
		double x1,y1,z1;
		int xm1,ym1,zm1;
		double B,a2,b2,c2;
		double x2,y2,z2;
		int xm2,ym2,zm2;
		double C,a3,b3,c3;
		double x3,y3,z3;
		a1=g1.a,a2=g2.a;
		b1=g1.b,b2=g2.b;
		c1=g1.c,c2=g2.c;
		A=g1.A,x1=g1.x1,y1=g1.y1,z1=g1.z1;
		B=g2.A,x2=g2.x1,y2=g2.y1,z2=g2.z1;
		xm1=g1.xm,ym1=g1.ym,zm1=g1.zm;
		xm2=g2.xm,ym2=g2.ym,zm2=g2.zm;

		C=j.A,a3=j.a,b3=j.b,c3=j.c,x3=j.x1,y3=j.y1,z3=j.z1;

		double Q1,Q2,Q3,Q;
		Q1=pi/pow(a1*a2+a1*a3+a2*a3,0.5)*exp(-a1*a2*a3*pow(x1-x2,2)/(a1*a2+a1*a3+a2*a3));
		Q2=pi/pow(b1*b2+b1*b3+b2*b3,0.5)*exp(-b1*b2*b3*pow(y1-y2,2)/(b1*b2+b1*b3+b2*b3));
		Q3=pi/pow(c1*c2+c1*c3+c2*c3,0.5)*exp(-c1*c2*c3*pow(z1-z2,2)/(c1*c2+c1*c3+c2*c3));
		Q=A*B*C*Q1*Q2*Q3;

		return Q;
	}
	ostream& operator<<(ostream& out, gaussjifen& gs)	//输出高斯函数的组成信息
	{
		out<<gs.A<<'\t'<<gs.a<<'\t'<<gs.b<<'\t'<<gs.c<<'\t'<<gs.x1<<'\t'<<gs.y1<<'\t'<<gs.z1<<'\t'<<gs.xm<<'\t'<<gs.ym<<'\t'<<gs.zm<<endl;
		return out;
	}
	ostream& operator<<(ostream& out, gaussgroup& gs)	//输出高斯函数的组成信息
	{
		for (int i=0;i<gs.n;i++)
		{
			out<<gs.G[i].A<<'\t'<<gs.G[i].a<<'\t'<<gs.G[i].b<<'\t'<<gs.G[i].c<<'\t'<<gs.G[i].x1<<'\t'<<gs.G[i].y1<<'\t'<<gs.G[i].z1<<'\t'<<gs.G[i].xm<<'\t'<<gs.G[i].ym<<'\t'<<gs.G[i].zm<<endl;
		}
		return out;
	}