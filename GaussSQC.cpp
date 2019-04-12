#include<iostream>
#include<ctime>
#include<omp.h>
#include<fstream>
#include"normalization.h"
#include<math.h>
using namespace std;

double b=1E-05;	//截断精度
int core;	//设置计算线程数
int n=3;	//原子个数
char *atn;	//元素名
double *qm;	//原子核电荷数

int g=10;	//高斯函数个数（不计入加1常数）

int em=1;	//元素种类数
int *m;	//同一元素基组轨道总数(不需要输入，以后需要计算出来)
char *en;	//元素名
int *om;	//基组轨道类型数
char **on;	//基组轨道名
int **og;	//基组某一轨道包含高斯函数的个数

int main()
{
	clock_t t1,t2;
	t1=clock();

	ifstream infile("infile.txt",ios::in);	//对截断、原子个数、元素名、核电荷数、原子坐标的接口
	infile>>b;
	infile>>core;
	infile>>n;
	atn=new char[n];	//元素名
	qm=new double[n];	//核电荷数
	double *po=new double[3*n];	//原子坐标
	for (int i=0;i<n;i++)
	{
		infile>>atn[i];
		infile>>qm[i];
		infile>>po[i*3];
		infile>>po[i*3+1];
		infile>>po[i*3+2];
	}
	double *xc=new double[n];	//原子坐标的x分量
	double *yc=new double[n];	//原子坐标的y分量
	double *zc=new double[n];	//原子坐标的z分量
	for (int i=0;i<n;i++)
	{
		xc[i]=po[i*3+0];
	}
	for (int i=0;i<n;i++)
	{
		yc[i]=po[i*3+1];
	}
	for (int i=0;i<n;i++)
	{
		zc[i]=po[i*3+2];
	}
	cout<<"读取infile完毕"<<endl;
	infile.close();

	ifstream infile2("r_gauss.txt",ios::in);	//对1/r高斯展开的接口
	double jifen0;
	infile2>>g;
	infile2>>jifen0;	//保存展开的常数项
	double *K=new double[g];	//保存展开的系数
	double *w=new double[g];	//保存展开的指数
	for (int i=0;i<g;i++)
	{
		infile2>>K[i];
	}
	for (int i=0;i<g;i++)
	{
		infile2>>w[i];
	}
	infile2.close();
	cout<<"读取r_gauss完毕"<<endl;

	ifstream infile3("set_basis.txt",ios::in);	//对基组系数的接口

	infile3>>em;	//元素种类数

	m=new int[em];	//保存同一元素的基组轨道总数
	en=new char[em];	//保存每个元素的元素名
	om=new int[em];	//保存每个元素的轨道类型数
	on=new char*[em];	//保存轨道名
	og=new int*[em];	//保存轨道包含高斯函数的个数
	double ***ai=new double**[em];	//保存读取的指数系数
	double ***ci=new double**[em];	//保存读取的初级系数
	double ***Ai=new double**[em];	//保存修正后的最终系数

	for(int i=0;i<em;i++)	//对元素的循环
	{
		infile3>>en[i];
		infile3>>om[i];
		on[i]=new char[om[i]];
		og[i]=new int[om[i]];
		ai[i]=new double*[om[i]];
		ci[i]=new double*[om[i]];
		Ai[i]=new double*[om[i]];
		m[i]=0;
		for(int j=0;j<om[i];j++)	//对轨道的循环
		{
			infile3>>on[i][j];
			infile3>>og[i][j];
			switch(on[i][j])
			{
				case 'S':m[i]=m[i]+1;break;
				case 'X':m[i]=m[i]+1;break;
				case 'x':m[i]=m[i]+1;break;
				case 'P':m[i]=m[i]+3;break;
				case 'D':m[i]=m[i]+5;break;
				case 'F':m[i]=m[i]+10;break;
			}
			ai[i][j]=new double[og[i][j]];
			ci[i][j]=new double[og[i][j]];
			Ai[i][j]=new double[og[i][j]];
			for(int k=0;k<og[i][j];k++)	//对轨道高斯函数的循环
			{
				infile3>>ai[i][j][k];
			}
			for(int k=0;k<og[i][j];k++)	//对轨道高斯函数的循环
			{
				infile3>>ci[i][j][k];
			}
			for(int k=0;k<og[i][j];k++)	//对轨道高斯函数的循环
			{
				Ai[i][j][k]=ci[i][j][k]*pow(2*ai[i][j][k]/pi,0.75);	//对读取的基组系数做修正
			}
		}
	}
	infile3.close();
	cout<<"读取set_basis完毕"<<endl;

	int aom=0;	//保存系统所需要的总轨道数
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<em;j++)
		{
			if(atn[i]==en[j])
			{
				aom=aom+m[j];
			}
		}
	}
	cout<<"总轨道数为："<<aom<<endl;
	double *co=new double[aom*aom];	//总的正交归一系数矩阵
	gaussgroup **G;
	G=new gaussgroup*[aom];
	int *l;	//保存组成每个基函数的高斯函数个数
	l=new int[aom];
	int ao=0;	//用来遍历所有轨道
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<em;j++)	//j为元素指标
		{
			if(atn[i]==en[j])
			{
				for(int k=0;k<om[j];k++)	//k为轨道指标
				{
					switch(on[j][k])
					{
						case 'S':
							{
								G[ao]=new gaussgroup[og[j][k]];
								l[ao]=og[j][k];
								for(int p=0;p<og[j][k];p++)
								{
									G[ao][p]=gaussgroup(Ai[j][k][p],ai[j][k][p],ai[j][k][p],ai[j][k][p],xc[i],yc[i],zc[i],0,0,0);
								}
								ao++;
								break;
							}
						case 'X':
							{
								G[ao]=new gaussgroup[og[j][k]];
								l[ao]=og[j][k];
								for(int p=0;p<og[j][k];p++)
								{
									G[ao][p]=gaussgroup(Ai[j][k][p],ai[j][k][p],ai[j][k][p],ai[j][k][p],xc[i],yc[i],zc[i],1,0,0);
								}
								ao++;
								break;
							}
						case 'x':
							{
								G[ao]=new gaussgroup[og[j][k]];
								l[ao]=og[j][k];
								for(int p=0;p<og[j][k];p++)
								{
									G[ao][p]=gaussgroup(Ai[j][k][p],ai[j][k][p],ai[j][k][p],ai[j][k][p],xc[i],yc[i],zc[i],0,1,0);
								}
								ao++;
								break;
							}
						case 'P':
							{
								G[ao]=new gaussgroup[og[j][k]];
								l[ao]=og[j][k];
								for(int p=0;p<og[j][k];p++)
								{
									G[ao][p]=gaussgroup(Ai[j][k][p],ai[j][k][p],ai[j][k][p],ai[j][k][p],xc[i],yc[i],zc[i],1,0,0);
								}
								ao++;
								G[ao]=new gaussgroup[og[j][k]];
								l[ao]=og[j][k];
								for(int p=0;p<og[j][k];p++)
								{
									G[ao][p]=gaussgroup(Ai[j][k][p],ai[j][k][p],ai[j][k][p],ai[j][k][p],xc[i],yc[i],zc[i],0,1,0);
								}
								ao++;
								G[ao]=new gaussgroup[og[j][k]];
								l[ao]=og[j][k];
								for(int p=0;p<og[j][k];p++)
								{
									G[ao][p]=gaussgroup(Ai[j][k][p],ai[j][k][p],ai[j][k][p],ai[j][k][p],xc[i],yc[i],zc[i],0,0,1);
								}
								ao++;
								break;
							}
						case 'D':
							{
								G[ao]=new gaussgroup[og[j][k]*2];
								l[ao]=og[j][k]*2;
								for(int p=0;p<og[j][k];p++)
								{
									G[ao][p]=gaussgroup(Ai[j][k][p],ai[j][k][p],ai[j][k][p],ai[j][k][p],xc[i],yc[i],zc[i],2,0,0);
								}
								for(int p=og[j][k];p<og[j][k]*2;p++)
								{
									G[ao][p]=gaussgroup(-Ai[j][k][p-og[j][k]],ai[j][k][p-og[j][k]],ai[j][k][p-og[j][k]],ai[j][k][p-og[j][k]],xc[i],yc[i],zc[i],0,2,0);
								}
								ao++;
								G[ao]=new gaussgroup[og[j][k]*3];
								l[ao]=og[j][k]*3;
								for(int p=0;p<og[j][k];p++)
								{
									G[ao][p]=gaussgroup(2*Ai[j][k][p],ai[j][k][p],ai[j][k][p],ai[j][k][p],xc[i],yc[i],zc[i],0,0,2);
								}
								for(int p=og[j][k];p<og[j][k]*2;p++)
								{
									G[ao][p]=gaussgroup(-Ai[j][k][p-og[j][k]],ai[j][k][p-og[j][k]],ai[j][k][p-og[j][k]],ai[j][k][p-og[j][k]],xc[i],yc[i],zc[i],2,0,0);
								}
								for(int p=og[j][k]*2;p<og[j][k]*3;p++)
								{
									G[ao][p]=gaussgroup(-Ai[j][k][p-og[j][k]*2],ai[j][k][p-og[j][k]*2],ai[j][k][p-og[j][k]*2],ai[j][k][p-og[j][k]*2],xc[i],yc[i],zc[i],0,2,0);
								}
								ao++;
								G[ao]=new gaussgroup[og[j][k]];
								l[ao]=og[j][k];
								for(int p=0;p<og[j][k];p++)
								{
									G[ao][p]=gaussgroup(Ai[j][k][p],ai[j][k][p],ai[j][k][p],ai[j][k][p],xc[i],yc[i],zc[i],1,1,0);
								}
								ao++;
								G[ao]=new gaussgroup[og[j][k]];
								l[ao]=og[j][k];
								for(int p=0;p<og[j][k];p++)
								{
									G[ao][p]=gaussgroup(Ai[j][k][p],ai[j][k][p],ai[j][k][p],ai[j][k][p],xc[i],yc[i],zc[i],1,0,1);
								}
								ao++;
								G[ao]=new gaussgroup[og[j][k]];
								l[ao]=og[j][k];
								for(int p=0;p<og[j][k];p++)
								{
									G[ao][p]=gaussgroup(Ai[j][k][p],ai[j][k][p],ai[j][k][p],ai[j][k][p],xc[i],yc[i],zc[i],0,1,1);
								}
								ao++;
								break;
							}
						case 'F':
							{
								G[ao]=new gaussgroup[og[j][k]];
								l[ao]=og[j][k];
								for(int p=0;p<og[j][k];p++)
								{
									G[ao][p]=gaussgroup(Ai[j][k][p],ai[j][k][p],ai[j][k][p],ai[j][k][p],xc[i],yc[i],zc[i],3,0,0);
								}
								ao++;
								G[ao]=new gaussgroup[og[j][k]];
								l[ao]=og[j][k];
								for(int p=0;p<og[j][k];p++)
								{
									G[ao][p]=gaussgroup(Ai[j][k][p],ai[j][k][p],ai[j][k][p],ai[j][k][p],xc[i],yc[i],zc[i],0,3,0);
								}
								ao++;
								G[ao]=new gaussgroup[og[j][k]];
								l[ao]=og[j][k];
								for(int p=0;p<og[j][k];p++)
								{
									G[ao][p]=gaussgroup(Ai[j][k][p],ai[j][k][p],ai[j][k][p],ai[j][k][p],xc[i],yc[i],zc[i],0,0,3);
								}
								ao++;
								G[ao]=new gaussgroup[og[j][k]];
								l[ao]=og[j][k];
								for(int p=0;p<og[j][k];p++)
								{
									G[ao][p]=gaussgroup(Ai[j][k][p],ai[j][k][p],ai[j][k][p],ai[j][k][p],xc[i],yc[i],zc[i],2,1,0);
								}
								ao++;
								G[ao]=new gaussgroup[og[j][k]];
								l[ao]=og[j][k];
								for(int p=0;p<og[j][k];p++)
								{
									G[ao][p]=gaussgroup(Ai[j][k][p],ai[j][k][p],ai[j][k][p],ai[j][k][p],xc[i],yc[i],zc[i],1,2,0);
								}
								ao++;
								G[ao]=new gaussgroup[og[j][k]];
								l[ao]=og[j][k];
								for(int p=0;p<og[j][k];p++)
								{
									G[ao][p]=gaussgroup(Ai[j][k][p],ai[j][k][p],ai[j][k][p],ai[j][k][p],xc[i],yc[i],zc[i],2,0,1);
								}
								ao++;
								G[ao]=new gaussgroup[og[j][k]];
								l[ao]=og[j][k];
								for(int p=0;p<og[j][k];p++)
								{
									G[ao][p]=gaussgroup(Ai[j][k][p],ai[j][k][p],ai[j][k][p],ai[j][k][p],xc[i],yc[i],zc[i],1,0,2);
								}
								ao++;
								G[ao]=new gaussgroup[og[j][k]];
								l[ao]=og[j][k];
								for(int p=0;p<og[j][k];p++)
								{
									G[ao][p]=gaussgroup(Ai[j][k][p],ai[j][k][p],ai[j][k][p],ai[j][k][p],xc[i],yc[i],zc[i],0,2,1);
								}
								ao++;
								G[ao]=new gaussgroup[og[j][k]];
								l[ao]=og[j][k];
								for(int p=0;p<og[j][k];p++)
								{
									G[ao][p]=gaussgroup(Ai[j][k][p],ai[j][k][p],ai[j][k][p],ai[j][k][p],xc[i],yc[i],zc[i],0,1,2);
								}
								ao++;
								G[ao]=new gaussgroup[og[j][k]];
								l[ao]=og[j][k];
								for(int p=0;p<og[j][k];p++)
								{
									G[ao][p]=gaussgroup(Ai[j][k][p],ai[j][k][p],ai[j][k][p],ai[j][k][p],xc[i],yc[i],zc[i],1,1,1);
								}
								ao++;
								break;
							}
					}
				}
			}
		}
	}
	cout<<"总遍历的轨道数为："<<ao<<endl;

	co=normalization(G,l,aom);	//求正交归一化系数
	for(int i=0;i<aom;i++)
	{
		for(int j=i+1;j<aom;j++)
		{
			co[i*aom+j]=0;
		}
	}
	for(int i=0;i<aom;i++)
	{
		for(int j=0;j<aom;j++)
		{
			cout<<"第"<<i<<"行第"<<j<<"列正交归一系数为："<<co[i*aom+j]<<endl;
		}
	}


	int gm=0;	//保存总的高斯函数个数
	for(int i=0;i<aom;i++)
	{
		gm=gm+l[i];
	}
	cout<<"总高斯函数个数为："<<gm<<endl;
	gaussgroup **A=new gaussgroup*[aom];
	for(int i=0;i<aom;i++)
	{
		A[i]=new gaussgroup[gm];
	}
	for(int i=0;i<aom;i++)	//遍历所有的正交归一轨道
	{
		int ph=0;	//遍历phi轨道的所有高斯函数
		for(int j=0;j<aom;j++)	//遍历所有的基函数轨道
		{
			for(int k=0;k<l[j];k++)	//遍历所有的高斯函数
			{
				A[i][ph]=G[j][k]*co[i*aom+j];
				ph++;
			}
		}
		cout<<"phi轨道的高斯函数个数为："<<ph<<endl;
	}
	//计算qm是否正确
	for(int i=0;i<n;i++)
	{
		cout<<"qm为"<<qm[i]<<endl;
	}

	ofstream outfile3("Fbase.txt",ios::out|ios::app);
	outfile3<<"提示，输出顺序为：总轨道数，每个轨道高斯函数个数，A,a,b,c,x1,y1,z1,xm,ym,zm"<<endl;
	outfile3<<aom<<endl;
	outfile3<<gm<<endl;
	for(int i=0;i<aom;i++)
	{
		for(int j=0;j<gm;j++)
		outfile3<<A[i][j];
		outfile3<<"*******************************************************************"<<endl;
	}
	outfile3.close();

	double H0=0;	//原子核间势能
	for (int i=0;i<n-1;i++)
	{
		for (int j=i+1;j<n;j++)
		{
			H0=H0+(qm[i]*qm[j])/sqrt(pow(xc[i]-xc[j],2)+pow(yc[i]-yc[j],2)+pow(zc[i]-zc[j],2));
		}
	}
	cout<<"原子核间势能："<<H0<<endl;

	ofstream outfile0("Hv.txt",ios::out|ios::app);
	outfile0<<"nuclear distance"<<endl;
	for (int i=0;i<n-1;i++)
	{
		outfile0<<sqrt(pow(xc[i]-xc[i+1],2)+pow(yc[i]-yc[i+1],2)+pow(zc[i]-zc[i+1],2))<<','<<'\t';
	}
	outfile0<<endl<<endl;
	outfile0<<"nuclear potential energy"<<endl;
	outfile0<<H0<<','<<endl;
	outfile0<<endl;

	t2=clock();
	cout<<"程序运行时间："<<(t2-t1)/CLK_TCK<<'s'<<endl;

	double *sum=new double[aom*aom];
	#pragma omp parallel for num_threads(core)
	for (int s1=0;s1<aom;s1++)
	{
		for (int s2=0;s2<aom;s2++)
		{
			double gycs=0;	//定义归一测试变量
			sum[s1*aom+s2]=0;
			for (int i=0;i<gm;i++)		//高斯轨道数
			{
				for (int j=0;j<gm;j++)		//高斯轨道数
				{
					gaussgroup guiyi1=A[s1][i];	
					gaussgroup guiyi2=A[s2][j];
					gycs=gycs+(guiyi1*guiyi2).jifen();
					for (int c1=0;c1<g;c1++)	//高斯展开函数个数
					{
						for (int c2=0;c2<n;c2++)		//产生势能的原子个数
						{
							if (c1<g-1)
							{
								gaussgroup jifen1(K[c1]/w[c1]/sqrt(pi/2),2/w[c1]/w[c1],2/w[c1]/w[c1],2/w[c1]/w[c1],xc[c2],yc[c2],zc[c2],0,0,0);
								sum[s1*aom+s2]=sum[s1*aom+s2]-qm[c2]*(guiyi1*jifen1*guiyi2).jifen();
							}
							else
							{
								sum[s1*aom+s2]=sum[s1*aom+s2]-qm[c2]*jifen0*(guiyi1*guiyi2).jifen();
							}
						}
					}
				}
			}
			cout<<"行"<<s1<<"列"<<s2<<":"<<sum[s1*aom+s2]<<"归一测试："<<gycs<<endl;
			gycs=0;
		}
	}

	ofstream outfile("Hv.txt",ios::out|ios::app);
	outfile<<"potential energy matrix"<<endl;
	for (int i=0;i<aom;i++)
	{
		for (int j=0;j<aom;j++)
		{
			if ((abs(sum[i*aom+j])>b))
			outfile<<sum[i*aom+j]<<','<<'\t';
			else
			outfile<<0<<','<<'\t';
		}
		outfile<<endl;
	}
	outfile<<endl;

	t2=clock();
	cout<<"程序运行时间："<<(t2-t1)/CLK_TCK<<'s'<<endl;
	
	double *sum2=new double[aom*aom];
	#pragma omp parallel for num_threads(core)
	for (int s1=0;s1<aom;s1++)
	{
		for (int s2=0;s2<aom;s2++)
		{
			double gycs=0;	//定义归一测试变量
			sum2[s1*aom+s2]=0;
			for (int i=0;i<gm;i++)
			{
				for (int j=0;j<gm;j++)
				{
					gaussgroup guiyi1=A[s1][i];	
					gaussgroup guiyi2=A[s2][j];
					gycs=gycs+(guiyi1*guiyi2).jifen();	
					sum2[s1*aom+s2]=sum2[s1*aom+s2]+guiyi1/guiyi2;
				}
			}
			sum2[s1*aom+s2]=-sum2[s1*aom+s2]/2;
			cout<<"行"<<s1<<"列"<<s2<<":"<<sum2[s1*aom+s2]<<"归一测试："<<gycs<<endl;
			gycs=0;
		}
	}
	

	outfile<<"kinetic energy matrix"<<endl;
	for (int i=0;i<aom;i++)
	{
		for (int j=0;j<aom;j++)
		{
			if ((abs(sum2[i*aom+j])>b))
			outfile<<sum2[i*aom+j]<<','<<'\t';
			else
			outfile<<0<<','<<'\t';
		}
		outfile<<endl;
	}
	outfile<<endl;
	outfile<<"*******************************************************************"<<endl;
	outfile.close();

	t2=clock();
	cout<<"程序运行时间："<<(t2-t1)/CLK_TCK<<'s'<<endl;
	

	double *sum3=new double[aom*aom*aom*aom];
	int dl=aom*(aom+1)/2;
	double dcore=core;
	#pragma omp parallel for num_threads(core)
	for (int di=1;di<=core;di++)
	{
		int dm=ceil(dl*pow((di-1)/dcore,0.5))+1;
		for(;dm<=ceil(dl*pow(di/dcore,0.5))&&dm>ceil(dl*pow((di-1)/dcore,0.5));dm++)
		{
			int s1=0;
			int s3=0;
			int sm=1;
			int ddm=dm;
			while(ddm>sm)
			{
				ddm=ddm-sm;
				sm++;
				s1++;
			}
			s3=ddm-1;
			for (int s2=0;s2<=s1;s2++)	//考虑对称性使s2<=s1
			{
				for (int s4=0;s3!=s1&&s4<=s3||s3==s1&&s4<=s2;s4++)	//考虑对称性使s3!=s1&&s4<=s3或s3==s1&&s4<=s2
				{
					double gycs=0;	//定义归一测试变量
					sum3[(s1*aom+s3)*aom*aom+s2*aom+s4]=0;
					/*if (((s1==s3&&s2==s4)||(s1==s2&&s3==s4))==0)
						continue;	//调整计算的双电子相互作用系数*/
					for (int i=0;i<gm;i++)
					{
						for (int j=0;j<gm;j++)
						{
							for (int k=0;k<gm;k++)
							{
								for (int l=0;l<gm;l++)
								{
									gaussgroup guiyi1=A[s1][i];	//参数改为1或2改变归一化波函数
									gaussgroup guiyi2=A[s2][j];
									gaussgroup guiyi3=A[s3][k];
									gaussgroup guiyi4=A[s4][l];
									gycs=gycs+(guiyi1*guiyi2).jifen()*(guiyi3*guiyi4).jifen();
									for (int c1=0;c1<g;c1++)
									{
										if (c1<g-1)
										{
											gaussgroup jifen1(K[c1]/w[c1]/sqrt(pi/2),2/w[c1]/w[c1],2/w[c1]/w[c1],2/w[c1]/w[c1],0,0,0,0,0,0);
											sum3[(s1*aom+s3)*aom*aom+s2*aom+s4]=sum3[(s1*aom+s3)*aom*aom+s2*aom+s4]+dgjifen(jifen1,(guiyi1*guiyi2),(guiyi3*guiyi4));	//对所有的以x1为变量的轨道的乘积做%运算
											sum3[(s2*aom+s3)*aom*aom+s1*aom+s4]=sum3[(s1*aom+s3)*aom*aom+s2*aom+s4];
											sum3[(s1*aom+s4)*aom*aom+s2*aom+s3]=sum3[(s1*aom+s3)*aom*aom+s2*aom+s4];
											sum3[(s2*aom+s4)*aom*aom+s1*aom+s3]=sum3[(s1*aom+s3)*aom*aom+s2*aom+s4];
											sum3[(s3*aom+s1)*aom*aom+s4*aom+s2]=sum3[(s1*aom+s3)*aom*aom+s2*aom+s4];
											sum3[(s4*aom+s1)*aom*aom+s3*aom+s2]=sum3[(s1*aom+s3)*aom*aom+s2*aom+s4];
											sum3[(s3*aom+s2)*aom*aom+s4*aom+s1]=sum3[(s1*aom+s3)*aom*aom+s2*aom+s4];
											sum3[(s4*aom+s2)*aom*aom+s3*aom+s1]=sum3[(s1*aom+s3)*aom*aom+s2*aom+s4];	//对称性填充
										}
										else
										{
											sum3[(s1*aom+s3)*aom*aom+s2*aom+s4]=sum3[(s1*aom+s3)*aom*aom+s2*aom+s4]+jifen0*(guiyi1*guiyi2).jifen()*(guiyi3*guiyi4).jifen();
											sum3[(s2*aom+s3)*aom*aom+s1*aom+s4]=sum3[(s1*aom+s3)*aom*aom+s2*aom+s4];
											sum3[(s1*aom+s4)*aom*aom+s2*aom+s3]=sum3[(s1*aom+s3)*aom*aom+s2*aom+s4];
											sum3[(s2*aom+s4)*aom*aom+s1*aom+s3]=sum3[(s1*aom+s3)*aom*aom+s2*aom+s4];
											sum3[(s3*aom+s1)*aom*aom+s4*aom+s2]=sum3[(s1*aom+s3)*aom*aom+s2*aom+s4];
											sum3[(s4*aom+s1)*aom*aom+s3*aom+s2]=sum3[(s1*aom+s3)*aom*aom+s2*aom+s4];
											sum3[(s3*aom+s2)*aom*aom+s4*aom+s1]=sum3[(s1*aom+s3)*aom*aom+s2*aom+s4];
											sum3[(s4*aom+s2)*aom*aom+s3*aom+s1]=sum3[(s1*aom+s3)*aom*aom+s2*aom+s4];
										}
									}
								}
							}
						}
					}
					cout<<"行"<<s1*aom+s3<<"列"<<s2*aom+s4<<":"<<sum3[(s1*aom+s3)*aom*aom+s2*aom+s4]<<"归一测试："<<gycs<<endl;
					gycs=0;
				}
			}
		}
	}
	

	ofstream outfile2("Ht.txt",ios::out|ios::app);
	outfile2<<"potential energy of electron"<<endl;
	for (int i=0;i<aom*aom;i++)
	{
		for (int j=0;j<aom*aom;j++)
		{
			if (abs(sum3[i*aom*aom+j])>b)
			outfile2<<sum3[i*aom*aom+j]<<','<<'\t';
			else
			outfile2<<0<<','<<'\t';
		}
		outfile2<<endl;
	}
	outfile2<<endl;
	outfile2<<"*******************************************************************"<<endl;
	outfile2.close();

	delete []atn;delete []qm;
	delete []po;delete []zc;delete []co;
	delete []K;delete []w;
	delete []m;delete []en;delete []om;delete []on;delete []og;
	delete []ai;delete []ci;delete []Ai;
	for(int i=0;i<aom;i++)
	{
		delete []G[i];
		delete []A[i];
	}
	delete []l;
	delete []A;delete []G;
	delete []sum;delete []sum2;delete []sum3;


	t2=clock();
	cout<<"程序运行时间："<<(t2-t1)/CLK_TCK<<'s'<<endl;

	return 0;
}