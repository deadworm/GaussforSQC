#include<iostream>
#include<ctime>
#include<omp.h>
#include<fstream>
#include"normalization.h"
#include<math.h>
using namespace std;

double b=1E-05;	//�ضϾ���
int core;	//���ü����߳���
int n=3;	//ԭ�Ӹ���
char *atn;	//Ԫ����
double *qm;	//ԭ�Ӻ˵����

int g=10;	//��˹�����������������1������

int em=1;	//Ԫ��������
int *m;	//ͬһԪ�ػ���������(����Ҫ���룬�Ժ���Ҫ�������)
char *en;	//Ԫ����
int *om;	//������������
char **on;	//��������
int **og;	//����ĳһ���������˹�����ĸ���

int main()
{
	clock_t t1,t2;
	t1=clock();

	ifstream infile("infile.txt",ios::in);	//�Խضϡ�ԭ�Ӹ�����Ԫ�������˵������ԭ������Ľӿ�
	infile>>b;
	infile>>core;
	infile>>n;
	atn=new char[n];	//Ԫ����
	qm=new double[n];	//�˵����
	double *po=new double[3*n];	//ԭ������
	for (int i=0;i<n;i++)
	{
		infile>>atn[i];
		infile>>qm[i];
		infile>>po[i*3];
		infile>>po[i*3+1];
		infile>>po[i*3+2];
	}
	double *xc=new double[n];	//ԭ�������x����
	double *yc=new double[n];	//ԭ�������y����
	double *zc=new double[n];	//ԭ�������z����
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
	cout<<"��ȡinfile���"<<endl;
	infile.close();

	ifstream infile2("r_gauss.txt",ios::in);	//��1/r��˹չ���Ľӿ�
	double jifen0;
	infile2>>g;
	infile2>>jifen0;	//����չ���ĳ�����
	double *K=new double[g];	//����չ����ϵ��
	double *w=new double[g];	//����չ����ָ��
	for (int i=0;i<g;i++)
	{
		infile2>>K[i];
	}
	for (int i=0;i<g;i++)
	{
		infile2>>w[i];
	}
	infile2.close();
	cout<<"��ȡr_gauss���"<<endl;

	ifstream infile3("set_basis.txt",ios::in);	//�Ի���ϵ���Ľӿ�

	infile3>>em;	//Ԫ��������

	m=new int[em];	//����ͬһԪ�صĻ���������
	en=new char[em];	//����ÿ��Ԫ�ص�Ԫ����
	om=new int[em];	//����ÿ��Ԫ�صĹ��������
	on=new char*[em];	//��������
	og=new int*[em];	//������������˹�����ĸ���
	double ***ai=new double**[em];	//�����ȡ��ָ��ϵ��
	double ***ci=new double**[em];	//�����ȡ�ĳ���ϵ��
	double ***Ai=new double**[em];	//���������������ϵ��

	for(int i=0;i<em;i++)	//��Ԫ�ص�ѭ��
	{
		infile3>>en[i];
		infile3>>om[i];
		on[i]=new char[om[i]];
		og[i]=new int[om[i]];
		ai[i]=new double*[om[i]];
		ci[i]=new double*[om[i]];
		Ai[i]=new double*[om[i]];
		m[i]=0;
		for(int j=0;j<om[i];j++)	//�Թ����ѭ��
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
			for(int k=0;k<og[i][j];k++)	//�Թ����˹������ѭ��
			{
				infile3>>ai[i][j][k];
			}
			for(int k=0;k<og[i][j];k++)	//�Թ����˹������ѭ��
			{
				infile3>>ci[i][j][k];
			}
			for(int k=0;k<og[i][j];k++)	//�Թ����˹������ѭ��
			{
				Ai[i][j][k]=ci[i][j][k]*pow(2*ai[i][j][k]/pi,0.75);	//�Զ�ȡ�Ļ���ϵ��������
			}
		}
	}
	infile3.close();
	cout<<"��ȡset_basis���"<<endl;

	int aom=0;	//����ϵͳ����Ҫ���ܹ����
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
	cout<<"�ܹ����Ϊ��"<<aom<<endl;
	double *co=new double[aom*aom];	//�ܵ�������һϵ������
	gaussgroup **G;
	G=new gaussgroup*[aom];
	int *l;	//�������ÿ���������ĸ�˹��������
	l=new int[aom];
	int ao=0;	//�����������й��
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<em;j++)	//jΪԪ��ָ��
		{
			if(atn[i]==en[j])
			{
				for(int k=0;k<om[j];k++)	//kΪ���ָ��
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
	cout<<"�ܱ����Ĺ����Ϊ��"<<ao<<endl;

	co=normalization(G,l,aom);	//��������һ��ϵ��
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
			cout<<"��"<<i<<"�е�"<<j<<"��������һϵ��Ϊ��"<<co[i*aom+j]<<endl;
		}
	}


	int gm=0;	//�����ܵĸ�˹��������
	for(int i=0;i<aom;i++)
	{
		gm=gm+l[i];
	}
	cout<<"�ܸ�˹��������Ϊ��"<<gm<<endl;
	gaussgroup **A=new gaussgroup*[aom];
	for(int i=0;i<aom;i++)
	{
		A[i]=new gaussgroup[gm];
	}
	for(int i=0;i<aom;i++)	//�������е�������һ���
	{
		int ph=0;	//����phi��������и�˹����
		for(int j=0;j<aom;j++)	//�������еĻ��������
		{
			for(int k=0;k<l[j];k++)	//�������еĸ�˹����
			{
				A[i][ph]=G[j][k]*co[i*aom+j];
				ph++;
			}
		}
		cout<<"phi����ĸ�˹��������Ϊ��"<<ph<<endl;
	}
	//����qm�Ƿ���ȷ
	for(int i=0;i<n;i++)
	{
		cout<<"qmΪ"<<qm[i]<<endl;
	}

	ofstream outfile3("Fbase.txt",ios::out|ios::app);
	outfile3<<"��ʾ�����˳��Ϊ���ܹ������ÿ�������˹����������A,a,b,c,x1,y1,z1,xm,ym,zm"<<endl;
	outfile3<<aom<<endl;
	outfile3<<gm<<endl;
	for(int i=0;i<aom;i++)
	{
		for(int j=0;j<gm;j++)
		outfile3<<A[i][j];
		outfile3<<"*******************************************************************"<<endl;
	}
	outfile3.close();

	double H0=0;	//ԭ�Ӻ˼�����
	for (int i=0;i<n-1;i++)
	{
		for (int j=i+1;j<n;j++)
		{
			H0=H0+(qm[i]*qm[j])/sqrt(pow(xc[i]-xc[j],2)+pow(yc[i]-yc[j],2)+pow(zc[i]-zc[j],2));
		}
	}
	cout<<"ԭ�Ӻ˼����ܣ�"<<H0<<endl;

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
	cout<<"��������ʱ�䣺"<<(t2-t1)/CLK_TCK<<'s'<<endl;

	double *sum=new double[aom*aom];
	#pragma omp parallel for num_threads(core)
	for (int s1=0;s1<aom;s1++)
	{
		for (int s2=0;s2<aom;s2++)
		{
			double gycs=0;	//�����һ���Ա���
			sum[s1*aom+s2]=0;
			for (int i=0;i<gm;i++)		//��˹�����
			{
				for (int j=0;j<gm;j++)		//��˹�����
				{
					gaussgroup guiyi1=A[s1][i];	
					gaussgroup guiyi2=A[s2][j];
					gycs=gycs+(guiyi1*guiyi2).jifen();
					for (int c1=0;c1<g;c1++)	//��˹չ����������
					{
						for (int c2=0;c2<n;c2++)		//�������ܵ�ԭ�Ӹ���
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
			cout<<"��"<<s1<<"��"<<s2<<":"<<sum[s1*aom+s2]<<"��һ���ԣ�"<<gycs<<endl;
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
	cout<<"��������ʱ�䣺"<<(t2-t1)/CLK_TCK<<'s'<<endl;
	
	double *sum2=new double[aom*aom];
	#pragma omp parallel for num_threads(core)
	for (int s1=0;s1<aom;s1++)
	{
		for (int s2=0;s2<aom;s2++)
		{
			double gycs=0;	//�����һ���Ա���
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
			cout<<"��"<<s1<<"��"<<s2<<":"<<sum2[s1*aom+s2]<<"��һ���ԣ�"<<gycs<<endl;
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
	cout<<"��������ʱ�䣺"<<(t2-t1)/CLK_TCK<<'s'<<endl;
	

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
			for (int s2=0;s2<=s1;s2++)	//���ǶԳ���ʹs2<=s1
			{
				for (int s4=0;s3!=s1&&s4<=s3||s3==s1&&s4<=s2;s4++)	//���ǶԳ���ʹs3!=s1&&s4<=s3��s3==s1&&s4<=s2
				{
					double gycs=0;	//�����һ���Ա���
					sum3[(s1*aom+s3)*aom*aom+s2*aom+s4]=0;
					/*if (((s1==s3&&s2==s4)||(s1==s2&&s3==s4))==0)
						continue;	//���������˫�����໥����ϵ��*/
					for (int i=0;i<gm;i++)
					{
						for (int j=0;j<gm;j++)
						{
							for (int k=0;k<gm;k++)
							{
								for (int l=0;l<gm;l++)
								{
									gaussgroup guiyi1=A[s1][i];	//������Ϊ1��2�ı��һ��������
									gaussgroup guiyi2=A[s2][j];
									gaussgroup guiyi3=A[s3][k];
									gaussgroup guiyi4=A[s4][l];
									gycs=gycs+(guiyi1*guiyi2).jifen()*(guiyi3*guiyi4).jifen();
									for (int c1=0;c1<g;c1++)
									{
										if (c1<g-1)
										{
											gaussgroup jifen1(K[c1]/w[c1]/sqrt(pi/2),2/w[c1]/w[c1],2/w[c1]/w[c1],2/w[c1]/w[c1],0,0,0,0,0,0);
											sum3[(s1*aom+s3)*aom*aom+s2*aom+s4]=sum3[(s1*aom+s3)*aom*aom+s2*aom+s4]+dgjifen(jifen1,(guiyi1*guiyi2),(guiyi3*guiyi4));	//�����е���x1Ϊ�����Ĺ���ĳ˻���%����
											sum3[(s2*aom+s3)*aom*aom+s1*aom+s4]=sum3[(s1*aom+s3)*aom*aom+s2*aom+s4];
											sum3[(s1*aom+s4)*aom*aom+s2*aom+s3]=sum3[(s1*aom+s3)*aom*aom+s2*aom+s4];
											sum3[(s2*aom+s4)*aom*aom+s1*aom+s3]=sum3[(s1*aom+s3)*aom*aom+s2*aom+s4];
											sum3[(s3*aom+s1)*aom*aom+s4*aom+s2]=sum3[(s1*aom+s3)*aom*aom+s2*aom+s4];
											sum3[(s4*aom+s1)*aom*aom+s3*aom+s2]=sum3[(s1*aom+s3)*aom*aom+s2*aom+s4];
											sum3[(s3*aom+s2)*aom*aom+s4*aom+s1]=sum3[(s1*aom+s3)*aom*aom+s2*aom+s4];
											sum3[(s4*aom+s2)*aom*aom+s3*aom+s1]=sum3[(s1*aom+s3)*aom*aom+s2*aom+s4];	//�Գ������
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
					cout<<"��"<<s1*aom+s3<<"��"<<s2*aom+s4<<":"<<sum3[(s1*aom+s3)*aom*aom+s2*aom+s4]<<"��һ���ԣ�"<<gycs<<endl;
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
	cout<<"��������ʱ�䣺"<<(t2-t1)/CLK_TCK<<'s'<<endl;

	return 0;
}