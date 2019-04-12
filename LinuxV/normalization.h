#include<iostream>
#include<fstream>
#include"gaussclass.h"
#include<math.h>
using namespace std;

double norm(gaussgroup *A, gaussgroup *B, int la, int lb)
{
	double sum=0;
	for(int i=0;i<la;i++)
	{
		for(int j=0;j<lb;j++)
		{
			sum=sum+(A[i]*B[j]).jifen();
		}
	}
	return sum;
}

double* normalization(gaussgroup **A, int *l, int len)	//其中len表示正交归一化的基函数个数，l表示每个基函数包含的高斯函数个数
{
	double *nz=new double[len*len];	//记录归一化系数
	double *jf=new double[len*len];	//记录之前归一化的函数与当前未归一化的函数的积分值

	nz[0]=1/sqrt(norm(A[0],A[0],l[0],l[0]));
	for(int i=1;i<len;i++)
	{
		for(int k1=0;k1<i;k1++)	//用来求之前归一化的函数与未归一化的函数的积分值，k1代表之前归一化的函数的指标
		{
			double sum=0;
			for(int k=0;k<=k1;k++)	//k代表归一化函数里基函数的指标
			{
				sum=sum+nz[k1*len+k]*norm(A[k],A[i],l[k],l[i]);
			}
			jf[i*len+k1]=sum;
		}
		for(int j=0;j<i;j++)	//用来求正交化的系数，j代表完全拆分后，系数不为1的基函数的指标
		{
			double sum=0;
			for(int k=j;k<i;k++)	//k代表之前归一化的函数的指标
			{
				sum=sum-nz[k*len+j]*jf[i*len+k];
			}
			nz[i*len+j]=sum;
		}
		nz[i*len+i]=1;
		double sum=0;
		for(int m=0;m<=i;m++)	//归一化
		{
			for(int n=0;n<=i;n++)
			{
				sum=sum+nz[i*len+m]*nz[i*len+n]*norm(A[m],A[n],l[m],l[n]);
			}
		}
		for(int j=0;j<=i;j++)
		{
			nz[i*len+j]=nz[i*len+j]/sqrt(sum);
		}
	}
	return nz;
}