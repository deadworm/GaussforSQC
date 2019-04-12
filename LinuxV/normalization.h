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

double* normalization(gaussgroup **A, int *l, int len)	//����len��ʾ������һ���Ļ�����������l��ʾÿ�������������ĸ�˹��������
{
	double *nz=new double[len*len];	//��¼��һ��ϵ��
	double *jf=new double[len*len];	//��¼֮ǰ��һ���ĺ����뵱ǰδ��һ���ĺ����Ļ���ֵ

	nz[0]=1/sqrt(norm(A[0],A[0],l[0],l[0]));
	for(int i=1;i<len;i++)
	{
		for(int k1=0;k1<i;k1++)	//������֮ǰ��һ���ĺ�����δ��һ���ĺ����Ļ���ֵ��k1����֮ǰ��һ���ĺ�����ָ��
		{
			double sum=0;
			for(int k=0;k<=k1;k++)	//k�����һ���������������ָ��
			{
				sum=sum+nz[k1*len+k]*norm(A[k],A[i],l[k],l[i]);
			}
			jf[i*len+k1]=sum;
		}
		for(int j=0;j<i;j++)	//��������������ϵ����j������ȫ��ֺ�ϵ����Ϊ1�Ļ�������ָ��
		{
			double sum=0;
			for(int k=j;k<i;k++)	//k����֮ǰ��һ���ĺ�����ָ��
			{
				sum=sum-nz[k*len+j]*jf[i*len+k];
			}
			nz[i*len+j]=sum;
		}
		nz[i*len+i]=1;
		double sum=0;
		for(int m=0;m<=i;m++)	//��һ��
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