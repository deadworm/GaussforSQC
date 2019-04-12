#include<iostream>
#include<fstream>
#include<math.h>
using namespace std;
inline int factor1(int m)
{
	int n=1;
	for(int i=m;i>0;i--)
	{
		n=n*i;
	}
	return n;
}
inline int factor2(int m)
{
	int n=1;
	for(int i=m;i>0;i=i-2)
	{
		n=n*i;
	}
	return n;
}