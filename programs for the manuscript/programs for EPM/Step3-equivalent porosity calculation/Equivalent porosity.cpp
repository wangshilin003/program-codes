# include<iostream> //头文件
# include<cmath>
#include<math.h>
# include<cstdlib>
# include<iomanip>
# include<fstream>
# include<sstream>
# include<string>
# include <stdio.h>
# include <stdlib.h>


using namespace std;
const int N=2;//number of the permeability data

double f(double x);
double k[200];
double Time[N],kLBM[N], k1LBM[N], k2LBM[N], k3LBM[N], Biomass[N], Biomass1[N], Biomass2[N], Biomass3[N];
double epsilon[N];

int main() 
{
 
	int i;
    for( i=0;i<=199;i++)
		k[i]=f(i*1.0/200);

   int p2;
	
    ifstream fin1( "KfromLBM.dat");
    for(p2=0;p2<N;p2++)
       fin1>> Time[p2]>>kLBM[p2]>>  Biomass[p2] >> Biomass1[p2] >> Biomass2[p2] >> Biomass3[p2];
	for (int p1 = 0; p1 < N; p1++)
		cout << Time[p1]<< "   "<<kLBM[p1]<<"     " <<Biomass[p1]<<"   " << Biomass1[p1] << "   " << Biomass2[p1] <<"   " << Biomass3[p1]<<endl;
	
 
   epsilon[0]=0;
  
   for(int j=0;j<N;j++)
   {
   int go=1;
    i=0;  
	do 
	{ i=i+1; 
	if (k[i]>kLBM[j]) {go=0; 
	             epsilon[j]=((k[i]-kLBM[j])*(i-1)*1.0/200+(kLBM[j]-k[i-1])*(i)*1.0/200)/(k[i]-k[i-1]);//线性插值
				}

	}
	while(go==1 && i<=199);
   }
	ofstream   out1( "Equivalent porosity.dat",ios::app);
    for( i=0;i<N;i++)
	out1<<Time[i]<<"   "<<epsilon[i] <<  endl;
	for (i = 0; i < N; i++)
		cout <<Time[i] << "   "<<epsilon[i] << endl;
 
  return 0;

}

double f(double x)
{  double d=1.95E-4,c=1.0/150;
   double k=pow(x,3) / pow(1-x,2)*d*d*c;
   return k;
}