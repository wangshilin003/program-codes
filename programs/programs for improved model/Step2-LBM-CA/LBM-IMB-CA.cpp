///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 
//  
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


///This computer code is for LBM-IMB-CA
# include<iostream>
# include<cmath>
# include<cstdlib>
# include<iomanip>
# include<fstream>
# include<sstream>
# include<string>

using namespace std;
const double Pi=3.14159;
const int Q=9;		//D2Q9 model
const int n=100;	//x direction
const int m=100;	//y direction

int e[Q][2]={{0,0},{1,0},{0,1},{-1,0},{0,-1},{1,1},{-1,1},{-1,-1},{1,-1}};
double w[Q]={4.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/36,1.0/36,1.0/36,1.0/36};

////i=1----n,j=1----m is a lattice node.i=1,j=1--m is the left boundary;i=n,j=1--m is the right boundary;
////j=1,i=1--n is the lower boundary;j=m£¬i=1--n is the upper boundary;
////i=0,i=n+1,j=0,j=m+1 are all virtual nodes.
double rho[n+2][m+2],u[n+2][m+2][2],u0[n+2][m+2][2],magU[n+2][m+2],
	   f[n+2][m+2][Q],fnew[n+2][m+2][Q],ftemp[n+2][m+2][Q];
double u_scale2[n+2][m+2][2];//Velocity field under nutrient scale 2

double temperature=25;//Celsius
const int gujiaN=39;//Number of grains

double gujia[gujiaN+1][2];//grain center position
double gujiara[gujiaN+1];//grain radius
double biora[gujiaN+1];//Biomass radius
int staticsolid[n+2][m+2];//Static solid , 0 table liquid, 1 table solid
int staticp[gujiaN+1][n+2][m+2];//Static solids
int statices00p[gujiaN+1][n+2][m+2];//Static solids
int statices00[n+2][m+2];//Pure grain static
double es[n+2][m+2],es00[n+2][m+2],Bs[n+2][m+2],es1[gujiaN+1][n+2][m+2], es2[gujiaN+1][n+2][m+2], vp[n+2][m+2][2];  //Percentage of solids, percentage of solids in the grain, speed of solid motion
double  M[n+2][m+2],Mold[n+2][m+2];//Density on solid microorganism + extracellular polymer space       Lattice unit
double excess[n+2][m+2];//Excess biomass
int i,j,k,p1,ip,jp,contin1,contin2,contin3,  njs2,njs3,totalnjs;// Large loop control variables


double dx=1,dy=1;//Lattice size
double dt=1;// Time step and relaxation time
double c=dx/dt; //1.0
double cs=c/pow(3.0,0.5);


/////////////////////////////////

double  rho0_phy;    //Fluid density "kg/m3", determined by temperature
double niu_phy ;//Kinematic viscosity  "m^2/s", determined by temperature

double u1_phy = 1.20e-4*0.9987;// The actual velocity "m/s" is obtained from the REV calculation
double inlet_concf_phy = 2.5e-2*0.778;//Left side inlet nutrient concentration "kg/m^3" 

double H_phy = 2.0e-3;//Height"m"
double L_phy = 2.0e-3;//Length"m"

double Sc=7;//Smithian number
double D_phy ;//Nutrient dispersion coefficientm^2/s
double Ki_phy ;//Half-saturation constants"kg/m^3"  
double kgmax_phy;// Maximum growth rate "s^-1" 
double densitymax_phy ;//Biomass density of the grid cell  "kg/m^3" 
double Y;//Yield of biomass on substrate, "kg organisms/kg nutrients " 
double ms_phy ;//Maintenance coefficient "kg nutrients/kg organisms s^-1 = s^-1 " 

///////////// -----------------  Flow scale1   ---------------///////////////
double Length_scale=L_phy/(n-1);//length scale  
double rho0=1.0;//Fluid Density in Lattice unit
double density_scale1;//Density scale, 
double niu=0.05;// Kinematic viscosity in lattice unit

double  Time_scale1;//Time Scale
double velocity_scale1 ;//Velocity Scale

double u1;//Left inlet Velocity in  lattice unit

double tauBGK=(3.0*niu+0.5);//The relaxation time of f. 
/////////////// --------------  Nutrient scale2   -------------- //////////////////
 double dt2=1;//Time step at scale 2
 double Time_scale2;//Nutrient time scale

 double concf_scale2;//Nutrient concentration scale

 double viscosity_scale2;//Kinematic viscosity at scale 2 
 double velocity_scale2;// Velocity scale 2
 double ki;//Half-saturation constants in Lattice units  
 
 double D2;//Nutrient diffusion coefficient £¨Lattice units£©
 double cf1;//Left inlet nutrient concentration£¨Lattice units£©
 double u1_scale2;//Left inlet nutrient velocity£¨Lattice units£©
 double tauc;//Nutrient relaxation time 


 /////////////---------------  Microbial growth scale 3   -------------///////////////
 double dt3=1;// Time step at scale 3
 double Time_scale3;//Microbial growth time scale 
 double density_scale3;///Density scale  
 double velocity_scale3;//Velocity scale 3
 double kgmax;//Maximum growth rate in  lattice units  
 
 double ms;//Maintenance coefficient£¨Lattice units£©
 double cxm;//Maximum concentration of biomass £¨Lattice units£©
 double rhop;//Biomass spatial density £¨Lattice units£©£¬ 1 spatial lattice dxdy has at most cxm amount of biomass in it 

////////////////////////////////////////


const double a=8;//Cylindrical radius 

double solidc=-5*cf1;//Solids concentration for displaying

double solidMicr=0;//Total solid microbial mass accumulation  in lattice units  

double u_gra[n + 2][m + 2][4];
double shear[n + 2][m + 2][4], shear_max[n + 2][m + 2], det_co[n + 2][m + 2], u_temp[n + 2][m + 2][2], l_stress[n + 2][m + 2][4];
int dir[4][2] = { {1,0},{-1,0},{0 - 1},{0,1} };
double miu_f;
double shear_critical_phy ;
double C_det ;
double shear_critical;
int spread;//spreading again?
////////////////////////////////////////

void F_shear();
void init();//Initialization functions
void gujiainit();// Initial function of grain position and radius
int poro(double x0,double y0,int gn); //Porous media generating function
void property(double t);//Determination of physical parameters and lattice speed etc. by temperature
int biocycle(double x0,double y0,int gn);// Generation of circular biomass
int opposite(int di);//For the opposite direction of the bounce boundary

void fboundary();//Boundary Conditions
void Macrovelcity();//Calculation of macro variables
double feq(int k,double rho,double u[2]);//Equilibrium function
void collesion();//Collision, note the scope of only fluid area collision, solid walls and objects are not involved in the collision
void streaming();//Migration
void gboundary();//Concentration boundary conditions
void MacroConcen();//Calculation of macroscopic variables in terms of concentration
double geq(int k,double con1,double u[2]); //Concentration distribution equilibrium function
void collesionforg();//Concentration distribution function collision
void streamingforg();//Concentration distribution function migration


void Solidnew1();// Solids Renewal
void report(int njs1);//Output Data
void output(int mjs);//Output results
void deposit();//Deposition
void spreading();//Distribution
void output2();
void output3( int njs1);	
void permeabilityK(int njs1);//Permeability
void  report_construction(int njs1);
void porosection(int njs);
double random(double start, double end);//Random number generation function

double Error();//Flow field error
double Errorcf();//Concentration error
///////////////////////////////////////////////////////////

double Sumes,Sumes0;//Total solids percentage
double changees;//es Variation
double Porositysolid0,Porosity0, Porosity, Sumsolid0, Sumsolid, Sumsolid00,changesolid;
double Bss(double e);
/////////////////////////////////////////////

double con[n+2][m+2],con0[n+2][m+2],g[n+2][m+2][Q],gnew[n+2][m+2][Q],gtemp[n+2][m+2][Q];//concentration field
double conshow[n+2][m+2];//For display purposes, the solids section is deliberately set to a negative concentration


int Neisolid;    int dire[9]; int addk; double Mwponint;  
int day;
int njs;
//int shortdays=1;//Number of days per microbial growth in a 2-stage cycle
//int longdays=2;//Total number of days of circulating microbial growth at level 1

double fiw(int i1, int j1, int d1, int d2, int d3,int nx,int ny); //Request fi, nutrient rebound
////////////////////////////////////////////////////////////////


int main()///Main Functions
{
    using namespace std;

int mmax,readdata;
    ///////////////////////////////////  Stage 1 without moving particles Flow field calculation    //////////////
   
     property(temperature); 
   	init(); 
	output(0);
///////////////////////////////////////////	
 //---------------------------- Flow calculation-----------------//
	               njs2=0;  contin3=1; 

					do
					{
						 njs2=njs2+1;
						 contin3=0;
            
						 collesion();
						 streaming();
						 fboundary();
						 Macrovelcity();   
						 
	                       for(i=0;i<=n+1;i++) 
				               for(j=0;j<=m+1;j++)
							   {		   
					             u_scale2[i][j][0]=u[i][j][0]*velocity_scale1/velocity_scale2; //Velocity converted to grid velocity under the concentration scale
						         u_scale2[i][j][1]=u[i][j][1]*velocity_scale1/velocity_scale2;
							   }
                        
					   	      collesionforg();
						      streamingforg();
						      gboundary();
						      MacroConcen(); 
						 if(njs2%1==0) {	
											 cout << " ---Don't go away, you need to input the number of days for microbial growth later ------ "<<endl;
											 cout<<"  Busy in calculating flow field    The "<< njs2<<"  loop  Error  " << Error()<<endl; 
											 cout<<endl;
										}
	   
						 if( njs2<30000) contin3=1;
                       
					}
					while(contin3==1);

	               
		 //----------------------------------------- Flow calculation End----------------------------------------------------------------------------//
 day=0,njs=0;
 permeabilityK(0);//Initial Permeability


 day = 10;
  mmax=njs+int(day*24*3600/Time_scale3);
  printf(" mmax=%d   \n",  mmax);

///////////
	
     while(njs<mmax) 
     {
       njs++;  
				     

                    	if(njs %int(0.5*3600/Time_scale3) ==0) permeabilityK(njs); //The permeability is written once every half hour	   
						Solidnew1(); //Solids renewed by precipitation or dissolution
						if(njs%1==0) {    cout<<"--------biomass ----  "<<solidMicr<<  "occupying the grids   "<<Sumes<<endl;
						         
									    	cout<<"----------  Microbial growth-------   "<<endl;
										    
										}
           

                        //-----------------Nutrients in the biomass are consumed once for biogenesis ---------//
			                
                         ////
    			        njs3=0;  contin2=1;  
				         do
						 {
						     njs3=njs3+1;
						     contin2=0;
                            
                                                     collesion(); 
						      streaming();
						      fboundary();
						      Macrovelcity(); 
                       
                             for(i=0;i<=n+1;i++) 
				               for(j=0;j<=m+1;j++)
							   {		   
					             u_scale2[i][j][0]=u[i][j][0]*velocity_scale1/velocity_scale2; 
						         u_scale2[i][j][1]=u[i][j][1]*velocity_scale1/velocity_scale2;
							   }
                        
					   	      collesionforg();
						      streamingforg();
						      gboundary();
						      MacroConcen();
	      
							  if(njs3%10==0) {   cout<<"--------biomass ----  "<<solidMicr<<  "occupying the grids   "<<Sumes<<endl;
		                                      cout << " Concentration   The  "<< njs3<<"   loop.  Errorcf  " << Errorcf()<<endl;									     
                                            
						                      cout<<endl; 
											}
		                    if(njs< int(24*3600/Time_scale3) && njs3<1000 && Errorcf()>8.0e-3) contin2=1;  
	                        if( njs>=int(24*3600/Time_scale3) && njs3<1000 && Errorcf()>8.0e-2) contin2=1;
						 }
	                     while(contin2==1);
               
				        //------------------------ Timely nutrient concentration calculation End---------------------------------------//
                         
						 if(njs%int(0.5*3600/Time_scale3) ==0) 
						 {
							 report(njs);  //Write data once every half hour
	                         output3(njs);
                           report_construction(njs);
						   porosection(njs);
						 }
						 if(njs%2==0) output2();

                       ///////				  
			   
        //----------------------------------------- Biological growth ends --------------------------------------------------------------------------//
        

		//---------------------------- Flow calculation-----------------  -------//
			
			  if(njs%int(1*3600/Time_scale3) ==0) // Flow field to be recalculated after 1 hour of biological growth
				{  njs2=0;  contin3=1; 

					do
					{
						 njs2=njs2+1;
						 contin3=0;
            
						 collesion();
						 streaming();
						 fboundary();
						 Macrovelcity();   
						 
	                       for(i=0;i<=n+1;i++) 
				               for(j=0;j<=m+1;j++)
							   {		   
					             u_scale2[i][j][0]=u[i][j][0]*velocity_scale1/velocity_scale2; 
						         u_scale2[i][j][1]=u[i][j][1]*velocity_scale1/velocity_scale2;
							   }
                        
					   	      collesionforg();
						      streamingforg();
						      gboundary();
						      MacroConcen(); 

						 if(njs2%1==0) {	
											 cout << "-----Loop again for flow. The   "<< njs2<<"  loop, Error  " << Error()<<endl;                         
											 cout<<endl;
										}
						 if (Error() == 100) break;
                         if( njs2<1000  ) contin3=1;
					}
				 	while(contin3==1);

	               	}//   if(njs%int(1*24*3600/Time_scale3) ==0) // The flow field has to be recalculated after one day of biological growth End

		        //--------------------------------------------------- Flow calculation End------------------------------------------------//
				

   }
  
 
 return 0;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

void property(double t)
{
	if (t >= 5 && t < 10) { niu_phy = ((10 - t) * 1.519e-6 + (t - 5) * 1.308e-6) / 5.0; }       //Interpolation calculation of kinematic viscosity coefficients
	if (t >= 10 && t < 15) { niu_phy = ((15 - t) * 1.308e-6 + (t - 10) * 1.14e-6) / 5.0; }
	if (t >= 15 && t < 20) { niu_phy = ((20 - t) * 1.14e-6 + (t - 15) * 1.007e-6) / 5.0; }
	if (t >= 20 && t < 25) { niu_phy = ((25 - t) * 1.007e-6 + (t - 20) * 0.897e-6) / 5.0; }
	if (t >= 25 && t < 30) { niu_phy = ((30 - t) * 0.897e-6 + (t - 25) * 0.804e-6) / 5.0; }

	if (t < 6) { rho0_phy = 1000; }
	if (t >= 6 && t < 8) { rho0_phy = ((8 - t) * 1000 + (t - 6) * 999.9) / 2.0; }
	if (t >= 8 && t < 10) { rho0_phy = ((10 - t) * 999.9 + (t - 8) * 999.7) / 2.0; }
	if (t >= 10 && t < 15) { rho0_phy = ((15 - t) * 999.7 + (t - 10) * 999.1) / 5.0; }
	if (t >= 15 && t < 20) { rho0_phy = ((20 - t) * 999.1 + (t - 15) * 998.2) / 5.0; }
	if (t >= 20 && t < 25) { rho0_phy = ((25 - t) * 998.2 + (t - 20) * 997.1) / 5.0; }
	if (t >= 25 && t < 30) { rho0_phy = ((30 - t) * 997.1 + (t - 25) * 995.7) / 5.0; }

	///////// --------Biological growth parameters --------------- //////////
	D_phy = niu_phy / Sc;
	Ki_phy = 7.64e-3;
	kgmax_phy = 1.91e-5;
	densitymax_phy = 50;
	Y = 0.53;
	ms_phy = 6.8e-6;


	////////// ---------scale1 --------------/////////////
	density_scale1 = rho0_phy / rho0;
	Time_scale1 = Length_scale * Length_scale / (niu_phy / niu);
	velocity_scale1 = Length_scale / Time_scale1;
	u1 = u1_phy / velocity_scale1;

	////////// ---------scale2 --------------/////////////
	Time_scale2 = 1 * Time_scale1;
	concf_scale2 = inlet_concf_phy;
	viscosity_scale2 = Length_scale * Length_scale / Time_scale2;

	velocity_scale2 = Length_scale / Time_scale2;
	ki = Ki_phy / concf_scale2;
	D2 = D_phy / viscosity_scale2;
	cf1 = inlet_concf_phy / concf_scale2;
	u1_scale2 = u1_phy / velocity_scale2;
	//tauc=3*D2*dt/dx/dx+0.5;

  ///////////// --------------- Microbial growth scale3 ------------- ///////////////
	dt3 = 1;
	Time_scale3 = 100;
	density_scale3 = densitymax_phy;

	velocity_scale3 = Length_scale / Time_scale3;
	kgmax = kgmax_phy * Time_scale3;

	ms = ms_phy * Time_scale3;
	cxm = densitymax_phy / density_scale3;
	rhop = cxm / dx / dy;
	/////
	miu_f = rho0 * niu;
	shear_critical_phy = 1.39;
	C_det = 5.0;
	shear_critical = shear_critical_phy / density_scale1 / velocity_scale1 / velocity_scale1;
}





//---------------------------------------------------------
void gujiainit()////// grain position and radius
{ 
   ////////////grain center coordinates (grain does not overlap) Random  ///////////
    int filled[n+2][m+2]; int c=1;
     for(i=1;i<=n;i++) 
	      for(j=1;j<=m;j++)
             filled[i][j]=0;//0 without grain, 1 with grain

   int p2;
    ifstream fin1( "grains for LBM.dat");
    for(p2=1;p2<=gujiaN;p2++)
       fin1>>gujia[p2][0]>>gujia[p2][1]>>gujiara[p2];

     for(p2=1;p2<=gujiaN;p2++)
        biora[p2]=gujiara[p2]+0.006;

}


void init()
{
	std::cout<<"tauf= "<<tauBGK<<endl;
    gujiainit();
   
	  /////////////////////////////////////////
	 
    
    int ee,rr;
  
  for(p1=1;p1<=gujiaN;p1++)
	for(i=1;i<=n;i++) 
	 for(j=1;j<=m;j++)
	 { 
		staticp[p1][i][j]=0;
                if(biocycle(i,j,p1)==1) { staticp[p1][i][j]=1; // the inside of biocycle is solid
		                           M[i][j]=cxm;	 
									}
		if(poro(i,j,p1)==1) { statices00p[p1][i][j]=1;//pure grain
		                           M[i][j]=0;	 // No organisms in the grain
									}	

	   es1[p1][i][j]=0; es2[p1][i][j]=0;  
       int nn=20;//number of sub-grids
       for(ee=1;ee<=nn;ee++)
	 	  for(rr=1;rr<=nn;rr++)
		  { if(poro(i-0.5+1.0/nn*ee, j-0.5+1.0/nn*rr,p1)==1) 			 
                 es1[p1][i][j] = es1[p1][i][j] + 1.0/nn/nn;
            
             if(biocycle(i-0.5+1.0/nn*ee, j-0.5+1.0/nn*rr,p1)==1) 			 
                 es2[p1][i][j] = es2[p1][i][j] + 1.0/nn/nn;
		  }

	 }

    /////////////////

      	for(i=1;i<=n;i++) 
	    	for(j=1;j<=m;j++)
            { es[i][j]=0; 
			  staticsolid[i][j]=0;
			  vp[i][j][0]=0;
			  vp[i][j][1]=0;
			}
    
  	   for(i=1;i<=n;i++) 
		for(j=1;j<=m;j++)
		{es00[i][j]=es1[1][i][j]; es[i][j]=es2[1][i][j]; 
		 for(p1=1;p1<=gujiaN;p1++)
		 {if(es1[p1][i][j]>es00[i][j]) es00[i][j]=es1[p1][i][j];  
        
          if(es2[p1][i][j]>es[i][j]) es[i][j]=es2[p1][i][j];  
		 }

          staticsolid[i][j]=staticp[1][i][j];  statices00[i][j]=statices00p[1][i][j]; 
         for(p1=1;p1<=gujiaN;p1++)
		 {if(staticp[p1][i][j]>staticsolid[i][j])  staticsolid[i][j]=staticp[p1][i][j];
          if(statices00p[p1][i][j]>statices00[i][j])  statices00[i][j]=statices00p[p1][i][j];
		 }
             
		}
	   for (i = 1; i <= n; i++)
		   for (j = 1; j <= m; j++)
		   {
			   M[i][j] = es[i][j] - es00[i][j];
		   }
/////////
    

      for(i=0;i<=n+1;i++) //
		for(j=0;j<=m+1;j++)
		{
		  u[i][j][0]=0;
		  u[i][j][1]=0;    
	      rho[i][j]=rho0; 
         
		  con[i][j]=0;//
		   con0[i][j]=0;
		}




		for(j=1;j<=m;j++) //Left boundary
		{	
	 	  u[1][j][0]=u1;
		  con[1][j]=cf1;
                }
	    for(i=0;i<=n+1;i++) 
		for(j=0;j<=m+1;j++)
		  for(k=0;k<Q;k++)
		  {
			  f[i][j][k]=feq(k,rho[i][j],u[i][j]);
              fnew[i][j][k]=feq(k,rho[i][j],u[i][j]);
			  ftemp[i][j][k]=feq(k,rho[i][j],u[i][j]);

              g[i][j][k]=geq(k,con[i][j],u_scale2[i][j]);
              gnew[i][j][k]=geq(k,con[i][j],u_scale2[i][j]);
			  gtemp[i][j][k]=geq(k,con[i][j],u_scale2[i][j]);

  		  }

	 for(i=0;i<=n+1;i++) 
		for(j=0;j<=m+1;j++)
         if(statices00[i][j]==1) 
		 { for(k=0;k<Q;k++)
			{ g[i][j][k]=0; gnew[i][j][k]=0; gtemp[i][j][k]=0;}
		 }



       for(i=0;i<=n+1;i++) // Show pure grain and biomass
		for(j=0;j<=m+1;j++)
		{	conshow[i][j]=0;
            if(es00[i][j]>0) conshow[i][j]=solidc;
			if(M[i][j]>0)   conshow[i][j]=1;

		}
 Sumes0 = 0,Sumsolid0=0;
	   for (i = 0; i <= n + 1; i++) // Show pure grain and biomass
		   for (j = 0; j <= m + 1; j++)
		   {
			   Sumes0 = Sumes0 + M[i][j];
			   Sumsolid0 = Sumsolid0 + es[i][j];
			   Sumsolid00 = Sumsolid00 + es00[i][j];
		   }
	   Porosity0 = 1 - Sumsolid0 / (n - 1) / (m - 1);
	   Porositysolid0 = 1 - Sumsolid00 / (n - 1) / (m - 1);

    double sumM=0; 
     for(i=1;i<=n;i++) 
	for(j=1;j<=m;j++)
	   sumM=sumM+M[i][j];
     
	sumM=sumM/(n*m)*density_scale3;
     ofstream   out1( "Initial biomass.dat ");
	out1<< sumM<<  endl;

	cout<<"initial  ok"<<endl;	 
}




double feq(int k,double rho,double u[2]) 
{	double eu,uv,feq;
	eu=(e[k][0]*u[0]+e[k][1]*u[1]);
	uv=(u[0]*u[0]+u[1]*u[1]);
	feq=w[k]*rho*(1+eu/cs/cs+eu*eu/2.0/pow(cs,4)-uv/2.0/cs/cs);
	return feq;
}

void fboundary()
{    
   
        for(j=1;j<=m;j++)
		{  	
		  double uw=u1;//Left boundary inlet velocity
                  double vw=0;
                  double rhow=1.0/(1-uw)*(f[1][j][0]+f[1][j][2]+f[1][j][4]+2*(f[1][j][3]+f[1][j][6]+f[1][j][7]));//Density
                     f[1][j][1]=f[1][j][3]+2*rhow*uw/3.0;
	                 f[1][j][5]=f[1][j][7]-(f[1][j][2]-f[1][j][4])/2.0+rhow*uw/6.0+rhow*vw/2.0;
	                 f[1][j][8]=f[1][j][6]+(f[1][j][2]-f[1][j][4])/2.0+rhow*uw/6.0-rhow*vw/2.0;
  

                     f[n][j][3]=f[n-1][j][3];//Right boundary Open Interface
	                 f[n][j][7]=f[n-1][j][7];
	                 f[n][j][6]=f[n-1][j][6];

		   }
        
          /////Upper and lower boundaries are periodic boundaries
              for(i=1;i<=n;i++)
                { f[i][m][4]=f[i][1][4];
                  f[i][m][7]=f[i][1][7];
                  f[i][m][8]=f[i][1][8];

                  f[i][1][5]=f[i][m][5]; 
                  f[i][1][2]=f[i][m][2];
                  f[i][1][6]=f[i][m][6];
                 }

        
}

void Macrovelcity()
{
    for(j=1;j<=m;j++)
         for(i=1;i<=n;i++)
	 {   rho[i][j]=0;
	      for(k=0;k<Q;k++)
                 rho[i][j]=rho[i][j]+f[i][j][k];
           }
  

      for(i=1;i<=n;i++)
	 	 for(j=1;j<=m;j++)
		 { 
		    u0[i][j][0]=u[i][j][0];
			u0[i][j][1]=u[i][j][1];
			 
			u[i][j][0]=0;
			u[i][j][1]=0;
			for(k=0;k<Q;k++)
			{  
                  u[i][j][0]=u[i][j][0]+e[k][0]*f[i][j][k];
		          u[i][j][1]=u[i][j][1]+e[k][1]*f[i][j][k];
			}
			if (es[i][j] == 1)
			{
				u[i][j][0] = 0;
				u[i][j][1] = 0;
			}
		       u[i][j][0]=u[i][j][0]/rho[i][j];
               u[i][j][1]=u[i][j][1]/rho[i][j];
                             
              magU[i][j]=pow(u[i][j][0]*u[i][j][0]+u[i][j][1]*u[i][j][1],0.5);
			  
		 }  
	  F_shear();
}
void F_shear()
{
	for (i = 1; i <= n; i++)
		for (j = 1; j <= m; j++)
		{
			u_temp[i][j][0] = u[i][j][0];
			u_temp[i][j][1] = u[i][j][1];
			if (es[i][j] > 0)
			{
				u_temp[i][j][0] = 0;
				u_temp[i][j][1] = 0;
			}

		}
	for (i = 1; i <= n; i++)
		for (j = 1; j <= m; j++)
			for (int k1 = 0; k1 <= 3; k1++)
			{
				if (es[i][j] == 0)
				{
					if (es[i + dir[k1][0]][j + dir[k1][1]] == 0)
						l_stress[i][j][k1] = 1;
					if (0 < es[i + dir[k1][0]][j + dir[k1][1]])
						l_stress[i][j][k1] = 1.5 - es[i][j];
				}
				if (es[i][j] > 0)
				{
					if (es[i + dir[k1][0]][j + dir[k1][1]] == 0)
						l_stress[i][j][k1] = 1.5 - es[i][j];
					if (es[i + dir[k1][0]][j + dir[k1][1]] > 0)
						l_stress[i][j][k1] = 10000;
				}
				if (k1 <= 1) u_gra[i][j][k1] = u_temp[i + dir[k1][0]][j + dir[k1][1]][1] - u_temp[i][j][1];
				if (k1 > 1)u_gra[i][j][k1] = u_temp[i + dir[k1][0]][j + dir[k1][1]][0] - u_temp[i][j][0];
				shear[i][j][k1] = miu_f * fabs(u_gra[i][j][k1]) / l_stress[i][j][k1];
			}
	for (i = 1; i <= n; i++)
		for (j = 1; j <= m; j++)
		{
			shear_max[i][j] = 0;
			for (int k1 = 0; k1 <= 3; k1++)
			{

				if (shear_max[i][j] < shear[i][j][k1])
					shear_max[i][j] = shear[i][j][k1];
			}
		}
	for (i = 1; i <= n; i++)
		for (j = 1; j <= m; j++)
		{
			det_co[i][j] = 0;
			if (M[i][j] > 0)
			{
				if (shear_max[i][j] < shear_critical)
					det_co[i][j] = 0;
				if (shear_critical <= shear_max[i][j] && shear_max[i][j] <= 2 * shear_critical)
					det_co[i][j] = exp(C_det * (shear_max[i][j] - 2 * shear_critical) / shear_critical);
				if (shear_max[i][j] > 2 * shear_critical)
					det_co[i][j] = 1;
			}
			if (M[i][j] == 0)
				det_co[i][j] = 0;
			M[i][j] = M[i][j] * (1 - det_co[i][j]);
		}
}

void collesion()  
{ 
          for(i=1;i<=n;i++)
	 	 for(j=1;j<=m;j++)		
		  { Bs[i][j]=Bss(es[i][j]);			 
               for(k=0;k<Q;k++)
				{ double omega=(f[i][j][opposite(k)]-feq(opposite(k),rho[i][j],u[i][j]))-(f[i][j][k]-feq(k,rho[i][j],vp[i][j]));		  
			      fnew[i][j][k]=f[i][j][k]-(1-Bs[i][j])*(dt/tauBGK)*(f[i][j][k]-feq(k,rho[i][j],u[i][j]))+Bs[i][j]*omega;		
				}
			
          }
}

void streaming() 
{ 
    for(i=1;i<=n;i++)
	 	 for(j=1;j<=m;j++)
		   	for(k=0;k<Q;k++)
			{
               int nexti=i+e[k][0];
			   int nextj=j+e[k][1];
			   ftemp[nexti][nextj][k]=fnew[i][j][k];
			}

	 for(i=0;i<=n+1;i++)
	 	 for(j=0;j<=m+1;j++)
		   	for(k=0;k<Q;k++)			
			   f[i][j][k]=ftemp[i][j][k];
				 

}



int opposite(int di)
{  int op;
  if(di==0) op=0;
  if(di==1) op=3;
  if(di==2) op=4;
  if(di==3) op=1;
  if(di==4) op=2;
  if(di==5) op=7;
  if(di==6) op=8;
  if(di==7) op=5;
  if(di==8) op=6;
  return op;
}


 

int poro(double x0,double y0,int gn)
{ int mate=0;					
  if ( (y0-gujia[gn][1])*(y0-gujia[gn][1])+(x0-gujia[gn][0])*(x0-gujia[gn][0])<=gujiara[gn]*gujiara[gn])// circle
		 mate=1;

 return mate;
}

int biocycle(double x0,double y0,int gn)
{
  int mate=0;					
  if ( (y0-gujia[gn][1])*(y0-gujia[gn][1])+(x0-gujia[gn][0])*(x0-gujia[gn][0])<=biora[gn]*biora[gn])// circle
	 
  mate=1;

 return mate;
}


void output(int mjs) //For Tecplot
{
	ostringstream name;
	name<<"Results for tecplot"<<mjs<<".dat";
	ofstream out(name.str().c_str());

	out<<"Title= \"FE-Volume Brick Data\""<<endl
	<<"VARIABLES=\"X\",\"Y\",\"U\",\"V\",\"rho\",\"magU\",\"es\",\"concen\",\"Mij\",\"conshow\"    "
		<<endl<<"ZONE Nodes="<<(n+1)*(m+1)<<","<<"Elements="<<n*m<<","<<"DATAPACKING=point"<<","<<"ZONETYPE=FEQUADRILATERAL"<<endl;
			for(j=0;j<=m;j++)
				for(i=0;i<=n;i++)			
					out<<double(i)<<" "<<double(j)<<" "<<u[i][j][0]<<" "<<
						u[i][j][1]<<"       "<<rho[i][j]<<"     "<<magU[i][j]<<"  "<<es[i][j]<<"  "<<con[i][j]<<"  "<<M[i][j]<<"  "<< conshow[i][j]<<  endl;

		       for(j=1;j<=m;j++)
		        	for(i=1;i<=n;i++)
                        out<<i+(j-1)*(n+1)<<" "<<i+j*(n+1)<<" "<<j*(n+1)+i+1<<" "<<(j-1)*(n+1)+i+1<<endl;

}



//////////////////////////////  Nutrients //////////////////////

double gfeq(int k,double con1,double u[2]) 
{
	double eu,uv,gfeq;
	eu=e[k][0]*u[0]+e[k][1]*u[1];
	uv=u[0]*u[0]+u[1]*u[1];
	gfeq=w[k]*con1*(1+3*eu/c+4.5*eu*eu/pow(c,2)-1.5*uv/c/c);
   
	return gfeq;
}



double sourgf(int i1, int j1, int k1, int nt) // construct the nutrient source term distribution function
{ double rs,qfso; double uu[2];
	
    uu[0]=u_scale2[i1][j1][0]; uu[1]=u_scale2[i1][j1][1]; 
	if(nt==1)
      rs=(kgmax/Y+ms)*M[i1][j1]*(con[i1][j1]/(ki+con[i1][j1]));
                                                            
    if(nt==0)
      rs=(kgmax/Y+ms)*M[i1][j1]*(con0[i1][j1]/(ki+con0[i1][j1]));

   double rs2=rs*density_scale3/Time_scale3/(concf_scale2/Time_scale2);
   qfso=-rs2;// macroscopic source term calculation

   

   double eu=e[k1][0]*uu[0]+e[k1][1]*uu[1]; 
  
   return w[k1]*qfso*(1+(tauc-0.5)/tauc*eu/cs/cs);
}



double geq(int k,double con1,double u[2]) 
{
	double eu,uv,geq;
	eu=(e[k][0]*u[0]+e[k][1]*u[1]);
	uv=(u[0]*u[0]+u[1]*u[1]);
   	geq=w[k]*con1*(1+3*eu/c+4.5*eu*eu/pow(c,4)-1.5*uv/c/c);
	return geq;
}

void gboundary() 
{ 
	      for(j=1;j<=m;j++)
		{  			 			
		     double cw=cf1;//left boundary fixed inlet concentration
               
                 g[1][j][1]=geq(1,cw,u[1][j])+geq(3,cw,u_scale2[1][j])-g[1][j][3];
	             g[1][j][5]=geq(5,cw,u[1][j])+geq(7,cw,u_scale2[1][j])-g[1][j][7];
	             g[1][j][8]=geq(8,cw,u[1][j])+geq(6,cw,u_scale2[1][j])-g[1][j][6];
  

                 g[n][j][3]=g[n-1][j][3];
	             g[n][j][7]=g[n-1][j][7];
	             g[n][j][6]=g[n-1][j][6];

                
		   }
        

        
              for(i=1;i<=n;i++)
                { g[i][m][4]=g[i][1][4];
                  g[i][m][7]=g[i][1][7];
                  g[i][m][8]=g[i][1][8];

                  g[i][1][5]=g[i][m][5]; 
                  g[i][1][2]=g[i][m][2];
                  g[i][1][6]=g[i][m][6];
                 }
   
    
      double gfTmp[Q]; int d1,d2,d3; int nx,ny; double fi;
        for(i=2;i<=n-1;i++)
         for(j=2;j<=m-1;j++)
            if( es00[i][j]>0 && es00[i][j]<1.0 ) // presence of grain
                  {
                
					for(k=0;k<Q;k++)
				        gfTmp[k]=g[i][j][k];
                    if(statices00[i+e[5][0]][j+e[5][1]]==0 && statices00[i+e[2][0]][j+e[2][1]]==0 && statices00[i+e[6][0]][j+e[6][1]]==0 ) 
					  { d1=5;d2=2;d3=6;// point to  the fluid
					    nx=0; ny=1;//normal direction                       
					  }
                    if(statices00[i+e[2][0]][j+e[2][1]]==0 && statices00[i+e[6][0]][j+e[6][1]]==0 && statices00[i+e[3][0]][j+e[3][1]]==0 )                
					  { d1=2;d2=6;d3=3;
					    nx=-1; ny=1;//normal direction                       
					  }

                   if(statices00[i+e[6][0]][j+e[6][1]]==0 && statices00[i+e[3][0]][j+e[3][1]]==0 && statices00[i+e[7][0]][j+e[7][1]]==0 )               
					  { d1=6;d2=3;d3=7;
					    nx=-1; ny=0;//normal direction                       
					  }
                   if(statices00[i+e[3][0]][j+e[3][1]]==0 && statices00[i+e[7][0]][j+e[7][1]]==0 && statices00[i+e[4][0]][j+e[4][1]]==0 )                
					  { d1=3;d2=7;d3=4;
					    nx=-1; ny=-1;//Normal direction                       
					  }
                    if(statices00[i+e[7][0]][j+e[7][1]]==0 && statices00[i+e[4][0]][j+e[4][1]]==0 && statices00[i+e[8][0]][j+e[8][1]]==0 )                
					  { d1=7;d2=4;d3=8;
					    nx=0; ny=-1;//Normal direction                       
					  }

                    if(statices00[i+e[4][0]][j+e[4][1]]==0 && statices00[i+e[8][0]][j+e[8][1]]==0 && statices00[i+e[1][0]][j+e[1][1]]==0 )               
					  { d1=4;d2=8;d3=1;
					    nx=1; ny=-1;//Normal direction                       
					  }

                   if( statices00[i+e[8][0]][j+e[8][1]]==0 && statices00[i+e[1][0]][j+e[1][1]]==0 && statices00[i+e[5][0]][j+e[5][1]]==0)               
					  { d1=8;d2=1;d3=5;
					    nx=1; ny=0;//normal direction                       
					  }
                   if( statices00[i+e[1][0]][j+e[1][1]]==0 && statices00[i+e[5][0]][j+e[5][1]]==0 && statices00[i+e[2][0]][j+e[2][1]]==0)                
					  { d1=1;d2=5;d3=2;
					    nx=1; ny=1;//normal direction                       
					  }


                   fi=fiw(i,j,d1,d2,d3,nx,ny);
						                       
                  
                   
                   g[i][j][d1]=-g[i][j][opposite(d1)]+2*w[d1]*fi;
				   g[i][j][d2]=-g[i][j][opposite(d2)]+2*w[d2]*fi;
				   g[i][j][d3]=-g[i][j][opposite(d3)]+2*w[d3]*fi;

                   for(k=0;k<Q;k++) 
                      g[i][j][k]=gfTmp[k];
                   }

  


}

void MacroConcen()
{
    
       for(j=1;j<=m;j++)
         for(i=1;i<=n;i++)
          {  
			 con0[i][j]=con[i][j];

			 con[i][j]=0;
	         for(k=0;k<Q;k++)
                 con[i][j]=con[i][j]+g[i][j][k];

           
		 }  
}



void collesionforg() 
{  
	double Dcf2;
      for(i=1;i<=n;i++)
	    for(j=1;j<=m;j++)
	
			{ Dcf2=(1-es[i][j])*D2+es[i][j]*0.5*D2;
               tauc=3*Dcf2*dt2/dx/dx+0.5;
			
			  for(k=0;k<Q;k++)				  	           
                      gnew[i][j][k]=g[i][j][k]-dt2/tauc*(1-es00[i][j])*(g[i][j][k]-geq(k,con[i][j],u_scale2[i][j]))
                                    +(1-es00[i][j])*dt2*sourgf(i,j,k,1)+0.5*dt2*dt2*(sourgf(i,j,k,1)-sourgf(i,j,k,0))/dt2;  
			 
                                                                    
			}	
}

void streamingforg() 
{  for(i=1;i<=n;i++)
	 for(j=1;j<=m;j++)
		for(k=0;k<Q;k++)
			{
               int nexti=i+e[k][0];
			   int nextj=j+e[k][1];
			   gtemp[nexti][nextj][k]=gnew[i][j][k];
			}

	 for(i=0;i<=n+1;i++)
	 	 for(j=0;j<=m+1;j++)
		   	for(k=0;k<Q;k++)			
			   g[i][j][k]=gtemp[i][j][k];
				 

}

/////////////////
void Solidnew1() // Solid update
{     

      deposit();// deposit

		 ////////////////
	  Sumes = 0; Sumsolid = 0;
	  for(i=1;i<=n;i++) 
	   for(j=1;j<=m;j++)
	   { 
	    conshow[i][j]=0; 
        if(M[i][j]>0) { conshow[i][j]=1; Sumes=Sumes+M[i][j];} 
		if(es00[i][j]>0) conshow[i][j]=solidc;
    	Sumsolid = Sumsolid + es[i][j];
       }
	      
     changees = (Sumes-Sumes0)/Sumes0; 
	 Porosity = 1 - Sumsolid / (m - 1) / (n - 1);
	 changesolid = (Sumsolid - Sumsolid0) / Sumsolid0;
}


void deposit() 
{
    double dM = 0;
   
   for(i=1;i<=n;i++)
	 for(j=1;j<=m;j++) 
	 { excess[i][j]=0; 
	   Mold[i][j]=M[i][j];
	 }
    spread=0; 
    for(i=2;i<=n-1;i++)
	 for(j=2;j<=m-1;j++)  
			{ double rs=(kgmax/Y+ms)*M[i][j]*(con[i][j]/(ki+con[i][j]));//ratio 3 under Corresponding Effect of Diffusive and Convective
                                                           
		        dM=dt3*(Y*(rs-ms*M[i][j]));
             
				
				solidMicr=solidMicr+dM;
			    if(solidMicr<0) solidMicr=0;
	                             					
                     if(dM>0) ///////// growth becomes larger ///////
                             { 
                                
								 if((dM+Mold[i][j])/rhop+es00[i][j]<1.0 ) {M[i][j]=M[i][j]+dM;													   
                                                                            es[i][j]=es00[i][j]+M[i][j]/rhop;// solids share = grain share + cured biomass volume share ,now es00[i][j]=0,no grain 
													  
																			}
                                 if((dM+Mold[i][j])/rhop+es00[i][j]>1.0 ) 
									{ 
									  M[i][j]=(1-es00[i][j])*rhop; 
								      excess[i][j] =dM-(M[i][j]-Mold[i][j]);
                                      spread=1;//to allocate the excess amount, to allocate
									 
									  es[i][j]=1;						
								      staticsolid[i][j]=1;
                                      u[i][j][0]=0; u[i][j][1]=0;
									}  																																				
							 }
                   
				}
           
      //////// ////

        
	 while (spread==1) 
	     spreading(); 
			
}


void spreading()
{ 
     
            spread=0;
            for(i=2;i<=n-1;i++)
	          for(j=2;j<=m-1;j++)  
                if(excess[i][j]>0) 
				{
					 Neisolid=0; 
                      for(k=1;k<=8;k++)																			
						if (es[i+e[k][0]][j+e[k][1]]<1 )
							{ Neisolid=Neisolid+1;
							  dire[Neisolid]=k;
							}
                     	
					if( Neisolid>0) 						
					{ for(addk=1;addk<=Neisolid;addk++)
                                  {
                                     
								        int idx=i+e[dire[addk]][0]; int jdy=j+e[dire[addk]][1];
										 Mwponint=M[idx][jdy];
										 if(Mwponint/rhop+es00[idx][jdy]+excess[i][j]/double(Neisolid)/rhop <1.0)
										 {
										   M[idx][jdy]=M[idx][jdy]+excess[i][j]/double(Neisolid);
                                           es[idx][jdy]=es00[idx][jdy]+M[idx][jdy]/rhop; 
                                         
										 }
									    if(Mwponint/rhop+es00[idx][jdy]+excess[i][j]/double(Neisolid)/rhop >=1.0)
											{M[idx][jdy]=(1-es00[idx][jdy])*rhop;
											 es[idx][jdy]=1;
											 staticsolid[idx][jdy]=1;
										     u[idx][jdy][0]=0; u[idx][jdy][1]=0;
											 spread=1;//show that there is still redistribution
                                             excess[idx][jdy] =excess[idx][jdy]+(Mwponint/rhop+es00[idx][jdy]+excess[i][j]/double(Neisolid)/rhop -1.0)*rhop;
											}
 
									 } // for(addk=1;addk<=Neisolid;addk++)
					}
					else 
					{ int dir2=rand()% 8 + 1;
					   excess[i+e[dir2][0]][j+e[dir2][1]] =excess[i+e[dir2][0]][j+e[dir2][1]]+excess[i][j];
                       spread=1;//allocation switch is on, also have to allocate again 
					}

                 excess[i][j]=0;
				}//if(excess[i][j]>0) //there is excess
  
}



double Bss(double e)
{return e*(tauBGK-0.5)/((1-e)+(tauBGK-0.5));
}


double random(double start, double end)
{   
	return start+(end-start)*rand()/(RAND_MAX + 1.0);
}

void  report(int njs) 
{
  ofstream   out1( "Time/s     biomass     occupying grids .dat ",ios::app);
  out1<<njs*Time_scale3/3600<<"     " <<solidMicr<<  "   "<<Sumes<<endl;			 
}
void  report_construction(int njs)
{
	ofstream   out1("changes in grains structure ", ios::app);
	

	out1 << njs*Time_scale3/3600<< "   " << Sumsolid0 << "   " << Sumsolid << "   "  << Porositysolid0 << "   " << Porosity0 << "   " << Porosity << "¡¡¡¡" << changesolid << endl;
}
void porosection(int njs)//porosity of different section
{
	double sumes1 = 0, sumes2 = 0, sumes3 = 0;//solid volume of  different section 
	double porosec_1 = 0, porosec_2 = 0, porosec_3 = 0;//porosity of different section 
	for (i = 1; i <= 25; i++)
		for (j = 1; j <= m; j++)
			sumes1 = sumes1 + es[i][j];// First segment biomass kg/m^3
	for (i = 26; i <= 50; i++)
		for (j = 1; j <= m; j++)
			sumes2 = sumes2 +es[i][j];// Second segment biomass kg/m^3
	for (i = 51; i <= 100; i++)
		for (j = 1; j <= m; j++)
			sumes3 = sumes3 + es[i][j];// Third segment biomass kg/m^3
	porosec_1 = sumes1 / 25 / m;
	porosec_2 = sumes2 / 25 / m;
	porosec_3 = sumes3 / 50 / m;
	
	ofstream   out1("porosity of different section.dat ", ios::app);
	out1 << njs * Time_scale3 / 3600 << "          " << sumes1 << "    " << sumes2 << "     " << sumes3 <<"    " << porosec_1 << "    " << porosec_2 << "     " << porosec_3 << endl;


}

double Error() //flow field error in a step
{
      double error;
	  double temp1=0;
	  double temp2=0;

	  double errSum=0;
	  int nSum=0;
      
            			  
	  /////
		   for(i=1;i<=n;i++)
				for(j=1;j<=m;j++)
				{
				 temp1= ((u[i][j][0]-u0[i][j][0])*(u[i][j][0]-u0[i][j][0])					
					         +(u[i][j][1]-u0[i][j][1])*(u[i][j][1]-u0[i][j][1]));
					 
				 temp2=(u[i][j][0]*u[i][j][0]+u[i][j][1]*u[i][j][1]);	
			
				temp1=sqrt(temp1);
				temp2=sqrt(temp2);
                if(temp2>1e-10)   //Calculation error after flow field has flow rate 
					{errSum=errSum+temp1/(temp2+1e-30);
					 nSum=nSum+1;					
					}
				}

         if(nSum>0)  error=errSum/nSum; //Average relative error
         else error=100;
  return  error;

}


double Errorcf() //nutrient field error in a step
{

      double error;
	  double temp1=0;
	  double temp2=0;

	  double errSum=0;
	  int nSum=0;
      
            			  
	  /////
		   for(i=1;i<=n;i++)
				for(j=1;j<=m;j++)
				{
				 temp1=fabs(con[i][j]-con0[i][j]);
					 
				 temp2=con[i][j];	
			
                  if(temp2>1e-10*cf1)   //Calculation error after flow field has flow rate 
					{errSum=errSum+temp1/(temp2+1e-30);
					 nSum=nSum+1;					
					}
				}

         if(nSum>0)  error=errSum/nSum; //Average relative error
         else error=100;
  return  error;


}




void output2() //Output--------------------------------------------------------------------------------------------------------------------------Function22
{
	ostringstream name;
	name<<"20220310 ÐÂÕõÔú³¡Á¿.dat";
	ofstream out(name.str().c_str());

	out<<"Title= \"FE-Volume Brick Data\""<<endl
	<<"VARIABLES=\"X\",\"Y\",\"U\",\"V\",\"rho\",\"magU\",\"es\",\"concen\",\"Mij\",\"conshow\"    "
		<<endl<<"ZONE Nodes="<<(n+1)*(m+1)<<","<<"Elements="<<n*m<<","<<"DATAPACKING=point"<<","<<"ZONETYPE=FEQUADRILATERAL"<<endl;//Two-dimensional with this type ZONETYPE=FEQUADRILATERAL; three-dimensional with this type ZONETYPE=FEBRICK
			for(j=0;j<=m;j++)
				for(i=0;i<=n;i++)			
					out<<double(i)<<" "<<double(j)<<" "<<u[i][j][0]<<" "<<
						u[i][j][1]<<"       "<<rho[i][j]<<"     "<<magU[i][j]<<"  "<<es[i][j]<<"  "<<con[i][j]<<"  "<<M[i][j]<<"  "<< conshow[i][j]<<  endl;

		       for(j=1;j<=m;j++)
		        	for(i=1;i<=n;i++)
                        out<<i+(j-1)*(n+1)<<" "<<i+j*(n+1)<<" "<<j*(n+1)+i+1<<" "<<(j-1)*(n+1)+i+1<<endl;

}

double fiw(int i1, int j1,int d1, int d2, int d3,int nx,int ny)
{
	double Ax=2*w[d1]*e[d1][0]+2*w[d2]*e[d2][0]+2*w[d3]*e[d3][0];
	double Ay=2*w[d1]*e[d1][1]+2*w[d2]*e[d2][1]+2*w[d3]*e[d3][1];
	double Bx=-(g[i1][j1][opposite(d1)]*e[d1][0]+g[i1][j1][opposite(d2)]*e[d2][0]+g[i1][j1][opposite(d3)]*e[d3][0]);
   	double By=-(g[i1][j1][opposite(d1)]*e[d1][1]+g[i1][j1][opposite(d2)]*e[d2][1]+g[i1][j1][opposite(d3)]*e[d3][1]);

	double Cx=0;double Cy=0; int k2;
	for(k2=1;k2<=8;k2++)
	  if(k2!=d1 && k2!=d2 && k2!=d3)
		{Cx=Cx+g[i1][j1][k2]*e[k2][0];
	     Cy=Cy+g[i1][j1][k2]*e[k2][1];
		}
	

    double fi1=(-(Cx+Bx)*nx-(Cy+By)*ny)/(Ax*nx+Ay*ny);
	return fi1;

}
void output3(int njs) //Output--------------------------------------------------------------------------------------------------------------------------Function22
{
	ostringstream name;
	name << "Fields after microbial growth " << njs << ".dat";
	ofstream out(name.str().c_str());

	out << "Title= \"FE-Volume Brick Data\"" << endl
		<< "VARIABLES=\"X\",\"Y\",\"U\",\"V\",\"rho\",\"magU\",\"u  m/s\"\"es\",\"concen\",\"C  mg/L\",\"Mij\",\"conshow\",\"det\",\"shear\"   "
		<< endl << "ZONE Nodes=" << (n + 1) * (m + 1) << "," << "Elements=" << n * m << "," << "DATAPACKING=point" << "," << "ZONETYPE=FEQUADRILATERAL" << endl;//Two-dimensional with this type ZONETYPE=FEQUADRILATERAL; three-dimensional with this type ZONETYPE=FEBRICK
	for (j = 0; j <= m; j++)
		for (i = 0; i <= n; i++)
			out << double(i) << " " << double(j) << " " << u[i][j][0] << " " <<
			u[i][j][1] << "       " << rho[i][j] << "     " << magU[i][j] << "  " << magU[i][j] * velocity_scale1 << "  " << es[i][j] << "  " << con[i][j] << "  " << con[i][j] * concf_scale2 << "  " << M[i][j] << "  " << conshow[i][j] << "  " << det_co[i][j] << "  " << shear_max[i][j] <<  endl;

	for (j = 1; j <= m; j++)
		for (i = 1; i <= n; i++)
			out << i + (j - 1) * (n + 1) << " " << i + j * (n + 1) << " " << j * (n + 1) + i + 1 << " " << (j - 1) * (n + 1) + i + 1 << endl;

}

void permeabilityK(int njs)//Permeability
{
   double pin1,pout1;//Import and export pressure
   double spin=0.0; 
   double spout=0.0;
   double uout = 0.0;
	  for(int j=1;j<=m;j++)
	  {	spin=spin+rho[1][j]*cs*cs;
        spout=spout+rho[n][j]*cs*cs;
		uout = uout + u[n][j][0];
	  }
	  pin1 = spin / m; pout1 = spout / m; uout = uout / m;
	 
	double K1=(niu*rho0*(n-1)*uout)/(pin1-pout1);


    double sumM=0,sum1=0,sum2=0,sum3=0;//LBM Calculate the total biomass in the domain  
      for (i = 1; i <= 25; i++)
		for (j = 1; j <= m; j++)
			sum1 = sum1 + M[i][j];// First segment biomass kg/m^3
	   for (i =26; i <= 50; i++)
		for (j = 1; j <= m; j++)
			sum2 = sum2 + M[i][j];// Second segment biomass kg/m^3
       for (i =51; i <= 100; i++)
		for (j = 1; j <= m; j++)
			sum3 = sum3 + M[i][j];// Third segment biomass kg/m^3
       for(i=1;i<=n;i++) 
    	for(j=1;j<=m;j++)
	   sumM=sumM+M[i][j];// biomass kg/m^3
	sum1=sum1/(25*m)*density_scale3;
	sum2=sum2/(25*m)*density_scale3;
    sum3=sum3/(50*m)*density_scale3;
    sumM=sumM/(n*m)*density_scale3;
    ofstream   out1( "Record every half hour   Permeability   biomass.dat ",ios::app);
   	out1<<njs*Time_scale3/3600<<"          "<<  K1*Length_scale*Length_scale<<"       " << 	sum1 << "    " << 	sum2 << "     " << sum3 << "      "<< sumM<<  endl;
        	      	

}
