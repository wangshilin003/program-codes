#include <stdio.h>
#include <stdlib.h>
#include <math.h>
# include<iostream> 
# include<cmath>
# include<cstdlib>
# include<iomanip>
# include<fstream>
# include<sstream>
# include<string>
# include<ctime>
#include <vector>
using namespace std;


#define L_phy 0.44//length of computational zone 0.44 m
#define H_phy 0.1//height of computational zone 0.1 m

#define L 0.44//1.0e-6//length in dugks   
#define H 0.1
//#define H H_phy/L_phy*L
#define M0 50// grids in  y 
#define N0 220 //grids in  x 
#define M (M0+2)
#define N (N0+2)
#define M1 (M+1)
#define N1 (N+1)
#define Qx 3
#define Qy 3
#define imin 1
#define imax (N-2)
#define jmin 1
#define jmax (M-2)

#define CFL 0.6
#define rho0 1.0

#define epsil0 0.734716//porosity
#define epsilon1 0.734716

#define nuliquid_phy 1.0e-6//The kinetic viscosity m2/s
#define u_inlet_phy 1.2e-4//inlet velocity m/s
#define c_oninlet_phy 2.5e-4//inlet concentration kg/m3
#define D_phy 1.28e-9//diffusion coefficient of  concentration  m2/s

#define Gx 0
#define Gy 0 
#define a 2.0

void gks_ini(void);
void SlopeY(int i);
void SlopeX(int j);
void InterpX(int j);
void InterpY(int i);

void gInterpX(int j);
void gInterpY(int i);

void boundary(void);
void gboundary(void);

void Evol(void);
void gEvol(void);
double feq(int kx, int ky, double Ux, double Uy, double RHO, double Eps);
//double geq(int kx, int ky, double Ux, double Uy, double CON);
double geq(int dcon, double Ux, double Uy, double CON);
void decpoints_init();
 
void output_Aver_UandF();
void output_decpoint_UandF();


void output(int n, double m);
void output_u_con();//output dimensionless velocity and concentration in microbial growth zone
void permeabilityK();

double Timescale;
double Lengthscale;
double nuscale;
double velocity_scale;
double concf_scale;

double dx,dy, dt,RT,w,wb,tau,nu,Cs,wg,wbg,D,gtau;
double f[M][N][Qx][Qy]; // f~ at the cell center
double f_plus[M][N][Qx][Qy]; // f_bar^+ at the cell center
double rho[M][N],ux[M][N],uy[M][N], ux0[M][N], uy0[M][N]; //cell center
double xc[N],yc[M],ti[M1];
double xf_face[N1][Qx][Qy], yf_face[M1][Qx][Qy];      // cell interface
double p[M][N];//press

double dp;
double epsil[M][N][Qx][Qy];
double g[M][N][5]; // f~ at the cell center
double g_plus[M][N][5]; // f_bar^+ at the cell center
double con[M][N]; //cell center
double xg_face[N1][5], yg_face[M1][5];      // cell interface

double epsilon[M][N],K[M][N],Fe[M][N],Fx[M][N],Fy[M][N]; //porosity£¬permeability¡¢structure parameters¡¢total force

double ex[Qx]={0., 1., -1.};
double ey[Qy]={0., 1., -1.};  
int re[Qx]={0,2,1}; 
int step;
int contin2;
double Timescale2;

double tpx[Qx]={2.0/3.0, 1.0/6.0, 1.0/6.0};
double tpy[Qy]={2.0/3.0, 1.0/6.0, 1.0/6.0}; 
double wgg[5]={1.0/3.0, 1.0/6.0, 1.0/6.0, 1.0/6.0, 1.0/6.0};

double u_inlet;
double con_inlet;
double yfFxold[M][N], yfFyold[M][N],xfFxold[M][N], xfFyold[M][N]; 
double ErrorU();
int sjp[M][N];
int sjcs[M][N];

////////////////////////////
double K_init_phy = 1.18129E-9;
double K_init,K_min;
int main()
{
    srand(time(0));
    int m, readdata, mmax;
    double err, u_old;



    RT = 1.0 / 3;
    Cs = pow(RT / sqrt(3 * RT), 0.5);

    for (m = 0; m < Qx; m++) { ex[m] = sqrt(3 * RT) * ex[m]; ey[m] = sqrt(3 * RT) * ey[m]; }

    Lengthscale = L_phy / L;
    nu = 1.0e-10;
    nuscale = nuliquid_phy / nu;

    velocity_scale = nuscale/Lengthscale;
    Timescale = Lengthscale / velocity_scale;
    Timescale2 = 10.0;
 
    concf_scale = c_oninlet_phy;


    u_inlet = u_inlet_phy / velocity_scale;
    
    con_inlet = c_oninlet_phy / concf_scale;
    D = D_phy / nuscale;

    dp = 1.75e-4 / Lengthscale; //dimensionless diameter
    K_init = K_init_phy / Lengthscale / Lengthscale;
    K_min = K_init * 5.0e-5;
    
    dx = L / (imax - imin + 1);
    dy = H / (jmax - jmin + 1);
    //////////

    ////////////
    xc[1] = 0 + dx / 2.0;
    for (m = 2; m <= imax; m++)
        xc[m] = xc[1] + dx * (m - 1);
    xc[0] = -xc[1];
    xc[imax + 1] = xc[imax] + dx;

    yc[1] = 0 + dy / 2.0;
    for (m = 2; m <= jmax; m++)
        yc[m] = yc[1] + dy * (m - 1);
    yc[0] = -yc[1];
    yc[jmax + 1] = yc[jmax] + dy;

    dt = CFL * dy / sqrt(6 * RT);     //time step

    tau = nu / Cs / Cs+0.001;      //Relaxation coefficient
    gtau = D / Cs / Cs+0.001 ;      //Relaxation coefficient for concentration

    w = 1.5 * dt / (2 * tau + dt);      //coefficient for computing f+ with f
    wb = dt / (2 * tau + 0.5 * dt);     //for interface

    wg = 1.5 * dt  / (2 * gtau + dt);    //coefficient for computing g+ with g
    wbg = dt  / (2 * gtau + 0.5 * dt );   //for interface


    printf("dt=%e  tau=%e  w=%e\n", dt, tau, w);
    gks_ini();  //initialization

    u_old = ux[M / 2 + 1][N / 2 + 1];//the center velocity 
    m = 0;
    int contin = 0;
    while (contin == 0)
    {
        m++;   // step (the current calculation steps from starting)
        Evol();
        gEvol();

        if (m % 1 == 0)
        {
            err = ErrorU();
            u_old = ux[M / 2][N / 2];
            printf("err=%e  u=%e   m=%d\n", err, u_old, m);
        }
        if (err < 1.0e-8)contin = 1;
    }
    output(0, m);
    permeabilityK();
    mmax = 100000;
    int m1 = 0;
    step = 0;
    while (m1 <= mmax)
    {


        Evol();
        gEvol();
        if (m1 % 1 == 0)
        {
            err = ErrorU();
            u_old = ux[M / 2][N / 2];
            printf("err=%e  u=%e   m1=%d\n", err, u_old, m1);
        }
        m = 0;
        contin2 == 0;
        while (contin2 == 0)
        {
            m++;   // step (the current calculation steps from starting)
            Evol();
            gEvol();

            if (m % 1 == 0)
            {
                err = ErrorU();
                u_old = ux[M / 2][N / 2];
                printf("err=%e  u=%e   m1=%d     m=%d\n", err, u_old, m1, m);
            }
            if (err < 1.0e-6)contin2 = 1;
        }

        step = m1;

        if (int(m1 * Timescale2) % 1800 == 0)
        {
            output(m1, step * Timescale2 / 3600.0);
            permeabilityK();
        }
        m1++;
    }
    output(m1, step * Timescale2 / 3600.0);
    permeabilityK();
    return 0;
}

void gks_ini()
{
    int i, j, kx, ky, dcon;
    for (j = jmin - 1; j <= jmax + 1; j++)  for (i = imin - 1; i <= imax + 1; i++)
    {
        epsilon[j][i] = epsil0;  sjcs[j][i] = 0;
        if (i > 10 && i < 12)
        {
            sjp[j][i] = 1;
        }
    }
    for (j = jmin - 1; j <= jmax + 1; j++) for (i = imin - 1; i <= imax + 1; i++)  for (kx = 0; kx < Qx; kx++) for (ky = 0; ky < Qy; ky++)
    {
        epsilon[j][i] = epsil0;
        Fe[j][i] = 1.75 / sqrt(150 * epsilon[j][i] * epsilon[j][i] * epsilon[j][i]);
        K[j][i] = pow(epsilon[j][i], 3) * pow(dp, 2) / (150 * pow(1 - epsilon[j][i], 2));
        Fx[j][i] = -epsilon[j][i] * nu / K[j][i] * ux[j][i] - epsilon[j][i] * Fe[j][i] / sqrt(K[j][i]) * sqrt(ux[j][i] * ux[j][i] + uy[j][i] * uy[j][i]) * ux[j][i] + epsilon[j][i] * Gx;//Resistance in x-direction
        Fy[j][i] = -epsilon[j][i] * nu / K[j][i] * uy[j][i] - epsilon[j][i] * Fe[j][i] / sqrt(K[j][i]) * sqrt(ux[j][i] * ux[j][i] + uy[j][i] * uy[j][i]) * uy[j][i] + epsilon[j][i] * Gy;//Resistance in y-direction
    }


    for (j = jmin - 1; j <= jmax + 1; j++) for (i = imin - 1; i <= imax + 1; i++)
    {
        ux[j][i] = uy[j][i] = 0.0; rho[j][i] = rho0; ux[j][imin - 1] = u_inlet; con[j][imin - 1] = con_inlet;
        for (kx = 0; kx < Qx; kx++) for (ky = 0; ky < Qy; ky++)
        {
            f[j][i][kx][ky] = feq(kx, ky, ux[j][i], uy[j][i], rho[j][i], epsilon[j][i]);

        }

        for (dcon = 0; dcon <= 4; dcon++)
        {
            g[j][i][dcon] = geq(dcon, ux[j][i], uy[j][i], con[j][i]);
        }
    }


}

double feq(int kx, int ky, double Ux, double Uy, double RHO , double Eps)
{
  double uv,eu,x;

  eu=(ex[kx]*Ux+ey[ky]*Uy)/RT;
  uv=(Ux*Ux+Uy*Uy)/RT;
  x=tpx[kx]*tpy[ky]*RHO*(1.0+eu+0.5*eu*eu/Eps-0.5*uv/Eps);

  return x;
}




void boundary()
{
  int i,j, kx, ky;
  double epsL, epsR;

 //////// inlet boundary at left and outlet boundary at right
  
  epsL=(xc[imin]-xc[imin-1])/(xc[imin+1]-xc[imin]);  
  epsR=(xc[imax+1]-xc[imax])/(xc[imax]-xc[imax-1]);
  for(j=jmin;j<=jmax;j++) for(ky=0;ky<Qy;ky++) for(kx=0;kx<Qx;kx++) 
  {
	f_plus[j][imin-1][kx][ky]=(1+epsL)*f_plus[j][imin][kx][ky]-epsL*f_plus[j][imin+1][kx][ky];
    f_plus[j][imax+1][kx][ky]=(1+epsR)*f_plus[j][imax][kx][ky]-epsR*f_plus[j][imax-1][kx][ky];
  }
  
  
  for(j=jmin;j<=jmax;j++) // left 
  { double rho_b;
    rho_b=0.0;
    for(kx=0;kx<Qx;kx++)  //veloctiy inlet
    {
      if(ex[kx]==0) for(ky=0;ky<Qy;ky++) rho_b+=0.5*(f_plus[j][imin-1][kx][ky]+f_plus[j][imin][kx][ky]);
      else if(ex[kx]<0) for(ky=0;ky<Qy;ky++) rho_b+=(f_plus[j][imin-1][kx][ky]+f_plus[j][imin][kx][ky]); 
    }
    for(kx=0;kx<Qx;kx++) 
    {
      if(ex[kx]>0) for(ky=0;ky<Qy;ky++)
      f_plus[j][imin-1][kx][ky]=f_plus[j][imin-1][re[kx]][re[ky]] +4*rho_b*tpx[kx]*tpy[ky]*(ex[kx]*u_inlet)/RT
	                                   +f_plus[j][imin][re[kx]][re[ky]]-f_plus[j][imin][kx][ky];      
    } 
  }
  
 

  //////
  for(j=jmin;j<=jmax;j++) //right  
  {
    for(kx=0;kx<Qx;kx++) //outflow
    {
      if(ex[kx]<0) for(ky=0;ky<Qy;ky++)
      f_plus[j][imax+1][kx][ky]=f_plus[j][imax][kx][ky];
     
    }
  }


   // top & bottom
  epsL=(yc[jmin]-yc[jmin-1])/(yc[jmin+1]-yc[jmin]);
  epsR=(yc[jmax+1]-yc[jmax])/(yc[jmax]-yc[jmax-1]);
  for(i=imin-1;i<=imax+1;i++) for(ky=0;ky<Qy;ky++) for(kx=0;kx<Qx;kx++)
  {
    f_plus[jmin-1][i][kx][ky]=(1+epsL)*f_plus[jmin][i][kx][ky]-epsL*f_plus[jmin+1][i][kx][ky];
    f_plus[jmax+1][i][kx][ky]=(1+epsR)*f_plus[jmax][i][kx][ky]-epsR*f_plus[jmax-1][i][kx][ky];
  }

  for(i=imin-1;i<=imax+1;i++) //bottom 
  {
      for (ky = 0; ky < Qy; ky++)
      {
          if (ey[ky] > 0) for (kx = 0; kx < Qx; kx++)
              f_plus[jmin - 1][i][kx][ky] = f_plus[jmax + 1][i][kx][ky];//periodic boundary
              //f_plus[jmin - 1][i][kx][ky] = f_plus[jmin - 1][i][re[kx]][re[ky]] + f_plus[jmin][i][re[kx]][re[ky]] - f_plus[jmin][i][kx][ky];//bounce back
      }
  }

  for(i=imin-1;i<=imax+1;i++) // top
  {
   
      for (ky = 0; ky < Qy; ky++)
      {
          if (ey[ky] < 0) for (kx = 0; kx < Qx; kx++) // 
              f_plus[jmax + 1][i][kx][ky] = f_plus[jmin - 1][i][kx][ky];//periodic boundary
              //f_plus[jmax + 1][i][kx][ky] = f_plus[jmax + 1][i][re[kx]][re[ky]] + f_plus[jmax][i][re[kx]][re[ky]] - f_plus[jmax][i][kx][ky];//bounce back                 
      }
  }
  
}
///////////////////////////


void InterpX(int j)   //    f at cell interface: X-direction
{
  int i, kx, ky, iL,jL,jR;
  double x, y, fc, dfx, dfy, ux_face, uy_face, rho_face;
  double hR,hL,AL,AR,AC;
  double FM,s;
  jL=j-1;  jR=j+1;
  hL=yc[j]-yc[jL];
  hR=yc[jR]-yc[j];
  AC=(hR-hL)/(hL*hR); AL=hR/(hL*(hR+hL)); AR=hL/(hR*(hR+hL)); 
  for(i=imin;i<=imax+1;i++) // inner nodes
  {
   iL=i-1;
   dx=xc[i]-xc[iL];
	for(ky=0;ky<Qy;ky++) for(kx=0;kx<Qx;kx++)
    {
      fc=0.5*(f_plus[j][i][kx][ky]+f_plus[j][iL][kx][ky]);     
      dfx=(f_plus[j][i][kx][ky]-f_plus[j][iL][kx][ky])/dx;      
      dfy=0.5*(AC*(f_plus[j][iL][kx][ky]+f_plus[j][i][kx][ky])   
               + AR*(f_plus[jR][iL][kx][ky]+f_plus[jR][i][kx][ky])
               - AL*(f_plus[jL][iL][kx][ky]+f_plus[jL][i][kx][ky]));
      x=0.5*ex[kx]*dt; y=0.5*ey[ky]*dt;//half time step  
      xf_face[i][kx][ky]=fc-x*dfx-y*dfy;  
    }
  }


  ///// Macro velocity of interface in x-direction and update of xf_face
   
  for(i=imin;i<=imax+1;i++)  
  { 
    ux_face=uy_face=rho_face=0.0;
    for(kx=0;kx<Qx;kx++) for(ky=0;ky<Qy;ky++)
    {
      rho_face+=xf_face[i][kx][ky];
      ux_face+=ex[kx]*xf_face[i][kx][ky];
      uy_face+=ey[ky]*xf_face[i][kx][ky];
    }
    ux_face=(ux_face+0.5*rho_face*xfFxold[j][i]*dt)/rho_face;    
    uy_face=(uy_face+0.5*rho_face*xfFyold[j][i]*dt)/rho_face;
   
     double fepsi=epsilon[j][i];
	 double fFe=1.75/sqrt(150*fepsi*fepsi*fepsi);
	 double fK=pow(epsilon[j][i],3)*pow(dp,2)/(150*pow(1-epsilon[j][i],2));
	 double fFx=-fepsi*nu/fK*ux_face-fepsi*fFe/sqrt(fK)*sqrt(ux_face*ux_face+uy_face*uy_face)*ux_face+fepsi*Gx;
	 double fFy=-fepsi*nu/fK*uy_face-fepsi*fFe/sqrt(fK)*sqrt(ux_face*ux_face+uy_face*uy_face)*uy_face+fepsi*Gy;
	  xfFxold[j][i]=fFx; xfFyold[j][i]=fFy; 
	///
    for(kx=0;kx<Qx;kx++) for(ky=0;ky<Qy;ky++)
	{   FM=feq(kx,ky, ux_face, uy_face,rho_face,fepsi);	   
		s=(fFx*(ex[kx]-ux_face)+fFy*(ey[ky]-uy_face))/RT*FM;
		xf_face[i][kx][ky]=(1.0-0.5*wb)*xf_face[i][kx][ky]+0.5*wb*FM+tau*dt/2.0/(2*tau+dt/2.0)*s ;
	}
  }
}

void InterpY(int i)   //  f at cell interface
{
  int j, kx, ky, iL, iR, jL;
  double fc, x, y, ux_face, uy_face, rho_face, dfx, dfy;
  double hR,hL,AL,AR,AC;
  double FM,s;
  iL=i-1; iR=i+1;
  hL=xc[i]-xc[iL];hR=xc[iR]-xc[i];
  AC=(hR-hL)/(hL*hR); AL=hR/(hL*(hR+hL)); AR=hL/(hR*(hR+hL));

 //// Macro velocity of interface in y-direction and update of yf_face
  for(j=jmin;j<=jmax+1;j++)
  {
    jL=j-1;
    dy=yc[j]-yc[jL];
    for(ky=0;ky<Qy;ky++) for(kx=0;kx<Qx;kx++)
    {
      fc=0.5*(f_plus[j][i][kx][ky]+f_plus[jL][i][kx][ky]);
      dfy=(f_plus[j][i][kx][ky]-f_plus[jL][i][kx][ky])/dy;
      dfx=0.5*(AC*(f_plus[jL][i][kx][ky]+f_plus[j][i][kx][ky])
               + AR*(f_plus[jL][iR][kx][ky]+f_plus[j][iR][kx][ky])
               - AL*(f_plus[jL][iL][kx][ky]+f_plus[j][iL][kx][ky]));
      y=0.5*ey[ky]*dt; x=0.5*ex[kx]*dt;//half time step
      yf_face[j][kx][ky]=fc-x*dfx-y*dfy;
    }
  }

//
  for(j=jmin;j<=jmax+1;j++)
  {
    ux_face=uy_face=rho_face=0.0;
    for(kx=0;kx<Qx;kx++) for(ky=0;ky<Qy;ky++)
    {
      rho_face+=yf_face[j][kx][ky];
      ux_face+=ex[kx]*yf_face[j][kx][ky];
      uy_face+=ey[ky]*yf_face[j][kx][ky];
    }
    ux_face=(ux_face+0.5*rho_face*yfFxold[j][i]*dt)/rho_face;    uy_face=(uy_face+0.5*rho_face*yfFyold[j][i]*dt)/rho_face;

    ///
      double fepsi=epsilon[j][i];
	 double fFe=1.75/sqrt(150*fepsi*fepsi*fepsi);
	 double fK=pow(epsilon[j][i],3)*pow(dp,2)/(150*pow(1-epsilon[j][i],2));
	 double fFx=-fepsi*nu/fK*ux_face-fepsi*fFe/sqrt(fK)*sqrt(ux_face*ux_face+uy_face*uy_face)*ux_face+fepsi*Gx;
	 double fFy=-fepsi*nu/fK*uy_face-fepsi*fFe/sqrt(fK)*sqrt(ux_face*ux_face+uy_face*uy_face)*uy_face+fepsi*Gy;
	  yfFxold[j][i]=fFx; yfFyold[j][i]=fFy; 
	///
    for(kx=0;kx<Qx;kx++) for(ky=0;ky<Qy;ky++)
	{FM=feq(kx,ky, ux_face, uy_face,rho_face,fepsi);
	   	s=(fFx*(ex[kx]-ux_face)+fFy*(ey[ky]-uy_face))/RT*FM;
		yf_face[j][kx][ky]=(1.0-0.5*wb)*yf_face[j][kx][ky]+0.5*wb*FM+tau*dt/2.0/(2*tau+dt/2.0)*s ; 
	}
  }
}

void Evol()
{
    int i, j, kx, ky;
    double FM, s;
    double t1;
 
    for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            if (sjp[j][i] == 1 && contin2 == 1)
            {
                sjcs[j][i] = sjcs[j][i] + 1;
                t1 = sjcs[j][i] * Timescale2;
                if (K[j][i] >K_min)
                {
                    epsilon[j][i] = 0.734716 - 0.72593 / (1 + 3.49815E7 * exp(-0.18353 * t1 / 3600));// Equivalent porosity with time 
                    K[j][i] = pow(epsilon[j][i], 3) * pow(dp, 2) / (150 * pow(1 - epsilon[j][i], 2));                  
                }
                if (K[j][i] <= K_min)
                {
                    K[j][i] = K_min;
                    sjp[j][i] = 3;
                }
            }
        }
    //source
    for (j = jmin; j <= jmax; j++) for (i = imin; i <= imax; i++)
    {
       
        Fe[j][i] = 1.75 / sqrt(150 * epsilon[j][i] * epsilon[j][i] * epsilon[j][i]);       
        Fx[j][i] = -epsilon[j][i] * nu / K[j][i] * ux[j][i] - epsilon[j][i] * Fe[j][i] / sqrt(K[j][i]) * sqrt(ux[j][i] * ux[j][i] + uy[j][i] * uy[j][i]) * ux[j][i] + epsilon[j][i] * Gx;
        Fy[j][i] = -epsilon[j][i] * nu / K[j][i] * uy[j][i] - epsilon[j][i] * Fe[j][i] / sqrt(K[j][i]) * sqrt(ux[j][i] * ux[j][i] + uy[j][i] * uy[j][i]) * uy[j][i] + epsilon[j][i] * Gy;
    }

    // f_plus in each cell
    for (j = jmin; j <= jmax; j++) for (i = imin; i <= imax; i++) for (kx = 0; kx < Qx; kx++) for (ky = 0; ky < Qy; ky++)
    {
        FM = feq(kx, ky, ux[j][i], uy[j][i], rho[j][i], epsilon[j][i]);
        s = (Fx[j][i] * (ex[kx] - ux[j][i]) + Fy[j][i] * (ey[ky] - uy[j][i])) / RT * FM;
        f_plus[j][i][kx][ky] = f[j][i][kx][ky] - w * (f[j][i][kx][ky] - FM) + 3 * tau * dt / 2.0 / (2 * tau + dt) * s;
    }

    boundary();

    //update f: X-direction
    for (j = jmin; j <= jmax; j++)
    {
        InterpX(j);
        for (i = imin; i <= imax; i++)
        {
            dx = (xc[i + 1] - xc[i - 1]) * 0.5;
            for (kx = 0; kx < Qx; kx++) for (ky = 0; ky < Qy; ky++)
            {
                f[j][i][kx][ky] = (4.0 * f_plus[j][i][kx][ky] - f[j][i][kx][ky]) / 3.0
                    + ex[kx] * dt / dx * (xf_face[i][kx][ky] - xf_face[i + 1][kx][ky]);
            }
        }
    }

    //update f: Y-direction
    for (i = imin; i <= imax; i++)
    {
        InterpY(i);
        for (j = jmin; j <= jmax; j++)
        {
            dy = (yc[j + 1] - yc[j - 1]) * 0.5;
            for (kx = 0; kx < Qx; kx++) for (ky = 0; ky < Qy; ky++)
            {
                f[j][i][kx][ky] = f[j][i][kx][ky]
                    + ey[ky] * dt / dy * (yf_face[j][kx][ky] - yf_face[j + 1][kx][ky]);
            }
        }
    }

    //    update macroscopic variables in each cell
    for (j = jmin; j <= jmax; j++) for (i = imin; i <= imax; i++)
    {
        ux0[j][i] = ux[j][i];
        uy0[j][i] = uy[j][i];
        rho[j][i] = ux[j][i] = uy[j][i] = 0.0;
        for (kx = 0; kx < Qx; kx++) for (ky = 0; ky < Qy; ky++)
        {
            rho[j][i] += f[j][i][kx][ky];
            ux[j][i] += ex[kx] * f[j][i][kx][ky];
            uy[j][i] += ey[ky] * f[j][i][kx][ky];
        }
        ux[j][i] = (ux[j][i] + 0.5 * rho[j][i] * Fx[j][i] * dt) / rho[j][i];
        uy[j][i] = (uy[j][i] + 0.5 * rho[j][i] * Fy[j][i] * dt) / rho[j][i];
        p[j][i] = rho[j][i] * Cs * Cs;
    }
}
double ErrorU() 
{
    double error;
    double temp1 = 0;
    double temp2 = 0;
    double magU[M][N], magU0[M][N];
    int j, i;
    double errSum = 0;
    int nSum = 0;


    /////
    for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++)
        {
            magU[j][i] = pow(ux[j][i] * ux[j][i] + uy[j][i] * uy[j][i], 0.5);
            magU0[j][i] = pow(ux0[j][i] * ux0[j][i] + uy0[j][i] * uy0[j][i], 0.5);

            temp1 = temp1 + pow(magU[j][i] - magU0[j][i], 2);
            temp2 = temp2 + pow(magU[j][i], 2);
        }
    temp1 = sqrt(temp1);
    temp2 = sqrt(temp2);
    error = temp1 / temp2; 
    return  error;

}

/////////////////////////////////////  computer concentration //////////////////


double geq(int dcon, double Ux, double Uy, double CON)
{
  double eu,x;
  int exdcon[5]={0,1,0,-1,0};
  int eydcon[5]={0,0,1,0,-1};
  eu=(exdcon[dcon]*Ux+eydcon[dcon]*Uy)/RT;
 
  x=wgg[dcon]*CON*(1.0+eu);
  return x;
}


void gEvol()
{
  int i,j, dcon;
  double FM,s;

  // g_plus in each cell
  for(j=jmin;j<=jmax;j++) for(i=imin;i<=imax;i++) for(dcon=0;dcon<=4;dcon++) 
  {
	FM=geq(dcon, ux[j][i], uy[j][i],con[j][i]);
    s=0;  //source
    g_plus[j][i][dcon]=g[j][i][dcon]-wg*(g[j][i][dcon]-FM)+3*gtau*dt/2.0/(2*gtau+dt)*s;  
  }

  gboundary(); 

  //update g: X-direction
  for(j=jmin;j<=jmax;j++)
  {
	 gInterpX(j);
    for(i=imin;i<=imax;i++)
    {
      dx=(xc[i+1]-xc[i-1])*0.5;
      
     g[j][i][0]=(4.0*g_plus[j][i][0]-g[j][i][0])/3.0+0*dt/dx*(xg_face[i][0]-xg_face[i+1][0]);   
     g[j][i][1]=(4.0*g_plus[j][i][1]-g[j][i][1])/3.0+1*dt/dx*(xg_face[i][1]-xg_face[i+1][1]);
     g[j][i][2]=(4.0*g_plus[j][i][2]-g[j][i][2])/3.0+0*dt/dx*(xg_face[i][2]-xg_face[i+1][2]);  
     g[j][i][3]=(4.0*g_plus[j][i][3]-g[j][i][3])/3.0-1*dt/dx*(xg_face[i][3]-xg_face[i+1][3]);
     g[j][i][4]=(4.0*g_plus[j][i][4]-g[j][i][4])/3.0+0*dt/dx*(xg_face[i][4]-xg_face[i+1][4]);

    }
  }

  //update f: Y-direction
  for(i=imin;i<=imax;i++)
  {
	gInterpY(i);
    for(j=jmin;j<=jmax;j++)
    {
      dy=(yc[j+1]-yc[j-1])*0.5;
      g[j][i][0] =g[j][i][0]  +0*dt/dy*(yg_face[j][0]-yg_face[j+1][0]);                  
		     
      g[j][i][1]=g[j][i][1]+0*dt/dy*(yg_face[j][1]-yg_face[j+1][1]);                     
		    
     g[j][i][2]=g[j][i][2]+1.0*dt/dy*(yg_face[j][2]-yg_face[j+1][2]);                     
		     
     g[j][i][3]=g[j][i][3] +0*dt/dy*(yg_face[j][3]-yg_face[j+1][3]);                    
		    
     g[j][i][4]=g[j][i][4] +(-1)*dt/dy*(yg_face[j][4]-yg_face[j+1][4]);                    
   		    
    }
  }

  /////   update macroscopic variables in each cell
  for(j=jmin;j<=jmax;j++) for(i=imin;i<=imax;i++)
  {
    con[j][i]=0.0;
    for(dcon=0;dcon<=4;dcon++)
    {
      con[j][i]+=g[j][i][dcon];
     
    }   
   // con[j][i] = con[j][i] / rho[j][i];
  }
}

/////////////////////////////////////////////////////////
void gboundary()  
{
  int i,j, dcon;
  double epsL, epsR;

 ////////  left & right boundary
  
  epsL=(xc[imin]-xc[imin-1])/(xc[imin+1]-xc[imin]);  
  epsR=(xc[imax+1]-xc[imax])/(xc[imax]-xc[imax-1]);
  for(j=jmin;j<=jmax;j++) for(dcon=0;dcon<=4;dcon++) 
  {
    g_plus[j][imin-1][dcon]=(1+epsL)*g_plus[j][imin][dcon]-epsL*g_plus[j][imin+1][dcon];
    g_plus[j][imax+1][dcon]=(1+epsR)*g_plus[j][imax][dcon]-epsR*g_plus[j][imax-1][dcon];
  }
  
 
  for(j=jmin;j<=jmax;j++) //left 
  { 
   g_plus[j][imin - 1][1] = geq(1, u_inlet, 0, con_inlet) + geq(3, u_inlet, 0, con_inlet) - g_plus[j][imin - 1][3];
  }
 /////


  //////
  for(j=jmin;j<=jmax;j++) //right 
  {
   g_plus[j][imax+1][3]=g_plus[j][imax][3];

  }

   

   // top & bottom
  epsL=(yc[jmin]-yc[jmin-1])/(yc[jmin+1]-yc[jmin]);
  epsR=(yc[jmax+1]-yc[jmax])/(yc[jmax]-yc[jmax-1]);
  for(i=imin-1;i<=imax+1;i++) for(dcon=0;dcon<=4;dcon++) 
  {
    g_plus[jmin-1][i][dcon]=(1+epsL)*g_plus[jmin][i][dcon]-epsL*g_plus[jmin+1][i][dcon];
    g_plus[jmax+1][i][dcon]=(1+epsR)*g_plus[jmax][i][dcon]-epsR*g_plus[jmax-1][i][dcon];
  }

  for(i=imin-1;i<=imax+1;i++) //bottom
  {
  
    g_plus[jmin - 1][i][2] = g_plus[jmax + 1][i][2];
   
  }

  for(i=imin-1;i<=imax+1;i++) // top
  {  
    
    g_plus[jmax + 1][i][4] = g_plus[jmin - 1][i][4];
  }

}

/////////
void gInterpX(int j)   
{
  int i, dcon,kx,ky, iL,jL,jR;
  double x, y, fc, dfx, dfy, ux_face, uy_face, rho_face, con_face;
  double hR,hL,AL,AR,AC;
  double FM,s;
  int exdcon[5]={0,1,0,-1,0};
  int eydcon[5]={0,0,1,0,-1};
  jL=j-1;  jR=j+1;
  hL=yc[j]-yc[jL];
  hR=yc[jR]-yc[j];
  AC=(hR-hL)/(hL*hR); AL=hR/(hL*(hR+hL)); AR=hL/(hR*(hR+hL)); 
  for(i=imin;i<=imax+1;i++) // inner nodes
  {
   iL=i-1;
   dx=xc[i]-xc[iL];
	for(dcon=0;dcon<=4;dcon++) 
    {
      fc=0.5*(g_plus[j][i][dcon]+g_plus[j][iL][dcon]);      
      dfx=(g_plus[j][i][dcon]-g_plus[j][iL][dcon])/dx;     
      dfy=0.5*(AC*(g_plus[j][iL][dcon]+g_plus[j][i][dcon])   
               + AR*(g_plus[jR][iL][dcon]+g_plus[jR][i][dcon])
               - AL*(g_plus[jL][iL][dcon]+g_plus[jL][i][dcon]));
      x=0.5*exdcon[dcon]*dt; y=0.5*eydcon[dcon]*dt;//half time step 
      xg_face[i][dcon]=fc-x*dfx-y*dfy;  
    }
  }

 
  for(i=imin;i<=imax+1;i++)  
  {
     ux_face=uy_face=rho_face=con_face=0.0;
		
    for(dcon=0;dcon<=4;dcon++) 
    {
      con_face+=xg_face[i][dcon];
      
    }

   for(kx=0;kx<Qx;kx++) for(ky=0;ky<Qy;ky++)
    {
      
      rho_face+=xf_face[i][kx][ky];
      ux_face+=ex[kx]*xf_face[i][kx][ky];
      uy_face+=ey[ky]*xf_face[i][kx][ky];
    }
    ux_face=(ux_face+0.5*rho_face*Fx[j][i]*dt)/rho_face;    uy_face=(uy_face+0.5*rho_face*Fy[j][i]*dt)/rho_face;

    for(dcon=0;dcon<=4;dcon++) 
	{   FM=geq(dcon, ux_face, uy_face,con_face);
    	    s=0;  
	    xg_face[i][dcon]=(1.0-0.5*wbg)*xg_face[i][dcon]+0.5*wbg*FM+gtau*dt/2.0/(2*gtau+dt/2.0)*s ;
	}
  }
}

void gInterpY(int i)   //    f at cell interface
{
  int j,dcon, kx,ky, iL, iR, jL;
  double fc, x, y, ux_face, uy_face, rho_face, con_face, dfx, dfy;
  double hR,hL,AL,AR,AC;
  double FM,s;
  int exdcon[5]={0,1,0,-1,0};
  int eydcon[5]={0,0,1,0,-1};
  iL=i-1; iR=i+1;
  hL=xc[i]-xc[iL];hR=xc[iR]-xc[i];
  AC=(hR-hL)/(hL*hR); AL=hR/(hL*(hR+hL)); AR=hL/(hR*(hR+hL));

 
  for(j=jmin;j<=jmax+1;j++)
  {
    jL=j-1;
    dy=yc[j]-yc[jL];
    for(dcon=0;dcon<=4;dcon++) 
    {
      fc=0.5*(g_plus[j][i][dcon]+g_plus[jL][i][dcon]);
      dfy=(g_plus[j][i][dcon]-g_plus[jL][i][dcon])/dy;
      dfx=0.5*(AC*(g_plus[jL][i][dcon]+g_plus[j][i][dcon])
               + AR*(g_plus[jL][iR][dcon]+g_plus[j][iR][dcon])
               - AL*(g_plus[jL][iL][dcon]+g_plus[j][iL][dcon]));
      y=0.5*eydcon[dcon]*dt; x=0.5*exdcon[dcon]*dt;//half time step
      yg_face[j][dcon]=fc-x*dfx-y*dfy;
    }
  }

//
  for(j=jmin;j<=jmax+1;j++)
  {
    ux_face=uy_face=rho_face=con_face=0.0;
    for(dcon=0;dcon<=4;dcon++) 
    {
      con_face+=yg_face[j][dcon];
    }
    for(kx=0;kx<Qx;kx++) for(ky=0;ky<Qy;ky++)
    {
      rho_face+=yf_face[j][kx][ky];
      ux_face+=ex[kx]*yf_face[j][kx][ky];
      uy_face+=ey[ky]*yf_face[j][kx][ky];
    }

    ux_face=(ux_face+0.5*rho_face*Fx[j][i]*dt)/rho_face;    uy_face=(uy_face+0.5*rho_face*Fy[j][i]*dt)/rho_face;

   for(dcon=0;dcon<=4;dcon++)
	{FM=geq(dcon, ux_face, uy_face,con_face);
	   s=0;
		yg_face[j][dcon]=(1.0-0.5*wbg)*yg_face[j][dcon]+0.5*wbg*FM+gtau*dt/2.0/(2*gtau+dt/2.0)*s ; 
	}
  }
}


void output(int n, double m)
{
    int i, j;
    ostringstream name;
    name << n << "fields parameters " << m << ".dat";
    ofstream out(name.str().c_str());

    out << "Title= \"FE-Volume Brick Data\"" << endl
        << "VARIABLES=\"X\",\"Y\",\"U\",\"V\",\"con\",\"magU\",\"Poro\"  "
        << endl << "ZONE Nodes=" << (imax + 1) * (jmax + 1) << "," << "Elements=" << imax * jmax << "," << "DATAPACKING=point" << "," << "ZONETYPE=FEQUADRILATERAL" << endl;
    for (j = 0; j <= jmax; j++)
        for (i = 0; i <= imax; i++)
            out << xc[i] << " "
            << yc[j] << " " << ux[j][i] << " " <<
            uy[j][i] << "       " << con[j][i] << "    " << pow(ux[j][i] * ux[j][i] + uy[j][i] * uy[j][i], 0.5) << "    " << epsilon[j][i] << endl;


    for (j = 1; j <= jmax; j++)
        for (i = 1; i <= imax; i++)
            out << i + (j - 1) * (imax + 1) << " " << i + j * (imax + 1) << " " << j * (imax + 1) + i + 1 << " " << (j - 1) * (imax + 1) + i + 1 << endl;


}


void permeabilityK()
{
    double pin1, pout1;
    double spin = 0.0;
    double spout = 0.0;
    double u_out = 0.0;
    for (int j = jmin; j <= jmax; j++)
    {
        spin = spin + p[j][imin];
        spout = spout + p[j][imax];
        u_out = u_out + ux[j][imax];
    }
    pin1 = spin / jmax; pout1 = spout / jmax; u_out = u_out / jmax;
    double K1 = (nu * rho0 * (xc[imax] - xc[imin]) * u_out) / (pin1 - pout1);

    ofstream   out1("permeability with time .dat",ios::app);//(time/h- permeability/m2 )
    {out1 << step * Timescale2 / 3600.0 << "      " << K1 * Lengthscale * Lengthscale << endl; }

}
 