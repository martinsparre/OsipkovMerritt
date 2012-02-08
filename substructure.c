#include "gadget2conv.h"

#include <gsl/gsl_sf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

//Constants
#define G (1.0) // Gravitational constant
#define PI	(3.141592653589793238462643)
#define MVIR 1.0 //virial mass of the main cluster
#define RVIR 50.0
#define RS 1.0
#define RHO0 (1.0/2.0/PI)
#define R_AN 1.20 //ANISOTROPY RADIUS

//DENSITY PROFILE: (inner slope, outer slope, transition...)
#define ALPHA (1.0)
#define BETA (4.0)
#define GAMMA (1.0)


//Number of gridpoints:
#define NG	36000 //Number of Gridpoints
#define N_EDD 24000


double Mass; //Mass of each particle

//Structs:
struct Params
{
	double Rho0;
	double scale_length;
	double conc;
	double rvir;
	double Mvir;
	double r_an;
	double rdecay;
	double rmin;
	double rmax;
	int Ndm;
	double x[3];
	double v[3];
	double E[N_EDD];
	double DF[N_EDD];
};

struct Grid
{
	double R;
	double Rho;
	double Mass;
	double Sigma2r;
	double Phi;
	double MaxProb;
};

//gsl variables (used to generate random numbers...):
const gsl_rng_type * type;
	gsl_rng * r;
	gsl_rng * r1;
	
//Functions:
double GetMass(double Mmin, double Mmax);
void CreateCluster(struct Particle *P, struct Grid *Gr,struct Params * Cl);
void FindPosAndVel(double *Pos, double *Vel, struct Grid * Gr, struct Params * Cl);
void CreateGrid(struct Grid * Gr,  struct Params * Cl );
void Eddington( struct Grid *Gr, struct Params *Cl);
void CreateParticles(struct Particle *P,  struct Grid *Gr,  struct Params *Cl);
double dnu2dpsi2(double psi, struct Grid *Gr,struct Params *Cl);
double GetDF(double E, struct Params *Cl);
double GetDistFromPsi (double psi, struct Grid *Gr);
double rho(double r, struct Params *Cl);
double drhodr(double r, struct Params *Cl);
double d2rhodr2(double r, struct Params *Cl);
double uniform(void);
double Gaussian(double Sigma);
double dmax(double x, double y);

double drho_Q(double r, struct Params *Cl);
double d2rho_Q(double r, struct Params *Cl);


int main(int argc, char **argv)//first argument: filename. Second: Number of particles
{
	int i,Ndm;

	struct Particle *P;	
	struct Params *Cl;
	struct Grid *Gr;
	
	type = gsl_rng_default;
	r = gsl_rng_alloc (type);
	r1 = gsl_rng_alloc (type);
		
	Ndm = atoi(argv[2]);
	
	//Main structure:
	P = (struct Particle *) malloc(Ndm*sizeof(struct Particle));
	Gr = (struct Grid *) malloc(NG*sizeof(struct Grid));
	Cl = (struct Params *) malloc(sizeof(struct Params));

	//Define the main cluster
	Cl->Rho0 = 1.0/2.0/PI;
	Cl->scale_length=1.0;
	
	
	Cl->conc=RVIR/RS;
	Cl->rvir=RVIR;
	Cl->Mvir=MVIR;
	Cl->rdecay=0.1*RVIR;
	Cl->rmin=1.0e-6*RS;
	Cl->rmax=RVIR+3.0*Cl->rdecay;
	Cl->r_an = R_AN * RS;

	Cl->x[0] = Cl->x[1] = Cl->x[2] = 0.0;
	Cl->v[0] = Cl->v[1] = Cl->v[2] = 0.0;
	Cl->Ndm=Ndm;
		
	//First find the particle mass and calculate Eddingtons formula for the main cluster:
	CreateGrid(Gr,Cl );
	Mass = Gr[NG-1].Mass/Ndm;
	Eddington(Gr,Cl);//Eddingtons formula is now available via Cl->E[i] & Cl->DF[i]
	
	CreateParticles(P,Gr,Cl);

	for(i=0;i<Ndm;i++)//Particle ID's
	{
		P[i].ID=i+1;
	}


	WriteGadget2(P, Ndm,argv[1]);

	return 0;
}

double GetMass(double Mmin, double Mmax)
{
	double M,u;
	
	if(Mmax<Mmin)
	printf("Something is wrong in GetMass\t%G\t%G\n",Mmin,Mmax);
	
	//accept reject...
	do
	{
		u=uniform()*pow(Mmin,-1.9);
		M = uniform()*(Mmax-Mmin)+Mmin;
	}
	while(u>pow(M,-1.9));

	//printf("GetMass:\t%G\t%G\t%G\n",M,Mmin,Mmax);	
	if(M<Mmin || M>Mmax)
	printf("Something is wrong in GetMass\t%G\t%G\t%G\n",M,Mmin,Mmax);
	
	return M;
}

void CreateCluster(struct Particle *P,  struct Grid *Gr,  struct Params *Cl)
{
	int i;
	
	CreateGrid(Gr,Cl );
	Eddington(Gr,Cl);
	CreateParticles(P, Gr,Cl);
	
	for(i=0;i<Cl->Ndm;i++);
	{
		P[i].Pos[0]+=Cl->x[0];
		P[i].Pos[1]+=Cl->x[1];
		P[i].Pos[2]+=Cl->x[2];
		
		P[i].Vel[0]+=Cl->v[0];
		P[i].Vel[1]+=Cl->v[1];
		P[i].Vel[2]+=Cl->v[2];
	}
}


void FindPosAndVel(double *Pos, double *Vel, struct Grid * Gr, struct Params * Cl)
{
	int Ndm;
	struct Particle *Pi;	
	
	Pi = (struct Particle *) malloc(sizeof(struct Particle));

	/*
	 * This is very nasty... Cl->Ndm is set to 1 so only 1 particles is created
	 * when CreateParticles() is called. Cl->Ndm gets its correct original value before the function ends
	 */
	Ndm = Cl->Ndm;
	Cl->Ndm = 1;
	
	//create a single particle and copy its position and velocity to Pos and Vel
	CreateParticles(Pi, Gr,Cl); 
	Pos[0] = Pi[0].Pos[0];
	Pos[1] = Pi[0].Pos[1];
	Pos[2] = Pi[0].Pos[2];

	Vel[0] = Pi[0].Vel[0];
	Vel[1] = Pi[0].Vel[1];
	Vel[2] = Pi[0].Vel[2];
	
	Cl->Ndm = Ndm;	
	free(Pi);
}


void CreateGrid(struct Grid * Gr, struct Params * Cl )
{
	int i;
	double logdr,Sum;
	
	logdr = (log10(Cl->rmax)-log10(Cl->rmin)) / (1.0*NG);
	for(i=0;i<NG;i++)
	{
			Gr[i].R = pow(10.,log10(Cl->rmin)+i*logdr);
			Gr[i].Rho = rho(Gr[i].R, Cl);
	}
	
	for(i=0;i<NG;i++)
	{
		if(i==0)
		Gr[i].Mass = Gr[i].Rho * 4./3. * PI * Gr[i].R * Gr[i].R * Gr[i].R;
		else
		Gr[i].Mass = Gr[i-1].Mass + 4.*PI * 0.5*( Gr[i-1].R * Gr[i-1].R + Gr[i].R * Gr[i].R ) * (Gr[i].R-Gr[i-1].R) * 0.5*(Gr[i-1].Rho+Gr[i].Rho);
	}

	for(i=NG-1;i>=0;i--)
	{
		if(i==NG-1)
			Gr[i].Phi = 0.0;		
		else
			Gr[i].Phi = Gr[i+1].Phi - G * 0.5 *( Gr[i].Mass / (Gr[i].R*Gr[i].R) + Gr[i+1].Mass / (Gr[i+1].R*Gr[i+1].R) ) * (Gr[i+1].R - Gr[i].R);
	}	
	

	Sum = 0.0;
	for(i=NG-1;i>=0;i--)
	{
		if(i==NG-1)
		{
			// Do nothing
		}
		else
		{
			Sum += G * Gr[i].Mass * Gr[i].Rho / (Gr[i].R*Gr[i].R) *  (Gr[i+1].R - Gr[i].R);
		}
		
		Gr[i].Sigma2r = 1./ Gr[i].Rho * Sum;
	}
}


void Eddington( struct Grid *Gr, struct Params *Cl)
{
	int i,j,N1;
	double midp,umin,du,uj;
	double psi, dens,Pmax,vmax,v;

	for(i=0;i<N_EDD;i++)
	{
		Cl->E[i]=((double) i)/((double) N_EDD)*(-Gr[0].Phi);
	}
	
	Cl->DF[0] = 0.0;
	
	
	for(i=1;i<N_EDD;i++)
	{
		umin = atan( 1. / sqrt(Cl->E[i]) );
		
		N1 =300;
		midp=0.0; //midpoint integration

		du = (PI/2.0 - umin) / ( (double) N1 );
		
		for(j=0;j<N1;j++)
		{
			uj = umin + (j+0.5)*du;
			midp += 2./pow( sin(uj) ,2) * dnu2dpsi2( Cl->E[i] - 1./pow(tan(uj),2) ,Gr,Cl);
		}
			
		midp /= sqrt(8.)*PI*PI/du;
		
		//midp+= - 1./(sqrt(8.)*PI*PI) * 1./sqrt(Cl->E[i]) * drhodr(Gr[NG-1].R,Cl) * Gr[NG-1].R*Gr[NG-1].R / G / Gr[NG-1].Mass;

		Cl->DF[i] = midp;

	}

	
	double theta,a, r_an=Cl->r_an;
	for(i=0;i<NG;i+=200)
	{		
			psi = - Gr[i].Phi;
			dens = Gr[i].Rho;

			vmax = sqrt(2.*psi);
		
		
		Pmax=0;
		
		
		if(Gr[i].R > 0.1*r_an)
		{
			for(a=0.0;a<psi; a=a+psi/500.)
			{
				for(theta=0.0;theta<=PI/2.0;theta=theta+0.005*PI/2.0)
				{
					v=sqrt(
						   2.0*(psi-a) / (1.+Gr[i].R*Gr[i].R/r_an/r_an*sin(theta)*sin(theta))
						   );
				
					if(    2.*PI*v*v*sin(theta) * GetDF(psi-0.5*v*v*(1.+Gr[i].R*Gr[i].R/r_an/r_an*sin(theta)*sin(theta)), Cl)  /dens  > Pmax )
					{
						Pmax = 2.*PI*v*v*sin(theta) * GetDF(psi-0.5*v*v*(1.+Gr[i].R*Gr[i].R/r_an/r_an*sin(theta)*sin(theta)), Cl)  /dens;
//						Gr[i].v_mp = v;
//						Gr[i].theta_mp = theta;
					}
				}
			}
		}
		
		
		theta=PI/2.0;
		
		for(v=0.0;v<vmax;v=v+0.00001*vmax)
		{
			if(    2.*PI*v*v*sin(theta) * GetDF(psi-0.5*v*v*(1.+Gr[i].R*Gr[i].R/r_an/r_an*sin(theta)*sin(theta)),Cl)  /dens  > Pmax )
			{
				Pmax = 2.*PI*v*v*sin(theta) * GetDF(psi-0.5*v*v*(1.+Gr[i].R*Gr[i].R/r_an/r_an*sin(theta)*sin(theta)),Cl)  /dens;
//				Gr[i].v_mp = v;
//				Gr[i].theta_mp = theta;
			}
			
			if(psi-0.5*v*v*(1.+Gr[i].R*Gr[i].R/r_an/r_an*sin(theta)*sin(theta))<0.0)
				break;
		}

		Gr[i].MaxProb = 3.0*Pmax;
		
		for(j=i+1;j<i+200 && j<NG;j++)
		{
			Gr[j].MaxProb=Gr[i].MaxProb;
		}
	}
		
}


void CreateParticles(struct Particle *P, struct Grid *Gr, struct Params *Cl)
{
	int i,BinNo, Nmax, Nmin, Ndm;
	double random,Frac,TotalMass,Ri,Theta,Phi, v, uniform1, Psi;
	Ndm = Cl->Ndm;
	
	TotalMass = Gr[NG-1].Mass;
	
	for(i=0;i<Ndm;i++)
	{
		random = uniform()*TotalMass;
		
		while((random > Gr[NG-2].Mass) || (random < Gr[1].Mass))		
			random = uniform()*TotalMass;
		
		Nmax = NG-1;
		Nmin = 0;
		
		do//midpoint algorithm that finds the Bin number...
		{		
			BinNo=Nmin+(Nmax-Nmin)/2;	
			
			if(random<=Gr[BinNo].Mass)
				Nmax = BinNo;
			
			if (random>=Gr[BinNo].Mass)
				Nmin = BinNo;
		}
		while(  Nmax>Nmin+1  );
		
		BinNo = Nmin;

		Frac = (random - Gr[BinNo].Mass) / (Gr[BinNo+1].Mass - Gr[BinNo].Mass);
		
		Ri = Gr[BinNo].R + Frac * (Gr[BinNo+1].R - Gr[BinNo].R);
		
		//Sigma2i = Gr[BinNo].Sigma2r + Frac * (Gr[BinNo+1].Sigma2r - Gr[BinNo].Sigma2r);
		
		Theta = acos( 2.0 * uniform() - 1.);
		Phi = 2. * PI * uniform();
		
		P[i].Pos[0] = Ri * sin(Theta) * cos(Phi);
		P[i].Pos[1] = Ri * sin(Theta) * sin(Phi);
		P[i].Pos[2] = Ri * cos(Theta);
		
		Psi = - (Gr[BinNo].Phi + Frac * (Gr[BinNo+1].Phi - Gr[BinNo].Phi) );
		//dens = Gr[BinNo].Rho + Frac * (Gr[BinNo+1].Rho - Gr[BinNo].Rho) ;
	
	
		double r_an = Cl->r_an;
		double eta, ps, vr,vtheta, vphi;

		do{
			
			eta = uniform()*PI;
			v = uniform() *sqrt(2.*Psi);
			uniform1 = uniform() * Gr[BinNo].MaxProb;
				

			if(Psi-0.5*v*v*(   1.+pow( Ri/r_an *sin(eta),2)   ) < 0.0)
			{
				eta = uniform()*PI;
				v = uniform() * sqrt(2.*Psi);
				uniform1 = uniform() * Gr[BinNo].MaxProb;

			} 
		}
		while(uniform1 > 2.*PI*v*v*sin(eta) * 
			  GetDF(
					Psi-0.5*v*v*
						(   1.+pow( Ri/r_an *sin(eta),2)   )
					,Cl)/ rho(Ri,Cl)
			  );

		ps = 2. * PI * uniform();		

		vr = v * cos(eta);
		vtheta = v * sin(eta) * cos(ps);
		vphi = v * sin(eta) * sin(ps);

		
		P[i].Vel[0] = vr * sin(Theta) * cos(Phi) + vtheta * cos(Theta) * cos(Phi) - vphi * sin(Phi);
		P[i].Vel[1] = vr * sin(Theta) * sin(Phi) + vtheta * cos(Theta) * sin(Phi) + vphi * cos(Phi);
		P[i].Vel[2] = vr * cos(Theta) - vtheta * sin(Theta);
				
		P[i].ID = i;
		P[i].Mass = Mass;
	}
}




double dnu2dpsi2(double psi, struct Grid *Gr,struct Params *Cl)
{
	double r, dpsi, d2psi, temp, Frac, MassRa;

	int Nmax,Nmin,BinNo;
	
	Nmax = NG-1;
	Nmin = 0;
	
	do
	{		
		BinNo=Nmin+(Nmax-Nmin)/2;	
		
		if(-psi<=Gr[BinNo].Phi)
			Nmax = BinNo;
		
		if (-psi>=Gr[BinNo].Phi)
			Nmin = BinNo;
		
	}
	while(  Nmax>Nmin+1  );
	
	BinNo = Nmin;
	
	Frac = (-psi - Gr[BinNo].Phi) / (Gr[BinNo+1].Phi - Gr[BinNo].Phi);
	
	r = Gr[BinNo].R + Frac * (Gr[BinNo+1].R - Gr[BinNo].R);


	MassRa = Gr[BinNo].Mass + Frac * (Gr[BinNo+1].Mass - Gr[BinNo].Mass);

	dpsi = - G * MassRa / (r*r);
	d2psi = 2.*G*MassRa / (r*r*r) - 4. * PI * G * rho(r,Cl);
	

	temp = d2rho_Q(r,Cl)- drho_Q(r,Cl) / dpsi * d2psi;	

	temp /=  dpsi * dpsi;
	
	return temp;
	
}




double GetDF(double E, struct Params *Cl)
{

	double Nmax,Nmin,Frac;
	int BinNo;
	
	if(E<=0.0)
		return 0.0;
	else
	{
		Nmax = N_EDD-1;
		Nmin = 0;
		
		do//midpoint method
		{		
			BinNo=Nmin+(Nmax-Nmin)/2;	
			
			if(E<=Cl->E[BinNo])
				Nmax = BinNo;
			
			if (E>=Cl->E[BinNo])
				Nmin = BinNo;
			
		}
		while(  Nmax>Nmin+1  );
		
		BinNo = Nmin;
		
		if(BinNo==N_EDD || BinNo==N_EDD-1)
			BinNo=N_EDD-2;
	
		Frac = (E - Cl->E[BinNo]) / (Cl->E[BinNo+1] - Cl->E[BinNo] );	
	
		return Cl->DF[BinNo]+ Frac * (Cl->DF[BinNo+1] - Cl->DF[BinNo]);
	}
}

double GetDistFromPsi (double psi, struct Grid *Gr)
{
	double Nmax,Nmin,r,Frac;
	int BinNo;
	
	Nmax = NG-1;
	Nmin = 0;
	
	do
	{		
		BinNo=Nmin+(Nmax-Nmin)/2;	
		
		if(-psi<=Gr[BinNo].Phi)
			Nmax = BinNo;
		
		if (-psi>=Gr[BinNo].Phi)
			Nmin = BinNo;
		
	}
	while(  Nmax>Nmin+1  );
	
	BinNo = Nmin;
	
	Frac = (-psi - Gr[BinNo].Phi) / (Gr[BinNo+1].Phi - Gr[BinNo].Phi);
	
	r = Gr[BinNo].R + Frac * (Gr[BinNo+1].R - Gr[BinNo].R);

	return r;
}

double drho_Q(double r, struct Params *Cl)
{
	double r_an = Cl->r_an;
	return 2.*r*rho(r,Cl)/r_an/r_an + (1.+r*r/r_an/r_an) * drhodr(r,Cl);
}


double d2rho_Q(double r, struct Params *Cl)
{
	double r_an = Cl->r_an;
	return 2.*rho(r, Cl)/r_an/r_an + 4.*r*drhodr(r,Cl)/r_an/r_an + (1.+r*r/r_an/r_an) * d2rhodr2(r,Cl);
}


double rho(double r, struct Params *Cl)
{
	double temp, a = Cl->scale_length, C=Cl->conc, alpha=ALPHA, beta=BETA, gamma=GAMMA,Rho0=Cl->Rho0, Rvir=Cl->rvir, Rdecay=Cl->rdecay;
	
	if(r<=Cl->rvir)
		temp = Rho0 / pow(r / a, alpha) / pow(1. + pow(r / a, gamma), (beta - alpha) / gamma);
	else
		temp = Rho0 / pow(C, gamma) / pow(1.0 + pow(C, alpha), (beta - 1.0 * gamma) / alpha) * pow(r / Rvir, (-1.0 * gamma - 1.0 * beta * pow(C, alpha)) / (1.0 + pow(C, alpha)) + Rvir / Rdecay) * exp(-1.0 * (r - 1.0 * Rvir) / Rdecay);
	
	return temp;
}

double drhodr(double r, struct Params *Cl)
{
	double temp, a = Cl->scale_length, C=Cl->conc, alpha=ALPHA, beta=BETA, gamma=GAMMA,Rho0=Cl->Rho0, Rvir=Cl->rvir, Rdecay=Cl->rdecay;
	
	if(r<=Cl->rvir)
		temp = - Rho0 / pow(r / a, alpha) / pow(1. + pow(r / a, gamma), (beta - 1. * alpha) / gamma) * alpha / r - 1. * Rho0 / pow(r / a, alpha) / pow(1. + pow(r / a, gamma), (beta - 1. * alpha) / gamma) * (beta - 1. * alpha) * pow(r / a, gamma) / r / (1. + pow(r / a, gamma));
	else
		temp = Rho0 / pow(C, gamma) / pow(1.0 + pow(C, alpha), (beta - 1.0 * gamma) / alpha) * pow(r / Rvir, (-1.0 * gamma - 1.0 * beta * pow(C, alpha)) / (1.0 + pow(C, alpha)) + Rvir / Rdecay) * ((-1.0 * gamma - 1.0 * beta * pow(C, alpha)) / (1.0 + pow(C, alpha)) + Rvir / Rdecay) / r * exp(-1.0 * (r - 1.0 * Rvir) / Rdecay) - 1.0 * Rho0 / pow(C, gamma) / pow(1.0 + pow(C, alpha), (beta - 1.0 * gamma) / alpha) * pow(r / Rvir, (-1.0 * gamma - 1.0 * beta * pow(C, alpha)) / (1.0 + pow(C, alpha)) + Rvir / Rdecay) / Rdecay * exp(-1.0 * (r - 1.0 * Rvir) / Rdecay);
	
	return temp;
}

double d2rhodr2(double r, struct Params *Cl)
{
	double temp, a = Cl->scale_length, C=Cl->conc, alpha=ALPHA, beta=BETA, gamma=GAMMA,Rho0=Cl->Rho0, Rvir=Cl->rvir, Rdecay=Cl->rdecay;

	if(r<=Cl->rvir)
		temp = Rho0 / pow(r / a, alpha) / pow(1. + pow(r / a, gamma), (beta - 1. * alpha) / gamma) * alpha * alpha * pow(r, -2.) + 2. * Rho0 / pow(r / a, alpha) / pow(1. + pow(r / a, gamma), (beta - 1. * alpha) / gamma) * alpha * pow(r, -2.) * (beta - 1. * alpha) * pow(r / a, gamma) / (1. + pow(r / a, gamma)) + Rho0 / pow(r / a, alpha) / pow(1. + pow(r / a, gamma), (beta - 1. * alpha) / gamma) * alpha * pow(r, -2.) + Rho0 / pow(r / a, alpha) / pow(1. + pow(r / a, gamma), (beta - 1. * alpha) / gamma) * pow(beta - 1. * alpha, 2) * pow(pow(r / a, gamma), 2) * pow(r, -2) * pow(1. + pow(r / a, gamma), -2) - 1. * Rho0 / pow(r / a, alpha) / pow(0.1e1 + pow(r / a, gamma), (beta - 1. * alpha) / gamma) * (beta - 1. * alpha) * pow(r / a, gamma) * gamma * pow(r, -2) / (1. + pow(r / a, gamma)) + Rho0 / pow(r / a, alpha) / pow(1. + pow(r / a, gamma), (beta - 1. * alpha) / gamma) * (beta - 1. * alpha) * pow(r / a, gamma) * pow(r, -2.) / (1. + pow(r / a, gamma)) + Rho0 / pow(r / a, alpha) / pow(1. + pow(r / a, gamma), (beta - 1. * alpha) / gamma) * (beta - 1. * alpha) * pow(pow(r / a, gamma), 2) * pow(r, -2) * pow(0.1e1 + pow(r / a, gamma), -2) * gamma;
	else
		temp = Rho0 / pow(C, gamma) / pow(1.0 + pow(C, alpha), (beta - 1.0 * gamma) / alpha) * pow(r / Rvir, (-1.0 * gamma - 1.0 * beta * pow(C, alpha)) / (1.0 + pow(C, alpha)) + Rvir / Rdecay) * pow((-1.0 * gamma - 1.0 * beta * pow(C, alpha)) / (1.0 + pow(C, alpha)) + Rvir / Rdecay, 2.0) * pow(r, -2.0) * exp(-1.0 * (r - 1.0 * Rvir) / Rdecay) - 1.0 * Rho0 / pow(C, gamma) / pow(1.0 + pow(C, alpha), (beta - 1.0 * gamma) / alpha) * pow(r / Rvir, (-1.0 * gamma - 1.0 * beta * pow(C, alpha)) / (1.0 + pow(C, alpha)) + Rvir / Rdecay) * ((-1.0 * gamma - 1.0 * beta * pow(C, alpha)) / (1.0 + pow(C, alpha)) + Rvir / Rdecay) * pow(r, -2.0) * exp(-1.0 * (r - 1.0 * Rvir) / Rdecay) - 2.0 * Rho0 / pow(C, gamma) / pow(1.0 + pow(C, alpha), (beta - 1.0 * gamma) / alpha) * pow(r / Rvir, (-1.0 * gamma - 1.0 * beta * pow(C, alpha)) / (1.0 + pow(C, alpha)) + Rvir / Rdecay) * ((-1.0 * gamma - 1.0 * beta * pow(C, alpha)) / (1.0 + pow(C, alpha)) + Rvir / Rdecay) / r / Rdecay * exp(-1.0 * (r - 1.0 * Rvir) / Rdecay) + Rho0 / pow(C, gamma) / pow(1.0 + pow(C, alpha), (beta - 1.0 * gamma) / alpha) * pow(r / Rvir, (-1.0 * gamma - 1.0 * beta * pow(C, alpha)) / (1.0 + pow(C, alpha)) + Rvir / Rdecay) * pow(Rdecay, -2.0) * exp(-1.0 * (r - 1.0 * Rvir) / Rdecay);
	
	return temp;
}



double uniform(void) //returns number between 0 and 1
{
	double val=1.;
	
	val = gsl_rng_uniform(r1);
	return val;
}

double Gaussian(double Sigma)//returns a normal distributed number
{
	return gsl_ran_gaussian(r, Sigma);		
}

double dmax(double x, double y)
{
	if(x > y)
		return x;
	else
		return y;
}
