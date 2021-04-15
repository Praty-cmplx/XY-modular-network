#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h>
#include<time.h>

#include "MT19937-64.c"
#define PI 3.1415926
#define size 512 // System size


double spin[size]; //spin state
int V[size][size]; //adjecency matrix
double p_add; //Probability of adding a spin to a cluster
#include "Modular_network_generator.c" //Adjacency matrix

int nm=16, av_k=14; //nm is number of modules, av_k is average degree
float r=0.001, C; //r is the ratio of inter modular edges and intra modular edges, C is temperature kT/J

void network_generator(int a, int b, float c); 
void quench();
double energy();
double flip(double a, double b);
void Wolff_run();
double mag();
double mod_mag();

void main()
{	
	double C_min=5.1; // Temperaure space start
	double C_max=6.9; // Temperaure space end
	double C_res=0.1; // step size
	
	double Mag, ModMag, NRG;
	printf("XY Simulation - Wolff - r = %f \n", r);
	int eqb_runs=51200, avg_runs=5120; //avg_runs steps are averaged and recorded
	int networks=10, conditions=10, mod_size=size/nm; //number of networks, initial conditions

	int seed = (unsigned) time(NULL);
	init_genrand64(seed);

	FILE *fp1, *fp2, *fp3, *flog;
	char fname1[50];
	char fname2[50];
	char fname3[50];
	double total_time;
	clock_t start, end;
	start = clock();

	flog=fopen("Network_MCsteps_record.dat", "w");
	fprintf(flog, "r = %f \nsize = %d \naverage degree = %d \nnumber of modules = %d \n", r, (int)size, av_k, nm);
	fprintf(flog, "eqbuilibration runs = %d \naverage runs = %d \ninitial conditions = %d \nno. of networks = %d \n", eqb_runs, avg_runs, conditions, networks);

	for (double t=C_min; t<C_max+C_res; t=t+C_res)
	{
		C=t;
		printf("C = %f \n", C);
		sprintf(fname1, "Global_r_%.5f_C_%.2f.dat", r, C);
		sprintf(fname2, "Modular_r_%.5f_C_%.2f.dat", r, C);
		sprintf(fname3, "Energy_r_%.5f_C_%.2f.dat", r, C);
		fp1=fopen(fname1, "w");
		fp2=fopen(fname2, "w");
		fp3=fopen(fname3, "w");
		fprintf(fp1,"#Next line is the temperature \n%.2f \n", C);
		fprintf(fp2,"#Next line is the temperature \n%.2f \n", C);
		fprintf(fp3,"#Next line is the temperature \n%.2f \n", C);
		fprintf(fp1, "#Global magnetization\n");
		fprintf(fp2, "#Modular magnetization\n");
		fprintf(fp3, "#Energy\n");

		for (int i=0; i<networks; i++)
		{	
			network_generator(nm, av_k, r);
			printf("Network number - %d \n", i);
						
			for (int j=0; j<conditions; j++)
			{
				quench();
				for (int k1=0; k1<eqb_runs; k1++)
					Wolff_run();

				for (int k2=0; k2<avg_runs; k2++)
				{
					Wolff_run();
					Mag=mag();
					ModMag=mod_mag();
					NRG=energy();
					fprintf(fp1, "%f \n", Mag);
					fprintf(fp2, "%f \n", ModMag);
					fprintf(fp3, "%f \n", NRG);
				}
			}
		}

		fclose(fp1);
		fclose(fp2);
		fclose(fp3);
	}
	end=clock();
	total_time=((double) (end-start))/CLOCKS_PER_SEC;
	fprintf(flog, "Parallel : Time taken = %f \n", total_time);
	fclose(flog);
	printf("C : \n");
	printf("M_ar : \n");
}


void Wolff_run() //One step of Wolff algorithm for a given configuration of the system
{
	int node, current, nbr, x, stack[size];
	double ref, oldspin, newspin, temp[size];
	node = (int) floor(genrand64_real2()*size);
	if (node>size-1)
		printf("Error1\n");

	for (int i=0; i<size; i++)
		temp[i]=spin[i]; 
	
	stack[0]=node;
	temp[node]=-10;
	x=1;
	oldspin=spin[node];
	ref = genrand64_real2()*2.0*PI;
	newspin=spin[node]=flip(oldspin, ref);
	
	int h=0;
	while(x)
	{
		current=stack[x-1];
		x=x-1;
		for (int i=0; i<size; i++)
		{
			if (V[current][i]==1 && (temp[i]>-9)) 
			{
				nbr=i;
				double S = 2*cos(spin[nbr]-ref)*cos(spin[current]-ref)/C;
				p_add=1.0-exp(S);
				if ((S<0) && (genrand64_real1()<p_add))
				{
					x=x+1;
					stack[x-1]=nbr;
					spin[nbr]=flip(spin[nbr], ref);
					temp[nbr]=-10;
				}
			}
		}
	}
}

double flip(double x, double r) //Wolff algorithm flip
{
	double newspin, r_prime;
	r_prime = fmod(r+PI/2.0, 2*PI);
	newspin = fmod((x+2*(r_prime-x))+2*PI, 2*PI); 
	return(newspin);
}

double energy() //Computes energy of the system
{
	double H=0;
	for (int i=0; i<size; i++)
	{
		for (int j=i; j<size; j++)
		{
			if (V[i][j]==1)
			{
				H=H-cos(spin[i]-spin[j]);
			}
		}
	}
	return(H);
}

void quench() //Quenched state intializer
{
	for (int i=0; i<size; i++)
		spin[i]=genrand64_real2()*2.0*PI;
}

double mag() //Computes the absolute order parameter
{	
	double complex tmp=0;
	for (int i=0; i<size; i++)
	{
		tmp=tmp + cexp(I*spin[i]);
	}
	double avg = cabs(tmp)/size;
	return (avg);
}

double mod_mag() //Computes the absolute order parameter at modular level
{
	double mod_size = (double)size/(double)nm;
	int initial=0, final=mod_size;
	double avg=0;
	double complex tmp;
	for (int i=0; i<nm; i++)
	{	
		tmp=0;
		for (int j=initial; j<final; j++)
		{
			tmp = tmp + cexp(I*spin[j]);
		}
		avg = avg + cabs(tmp);
		initial=initial+mod_size;
		final = final+mod_size;
	}
	avg=avg/size;
	return (avg);
}
