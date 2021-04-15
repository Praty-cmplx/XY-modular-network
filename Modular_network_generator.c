#include <math.h>

void network_generator(int nm, int av_k, float r)
{
	//printf("net_seg_fault 1\n");
	int mod_size=(int)size/nm;
	float constant = (mod_size-1)+r*((nm-1)*mod_size);
	float rho1=av_k/constant;
	float rho2=av_k*r/constant;
	int initial=0; int final=mod_size;
	//Intra Modular Connections
	//printf("net_seg_fault 2\n");
	for (int k1=0; k1<size; k1++)
		for (int k2=0; k2<size; k2++)
			V[k1][k2]=0;
	//printf("net_seg_fault 3\n");
	for (int n=0; n<nm; n++)
	{
		for (int i=initial; i<final; i++)
		{
			for (int j=i+1; j<final; j++)
			{	
				if (genrand64_real1()<rho1)
				{
					V[i][j]=1;
					V[j][i]=1;
				}
			}	
		}
		initial=initial+mod_size;
		final=final+mod_size;
	}
	//printf("net_seg_fault 4\n");
	initial=0; final=mod_size;
	int cinitial=mod_size;
	int cfinal=size;
	for (int n=0; n<nm-1; n++)
	{
		//printf("net_seg_fault 5\n");
		for (int i=initial; i<final; i++)
		{
			//printf("net_seg_fault 6\n");
			for(int j=cinitial; j<cfinal; j++)
			{
				//printf("j = %d \n", j);
				//printf("net_seg_fault 7\n");
				if (genrand64_real1()<rho2)
				{
					//printf("net_seg_fault 8");
					V[i][j]=1;
					V[j][i]=1;
				}
			}
		}
		initial=initial+mod_size;
		final=final+mod_size;
		cinitial=cinitial+mod_size;
	}
	for (int k=0; k<size; k++)
		V[k][k]=0;
}

void network_printer()
{
	for (int i=0; i<size; i++)
	{	
		//printf("%d : ", i);
		for (int j=0; j<size; j++)
			printf("%d  ", V[i][j]);
		printf("\n");
			//if (V[i][j]==1 && i!=j)
			//	printf("%d  ", j);
		//printf("\n --------------------------------------\n");
	}
}