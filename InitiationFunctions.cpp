#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <iostream>
#include "Global.hpp"
#include "Functions.hpp"
#include "InitiationFunctions.hpp"


///////////////////////////////////
/////////Initiating matirces//////
/////////////////////////////////

void InitiateTypes(int types_mat[mating_types][G+1]){
    int i=0;
    int j=0;
    
    for (i=0; i<mating_types; i++)
        for (j=0; j<(G+1); j++)
            types_mat[i][j]=0;
}//InitiateTypes

//Initiate Nx2 matrix
void  InitiateX(int array[N][2], double initial_p){
	int i=0;
    int j=0;
    
    for(i=0; i<mating_types; i++)
        for(j=i*N/mating_types; j<(i+1)*N/mating_types; j++){
            array[j][0]=i;
            array[j][1]=0;
            //printf("j: %.d, i: %.d\n",j,i);
        }
}//InitiateX

void Introduce_mutant(int array[N][2], double initial_p){
    int i=0;
    
    for(i=0; i<N; i++)
        if(unifRand()<initial_p){
            array[i][1]=1;
        }
}

void InitiateP(double mat[runs][2]){
	int i=0;
	
	for (i=0; i<runs; i++){
		mat[i][0] = 0.0;
		mat[i][1] = 0.0;
	}
}//InitiateP

void IntiatepTraj(double mat[runs][maxSteps]){
	int i =0;
	int j=0;
	
	for (i=0; i<runs; i++) 
		for (j=0; j<maxSteps; j++) 
			mat[i][j]=0.0;
}//InitiatepTraj