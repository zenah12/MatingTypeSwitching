#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <iostream>
#include "Global.hpp"
#include "Functions.hpp"
#include "InitiationFunctions.hpp"
#include "WriteFiles.hpp"

//saving mycells.position file


//write pair stats file
void WriteTypes(int array[mating_types][G+1], int ind){
    char name[250];
    
    FILE *pfile=NULL;
    
    sprintf(name, "TYPES_N%.d_g%.d_runs%.d_pswitch%.3f0_pintro%.3f_Tmut%.d_types%.d_%.d_var%.2f_mean%.3f_pSize%.3f.txt.txt", N, g, runs, pswitch, p_initial, mut_time, mating_types, ind,size_mean,size_var, prob_size);
    
    pfile = fopen(name, "w+");
    
    int i=0;
    int j=0;

    for (i=0; i<(G+1); i++){
        for(j=0; j<mating_types; j++){
            fprintf(pfile, "%.d   ", array[j][i]+1);
        }
        fprintf(pfile, "\n");
    }
    
    fclose(pfile);
}//WritePfix


//write pair stats file
void Write_p(double array[G], int ind){
    char name[250];
    
    FILE *pfile=NULL;
    
    sprintf(name, "p_N%.d_g%.d_runs%.d_pswitch%.3f0_pintro%.3f_Tmut%.d_types%.d_%.d_var%.2f_mean%.3f_pSize%.3f.txt.txt", N, g, runs, pswitch, p_initial, mut_time, mating_types, ind,size_mean,size_var, prob_size);
    
    pfile = fopen(name, "w+");
    
    int i=0;
    
    for (i=0; i<(G); i++){
        fprintf(pfile, "%.5f\n", array[i]);
    }
    
    
    
    fclose(pfile);
}//WritePfix


//write pair stats file
void WriteRatioMated(double array[runs], int index){
    char name[250];
    
    FILE *pfile=NULL;
    
    sprintf(name, "RM_N%.d_g%.d_runs%.d_pswitch%.3f0_pintro%.3f_Tmut%.d_types%.d_N'%.2f_pchange%.3f_cost%.5f_exponent%.3f_hard.txt", N, g, runs, pswitch, p_initial, mut_time, mating_types, size_mean,prob_size,cost, exponent);
    
    
    pfile = fopen(name, "w+");
    
    int i=0;
    for (i=0; i<runs; i++)
        fprintf(pfile, "%.5f\n", array[i]);
    
    fclose(pfile);
}//WritePfix




//write pair stats file
void WritePfix(double array[runs][2], int index){
	char name[250];
		
	FILE *pfile=NULL;
	  
	sprintf(name, "Pfix_N%.d_g%.d_runs%.d_pswitch%.3f0_pintro%.3f_Tmut%.d_types%.d_N'%.2f_pchange%.3f_cost%.5f_exponent%.3f_hard.txt", N, g, runs, pswitch, p_initial, mut_time, mating_types, size_mean,prob_size,cost, exponent);
	
	
	pfile = fopen(name, "w+");
		
	int i=0;
	for (i=0; i<runs; i++)
        fprintf(pfile, "%.5f0	%.5f0\n", array[i][0], array[i][1]);
	
	fclose(pfile);
}//WritePfix

//write pair stats file
void WriteP_trajectory(double array[maxSteps]){
	char name[250];
    
	FILE *pfile=NULL;
    
	sprintf(name, "p_traj_N%.d_g%.d_runs%.d_pswitch%.3f0_types%.d.txt", N, g, runs, pswitch, mating_types);
	
	pfile = fopen(name, "w+");
    
	int i=0;

    while(array[i]!=100){
        fprintf(pfile, "%.5f0\n", array[i]);
        i++;
    }
	
	fclose(pfile);
}//WritePfix










