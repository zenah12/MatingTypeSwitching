#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <cmath>
#include <iostream>
#include "Global.hpp"
#include "Functions.hpp"
#include "InitiationFunctions.hpp"
#include "WriteFiles.hpp"

using namespace std;

//declare matrices
int types[mating_types][G+1];           //keeps track of mating type numbers over vegetative growth
int x[N][2];                            //x has the mating type and presence or absence of mating type gene for each cell
int x_double[2*N][2];                   //a matrix that can keep the genotypes of cells present in the population during vegetative growth
double p[runs][2];                      //p has the final frequency of the switching gene at the end of each reapeat (needed for prob of fixation)
double ratio_mated[runs];               //the ratio of mated cells in each run
double pTrajectory[runs][maxSteps];     //for each run, keep track of the frequency of the siwtching gene over time for each run
//declare parameters
int g=G;                                //number of vegetative rounds, g
double pswitch = 0.5;                   //switching rate, ps
int length = maxSteps*(g+1)+1;
double p_initial=0.05;                  //initial frequency of the mutant, qo
int mut_time = g-5;                     //time of intoducing the mutant
double exponent = 0.0;                  //exponent of cost function (k in the model)
double cost=0.0;                        //switching cost (c in the model)
int k=0;



//////////////////////////////////////////////////////////////////////////////
////////////////////////// M A I N   P R O G R A M //////////////////////////
////////////////////////////////////////////////////////////////////////////
int main(int argc, char ** argv) {
    
	//srand((unsigned)time(0));
	srand(time(NULL));

	int t=0;                //time step, generation
	int h=0;
	//int k=0;
    int i=0;
    int j=0;
    int count=0;
    
    InitiateTypes(types);
    InitiateP(p);
    IntiatepTraj(pTrajectory);
    
    for(h=0; h<runs; h++){//runs is the repeats used for each set of parameters
        InitiateX(x, p_initial);
        t=0;
        pTrajectory[h][t] = computeP(x);
        count=0;
        p_over_g[0] = computeP(x);
        compute_types(x, types, count);
        //have mut_time veg growth steps
        for(i=0; i<(mut_time-1); i++){
            grow(x,x_double);           //the population size doubles through mitosis
            back_to_N(x_double,x);      //the population size returns to carrying capacity
            count++;
            //compute_types(x, types, count);   //use when keeping track of the number of the frequency of mating types over time
            p_over_g[count] = computeP(x);
        }
        //introduce mutant
        Introduce_mutant(x, p_initial);
        //have another g-mut_time veg steps
        for(i=0; i<(g-mut_time); i++){
            grow(x,x_double);
            back_to_N(x_double,x);
            count++;
            //compute_types(x, types, count);
        }
        t++;

        pTrajectory[h][t] = computeP(x);

        //enter normal loop starting from mating
        while(t<maxSteps && pTrajectory[h][t]!=1 && pTrajectory[h][t]!=0){
            mate(x, ratio_mated, h);
            //mate_hard(x);  §§         //this function imposes the need for speedy mating and can be used instead of mate()
            for(i=0; i<g; i++){
                grow(x,x_double);
                back_to_N(x_double,x);
            }
            pTrajectory[h][t+1] = computeP(x);
            //compute_types(x, types, count);
            count++;
            t++;
        }
        p[h][0] = computeP(x); //p will be either 0 or 1 (fixed or lost)
        p[h][1] = t*1.0; //time to fixation or extinction
    }//for h
    
    //WriteRatioMated(ratio_mated,1);
    WritePfix(p, j);
    //Write_p(p_over_g, j);
    //WriteP_trajectory(p_traj_with_g);
    //WriteTypes(types, j);

	return 0;
}



