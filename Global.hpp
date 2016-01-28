#ifndef GLOBAL
#define GLOBAL //Global.hpp

//dimensions
#define N (5000)             //number of cells

//Define global fixed parameters
#define maxSteps (2500)     //maximum number of time steps per iteration
#define max (1000)
#define runs (10000)        //number of repeats for a single parameter set for averaging
#define repeats (1)         //number of repeats to average over pfixation
#define mating_types (10)   //the number of mating types
#define G (50)              //the number of vegetative growth rounds between sex

//define global parameters that can be modified in main
extern double pswitch;			//ps
extern int g;                  //vegetative generations
extern double p_initial;       //initial frequency of p
extern int mut_time;           //time at which mutation is intriduced
extern double cost;             //cost of switching
extern double exponent;         //exponent for cost function (k in the paper)
extern double prob_size;        //probability of having a size change in each run
extern double size_mean;        //variance in size sampling distribution
extern double size_var;         //variance in size sampling distribution

#endif
