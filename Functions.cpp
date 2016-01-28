#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <iostream>
#include "Global.hpp"
#include "Functions.hpp"
#include "InitiationFunctions.hpp"

//////////////
//FUNTIONS //
/////////////

//create random number from 0 to 1
double unifRand()
{
    return rand() / double(RAND_MAX);
}//unifRand

void copyAtoB(int a[N][2], int b[N][2]){
    int i=0;
    for (i=0; i<N; i++) {
        b[i][0]=a[i][0];
        b[i][1]=a[i][1];
    }
}//copyAtoB

//SampleWithReplacement
void SampleWithReplacement(int cells[N][2], int gens){
    int cellsTemp[N][2];
	int randomNum =0;
    int i=0;
    int j=0;

    InitiateX(cellsTemp,0);
    
    for(i=0; i<gens; i++){
        for (j=0; j<N; j++) {
            randomNum = rand()%N;
            cellsTemp[j][0] = cells[randomNum][0];
            cellsTemp[j][1] = cells[randomNum][1];
        }
        copyAtoB(cellsTemp, cells);		
    }
}//SampleWithReplacement

void switchType(int cells[N][2], double prob){
    int i=0;
    double sample = 0.0;
    
    for (i=0; i<N; i++)
		//if a cells has the switching gene
        if(cells[i][1]==1){
            sample = unifRand();
            if(sample<=prob){
                //switch with probability prob
                if(cells[i][0]==0) cells[i][0]=1;
                else    cells[i][0]=0;
            }
        }
}//switchType

//shuffle list; used for sampling randomization
void shuffleList(int list[N][2]){
	int i=0;
	int temp1=0;
	int temp2=0;
	
	for (i=0; i<N; i++) {
		int rval = rand()%N;
		temp1 = list[i][0];
		temp2 = list[i][1];
		list[i][0] = list[rval][0];
		list[i][1] = list[rval][1];
		list[rval][0] = temp1;
		list[rval][1] = temp2;
	}
}//shuffleList

//shuffle list of size 2N to randomize sampling
void shuffleList_2N(int list[2*N][2]){
	int i=0;
	int temp1=0;
	int temp2=0;
	
	for (i=0; i<2*N; i++) {
		int rval = rand()%(2*N);
		temp1 = list[i][0];
		temp2 = list[i][1];
		list[i][0] = list[rval][0];
		list[i][1] = list[rval][1];
		list[rval][0] = temp1;
		list[rval][1] = temp2;
	}
}//shuffleList_2N

//updates matric sum so the ith entry is equal to the number of
//individuals with the ith mating type allele
void sum_types_left(int cell_mat[N][2], int sum[mating_types]){
    int i=0;
    int j=0;
    
    for (i=0; i<mating_types; i++) {
        sum[i]=0;
    }
    
    for(i=0; i<N; i++){
        for(j=0; j<mating_types; j++){
            if(cell_mat[i][0]==j){
                sum[j]++;
            }
        }//for j
    }//for i

}//sum_types_left

//it returns the total number of distinct mating
//types remaining in the population
int count_types_left(int mat[mating_types]){
    int i=0;
    int count=0;
    
    for(i=0; i<mating_types; i++) if(mat[i]>0) count++;
    return(count);
}

//returns raraest mating type if whcih=1 and commonest
//mating type if which = 2. Only teo mating types
//should be left by now
int get_type(int type_sums[mating_types], int which){
    int i=0;
    int max_int=0;
    int min_int=0;
    int left[2];
    int left_index[2];
    left[0]=-1;
    left[1]=-2;
    left_index[0]=-1;
    left_index[1]=-2;
    int count=0;
    
    for(i=0; i<mating_types; i++){
        if(type_sums[i]>0){//this should be true only once
            left[count]=type_sums[i];
            left_index[count]=i;
            count++;
        }
    }
    //if at 50-50
    if(left[0]==left[1]){
        max_int = left_index[0];
        min_int = left_index[1];
    }
    else if(left[0]>left[1]){
        max_int = left_index[0];
        min_int = left_index[1];
    }
    else if(left[0]<left[1]){
        max_int = left_index[1];
        min_int = left_index[0];
    }
    if(which==1) return(max_int);
    else return(min_int);
}//get_type

void mate(int cells[N][2], double ratio[runs], int index){
    int i=0;
    int index1=0;
    int index2=0;
    int sums[mating_types];
    int count =0; //keeps trakc of mated cells and is used as index in cellsMATED
    double cellsMATED[N][2];
    //start by randomizing list
    shuffleList(cells);
    //make a matrix to store mated cells
    for(i=0; i<N;i++){
        cellsMATED[i][0]= -50.0;
        cellsMATED[i][1] = -50.0;
    }
    //save the total of each type in sums[]
    sum_types_left(cells, sums);
    //while more than one type present
    while(count_types_left(sums)>1){ //&& count<N/10.0){
        //sample two cells at random
        index1 = rand()%N;
        index2 = rand()%N;
        //make sure chosen cells haven't already mated, and are not the same type
        while(cells[index1][0]<0 || cells[index2][0]<0 || cells[index1][0]==cells[index2][0]){
            index1 = rand()%N;
            index2 = rand()%N;
        }//while
       //mate if different
        //add to mated list and recombine
        if(unifRand()<0.5)
        {
            cellsMATED[count][0] = cells[index1][0];
            cellsMATED[count][1] = cells[index1][1];
            cellsMATED[count+1][0] = cells[index2][0];
            cellsMATED[count+1][1] = cells[index2][1];
        }//if
        else
        {
            cellsMATED[count][0] = cells[index1][0];
            cellsMATED[count][1] = cells[index2][1];
            cellsMATED[count+1][0] = cells[index2][0];
            cellsMATED[count+1][1] = cells[index1][1];
        }//else
        count = count+2;
//        }//mate if different
        //set their mating type to -1 to indicated mated
        cells[index1][0]=-1.0;
        cells[index2][0]=-1.0;
        //get new sums
        sum_types_left(cells, sums);
    }//while
    //now count cells have mated and one type left
    int rval = 0;
    //now ready to sample back to N IF count>0
    if(count>0)
    {
        for(i=0; i<N; i++)
        {
            rval = rand()%(count);
            cells[i][0] = cellsMATED[rval][0];
            cells[i][1] = cellsMATED[rval][1];
        }//for
    }//if count>0
    
    sum_types_left(cells, sums);
    ratio[index] = (1.0*count)/(1.0*N);
}//mate

//finds the freuqncy of the switch gene in the population
double computeP(int cells[N][2]){
	int i=0;
	int sum=0;
	
	for(i=0; i<N;i++) sum=sum+cells[i][1];
	
	return((1.0*sum+0.0)/(1.0*N+0.0));
}//computeP

//this function imposes one veg growth step with the population size doubling and
//allows for switching depending on the switchinbg probability
void grow(int cells[N][2], int cells_2[2*N][2]){
    int i=0;
    int count =0;
    int type_sums[mating_types];
    double rarest=0;
    double switching_events=0;
    int rand_index=-1;
    
    for (i=0; i<N; i++) {
        //make one identical copy
        cells_2[count][0] = cells[i][0];
        cells_2[count][1] = cells[i][1];
        if(cells[i][1]==0){//if cell does not have the capacity to switch both dauughter cells are clonal
            cells_2[count+1][0] = cells[i][0];
            cells_2[count+1][1] = cells[i][1];
        }
        else{//else one daughter cell will switch with given probability
            if(unifRand()<pswitch){//switch with prob (pswitch)
                switching_events = switching_events + 1.0;
                //choose an index from 0 to mating_types -1
                rand_index = rand()%(mating_types);
                //if the type chosen to switch to is the same as own type, choose another type
                while(rand_index==cells[i][0]){
                    rand_index = rand()%(mating_types);
                }
                if(type_sums[cells[i][0]]>type_sums[rand_index]) rarest = rarest+1.0;
                //switch to randomly chosen type
                cells_2[count+1][0] = rand_index;
            }//if
            //else do not switch
            else    cells_2[count+1][0] = cells[i][0];
        }//else
        cells_2[count+1][1] = cells[i][1];
        count=count+2;
    }//for
}//grow

double switching_cost(double cost, double exponent){
    return(1.0 - cost*pow(pswitch, exponent));
}//switching_cost

void back_to_N(int cells2N[2*N][2], int c[N][2]){
    int i=0;
//    int rval=0;
//    int count=0;
    
    //shuffle 2N population and choose first N
    shuffleList_2N(cells2N);
    for (i=0; i<N; i++) {
        c[i][0] = cells2N[i][0];
        c[i][1] = cells2N[i][1];
    }
    //use the commented code below for costly swithcinbg
//    while(count<N){
//        rval = rand()%(2*N);
//        if(cells2N[rval][1] == 1){//if cell has switching gene
//            if(unifRand()<switching_cost(cost, exponent)){//accept given cell only with a probability
//                c[count][0] = cells2N[rval][0];
//                c[count][1] = cells2N[rval][1];
//                count++;
//            }
//        }//if
//        else{//if no switching gene accept with prob 1
//            c[count][0] = cells2N[rval][0];
//            c[count][1] = cells2N[rval][1];
//            count++;
//        }//else
//    }//while
}//back_to_N

void compute_types(int cells[N][2], int types[mating_types][G+1], int index){
    int i=0;
    
    for(i=0; i<N; i++){
        //add one to the types matrix for the mating type of the ith cell
        types[cells[i][0]][index]++;
       
    }
}//compute_types

//this function imposes a change in the population size from N to size and then back to N
void size_change(int cells[N][2]){
    int size=-1;
    while(size<0 || size<1){
        size = int(rand_Normal(N, size_var));
    }
    //printf("size: %.d\n",size);
    int tempCells[size][2];
    int rval=0;
    int i=0;
    
    //get temporary population
    for (i=0; i<size; i++) {
        rval = rand()%(N);
        tempCells[i][0] = cells[rval][0];
        tempCells[i][1] = cells[rval][1];
    }
    //sample back to N
    for(i=0; i<N; i++){
        rval = rand()%(size);
        cells[i][0] = tempCells[rval][0];
        cells[i][1] = tempCells[rval][1];
    }
}

//rand_Normal normal distribution genator
double rand_Normal (double mean, double sigma){
    double U1, U2, W, mult;
    static double X1, X2;
    static int call = 0;
    
    if (call == 1)
    {
        call = !call;
        return (mean + sigma * (double) X2);
    }
    
    do
    {
        U1 = -1 + ((double) rand () / RAND_MAX) * 2;
        U2 = -1 + ((double) rand () / RAND_MAX) * 2;
        W = pow (U1, 2) + pow (U2, 2);
    }
    while (W >= 1 || W == 0);
    
    mult = sqrt ((-2 * log (W)) / W);
    X1 = U1 * mult;
    X2 = U2 * mult;
    
    call = !call;
    
    return (mean + sigma * (double) X1);
}//rand_Normal

//this function can be used to impose selction
//for speedy mating
void mate_hard(int cells[N][2]){
    int i=0;
    int index1=0;
    int index2=0;
    int sums[mating_types];
    int count =0; //keeps trakc of mated cells and is used as index in cellsMATED
    double cellsMATED[N][2];
    int freq_ps=0;
    //start by randomizing list
    shuffleList(cells);
    //make a matrix to store mated cells
    for(i=0; i<N;i++){
        cellsMATED[i][0]= -50.0;
        cellsMATED[i][1] = -50.0;
    }
    //save the sums of 0,1,2 in sums[]
    sum_types_left(cells, sums);

    //for(i=0; i<mating_types; i++) printf("%.d type: %.d\n", i, sums[i]);
    for (i=0; i<N; i++) freq_ps = freq_ps + cells[i][1];
    //while more than one type present
    while(count_types_left(sums)>1){ //&& count<N/10.0){
        //sample two cells at random
        index1 = rand()%N;
        index2 = rand()%N;
        //make sure chosen cells haven't already mated, and are not the same type
        while(cells[index1][0]<0 || cells[index2][0]<0){
            index1 = rand()%N;
            index2 = rand()%N;
        }//while
        if(cells[index1][0]!=cells[index2][0]){//if not the same type allow them to mate
            //add to mated list and recombine
            if(unifRand()<0.5)
            {
                cellsMATED[count][0] = cells[index1][0];
                cellsMATED[count][1] = cells[index1][1];
                cellsMATED[count+1][0] = cells[index2][0];
                cellsMATED[count+1][1] = cells[index2][1];
            }
            else
            {
                cellsMATED[count][0] = cells[index1][0];
                cellsMATED[count][1] = cells[index2][1];
                cellsMATED[count+1][0] = cells[index2][0];
                cellsMATED[count+1][1] = cells[index1][1];
            }
            count = count+2;
        }//if
        //set their mating type to -1 to indicated mated or dead because of selection and failure to mate
        cells[index1][0]=-1.0;
        cells[index2][0]=-1.0;
        //get new sums
        sum_types_left(cells, sums);
    }//while
    //now count cells have mated and one type left
    int rval = 0;
    //now ready to sample back to N IF count>0
    if(count>0)
    {
        for(i=0; i<N; i++)
        {
            rval = rand()%(count);
            cells[i][0] = cellsMATED[rval][0];
            cells[i][1] = cellsMATED[rval][1];
            //printf("%.d. chosen cell: %.3f, %.3f\n", i, cells[i][0], cells[i][1]);
        }//for
    }//if count>0
    
    sum_types_left(cells, sums);
    freq_ps=0;
    for (i=0; i<N; i++) freq_ps = freq_ps + cells[i][1];
}//mate_hard







