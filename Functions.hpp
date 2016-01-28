double unifRand();

void copyAtoB(int a[N][2], int b[N][2]);

void SampleWithReplacement(int cells[N][2], int gens);

void switchType(int cells[N][2], double prob);

void shuffleList(int list[N][2]);

void shuffleList_2N(int list[2*N][2]);

void sum_types_left(int cell_mat[N][2], int sum[mating_types]);

int count_types_left(int mat[mating_types]);

int get_type(int type_sums[mating_types], int which);

void mate(int cells[N][2], double ratio[runs], int index);

double computeP(int cells[N][2]);

double switching_cost(double c, double k);

void grow(int cells[N][2], int cells_2[2*N][2]);

void back_to_N(int cells2N[2*N][2], int cells[N][2]);

void compute_types(int cells[N][2], int types[mating_types][G+1], int index);

void size_change(int cells[N][2]);

double rand_Normal (double mean, double sigma);

void mate_hard(int cells[N][2]);