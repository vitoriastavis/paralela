/*
	Trabalho 1 - OpenMP
	Vitoria Stavis de Araujo
	GRR20200243
	
	mpicc mpi_hello_world.c -o hello-world  
	mpirun -np 5 ./hello-world

	nao funciona pra mais 50000
	e nao ta indo mais threads

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"

#ifndef max
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#define STD_TAG 0

typedef unsigned short mtype;

/* Read sequence from a file to a char vector.
 Filename is passed as parameter */

char* read_seq(char *fname)
{
	//file pointer
	FILE *fseq = NULL;
	//sequence size
	int size = 0;
	//sequence pointer
	char *seq = NULL;
	//sequence index
	int i = 0;

	//open file
	fseq = fopen(fname, "rt");
	
	if (fseq == NULL )
	{
		printf("Error reading file %s\n", fname);
		exit(1);
	}

	//find out sequence size to allocate memory afterwards
	fseek(fseq, 0L, SEEK_END);
	size = ftell(fseq);
	rewind(fseq);

	//allocate memory (sequence)
	seq = (char *) calloc(size + 1, sizeof(char));
	if (seq == NULL )
	{
		printf("Erro allocating memory for sequence %s.\n", fname);
		exit(1);
	}

	//read sequence from file
	while (!feof(fseq))
	{
		seq[i] = fgetc(fseq);
		if ((seq[i] != '\n') && (seq[i] != EOF))
			i++;
	}
	//insert string terminator
	seq[i] = '\0';

	//close file
	fclose(fseq);

	//return sequence pointer
	return seq;
}

// print matrix
void printMatrix(char * seqA, char * seqB, int ** scoreMatrix,  int sizeA,
		 int sizeB) {

	int i, j;

	//print header
	printf("Score Matrix:\n");
	printf("========================================\n");

	//print LCS score matrix al with sequences

	printf("    ");
	printf("%5c   ", ' ');

	for (j = 0; j < sizeB; j++)
		printf("%5c   ", seqB[j]);

	printf("\n");

	for (i = 0; i < sizeA+1; i++)
	{
		if (i == 0)
			printf("    ");
		else
			printf("%c   ", seqA[i - 1]);
		for (j = 0; j < sizeB + 1; j++) {
			printf("%5d   ", scoreMatrix[i][j]);
		}
		printf("\n");
	}
	printf("========================================\n");
}

// return the position of a character in a sequence, -1 if not found
int get_idx(char* str, int len, char c)
{
	int i;
	
    for(i = 0; i < len; i++)
    {		
        if(str[i] == c)
        {
            return i;
        }
    }

    return -1; 
}

// return a alphabet sequence of the common characters between two sequences
char* get_seq_c(int* size_c, char* seq_a, char* seq_b, int size_a, int size_b)
{
	int i, index, size;

	char* seq_c = (char*) malloc((size_a + size_b) * sizeof(char));

	seq_c[0] = seq_a[0];

	size = 1;

	for(i = 1; i < size_a; i++)
	{
		for(index = 0; index < size; index++)
		{
			if (seq_c[index] == seq_a[i]) break;
		}		
		if (index == size)
		{
			seq_c[index] = seq_a[i];
			size++;
		}
	}

	for(i = 0; i < size_b; i++)
	{
		for(index = 0; index < size; index++)
		{
			if (seq_c[index] == seq_b[i]) break;
		}

		if (index == size)
		{
			seq_c[index] = seq_b[i];
			size++;
		}
	}

	seq_c = (char*) realloc(seq_c, size * sizeof(char));

	*size_c = size;

	return (char*) realloc(seq_c, (size + 1) * sizeof(char));	
}

// calculate p_matrix
void calc_p(mtype* p_matrix, char* seq_b, char* seq_c,  int size_b, int size_c, int rank, int chunk_size)
{	
	// printf("teste 0asd \n"); 

	char c_recv[chunk_size];
	// printf("teste p c recv \n");


    mtype p_recv[chunk_size*(size_b + 1)];


	// printf("teste p p recv \n");
	// printf("teste p 1 \n");
	//Scatter the char array chunks by sending each process a particular chunk
    MPI_Scatter(seq_c, chunk_size, MPI_CHAR, &c_recv, chunk_size, MPI_CHAR, 0, MPI_COMM_WORLD);
   	
	// printf("teste p 2 \n");
	//Scatter the char array chunks by sending each process a particular chunk
    MPI_Scatter(p_matrix, chunk_size * (size_b + 1), MPI_UNSIGNED_SHORT, &p_recv, chunk_size*(size_b+1), MPI_UNSIGNED_SHORT, 0, MPI_COMM_WORLD);
    
	// printf("teste p 3\n");
	// Broadcast the whole b  array to everybody
    MPI_Bcast(seq_b, size_b, MPI_CHAR, 0, MPI_COMM_WORLD);

    
	for(int i = 0; i < chunk_size; i++)
	{	
		for(int j = 1; j < size_b + 1; j++)
		{
			if (seq_b[j - 1] == c_recv[i]) 
			{		
				p_recv[(i * (size_b + 1)) + j] = j;		
			}
			else
			{				
				p_recv[(i * (size_b + 1)) + j] = p_recv[(i * (size_b + 1)) + j - 1];   
			}
		}
	}  
	// printf(" teste p final \n");
	//now gather all the calculated values of P matrix in process 0
    MPI_Gather(p_recv, chunk_size * (size_b + 1), MPI_UNSIGNED_SHORT, p_matrix, chunk_size*(size_b+1), MPI_UNSIGNED_SHORT, 0, MPI_COMM_WORLD);
}

// calculate s_matrix according to p_matrix
mtype lcs_mpi(int size_a, int size_b, int size_c, mtype** s_matrix, mtype* p_matrix, char* seq_a, char* seq_b, char* seq_c, int rank, int chunk_size)
{
	// m é do b e n é do a eu acho, e o size_c?
    MPI_Bcast(p_matrix, (size_c * (size_a + 1)), MPI_UNSIGNED_SHORT, 0, MPI_COMM_WORLD);
	
    for (int i = 1; i < size_b + 1; i++)
    {
		// a principio m é do b
		// mas no mpi ta usando m no for e A aqui embaixo
		int c = get_idx(seq_c, size_c, seq_a[i - 1]);

		// printf("acho q eh isso \n");
		mtype s_recv[chunk_size]; 
		// printf("acho q eh isso 2 \n");

		// Broadcast the  whole B  array to everybody
		// mas pq DP[i] e isso é o S?
		MPI_Scatter(s_matrix[i], chunk_size, MPI_UNSIGNED_SHORT, &s_recv, chunk_size, MPI_UNSIGNED_SHORT, 0, MPI_COMM_WORLD);
		// printf("teste lcs scatter s \n");
		int start_id = (rank * chunk_size);
		int end_id = (rank * chunk_size) + chunk_size;
		// printf("teste lcs starta nd end \n");

        int t, s;

        for (int j = start_id; j < end_id; j++)
		{		
			if (j == start_id && rank == 0) 
				j = j + 1;
		
			t = (0 - p_matrix[(c * (size_a + 1)) + j]) < 0;
            s = (0 - (s_matrix[i - 1][j] - (t * s_matrix[i-1][p_matrix[(c*(size_a+1))+j]-1])));
            s_recv[j - start_id] = ((t^1) || (s^0)) * (s_matrix[i-1][j]) + (!((t^1)||(s^0))) * (s_matrix[i-1][p_matrix[(c*(size_a+1))+j]-1]+1);
		}	

		// printf("teste lcs for \n");
		MPI_Allgather(s_recv, chunk_size, MPI_UNSIGNED_SHORT, s_matrix[i], chunk_size, MPI_UNSIGNED_SHORT, MPI_COMM_WORLD);
		// printf("teste lcs gather \n");
	}

	return s_matrix[size_b][size_a]; 
}


int main(int argc, char ** argv)
{
	
	/*
	if(argc <= 1)
	{
        printf("Error: No input file specified! Please specify the input file, and run again!\n");
        return 0;
    }
	*/
	
    MPI_Init(&argc, &argv);

	// Declare process-related vars
    //     // and initialize MPI
    int rank;
    int num_procs;
    int chunk_size_p, chunk_size_s;//chunk_size for P matrix and DP matrix

    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 		//grab this process's rank
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs); 	//grab the total num of processes

    double time_total, seq_time,  begin, end;	
	begin = MPI_Wtime();	

	// sequence pointers for both sequences
	char *seq_a, *seq_b, *seq_c;

	seq_c = NULL;	

	// sizes of both sequences
	int size_a, size_b, size_c;

	// score result
    mtype res_par;
	
	

	// read both sequences
	seq_a = read_seq("A50000.in");
	seq_b = read_seq("B50000.in");

	// find out sizes
	size_a = strlen(seq_a);
	size_b = strlen(seq_b);

	// alphabet sequence c
	seq_c = get_seq_c(&size_c, seq_a, seq_b, size_a, size_b);	

	chunk_size_p = (size_c / num_procs);
    chunk_size_s = ((size_b + 1) / num_procs);

	
	if(rank == 0)
	{
		printf("rank %d \n", rank);
		printf("chunk_p: %d chunk_dp: %d procs: %d\n", chunk_size_p, chunk_size_s, num_procs);
	}	

    // allocate and initiate LCS score matrix 	
	// printf("teste 1 \n");
	mtype** s_matrix;
	mtype* p_matrix;
	s_matrix = (mtype **)malloc((size_a+1) * sizeof(mtype *));
    for(int k = 0; k < size_a + 1; k++)
    {
        s_matrix[k] = (mtype *)calloc((size_b + 1), sizeof(mtype));
    }

    for(int k = 0;k < size_a + 1; k++)
    {
        for(int l=0;l<size_b+1;l++)
        {
            s_matrix[k][l]=0;
        }
    }
	// printf("teste 2 aloquei s \n");
	// allocate and initiate p_matrix in vector form
	p_matrix = (mtype *)malloc((size_c * (size_b + 1)) * sizeof(mtype));
	// printf("teste 3 aloquei p \n");

	// LCS parallel algorithm
	calc_p(p_matrix, seq_b, seq_c, size_b, size_c, rank, chunk_size_p);

	// printf("teste 4 calc p \n");

	res_par = lcs_mpi(size_a, size_b, size_c, s_matrix, p_matrix, seq_a, seq_b, seq_c, rank, chunk_size_s);

	// printf("teste 5 lcs mpi \n");
    end = MPI_Wtime();
	time_total = end - begin;	   
	   
    if (rank == 0)
    {
        printf("res par: %d, time par: %f \n", res_par, time_total);
    }

	
    //printf("%f   |   %f%% \n", time_total, seq*100/time_total);
	//printf("%f   |   %f%% \n", 1-(seq*1/time_total), seq*1/time_total);
	//printf("%f\n", time_total);
	//printf("%f\n", seq/time_total);

	free(p_matrix);
	free(s_matrix);

	// Shutdown MPI 
    MPI_Finalize();
    return 0;
}

