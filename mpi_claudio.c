/*
	Trabalho 1 - OpenMP
	Vitoria Stavis de seq_araujo
	GRR20200243
	
	mpicc mpi_hello_world.c -o hello-world  
	mpirun -np 5 ./hello-world

	versao claudio

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mpi.h"

// A utility function to find min of two integers 
int min(int a, int b) 
{ return (a < b)? a: b; } 
  
// A utility function to find min of three integers 
int min3(int a, int b, int c) 
{ return min(min(a, b), c);} 
  
// A utility function to find max of two integers 
int max(int a, int b) 
{ return (a > b)? a: b; } 

#define STD_TAG 0

typedef unsigned short mtype;

/* Read sequence from a file to a char vector.
 Filename is passed as parameter */

char* read_seq(char *fname)
{
	
	FILE *fseq = NULL;
	int size = 0;
	char *seq = NULL;

	int i = 0;
	
	fseq = fopen(fname, "rt");
	
	if (fseq == NULL )
	{
		printf("erro ao ler arquivo %s\n", fname);
		exit(1);
	}

	// encontrar tamanho da sequencia
	fseek(fseq, 0L, SEEK_END);
	size = ftell(fseq);
	rewind(fseq);

	// alocar memoria para sequencia

    // = (char *) malloc((size + 1) * sizeof(char));
	seq = (char *) calloc(size + 1, sizeof(char));
	if (seq == NULL )
	{
		printf("erro ao alocar memoria %s.\n", fname);
		exit(1);
	}

	// le a sequencia e adiciona o \0 no final
	while (!feof(fseq))
	{
		seq[i] = fgetc(fseq);
		if ((seq[i] != '\n') && (seq[i] != EOF))
			i++;
	}	
	seq[i] = '\0';

	fclose(fseq);

	return seq;
}

// print matrix
void printMatrix(char * seqseq_a, char * seqseq_b, int ** scoreMatrix,  int sizeseq_a,
		 int sizeseq_b) {

	int i, j;

	//print header
	printf("Score Matrix:\n");
	printf("========================================\n");

	//print Lseq_cS score matrix al with sequences

	printf("    ");
	printf("%5c   ", ' ');

	for (j = 0; j < sizeseq_b; j++)
		printf("%5c   ", seqseq_b[j]);

	printf("\n");

	for (i = 0; i < sizeseq_a+1; i++)
	{
		if (i == 0)
			printf("    ");
		else
			printf("%c   ", seqseq_a[i - 1]);
		for (j = 0; j < sizeseq_b + 1; j++) {
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

	*size_c = size+1;

	return (char*) realloc(seq_c, (size + 1) * sizeof(char));	
}

int lcs_parallel(char* s1, char *s2, int len_s1, int len_s2, int rank, int size) {
    int rows = len_s1 + 1;
    int cols = len_s2 + 1;

    int row, col;
    MPI_Status status;

    int dp[3][cols];

    for (int line=1; line<rows+cols; line++) {
        int curr_line = line % 3;
        int prev_line = (line-1) % 3;
        int prev_prev_line = (line-2) % 3;

        int start_col =  max(0, line-rows); 
        int count = min3(line, (cols-start_col), rows); 

        int start, end;
        if (count <= size) {
            start = rank;
            end = min(rank, count-1);
        } else {
            float block_len = (float)count / size;
            start = round(block_len*rank);
            end = round(block_len*(rank+1))-1;
        }
  
        for (int j=start; j<=end; j++) {
            row = min(rows, line)-j-1;
            col = start_col+j;

            if (row==0 || col==0) {
                dp[curr_line][col] = 0;
            }
            else if (s1[row - 1] == s2[col - 1]) {
                int upper_left = dp[prev_prev_line][col-1];
                dp[curr_line][col] = upper_left + 1;
            }
            else {
                int left = dp[prev_line][col-1];
                int up = dp[prev_line][col];
                dp[curr_line][col] = max(left, up);
            }

            if (j == start && rank > 0) {
                MPI_Send(&dp[curr_line][col], 1, MPI_INT, rank-1, 1, MPI_COMM_WORLD);
            }
            if (j == end && rank < size-1) {
                MPI_Send(&dp[curr_line][col], 1, MPI_INT, rank+1, 1, MPI_COMM_WORLD);
            }
        }

        int prev_index = start - 1;
        if (prev_index >= 0 && prev_index < count) {
            col = start_col+prev_index;
            MPI_Recv(&dp[curr_line][col], 1, MPI_INT, rank-1, 1, MPI_COMM_WORLD, &status);
        }

        int next_index = end + 1;
        if (next_index >= 0 && next_index < count) {
            col = start_col+next_index;
            MPI_Recv(&dp[curr_line][col], 1, MPI_INT, rank+1, 1, MPI_COMM_WORLD, &status);
        }
    }

    return dp[(rows+cols-1) % 3][cols-1];
}

void load_input(char **sa, int *sia, char **sb, int *sib, char **sc, int *sic, int rank) {
	char *seq_a, *seq_b, *seq_c;
	int size_a, size_b, size_c;	

	if(rank == 0)
	{
	    // ler sequencias a e b
        seq_a = read_seq("teste.in");
        seq_b = read_seq("teste2.in");
        size_a = strlen(seq_a);
        size_b = strlen(seq_b);

        // sequencia c
        seq_c = get_seq_c(&size_c, seq_a, seq_b, size_a, size_b);
        MPI_Bcast(&size_a, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&size_b, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&size_c, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
    else
    {
        MPI_Bcast(&size_a, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&size_b, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&size_c, 1, MPI_INT, 0, MPI_COMM_WORLD);
        seq_a = malloc(size_a+1);
        seq_b = malloc(size_b+1);
        seq_c = malloc(size_c+1);
    }

    MPI_Bcast(&seq_a[0], size_a, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(&seq_b[0], size_b, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(&seq_c[0], size_c, MPI_CHAR, 0, MPI_COMM_WORLD);

    *sa = seq_a;
    *sb = seq_b;
    *sc = seq_c;
	
    *sia = size_a;
    *sib = size_b;
    *sic = size_c;
}

int main(int argc, char ** argv)
{
	int rank;
    int num_procs;
    int chunk_size_p, chunk_size_s;		// chunk_size para matrizes P e S		
    mtype res_par;						// variavel do score

    MPI_Init(&argc, &argv);	
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 		// rank do processo
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs); 	// numero total de processos

	// variaveis para medir o tempo
    double time_total, seq_time, bs, begin, end;	
	begin = MPI_Wtime();
    
    char *seq_a, *seq_b, *seq_c;
	int size_a, size_b, size_c;	
    load_input(&seq_a, &size_a, &seq_b, &size_b, &seq_c, &size_c, rank);
	
	mtype** s_matrix;
	mtype* p_matrix;

	// alocar e iniciar matriz S 	
	s_matrix = (mtype **)malloc((size_a + 1) * sizeof(mtype *));
	for(int k = 0; k < size_a + 1; k++)
	{
		//s_matrix[k] = (mtype *) malloc((size_b + 1) * sizeof(mtype));
		s_matrix[k] = (mtype *)calloc((size_b + 1), sizeof(mtype));
	}

	for(int k = 0;k < size_a + 1; k++)
	{
		for(int l=0;l<size_b+1;l++)
		{
			s_matrix[k][l]=0;
		}
	}

	// alocar e iniciar a matriz p em forma de vetor
	p_matrix = (mtype *)malloc((size_c * (size_b + 1)) * sizeof(mtype));

	chunk_size_p = (size_c / num_procs);
	chunk_size_s = ((size_b + 1) / num_procs);	

	if(rank == 0)
	{
		bs = MPI_Wtime();
		seq_time = bs - begin;
		//printf("chunk_p: %d chunk_dp: %d procs: %d\n", chunk_size_p, chunk_size_s, num_procs);
	}	
    
	// LCS paralelo
	//calc_P(p_matrix, seq_b, size_b, seq_c, size_c, rank, chunk_size_p);
	//res_par = lcs_mpi(s_matrix, p_matrix, seq_a, seq_b, seq_c, size_b, size_a, size_c, rank, chunk_size_s);
    res_par = lcs_parallel(seq_a, seq_b, size_a, size_b, rank, num_procs);
      
	   
    if (rank == 0)
    {
		end = MPI_Wtime();
		time_total = end - begin;	
        //printf("res par: %d, time par: %f \n", res_par, time_total);
		printf("%f\n", (seq_time*100/time_total));
    }

	free(p_matrix);
	for(int l=0;l<size_b+1;l++)
    {
            free(s_matrix[l]);
    }
	free(s_matrix);

	// Shutdown MPI 
    MPI_Finalize();
    return 0;
}
