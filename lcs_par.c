/*
	Trabalho 1 - OpenMP
	Vitoria Stavis de Araujo
	GRR20200243
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#ifndef max
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

int n_threads = 0;

typedef unsigned short mtype;

/* Read sequence from a file to a char vector.
 Filename is passed as parameter */

char* read_seq(char *fname)
{
	//file pointer
	FILE *fseq = NULL;
	//sequence size
	mtype size = 0;
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

// matrix allocation
mtype ** allocate_matrix(int size_a, int size_b)
{
	int i;
	//Allocate memory for LCS score matrix
	mtype ** matrix = (mtype **) malloc((size_b + 1) * sizeof(mtype *));
	for (i = 0; i < (size_b + 1); i++)
		matrix[i] = (mtype *) calloc((size_a + 1), sizeof(mtype));
	return matrix;
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
void calc_p(mtype ** p_matrix, char* seq_a, char* seq_c,  int size_a, int size_c)
{	 	
    #pragma omp parallel for collapse(2) num_threads(n_threads) // UM POUCO MAIS RÁPIDO
	//#pragma omp parallel for num_threads(n_threads)
	for(int i = 0; i < size_c; i++)
	{
		//#pragma omp parallel for schedule(auto) num_threads(n_threads)
		for(int j = 2; j < size_a + 1; j++)
		{
			if (seq_a[j - 2] == seq_c[i]) 
			{		
				p_matrix[i][j] = j - 1;		
			}
			else
			{				
				p_matrix[i][j] = p_matrix[i][j - 1];   
			}
		}
	}   
}

// calculate s_matrix according to p_matrix
int lcs_parallel(int size_a, int size_b, int size_c, mtype ** s_matrix, mtype ** p_matrix, char* seq_a, char* seq_b, char* seq_c)
{
    int i, j, c;

	//#pragma omp parallel for num_threads(n_threads) RESULTADO ERRADO
    for (i = 1; i < size_b + 1; i++)
    {
		c = get_idx(seq_c, size_c, seq_b[i - 1]);

        #pragma omp parallel for schedule(auto) num_threads(n_threads)
		// #pragma omp parallel for schedule(static) num_threads(n_threads) private(j) DEMOROU MAIS
        for (j = 1; j < (size_a + 1); j++)
		{		
			if(seq_a[j - 1] == seq_b[i - 1])
			{
				s_matrix[i][j] = s_matrix[i - 1][j - 1] + 1;
			}
			else if(p_matrix[c][j] == 0)
			{
				s_matrix[i][j] = max(s_matrix[i - 1][j], 0);				
			}
			else
			{
				s_matrix[i][j] = max(s_matrix[i - 1][j], s_matrix[i - 1][p_matrix[c][j] - 1] + 1);
			}
		}	
	}

	return s_matrix[size_b][size_a]; 
}

int main(int argc, char ** argv)
{
    double time_total, seq_time,  begin, end;

	begin = omp_get_wtime();

	// sequence pointers for both sequences
	char *seq_a, *seq_b, *seq_c;

	seq_c = NULL;	

	// sizes of both sequences
	int size_a, size_b, size_c;

	// score result
    int res_par;

	n_threads = atoi(argv[1]);
	//printf("%d \n", n_threads);

	// read both sequences
	seq_a = read_seq("A50000.in");
	seq_b = read_seq("B50000.in");

	// find out sizes
	size_a = strlen(seq_a);
	size_b = strlen(seq_b);

	// alphabet sequence c
	seq_c = get_seq_c(&size_c, seq_a, seq_b, size_a, size_b);	

    // allocate and initiate LCS score matrix and p_matrix
	mtype ** s_matrix = allocate_matrix(size_a, size_b);	
	mtype ** p_matrix = allocate_matrix(size_a, size_c-1);

	// ending sequencial part	
	seq_time = omp_get_wtime() - begin;

	// LCS parallel algorithm
	calc_p(p_matrix, seq_a, seq_c, size_a, size_c);	  
	res_par = lcs_parallel(size_a, size_b, size_c, s_matrix, p_matrix, seq_a, seq_b, seq_c);
      
    end = omp_get_wtime();
	time_total = end - begin;	   
	   
	printf("res par: %d, time par: %f \n", res_par, time_total);
    //printf("%f   |   %f%% \n", time_total, seq*100/time_total);
	//printf("%f   |   %f%% \n", 1-(seq*1/time_total), seq*1/time_total);
	//printf("%f\n", time_total);
	//printf("%f\n", seq/time_total);
}

