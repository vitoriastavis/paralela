/*
	Trabalho 1 - OpenMP
	Vitoria Stavis de Araujo
	GRR20200243
	
	mpicc mpi_hello_world.c -o hello-world  
	mpirun -np 5 ./hello-world

	mpi3 mas dando certo só pra 1 ou 2 procs mais ou menos

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"

#ifndef max
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )

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


// se o tamanho da lcs for nao divisivel pelo numero de procs
// usra all gather v 
// ou criar colunas a mais com 
void printMatrixv(char * seqA, char * seqB, mtype * scoreMatrix,  int sizeA,
		 int sizeB)
        {

        int i, j;
    
        //print header
        printf("Score Matrix:\n");
        printf("======================================================================\n");

        printf("    ");
	    printf("%5c   ", ' ');

        for (j = 0; j < sizeB + 1; j++)
		    printf("%5c   ", seqB[j]);

        printf("\n");

        for(int i = 0; i < sizeA; i++)
        {
            if (i == 0)
			    printf("    ");
		    else
			    printf("%c   ", seqA[i - 1]);

            for(int j = 0; j < sizeB + 1; j++)
            {
                printf("%5d   ", scoreMatrix[(i*(sizeB+1))+j]);
            }
            printf("\n");
        }
        printf("===================================================================================\n\n\n");
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
	
	char c_recv[chunk_size];
    mtype p_recv[chunk_size*(size_b + 1)];
	
	//Scatter the char array chunks by sending each process a particular chunk
    //MPI_Scatter(seq_c, chunk_size, MPI_CHAR, &c_recv, chunk_size, MPI_CHAR, 0, MPI_COMM_WORLD);
   	
	//Scatter the char array chunks by sending each process a particular chunk
    //MPI_Scatter(p_matrix, chunk_size * (size_b + 1), MPI_UNSIGNED_SHORT, &p_recv, chunk_size*(size_b+1), MPI_UNSIGNED_SHORT, 0, MPI_COMM_WORLD);
    
	// Broadcast the whole b  array to everybody
    //MPI_Bcast(seq_b, size_b, MPI_CHAR, 0, MPI_COMM_WORLD);
    
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

	//now gather all the calculated values of P matrix in process 0
    MPI_Gather(p_recv, chunk_size * (size_b + 1), MPI_UNSIGNED_SHORT, p_matrix, chunk_size*(size_b+1), MPI_UNSIGNED_SHORT, 0, MPI_COMM_WORLD);
}

// calculate s_matrix according to p_matrix
mtype lcs_mpi(int size_a, int size_b, int size_c, mtype** s_matrix, mtype* p_matrix, char* seq_a, char* seq_b, char* seq_c, int rank, int chunk_size)
{	
    MPI_Bcast(p_matrix, (size_c * (size_a + 1)), MPI_UNSIGNED_SHORT, 0, MPI_COMM_WORLD);
	
    for (int i = 1; i < size_b + 1; i++)
    {
		// a principio m é do b
		// mas no mpi ta usando m no for e A aqui embaixo
		int c = get_idx(seq_c, size_c, seq_a[i - 1]);
		
		mtype s_recv[chunk_size]; 
		
		// Broadcast the  whole B  array to everybody	
		MPI_Scatter(s_matrix[i], chunk_size, MPI_UNSIGNED_SHORT, &s_recv, chunk_size, MPI_UNSIGNED_SHORT, 0, MPI_COMM_WORLD);
	
		int start_id = (rank * chunk_size);
		int end_id = (rank * chunk_size) + chunk_size;
		
        int t = 0;
		int s = 0;

        for (int j = start_id; j < end_id; j++)
		{		
			if (j == start_id && rank == 0) 
				j = j + 1;
		
			t = (0 - p_matrix[(c * (size_a + 1)) + j]) < 0;
            s = (0 - (s_matrix[i - 1][j] - (t * s_matrix[i-1][p_matrix[(c*(size_a+1))+j]-1])));
            s_recv[j - start_id] = ((t^1) || (s^0)) * (s_matrix[i-1][j]) + (!((t^1)||(s^0))) * (s_matrix[i-1][p_matrix[(c*(size_a+1))+j]-1]+1);
		}	
	
		MPI_Allgather(s_recv, chunk_size, MPI_UNSIGNED_SHORT, s_matrix[i], chunk_size, MPI_UNSIGNED_SHORT, MPI_COMM_WORLD);
	
	}

	printf("reank %d, number %d \n", rank, s_matrix[size_b][size_a]);
	return s_matrix[size_b][size_a]; 
}

void print_line(char * seqA, char * seqB, mtype * matrix,  int sizeA,
		 int sizeB)
{
    for(int i = 0; i < sizeB; i++)
        printf("%d ", matrix[i]);
	printf("\n");

}

void calc_P_matrix(mtype *P, char *b, int len_b, char *c, int len_c, int rank, int cs1, int cs2)
{
    char c_recv[cs1];
    mtype p_recv[cs1*(len_b+1)];

	
	MPI_Scatter(c, cs1, MPI_CHAR,&c_recv,cs1,MPI_CHAR, 0, MPI_COMM_WORLD); //apaguei n deu
	
	for(int i = 0; i < cs1; i++)
    {			
        for(int j = 2; j < len_b + 1 ; j++)
        {     			
            if(b[j-2] == c_recv[i]) 
            {
                p_recv[(i*(len_b+1))+j] = j-1;
            }
            else
            {
                p_recv[(i*(len_b+1))+j] = p_recv[(i*(len_b+1))+j-1];
            }        	        
        }

        MPI_Allgather(p_recv, cs1*(len_b+1), MPI_UNSIGNED_SHORT, P, cs1*(len_b+1), MPI_UNSIGNED_SHORT, MPI_COMM_WORLD);
    }
}


mtype lcs_v1(mtype **S, mtype *P, char *A, char *B, char *C, int m, int n, int u, int rank, int chunk_size)
{    
    for(int i = 1; i < m + 1; i++)
    {	   		
        int c_i = get_idx(C, u, A[i - 1]);
        
        short s_recv[chunk_size];      

        int start_id = (rank * chunk_size);
        int end_id = start_id + chunk_size;

		MPI_Scatter(S[i], chunk_size, MPI_UNSIGNED_SHORT,&s_recv,chunk_size,MPI_UNSIGNED_SHORT, 0, MPI_COMM_WORLD);

        for(int j = start_id; j < end_id; j++)
        {
           // if(j == start_id && rank == 0) j = j + 1;

            if(A[i-1] == B[j-1])
            {
                s_recv[j-start_id] = S[i-1][j-1] + 1;
            }
            else if(P[(c_i * (n + 1)) + j]==0)
            {
                s_recv[j-start_id] = max(S[i-1][j], 0);
            }
            else
            {
                s_recv[j-start_id] = max(S[i-1][j], S[i-1][P[(c_i*(n+1))+j]-1] + 1);
            }
        }
        
        MPI_Allgather(s_recv, chunk_size, MPI_UNSIGNED_SHORT,S[i], chunk_size, MPI_UNSIGNED_SHORT, MPI_COMM_WORLD);
    }

    return S[m][n-1];
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
    double time_total, seq_time,  begin, end;	
	begin = MPI_Wtime();	
	
	char *seq_a, *seq_b, *seq_c;
	seq_c = NULL;	

	// tamanho das sequencias
	int size_a, size_b, size_c;	

	// ler sequencias a e b
	seq_a = read_seq("A20000.in");
	seq_b = read_seq("B20000.in");

	size_a = strlen(seq_a);
	size_b = strlen(seq_b);

	// sequencia c
	seq_c = get_seq_c(&size_c, seq_a, seq_b, size_a, size_b);	

	chunk_size_p = (size_c / num_procs);
    chunk_size_s = ((size_b + 1) / num_procs);

	if(rank == 0)
	{
		printf("chunk_p: %d chunk_S: %d procs: %d\n", chunk_size_p, chunk_size_s, num_procs);
	}	

    // alocar e iniciar matriz S 	
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
            s_matrix[k][l] = 0;
        }
    }

	// alocar e iniciar a matriz p em forma de vetor
	p_matrix = (mtype *)malloc((size_c * (size_b + 1)) * sizeof(mtype));
	
	// LCS paralelo
	calc_P_matrix(p_matrix, seq_b, size_b, seq_c, size_c, rank, chunk_size_p, chunk_size_s);
	res_par = lcs_v1(s_matrix, p_matrix, seq_a, seq_b, seq_c, size_b, size_a, size_c, rank, chunk_size_s);
		
    end = MPI_Wtime();
	time_total = end - begin;	   
	   
    if (rank == 0)
    {
        printf("res par: %d, time par: %f \n", res_par, time_total);
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

