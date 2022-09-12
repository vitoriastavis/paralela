#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifndef max
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

typedef unsigned short mtype;

/* Read sequence from a file to a char vector.
 Filename is passed as parameter */

char* read_seq(char *fname)
{
	//file pointer
	FILE *fseq = NULL;
	//sequence size
	long size = 0;
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

void printMatrix(char * seqA, char * seqB, mtype ** scoreMatrix,  int sizeA,
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
mtype ** allocate_matrix(int sizeA, int sizeB)
{
	int i;
	//Allocate memory for LCS score matrix
	mtype ** scoreMatrix = (mtype **) malloc((sizeB + 1) * sizeof(mtype *));
	for (i = 0; i < (sizeB + 1); i++)
		scoreMatrix[i] = (mtype *) calloc((sizeA + 1), sizeof(mtype));
	return scoreMatrix;
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

// calculate s_matrix 
int lcs(mtype ** s_matrix, int size_a, int size_b, char* seq_a, char* seq_b)
{
	int i, j;

	for (i = 1; i < size_b + 1; i++)
	{
		for (j = 1; j < size_a + 1; j++)
		{
			if (seq_a[j - 1] == seq_b[i - 1])
			{			
				s_matrix[i][j] = s_matrix[i - 1][j - 1] + 1;
			}
			else
			{
				s_matrix[i][j] = max(s_matrix[i-1][j], s_matrix[i][j-1]);
			}
		}
	}
	return s_matrix[size_b][size_a];
}

int main(int argc, char ** argv)
{
	// time variables
    time_t start_time, stop_time;
	 
    double time_seq;
	start_time = time(NULL);

	// sequence pointers for both sequences
	char *seq_a, *seq_b, *seq_c;

	// initialize seq_c
	seq_c = NULL;

	// sizes of both sequences
	int size_a, size_b, size_c;

	// score result
	int res_seq; 

	// read both sequences
	seq_a = read_seq("teste.in");
	seq_b = read_seq("teste2.in");

	// find out sizes
	size_a = strlen(seq_a);
	size_b = strlen(seq_b);

	// alphabet sequence c
	seq_c = get_seq_c(&size_c, seq_a, seq_b, size_a, size_b);	

    // allocate and initiate LCS score matrix	
	mtype ** s_matrix = allocate_matrix(size_a, size_b);	
	//printMatrix(seq_a, seq_b, s_matrix, size_a, size_b);
	  
	// LCS sequential algorithm and time count
   
    res_seq = lcs(s_matrix, size_a, size_b, seq_a, seq_b);
    stop_time = time(NULL);
    time_seq = difftime(stop_time, start_time);	
	
    printf("res seq: %d, time seq: %f \n", res_seq, time_seq);   
	//printf("%f", time_seq);

}

