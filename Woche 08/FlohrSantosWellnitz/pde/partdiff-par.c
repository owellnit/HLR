/****************************************************************************/
/****************************************************************************/
/**                                                                        **/
/**                 TU München - Institut für Informatik                   **/
/**                                                                        **/
/** Copyright: Prof. Dr. Thomas Ludwig                                     **/
/**            Andreas C. Schmidt                                          **/
/**                                                                        **/
/** File:      partdiff-par.c                                              **/
/**                                                                        **/
/** Purpose:   Partial differential equation solver for Gauß-Seidel and   **/
/**            Jacobi method.                                              **/
/**                                                                        **/
/****************************************************************************/
/****************************************************************************/

/* ************************************************************************ */
/* Include standard header file.                                            */
/* ************************************************************************ */
#define _POSIX_C_SOURCE 200809L

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include <malloc.h>
#include <sys/time.h>
#include <mpi.h>

#include "partdiff-par.h"

const ROOT_RANK = 0;

const TAG_PREV_ROW = 1;
const TAG_NEXT_ROW = 2;

//Parameter für die calculation mit MPI
struct mpi_calc_arguments
{
    int rank;                   /*Aktueller rank*/
    int num_procs;              /*Anzahl der Gesamtprozesse*/
    MPI_Status status;          /*Der MPI-Status*/
    
    uint64_t absoluteStartRow;  /*Die Startreihe in der Gesamt-Matrix*/
    uint64_t matrixRows;        /*Anzahl der Teilmatrix-Reihen*/
    uint64_t matrixColumns;     /*Anzahl der Teilmatrix-Spalten*/
};

struct calculation_arguments
{
	uint64_t  N;              /* number of spaces between lines (lines=N+1)     */
	uint64_t  num_matrices;   /* number of matrices                             */
	double    h;              /* length of a space between two lines            */
	double    ***Matrix;      /* index matrix used for addressing M             */
	double    *M;             /* two matrices with real values                  */
};

struct calculation_results
{
	uint64_t  m;
	uint64_t  stat_iteration; /* number of current iteration                    */
	double    stat_precision; /* actual precision of all slaves in iteration    */
};

/* ************************************************************************ */
/* Global variables                                                         */
/* ************************************************************************ */

/* time measurement variables */
struct timeval start_time;       /* time when program started                      */
struct timeval comp_time;        /* time when calculation completed                */


/* ************************************************************************ */
/* initVariables: Initializes some global variables                         */
/* ************************************************************************ */
static
void
initVariables (struct calculation_arguments* arguments, struct calculation_results* results, struct options const* options)
{
	arguments->N = (options->interlines * 8) + 9 - 1;
	arguments->num_matrices = (options->method == METH_JACOBI) ? 2 : 1;
	arguments->h = 1.0 / arguments->N;

	results->m = 0;
	results->stat_iteration = 0;
	results->stat_precision = 0;
}

/* ************************************************************************ */
/* freeMatrices: frees memory for matrices                                  */
/* ************************************************************************ */
static
void
freeMatrices (struct calculation_arguments* arguments)
{
	uint64_t i;

	for (i = 0; i < arguments->num_matrices; i++)
	{
		free(arguments->Matrix[i]);
	}

	free(arguments->Matrix);
	free(arguments->M);
}

/* ************************************************************************ */
/* allocateMemory ()                                                        */
/* allocates memory and quits if there was a memory allocation problem      */
/* ************************************************************************ */
static
void*
allocateMemory (size_t size)
{
	void *p;

	if ((p = malloc(size)) == NULL)
	{
		printf("Speicherprobleme! (%" PRIu64 " Bytes angefordert)\n", size);
		exit(1);
	}

	return p;
}

/* ************************************************************************ */
/* allocateMatrices: allocates memory for matrices                          */
/* ************************************************************************ */
static
void
allocateMatrices (struct calculation_arguments* arguments)
{
	uint64_t i, j;

	uint64_t const N = arguments->N;

	arguments->M = allocateMemory(arguments->num_matrices * (N + 1) * (N + 1) * sizeof(double));
	arguments->Matrix = allocateMemory(arguments->num_matrices * sizeof(double**));

	for (i = 0; i < arguments->num_matrices; i++)
	{
		arguments->Matrix[i] = allocateMemory((N + 1) * sizeof(double*));

		for (j = 0; j <= N; j++)
		{
			arguments->Matrix[i][j] = arguments->M + (i * (N + 1) * (N + 1)) + (j * (N + 1));
		}
	}
}

/* ************************************************************************ */
/* allocateMatrices: allocates memory for matrices                          */
/* ************************************************************************ */
static
void
allocateMatrices_mpi (struct calculation_arguments* arguments, struct mpi_calc_arguments* mpiArgs)
{
    uint64_t i, j;
    
    //Reihen für Prozess ermitteln
    uint64_t rowsPerProcess = (arguments->N + 1) / mpiArgs->num_procs;
    
    //Anzahl restlicher Reihen
    int64_t remainRows = (arguments->N + 1) % mpiArgs->num_procs;
    
    //Zwei addieren, wegen letzte Reihe Vorgänger und erste Reihe Nachfolger
    mpiArgs->matrixRows  = rowsPerProcess + 2;
    
    mpiArgs->matrixColumns = arguments->N + 1;
    
    //Restliche Reihen aufteilen, wenn Prozessnummer kleiner Anzahl restlichen Reihen
    if (mpiArgs->rank < remainRows)
    {
        mpiArgs->matrixRows  += 1;
    }
    
    //Speicher für Matrixen reservieren
    arguments->M = allocateMemory(arguments->num_matrices * (mpiArgs->matrixRows  * mpiArgs->matrixColumns) * sizeof(double));
    arguments->Matrix = allocateMemory(arguments->num_matrices * sizeof(double**));
    
    
    for (i = 0; i < arguments->num_matrices; i++)
    {
        arguments->Matrix[i] = allocateMemory((mpiArgs->matrixRows) * sizeof(double*));
        
        for (j = 0; j < mpiArgs->matrixRows; j++)
        {
            arguments->Matrix[i][j] = arguments->M + (i * mpiArgs->matrixRows * mpiArgs->matrixColumns) + (j * mpiArgs->matrixColumns);
        }
    }
    
    //Startreihe ermitteln
    int actual_splitted_remain_rows = 0;
    //Wenn Rank keine Extra-Reihe bekommen hat, dann müssen alle Extra-Reihen hinzuaddiert werden
    if (mpiArgs->rank > remainRows)
    {
        actual_splitted_remain_rows = remainRows;
    }
    //Sonst müssen die bisherigen Extra-Reihen addiert werden, die bisher aufgeteilt wurden
    else
    {
        actual_splitted_remain_rows = mpiArgs->rank;
    }
    
    //Startreihe des Prozesses in der Gesamt-Matrix berechnen
    mpiArgs->absoluteStartRow = (rowsPerProcess * mpiArgs->rank) + actual_splitted_remain_rows;
    
    printf("Rank %d rows %lu N %lu from: %lu to: %lu\n", mpiArgs->rank, mpiArgs->matrixRows, arguments->N, mpiArgs->absoluteStartRow, mpiArgs->absoluteStartRow + mpiArgs->matrixRows - 3);
}

/* ************************************************************************ */
/* initMatrices: Initialize matrix/matrices and some global variables       */
/* ************************************************************************ */
static
void
initMatrices (struct calculation_arguments* arguments, struct options const* options)
{
	uint64_t g, i, j;                                /*  local variables for loops   */

	uint64_t const N = arguments->N;
	double const h = arguments->h;
	double*** Matrix = arguments->Matrix;

	/* initialize matrix/matrices with zeros */
	for (g = 0; g < arguments->num_matrices; g++)
	{
		for (i = 0; i <= N; i++)
		{
			for (j = 0; j <= N; j++)
			{
				Matrix[g][i][j] = 0.0;
			}
		}
	}

	/* initialize borders, depending on function (function 2: nothing to do) */
	if (options->inf_func == FUNC_F0)
	{
		for (g = 0; g < arguments->num_matrices; g++)
		{
			for (i = 0; i <= N; i++)
			{
				Matrix[g][i][0] = 1.0 - (h * i);
				Matrix[g][i][N] = h * i;
				Matrix[g][0][i] = 1.0 - (h * i);
				Matrix[g][N][i] = h * i;
			}

			Matrix[g][N][0] = 0.0;
			Matrix[g][0][N] = 0.0;
		}
	}
}

/* ************************************************************************ */
/* initMatrices: Initialize matrix/matrices and some global variables       */
/* ************************************************************************ */
static
void
initMatrices_mpi (struct calculation_arguments* arguments, struct options const* options, struct mpi_calc_arguments* mpiArgs)
{
    uint64_t g, i, j;                                /*  local variables for loops   */
    
    const uint64_t rows = mpiArgs->matrixRows;
    const uint64_t columns = mpiArgs->matrixColumns;
    
    double const h = arguments->h;
    double*** Matrix = arguments->Matrix;
    
    /* initialize matrix/matrices with zeros */
    for (g = 0; g < arguments->num_matrices; g++)
    {
        for (i = 1; i < (rows - 1); i++)
        {
            for (j = 0; j < columns; j++)
            {
                Matrix[g][i][j] = 0.0;
            }
        }
    }
    
    /* initialize borders, depending on function (function 2: nothing to do) */
    if (options->inf_func == FUNC_F0)
    {
        const uint64_t absoluteStartRow = mpiArgs->absoluteStartRow; // die erste Zeile wird überprungen
        
        for (g = 0; g < arguments->num_matrices; g++)
        {
            // linke und rechte Ränder
            for (uint64_t i = 1; i < (rows - 1); i++)
            {
                Matrix[g][i][0] = 1.0 - (h * (absoluteStartRow + i -1));
                Matrix[g][i][columns-1] = h * (absoluteStartRow + i -1);
            }
            
            if (mpiArgs->rank == ROOT_RANK)
            {
                // obere Kante initialisieren
                for (uint64_t i = 0; i < columns; i++)
                {
                    Matrix[g][1][i] = 1.0 - (h * i);
                }
                Matrix[g][1][columns - 1] = 0;
            }
            if (mpiArgs->rank == (mpiArgs->num_procs - 1))
            {
                // untere Kante initialisieren
                for (uint64_t i = 0; i < columns; i++)
                {
                    Matrix[g][rows - 2][i] = h * i;
                }
                Matrix[g][rows - 2][0] = 0;
            }
        }
    }
}

/* ************************************************************************ */
/* calculate: solves the equation                                           */
/* ************************************************************************ */
static
void
calculate (struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options)
{
	int i, j;                                   /* local variables for loops */
	int m1, m2;                                 /* used as indices for old and new matrices */
	double star;                                /* four times center value minus 4 neigh.b values */
	double residuum;                            /* residuum of current iteration */
	double maxresiduum;                         /* maximum residuum value of a slave in iteration */

	int const N = arguments->N;
	double const h = arguments->h;

	double pih = 0.0;
	double fpisin = 0.0;

	int term_iteration = options->term_iteration;

	/* initialize m1 and m2 depending on algorithm */
	if (options->method == METH_JACOBI)
	{
		m1 = 0;
		m2 = 1;
	}
	else
	{
		m1 = 0;
		m2 = 0;
	}

	if (options->inf_func == FUNC_FPISIN)
	{
		pih = PI * h;
		fpisin = 0.25 * TWO_PI_SQUARE * h * h;
	}

	while (term_iteration > 0)
	{
		double** Matrix_Out = arguments->Matrix[m1];
		double** Matrix_In  = arguments->Matrix[m2];

		maxresiduum = 0;

		/* over all rows */
		for (i = 1; i < N; i++)
		{
			double fpisin_i = 0.0;

			if (options->inf_func == FUNC_FPISIN)
			{
				fpisin_i = fpisin * sin(pih * (double)i);
			}

			/* over all columns */
			for (j = 1; j < N; j++)
			{
				star = 0.25 * (Matrix_In[i-1][j] + Matrix_In[i][j-1] + Matrix_In[i][j+1] + Matrix_In[i+1][j]);

				if (options->inf_func == FUNC_FPISIN)
				{
					star += fpisin_i * sin(pih * (double)j);
				}

				if (options->termination == TERM_PREC || term_iteration == 1)
				{
					residuum = Matrix_In[i][j] - star;
					residuum = (residuum < 0) ? -residuum : residuum;
					maxresiduum = (residuum < maxresiduum) ? maxresiduum : residuum;
				}

				Matrix_Out[i][j] = star;
			}
		}

		results->stat_iteration++;
		results->stat_precision = maxresiduum;

		/* exchange m1 and m2 */
		i = m1;
		m1 = m2;
		m2 = i;

		/* check for stopping calculation depending on termination method */
		if (options->termination == TERM_PREC)
		{
			if (maxresiduum < options->term_precision)
			{
				term_iteration = 0;
			}
		}
		else if (options->termination == TERM_ITER)
		{
			term_iteration--;
		}
	}

	results->m = m2;
}

/* ************************************************************************ */
/* calculate: solves the equation                                           */
/* ************************************************************************ */
static
void
calculate_mpi (struct calculation_arguments const* arguments, struct calculation_results *results, struct options const* options, struct mpi_calc_arguments* mpiArgs)
{
    int i, j;                                   /* local variables for loops  */
    int m1, m2;                                 /* used as indices for old and new matrices       */
    double star;                                /* four times center value minus 4 neigh.b values */
    double residuum;                            /* residuum of current iteration                  */
    double maxresiduum;                         /* maximum residuum value of a slave in iteration */
    double const h = arguments->h;
    uint64_t absoluteRow = mpiArgs->absoluteStartRow;
    if (mpiArgs->rank == ROOT_RANK)
    {
        absoluteRow = 1;
    }
    
    const int columns = mpiArgs->matrixColumns;
    const int startRow = (mpiArgs->rank == ROOT_RANK)? 2 : 1;
    const int stopRow = (mpiArgs->rank == (mpiArgs->num_procs -1))? mpiArgs->matrixRows -2 : mpiArgs->matrixRows -1;
    
    const int nextTarget = (mpiArgs->rank + 1) % mpiArgs->num_procs;
    const int prevTarget = ((mpiArgs->rank - 1) >= 0)? mpiArgs->rank -1 : mpiArgs->num_procs -1;
    
    printf("Rank %d start %d stop %d next %d prev %d\n", mpiArgs->rank, startRow, stopRow, nextTarget, prevTarget);
    
    double pih = 0.0;
    double fpisin = 0.0;
    
    int term_iteration = options->term_iteration;
    
    /* initialize m1 and m2 depending on algorithm */
    if (options->method == METH_JACOBI)
    {
        m1 = 0;
        m2 = 1;
    }
    else
    {
        m1 = 0;
        m2 = 0;
    }
    
    if (options->inf_func == FUNC_FPISIN)
    {
        pih = PI * h;
        fpisin = 0.25 * TWO_PI_SQUARE * h * h;
    }
    
    while (term_iteration > 0)
    {
        double** Matrix_Out = arguments->Matrix[m1];
        double** Matrix_In  = arguments->Matrix[m2];
        
        maxresiduum = 0;
        
        // Zeilen austauschen
        if (mpiArgs->rank == 0) // mit ROOT
        {
            //Wenn nicht Root-Rank, dann erste Reihe an Vorgänger senden
            if (mpiArgs->rank != ROOT_RANK)
            {
                MPI_Send(Matrix_In[1], mpiArgs->matrixColumns, MPI_DOUBLE, prevTarget, TAG_PREV_ROW, MPI_COMM_WORLD);
            }
            
            //Wenn nicht letzter Rank, dann Letzte Reihe an Nachfolger senden
            if (mpiArgs->rank != (mpiArgs->num_procs - 1))
            {
                MPI_Send(Matrix_In[stopRow - 1], mpiArgs->matrixColumns, MPI_DOUBLE, nextTarget, TAG_NEXT_ROW, MPI_COMM_WORLD);
            }
            
            //Wenn nicht Root-Rank, dann letzte Reihe vom Vorgänger empfangen
            if (mpiArgs->rank != ROOT_RANK)
            {
                MPI_Recv(Matrix_In[0], mpiArgs->matrixColumns, MPI_DOUBLE,prevTarget, TAG_NEXT_ROW, MPI_COMM_WORLD, &mpiArgs->status);
            }
            
            //Wenn nicht letzter Rank, dann erste Reihe vom Nachfolger empfangen
            if (mpiArgs->rank != (mpiArgs->num_procs - 1))
            {
                MPI_Recv(Matrix_In[stopRow], mpiArgs->matrixColumns, MPI_DOUBLE, nextTarget, TAG_PREV_ROW, MPI_COMM_WORLD, &mpiArgs->status);
            }
        }
        else // ohne ROOT
        {
            //Letzte Reihe vom Vorgänger Empfangen
            MPI_Recv(Matrix_In[0], mpiArgs->matrixColumns, MPI_DOUBLE,prevTarget, TAG_NEXT_ROW, MPI_COMM_WORLD, &mpiArgs->status);
            
            //Wenn nicht letzter Rank, erst Reihe vom Nachfolger empfangen
            if (mpiArgs->rank != (mpiArgs->num_procs - 1))
            {
                MPI_Recv(Matrix_In[stopRow], mpiArgs->matrixColumns, MPI_DOUBLE,nextTarget, TAG_PREV_ROW, MPI_COMM_WORLD, &mpiArgs->status);
            }
            
            //Erste Reihe an Vorgänger senden
            MPI_Send(Matrix_In[1], mpiArgs->matrixColumns, MPI_DOUBLE, prevTarget, TAG_PREV_ROW, MPI_COMM_WORLD);
            
            //Wenn nicht letzter Rank, letzte Reihe an Nachfolger senden
            if (mpiArgs->rank != (mpiArgs->num_procs - 1))
            {
                MPI_Send(Matrix_In[stopRow - 1], mpiArgs->matrixColumns, MPI_DOUBLE, nextTarget, TAG_NEXT_ROW, MPI_COMM_WORLD);
            }
        }
        
        /* over all rows */
        for (i = startRow; i < stopRow; i++)
        {
            double fpisin_i = 0.0;
            
            if (options->inf_func == FUNC_FPISIN)
            {
                fpisin_i = fpisin * sin(pih * (double)(absoluteRow + i - startRow));
            }
            
            /* over all columns */
            for (j = 1; j < (columns - 1); j++)
            {
                star = 0.25 * (Matrix_In[i-1][j] + Matrix_In[i][j-1] + Matrix_In[i][j+1] + Matrix_In[i+1][j]);
                
                if (options->inf_func == FUNC_FPISIN)
                {
                    star += fpisin_i * sin(pih * (double)j);
                }
                
                if (options->termination == TERM_PREC || term_iteration == 1)
                {
                    residuum = Matrix_In[i][j] - star;
                    residuum = (residuum < 0) ? -residuum : residuum;
                    maxresiduum = (residuum < maxresiduum) ? maxresiduum : residuum;
                }
                
                Matrix_Out[i][j] = star;
            }
        }
        
        results->stat_iteration++;
        
        double all_maxresiduum;
        MPI_Allreduce(&maxresiduum, &all_maxresiduum, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        
        maxresiduum = all_maxresiduum;
        results->stat_precision = maxresiduum;
        
        /* exchange m1 and m2 */
        i = m1;
        m1 = m2;
        m2 = i;
        
        /* check for stopping calculation, depending on termination method */
        if (options->termination == TERM_PREC)
        {
            if (maxresiduum < options->term_precision)
            {
                term_iteration = 0;
            }
        }
        else if (options->termination == TERM_ITER)
        {
            term_iteration--;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    results->m = m2;
}


/* ************************************************************************ */
/*  displayStatistics: displays some statistics about the calculation       */
/* ************************************************************************ */
static
void
displayStatistics (struct calculation_arguments const* arguments, struct calculation_results const* results, struct options const* options)
{
	int N = arguments->N;
	double time = (comp_time.tv_sec - start_time.tv_sec) + (comp_time.tv_usec - start_time.tv_usec) * 1e-6;

	printf("Berechnungszeit:    %f s \n", time);
	printf("Speicherbedarf:     %f MiB\n", (N + 1) * (N + 1) * sizeof(double) * arguments->num_matrices / 1024.0 / 1024.0);
	printf("Berechnungsmethode: ");

	if (options->method == METH_GAUSS_SEIDEL)
	{
		printf("Gauß-Seidel");
	}
	else if (options->method == METH_JACOBI)
	{
		printf("Jacobi");
	}

	printf("\n");
	printf("Interlines:         %" PRIu64 "\n",options->interlines);
	printf("Stoerfunktion:      ");

	if (options->inf_func == FUNC_F0)
	{
		printf("f(x,y) = 0");
	}
	else if (options->inf_func == FUNC_FPISIN)
	{
		printf("f(x,y) = 2pi^2*sin(pi*x)sin(pi*y)");
	}

	printf("\n");
	printf("Terminierung:       ");

	if (options->termination == TERM_PREC)
	{
		printf("Hinreichende Genaugkeit");
	}
	else if (options->termination == TERM_ITER)
	{
		printf("Anzahl der Iterationen");
	}

	printf("\n");
	printf("Anzahl Iterationen: %" PRIu64 "\n", results->stat_iteration);
	printf("Norm des Fehlers:   %e\n", results->stat_precision);
	printf("\n");
}

/****************************************************************************/
/** Beschreibung der Funktion DisplayMatrix:                               **/
/**                                                                        **/
/** Die Funktion DisplayMatrix gibt eine Matrix                            **/
/** in einer "ubersichtlichen Art und Weise auf die Standardausgabe aus.   **/
/**                                                                        **/
/** Die "Ubersichtlichkeit wird erreicht, indem nur ein Teil der Matrix    **/
/** ausgegeben wird. Aus der Matrix werden die Randzeilen/-spalten sowie   **/
/** sieben Zwischenzeilen ausgegeben.                                      **/
/****************************************************************************/
static
void
DisplayMatrix (struct calculation_arguments* arguments, struct calculation_results* results, struct options* options)
{
	int x, y;

	double** Matrix = arguments->Matrix[results->m];

	int const interlines = options->interlines;

	printf("Matrix:\n");

	for (y = 0; y < 9; y++)
	{
		for (x = 0; x < 9; x++)
		{
			printf ("%7.4f", Matrix[y * (interlines + 1)][x * (interlines + 1)]);
		}

		printf ("\n");
	}

	fflush (stdout);
}

/**
 * rank and size are the MPI rank and size, respectively.
 * from and to denote the global(!) range of lines that this process is responsible for.
 *
 * Example with 9 matrix lines and 4 processes:
 * - rank 0 is responsible for 1-2, rank 1 for 3-4, rank 2 for 5-6 and rank 3 for 7.
 *   Lines 0 and 8 are not included because they are not calculated.
 * - Each process stores two halo lines in its matrix (except for ranks 0 and 3 that only store one).
 * - For instance: Rank 2 has four lines 0-3 but only calculates 1-2 because 0 and 3 are halo lines for other processes. It is responsible for (global) lines 5-6.
 */
static
void
DisplayMatrix_mpi (struct calculation_arguments* arguments, struct calculation_results* results, struct options* options, int rank, int size, int from, int to)
{
    int const elements = 8 * options->interlines + 9;
    
    int x, y;
    double** Matrix = arguments->Matrix[results->m];
    MPI_Status status;
    
    /* first line belongs to rank 0 */
    if (rank == 0)
        from--;
    
    /* last line belongs to rank size - 1 */
    if (rank + 1 == size)
        to++;
    
    if (rank == 0)
        printf("Matrix:\n");
    
    for (y = 0; y < 9; y++)
    {
        int line = y * (options->interlines + 1);
        
        if (rank == 0)
        {
            /* check whether this line belongs to rank 0 */
            if (line < from || line > to)
            {
                /* use the tag to receive the lines in the correct order
                 * the line is stored in Matrix[0], because we do not need it anymore */
                MPI_Recv(Matrix[0], elements, MPI_DOUBLE, MPI_ANY_SOURCE, 42 + y, MPI_COMM_WORLD, &status);
            }
        }
        else
        {
            if (line >= from && line <= to)
            {
                /* if the line belongs to this process, send it to rank 0
                 * (line - from + 1) is used to calculate the correct local address */
                MPI_Send(Matrix[line - from + 1], elements, MPI_DOUBLE, 0, 42 + y, MPI_COMM_WORLD);
            }
        }
        
        if (rank == 0)
        {
            for (x = 0; x < 9; x++)
            {
                int col = x * (options->interlines + 1);
                
                if (line >= from && line <= to)
                {
                    /* this line belongs to rank 0 */
                    printf("%7.4f", Matrix[line][col]);
                }
                else
                {
                    /* this line belongs to another rank and was received above */
                    printf("%7.4f", Matrix[0][col]);
                }
            }
            
            printf("\n");
        }
    }
    
    fflush(stdout);
}


/* ************************************************************************ */
/*  main                                                                    */
/* ************************************************************************ */
int
main (int argc, char** argv)
{
	struct options options;
	struct calculation_arguments arguments;
	struct calculation_results results;
    
    struct mpi_calc_arguments mpiArgs;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiArgs.rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiArgs.num_procs);

	AskParams(&options, argc, argv);

	initVariables(&arguments, &results, &options);

    
    if (options.method == METH_JACOBI && mpiArgs.num_procs > 1)
    {
        allocateMatrices_mpi(&arguments, &mpiArgs);
        initMatrices_mpi(&arguments, &options, &mpiArgs);
        
        gettimeofday(&start_time, NULL);
        calculate_mpi(&arguments, &results, &options, &mpiArgs);
        gettimeofday(&comp_time, NULL);
        
        if (mpiArgs.rank == ROOT_RANK)
        {
            displayStatistics(&arguments, &results, &options);
        }
        
        DisplayMatrix_mpi(&arguments, &results, &options, mpiArgs.rank, mpiArgs.num_procs, mpiArgs.absoluteStartRow, mpiArgs.absoluteStartRow + mpiArgs.matrixRows -3);
    }
    else
    {
        allocateMatrices(&arguments);
        initMatrices(&arguments, &options);

        gettimeofday(&start_time, NULL);
        calculate(&arguments, &results, &options);
        gettimeofday(&comp_time, NULL);

        displayStatistics(&arguments, &results, &options);
        DisplayMatrix(&arguments, &results, &options);
        
    }
    
    freeMatrices(&arguments);
    MPI_Finalize();

	return 0;
}
