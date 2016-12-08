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
/** Purpose:   Partial differential equation solver for Gauss-Seidel and   **/
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

#define ROOT 0

#define TAG_SEND_UPPER_ROW 1
#define TAG_SEND_LOWER_ROW 2
#define TAG_SEND_RESIDUUM 3
#define TAG_SEND_TERMINATION 4

#define TERMINATION -5.0

struct comm_options
{
	int rank;
	int num_procs;
	MPI_Status status;

	uint64_t absoluteStartRow;
	uint64_t matrixRows;
	uint64_t matrixColumns;
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
		printf("Speicherprobleme! (%" PRIu64 " Bytes)\n", size);
		/* exit program */
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
allocateMatrices_mpi (struct calculation_arguments* arguments, struct comm_options* comm)
{
	uint64_t i, j;

	uint64_t chunkSize = (arguments->N + 1) / comm->num_procs;
	int64_t extra = (arguments->N + 1) % comm->num_procs;

	uint64_t rows = chunkSize + 2;
	uint64_t columns = arguments->N + 1;
	
	// Zwei extra Zeilen pro Prozess für die Daten die Empfangen werden müssen
	if (comm->rank < extra)
	{
		rows += 1;
	}

	arguments->M = allocateMemory(arguments->num_matrices * (rows * columns) * sizeof(double));
	arguments->Matrix = allocateMemory(arguments->num_matrices * sizeof(double**));

	for (i = 0; i < arguments->num_matrices; i++)
	{
		arguments->Matrix[i] = allocateMemory((rows) * sizeof(double*));

		for (j = 0; j < rows; j++)
		{
			arguments->Matrix[i][j] = arguments->M + (i * rows * columns) + (j * columns);
		}
	}

	comm->matrixRows = rows;
	comm->matrixColumns = columns;

	int prev_extras = 0;
	if (comm->rank > extra)
	{
		prev_extras = extra;
	}
	else
	{
		prev_extras = comm->rank;
	}
	comm->absoluteStartRow = (chunkSize * comm->rank) + prev_extras;

	// printf("Rank %d rows %lu N %lu from: %lu to: %lu\n", comm->rank, rows, arguments->N, comm->absoluteStartRow, comm->absoluteStartRow + comm->matrixRows -3);
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
		for (i = 1; i <= N; i++)
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
initMatrices_mpi (struct calculation_arguments* arguments, struct options const* options, struct comm_options* comm)
{
	uint64_t g, i, j;                                /*  local variables for loops   */

	const uint64_t rows = comm->matrixRows;
	const uint64_t columns = comm->matrixColumns;

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
		const uint64_t absoluteStartRow = comm->absoluteStartRow; // die erste Zeile wird überprungen

		for (g = 0; g < arguments->num_matrices; g++)
		{
			// linke und rechte Ränder
			for (uint64_t i = 0; i < rows; i++)
			{
				Matrix[g][i][0] = 1.0 - (h * (absoluteStartRow + i -1));
				Matrix[g][i][columns-1] = h * (absoluteStartRow + i -1);
			}

			if (comm->rank == ROOT)
			{
				// obere Kante initialisieren
				for (uint64_t i = 0; i < columns; i++)
				{
					Matrix[g][1][i] = 1.0 - (h * i);
				}
				Matrix[g][1][columns - 1] = 0;
			}
			if (comm->rank == (comm->num_procs - 1))
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
calculate (struct calculation_arguments const* arguments, struct calculation_results *results, struct options const* options)
{
	int i, j;                                   /* local variables for loops  */
	int m1, m2;                                 /* used as indices for old and new matrices       */
	double star;                                /* four times center value minus 4 neigh.b values */
	double residuum;                            /* residuum of current iteration                  */
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
	}

	results->m = m2;
}

/* ************************************************************************ */
/* calculate: solves the equation                                           */
/* ************************************************************************ */
static
void
calculate_mpi (struct calculation_arguments const* arguments, struct calculation_results *results, struct options const* options, struct comm_options* comm)
{
	int i, j;                                   /* local variables for loops  */
	int m1, m2;                                 /* used as indices for old and new matrices       */
	double star;                                /* four times center value minus 4 neigh.b values */
	double residuum;                            /* residuum of current iteration                  */
	double maxresiduum;                         /* maximum residuum value of a slave in iteration */
	double const h = arguments->h;
	uint64_t absoluteRow = comm->absoluteStartRow;
	if (comm->rank == ROOT)
	{
		absoluteRow = 1;
	}

	const int columns = comm->matrixColumns;
	const int startRow = (comm->rank == ROOT)? 2 : 1;
	const int stopRow = (comm->rank == (comm->num_procs -1))? comm->matrixRows -2 : comm->matrixRows -1;

	const int nextTarget = (comm->rank + 1) % comm->num_procs;
	const int prevTarget = ((comm->rank - 1) >= 0)? comm->rank -1 : comm->num_procs -1;

	//printf("Rank %d start %d stop %d next %d prev %d\n", comm->rank, startRow, stopRow, nextTarget, prevTarget);

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
		if ((comm->rank % 2) == 0) // mit ROOT
		{
			if (comm->rank != ROOT)
			{
				MPI_Send(Matrix_In[1], comm->matrixColumns, MPI_DOUBLE, prevTarget, TAG_SEND_UPPER_ROW, MPI_COMM_WORLD);	
			}
			if (comm->rank != (comm->num_procs - 1))
			{
				MPI_Send(Matrix_In[stopRow - 1], comm->matrixColumns, MPI_DOUBLE, nextTarget, TAG_SEND_LOWER_ROW, MPI_COMM_WORLD);
			}

			if (comm->rank != ROOT)
			{
				MPI_Recv(Matrix_In[0], comm->matrixColumns, MPI_DOUBLE,prevTarget, TAG_SEND_LOWER_ROW, MPI_COMM_WORLD, &comm->status);
			}
			if (comm->rank != (comm->num_procs - 1))
			{
				MPI_Recv(Matrix_In[stopRow], comm->matrixColumns, MPI_DOUBLE,nextTarget, TAG_SEND_UPPER_ROW, MPI_COMM_WORLD, &comm->status);
			}
		}
		else // ohne ROOT
		{
			MPI_Recv(Matrix_In[0], comm->matrixColumns, MPI_DOUBLE,prevTarget, TAG_SEND_LOWER_ROW, MPI_COMM_WORLD, &comm->status);
			if (comm->rank != (comm->num_procs - 1))
			{
				MPI_Recv(Matrix_In[stopRow], comm->matrixColumns, MPI_DOUBLE,nextTarget, TAG_SEND_UPPER_ROW, MPI_COMM_WORLD, &comm->status);
			}

			MPI_Send(Matrix_In[1], comm->matrixColumns, MPI_DOUBLE, prevTarget, TAG_SEND_UPPER_ROW, MPI_COMM_WORLD);	
			if (comm->rank != (comm->num_procs - 1))
			{
				MPI_Send(Matrix_In[stopRow -1], comm->matrixColumns, MPI_DOUBLE, nextTarget, TAG_SEND_LOWER_ROW, MPI_COMM_WORLD);
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
/* calculate: solves the equation                                           */
/* ************************************************************************ */
static
void
calculate_mpi2_gauss (struct calculation_arguments const* arguments, struct calculation_results *results, struct options const* options, struct comm_options* comm)
{
    int i, j;                                   /* local variables for loops  */
    double star;                                /* four times center value minus 4 neigh.b values */
    double residuum;                            /* residuum of current iteration                  */
    double maxresiduum, prevMaxResiduum;       /* maximum residuum value of a slave in iteration */
    double const h = arguments->h;
    uint64_t absoluteRow = comm->absoluteStartRow;
    if (comm->rank == ROOT)
    {
        absoluteRow = 1;
    }
    
    int stopSend = 0;
    int stopReceived = 0;

    const int columns = comm->matrixColumns;
    const int startRow = (comm->rank == ROOT)? 2 : 1;
    const int stopRow = (comm->rank == (comm->num_procs -1))? comm->matrixRows -2 : comm->matrixRows -1;
    
    const int nextTarget = (comm->rank + 1) % comm->num_procs;
    const int prevTarget = ((comm->rank - 1) >= 0)? comm->rank -1 : comm->num_procs -1;
    
    // printf("Rank %d start %d stop %d next %d prev %d\n", comm->rank, startRow, stopRow, nextTarget, prevTarget);
    
    double pih = 0.0;
    double fpisin = 0.0;
    int term_iteration = options->term_iteration;
    
	double** Matrix = arguments->Matrix[0];
    
    if (options->inf_func == FUNC_FPISIN)
    {
        pih = PI * h;
        fpisin = 0.25 * TWO_PI_SQUARE * h * h;
    }

    MPI_Request stopSignal;
    if (options->termination == TERM_PREC && comm->rank == ROOT)
    {
    	MPI_Irecv(&stopReceived, 1, MPI_INT, (comm->num_procs - 1), TAG_SEND_TERMINATION, MPI_COMM_WORLD, &stopSignal);
    }
    
	MPI_Request sendLowerToNext, sendUpperToPrev, receiveUpperFromNext, receiveLowerFromPrev;
	MPI_Request sendMaxResiduum, receiveMaxResiduum;

	sendLowerToNext = sendUpperToPrev = receiveUpperFromNext = MPI_REQUEST_NULL;
	sendMaxResiduum = receiveMaxResiduum = MPI_REQUEST_NULL;
    
    while (term_iteration > 0)
    {
		// Zeilen empfangen
        if (comm->rank != ROOT)
        {
        	// Untere Zeile des vorherigen prozesses empfangen
        	MPI_Irecv(Matrix[0], comm->matrixColumns, MPI_DOUBLE, prevTarget, TAG_SEND_LOWER_ROW, MPI_COMM_WORLD, &receiveLowerFromPrev);

        	// maxresiduum des vorherigen Prozesses empfangen
			if (options->termination == TERM_PREC)
			{
				MPI_Irecv(&prevMaxResiduum, 1, MPI_DOUBLE, prevTarget, TAG_SEND_RESIDUUM, MPI_COMM_WORLD, &receiveMaxResiduum);
			}

    		MPI_Wait(&sendMaxResiduum, &comm->status);
			MPI_Wait(&receiveLowerFromPrev, &comm->status);

			int itemCount;
			MPI_Get_count(&comm->status, MPI_DOUBLE, &itemCount);
			if (itemCount == 1 && Matrix[0][0] == TERMINATION)
			{
				stopReceived = 1;
			}
        }

        if (stopReceived)
		{
			// Kommunikation aufräumen
			MPI_Wait(&receiveUpperFromNext, &comm->status);
			// Um sicher zu stellen das es nicht zu einer race kondition kommt senden wir über TAG_SEND_LOWER_ROW nur einen Wert um zu beenden statt der Zeile
			double end = TERMINATION;
			MPI_Send(&end, 1, MPI_DOUBLE, nextTarget, TAG_SEND_LOWER_ROW, MPI_COMM_WORLD);

			results->stat_precision = maxresiduum;
			break;
		}
       
        maxresiduum = 0;      
        /* over all rows */
        for (i = startRow; i < stopRow; i++)
        {
            double fpisin_i = 0.0;
            
            if (options->inf_func == FUNC_FPISIN)
            {
                fpisin_i = fpisin * sin(pih * (double)(absoluteRow + i - startRow));
            }
            
            // auf die letzte Zeile warten
            if (i == (stopRow - 1) && comm->rank != (comm->num_procs - 1))
            {
            	MPI_Wait(&receiveUpperFromNext, &comm->status);
            }

            /* over all columns */
            for (j = 1; j < (columns - 1); j++)
            {
                star = 0.25 * (Matrix[i-1][j] + Matrix[i][j-1] + Matrix[i][j+1] + Matrix[i+1][j]);
                
                if (options->inf_func == FUNC_FPISIN)
                {
                    star += fpisin_i * sin(pih * (double)j);
                }
                
                if (options->termination == TERM_PREC || term_iteration == 1)
                {
                    residuum = Matrix[i][j] - star;
                    residuum = (residuum < 0) ? -residuum : residuum;
                    maxresiduum = (residuum < maxresiduum) ? maxresiduum : residuum;
                }
                
                Matrix[i][j] = star;
            }

            // die obere Zeile an den vorherigen Prozess senden
            if (i == startRow && comm->rank != ROOT)
            {
            	MPI_Isend(Matrix[1], comm->matrixColumns, MPI_DOUBLE, prevTarget, TAG_SEND_UPPER_ROW, MPI_COMM_WORLD, &sendUpperToPrev);
            }
        }
        
        results->stat_iteration++;

        if (comm->rank != (comm->num_procs - 1))
        {
        	// Untere Zeile an den nächsten Prozess senden
        	MPI_Isend(Matrix[stopRow - 1], comm->matrixColumns, MPI_DOUBLE, nextTarget, TAG_SEND_LOWER_ROW, MPI_COMM_WORLD, &sendLowerToNext);

        	// Obere Zeile des nächsten Prozesses empfangen
			MPI_Irecv(Matrix[stopRow], comm->matrixColumns, MPI_DOUBLE, nextTarget, TAG_SEND_UPPER_ROW, MPI_COMM_WORLD, &receiveUpperFromNext);
        }

        /* check for stopping calculation, depending on termination method */
        if (options->termination == TERM_PREC)
        {
        	if (comm->rank != ROOT)
        	{
        		MPI_Wait(&receiveMaxResiduum, &comm->status);
        		maxresiduum = (prevMaxResiduum < maxresiduum) ? maxresiduum : prevMaxResiduum;
        	}

        	if (comm->rank < (comm->num_procs - 1))
        	{
        		MPI_Isend(&maxresiduum, 1, MPI_DOUBLE, nextTarget, TAG_SEND_RESIDUUM, MPI_COMM_WORLD, &sendMaxResiduum);
        	}
        	else if (maxresiduum < options->term_precision && !stopSend) // der letzte Prozess mit ausreichender genauigkeit
        	{
        		stopSend = 1;
        		MPI_Send(&stopSend, 1, MPI_INT, ROOT, TAG_SEND_TERMINATION, MPI_COMM_WORLD);
        	}
        }
        else if (options->termination == TERM_ITER)
        {
            term_iteration--;
        }

        MPI_Wait(&sendLowerToNext, &comm->status);
        MPI_Wait(&sendUpperToPrev, &comm->status);
    }


    MPI_Barrier(MPI_COMM_WORLD);
    if (options->termination == TERM_ITER)
    {
    	double all_maxresiduum;
    	MPI_Allreduce(&maxresiduum, &all_maxresiduum, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    	results->stat_precision = all_maxresiduum;
    }
    
    MPI_Wait(&receiveUpperFromNext, &comm->status);
    results->m = 0;
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
		printf("Gauss-Seidel");
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
        MPI_Recv(Matrix[1], elements, MPI_DOUBLE, MPI_ANY_SOURCE, 42 + y, MPI_COMM_WORLD, &status);
      }
    }
    else
    {
      if (line >= from && line <= to)
      {
        /* if the line belongs to this process, send it to rank 0
         * (line - from + 1) is used to calculate the correct local address */
        MPI_Send(Matrix[(line - from) + 1], elements, MPI_DOUBLE, 0, 42 + y, MPI_COMM_WORLD);
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
          printf("%7.4f", Matrix[1 + line][col]);
        }
        else
        {
          /* this line belongs to another rank and was received above */
          printf("%7.4f", Matrix[1][col]);
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

	struct comm_options comm;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &comm.rank);
	MPI_Comm_size(MPI_COMM_WORLD, &comm.num_procs);

	/* get parameters */
	AskParams(&options, argc, argv, comm.rank);              /* ************************* */
	initVariables(&arguments, &results, &options);           /* ******************************************* */

	if (comm.num_procs > 1)
	{
		allocateMatrices_mpi(&arguments, &comm);
		initMatrices_mpi(&arguments, &options, &comm); 

		gettimeofday(&start_time, NULL);                   /*  start timer         */
		if (options.method == METH_JACOBI)
		{
			calculate_mpi(&arguments, &results, &options, &comm);                                      /*  solve the equation  */
		}
		else
		{
			calculate_mpi2_gauss(&arguments, &results, &options, &comm);   
		}
		gettimeofday(&comp_time, NULL);                   /*  stop timer          */

		if (comm.rank == (comm.num_procs - 1))
		{
			displayStatistics(&arguments, &results, &options);
		}

		MPI_Barrier(MPI_COMM_WORLD);
		DisplayMatrix_mpi(&arguments, &results, &options, comm.rank, comm.num_procs, comm.absoluteStartRow, comm.absoluteStartRow + comm.matrixRows -3);
		freeMatrices(&arguments);                                       /*  free memory     */
	}
	else
	{
		allocateMatrices(&arguments);        /*  get and initialize variables and matrices  */
		initMatrices(&arguments, &options);            /* ******************************************* */

		gettimeofday(&start_time, NULL);                   /*  start timer         */
		calculate(&arguments, &results, &options);                                      /*  solve the equation  */
		gettimeofday(&comp_time, NULL);                   /*  stop timer          */

		displayStatistics(&arguments, &results, &options);
		DisplayMatrix(&arguments, &results, &options);

		freeMatrices(&arguments);                                       /*  free memory     */
	}

	MPI_Finalize();

	return 0;
}
