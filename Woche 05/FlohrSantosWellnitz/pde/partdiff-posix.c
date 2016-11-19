/****************************************************************************/
/****************************************************************************/
/**                                                                        **/
/**                 TU München - Institut für Informatik                   **/
/**                                                                        **/
/** Copyright: Prof. Dr. Thomas Ludwig                                     **/
/**            Andreas C. Schmidt                                          **/
/**                                                                        **/
/** File:      partdiff-seq.c                                              **/
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
#include <pthread.h>

#include "partdiff-posix.h"

//Lock-Variable
pthread_mutex_t lock;

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

//Struct, um mit den einzelnen Threads Daten auszutauschen
struct calculate_thread_arguments
{
    double** 	matrixInput;           /* Matrix vor Berechnung */
    double** 	matrixOutput;          /* Matrix nach Berechnung */
    int 	spacesBetweenLines;    /* number of spaces between lines (lines=N+1)  */
    int 	firstRow;              /* Startreihe */
    int 	lastRow;               /* Letzte Reihe */
    double*   	maxresiduum;           /* Höchster Residuum Wert pro Thread */
    int* 	termIteration;         /* Anzahl der Schleifendurchläufe */
    uint64_t 	inferenceFunc;         /* Störfunktion */
    uint64_t 	termination;           /* Abbruchbedienung */
    double*	pih;                   /* Konstante für Störfunktion */
    double*	fpisin;                /* Konstante für Störfunktion */
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
/* calculate: Berechnung für einen Thread                                   */
/* ************************************************************************ */
void *doCalculate(void* args)
{
    //Cast der Übergebenen Argumente
    struct calculate_thread_arguments * arguments = (struct calculate_thread_arguments*) args;
    
    //Parameter aus Thread-Argumenten entnehmen
    double** Matrix_In = arguments->matrixInput;
    double** Matrix_Out = arguments->matrixOutput;
    int N = arguments->spacesBetweenLines;
    int lastRow = arguments->lastRow;
    double pih = *arguments->pih;
    double fpisin = *arguments->fpisin;
    
    double star;
    double residuum = 0;
    double* maxresiduum = arguments->maxresiduum;
    
    /* over all rows */
    for (int i = arguments->firstRow; i < lastRow; i++)
    {
        double fpisin_i = 0.0;
        
        if (arguments->inferenceFunc == FUNC_FPISIN)
        {
            fpisin_i = fpisin * sin(pih * (double)i);
        }
        
        /* over all columns */
        for (int j = 1; j < N; j++)
        {
            star = 0.25 * (Matrix_In[i-1][j] + Matrix_In[i][j-1] + Matrix_In[i][j+1] + Matrix_In[i+1][j]);
            
            if (arguments->inferenceFunc == FUNC_FPISIN)
            {
                star += fpisin_i * sin(pih * (double)j);
            }
            
            if (*arguments->termIteration == TERM_PREC || *arguments->termIteration == 1)
            {
                residuum = Matrix_In[i][j] - star;
                residuum = (residuum < 0) ? -residuum : residuum;
               
		//Prüfen, ob maxresiduum angepasst werden muss
		if(residuum > *maxresiduum)
		{
            		//maxresiduum für anderen Threads sperren
			pthread_mutex_lock(&lock);
            		
			//Prüfen, ob Anpassung nach dem Warten immernoch nötig ist
			if(residuum > *maxresiduum)
			{
				*maxresiduum = residuum;
			}
			pthread_mutex_unlock(&lock);
		}
		
        	//Maxresiduum für einzelne Thread-Versionen ermitelln (ohne Mutex)
		//maxresiduum = (residuum < maxresiduum) ? maxresiduum : residuum;
            }
            
            Matrix_Out[i][j] = star;
        }
    }
    
    arguments->maxresiduum = maxresiduum;
    
    return NULL;
}


/* ************************************************************************ */
/* calculate: solves the equation                                           */
/* ************************************************************************ */
static
void
calculate (struct calculation_arguments const* arguments, struct calculation_results* results, struct options const* options)
{
	unsigned int i;                                   /* local variables for loops */
	int m1, m2;                                 /* used as indices for old and new matrices */
	double maxresiduum;                         /* maximum residuum value of a slave in iteration */

	int const N = arguments->N;
	double const h = arguments->h;
    
	//Thread-Array mit Anzahl der eingegebenne Threads
    	pthread_t* threadArray = allocateMemory(sizeof(pthread_t) * options->number);
    	struct calculate_thread_arguments* thread_args =  allocateMemory(sizeof(struct calculate_thread_arguments) * options->number);

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
    
    	//Berechnung wie viele Reihen ein Thread berechnen muss
    	int rowsPerThread = (N - 1) / options->number;
    
    	//Berechnung der Restreihen, die nicht auf die Threads aufgeteilt werden konnten
    	int remainRows = (N - 1) % options->number;
    
    	//Start- und End-Reihen für die einzelnen Threads
    	int firstThreadRow = 1;
    	int lastThreadRow = 1;

    	//Zuweisung der festen Parameter für die einzelnen Threads
    	for(i = 0; i < options->number; i++)
    	{
        	//Berechnung der letzten Reihen, die in einem Thread berechnet werden muss
        	lastThreadRow += rowsPerThread;
        
        	//Aufteilen der Restreihen, bis keine mehr zu Verüfugung stehen
        	if (remainRows > 0)
        	{
                	remainRows = remainRows - 1;
                	lastThreadRow = lastThreadRow + 1;
        	}
        
        	//Parameter zuweisen
        	thread_args[i].firstRow = firstThreadRow;
        	thread_args[i].lastRow = lastThreadRow;
        	thread_args[i].inferenceFunc = options->inf_func;
        	thread_args[i].termIteration = &term_iteration;
        	thread_args[i].termination = options->termination;
        	thread_args[i].spacesBetweenLines = arguments->N;
        	thread_args[i].pih = &pih;
        	thread_args[i].fpisin = &fpisin;
        
        	//Berechnung der nächsten Startreihe
        	firstThreadRow = lastThreadRow;
    	}
	
    	while (term_iteration > 0)
    	{
        	//Zuweisung der variablen Parameter und starten der Threads
        	for(i = 0; i < options->number; i++)
        	{
                	thread_args[i].maxresiduum = &maxresiduum;
                	thread_args[i].matrixInput = arguments->Matrix[m2];
                	thread_args[i].matrixOutput = arguments->Matrix[m1];
            
               		pthread_create(&threadArray[i], NULL, &doCalculate, &thread_args[i]);
        	}
        
        	//maxresiduum = 0;
        	results->stat_iteration++;
        
        	//Alle Threads wieder zusammenführen und ggf. ohne Mutex: gemeinsames maxresiduum ermitteln
        	for(i = 0; i < options->number; i++)
        	{
            		if (pthread_join(threadArray[i], NULL))
            		{
                		printf("%s\n", "Join failed!");
                		return;
            		}
                
                    	//Auswerten der maxresiduum aller Threads, um das maxresiduum für alle zu ermitteln (ohne Mutex)
            		//maxresiduum = (thread_args[i].maxresiduum < maxresiduum) ? maxresiduum : thread_args[i].maxresiduum;
        	}
        
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
    
    	free(thread_args);
   	free(threadArray);

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

/* ************************************************************************ */
/*  main                                                                    */
/* ************************************************************************ */
int
main (int argc, char** argv)
{
	struct options options;
	struct calculation_arguments arguments;
	struct calculation_results results;

    	//Mutex initialisieren
	pthread_mutex_init (&lock, NULL);

	AskParams(&options, argc, argv);

	initVariables(&arguments, &results, &options);

	allocateMatrices(&arguments);
	initMatrices(&arguments, &options);

	gettimeofday(&start_time, NULL);
	calculate(&arguments, &results, &options);
	gettimeofday(&comp_time, NULL);

	displayStatistics(&arguments, &results, &options);
	DisplayMatrix(&arguments, &results, &options);

    	//Mutex freigeben
	pthread_mutex_destroy(&lock);
	freeMatrices(&arguments);

	return 0;
}
