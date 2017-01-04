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

//Konstanten
const int ROOT_RANK = 0;
const int TAG_PREVIOUS_ROW = 1;
const int TAG_NEXT_ROW = 2;

#define TAG_SEND_UPPER_ROW 1
#define TAG_SEND_LOWER_ROW 2
#define TAG_SEND_RESIDUUM 3
#define TAG_SEND_TERMINATION 4

//Parameter für die MPI-Kommunikation
struct mpi_calc_arguments
{
    //MPI-Parameter
    int rank;                       /*Aktueller rank*/
    int numberOfProcesses;          /*Anzahl der Gesamtprozesse*/
    MPI_Status status;              /*Der MPI-Status*/
    
    //Matrix-Paramter
    uint64_t startRowInTotalMatrix; /*Die Startreihe in der Gesamt-Matrix*/
    uint64_t matrixRows;            /*Anzahl der Teilmatrix-Reihen*/
    uint64_t matrixColumns;         /*Anzahl der Teilmatrix-Spalten*/
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
allocateMatricesMpi (struct calculation_arguments* arguments, struct mpi_calc_arguments* mpiArgs)
{
    uint64_t i, j;
    
    //Reihen für Prozess ermitteln
    const uint64_t rowsPerProcess = (arguments->N + 1) / mpiArgs->numberOfProcesses;
    
    //Anzahl restlicher Reihen
    const int64_t remainRows = (arguments->N + 1) % mpiArgs->numberOfProcesses;
    
    //Nur eine Reihe addieren, wenn erster oder letzter rank
    //Es wird nur Vor- oder Nachfolger benötigt
    if (mpiArgs->rank == ROOT_RANK || mpiArgs->rank == (mpiArgs->numberOfProcesses - 1))
    {
        mpiArgs->matrixRows  = rowsPerProcess + 1;
    }
    //Zwei addieren, wegen letzte Reihe Vorgänger und erste Reihe Nachfolger
    else
    {
        mpiArgs->matrixRows  = rowsPerProcess + 2;
    }
    
    //Anzahl Columns ermitteln
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
    
    int actual_splitted_remain_rows = 0;
    
    //Wenn Rank keine Extra-Reihe bekommen hat, dann müssen alle übergebliebenen Reihen hinzuaddiert werden
    if (mpiArgs->rank > remainRows)
    {
        actual_splitted_remain_rows = remainRows;
    }
    //Sonst müssen die bisherigen Extra-Reihen aufaddiert werden, die bisher aufgeteilt wurden (=Wert aktueller Rank)
    else
    {
        actual_splitted_remain_rows = mpiArgs->rank;
    }
    
    //Root-Rank startet bei 1, weil die erste Reihe der Matrix nicht berechnet wird
    if(mpiArgs->rank == ROOT_RANK)
    {
        mpiArgs->startRowInTotalMatrix = 1;
    }
    else
    {
        //Startreihe des Prozesses in der Gesamt-Matrix berechnen
        mpiArgs->startRowInTotalMatrix = (rowsPerProcess * mpiArgs->rank) + actual_splitted_remain_rows;
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
/* initMatrices: Initialize matrix/matrices and some global variables       */
/* ************************************************************************ */
static
void
initMatricesMpi (struct calculation_arguments* arguments, struct options const* options, struct mpi_calc_arguments* mpiArgs)
{
    uint64_t g, i, j;                                /*  local variables for loops   */
    
    const uint64_t rows = mpiArgs->matrixRows;
    const uint64_t columns = mpiArgs->matrixColumns;
    
    double const h = arguments->h;
    uint64_t const N = arguments->N;
    double*** Matrix = arguments->Matrix;
    
    /* initialize matrix/matrices with zeros */
    for (g = 0; g < arguments->num_matrices; g++)
    {
        for (i = 1; i < rows; i++)
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
        uint64_t startRowInTotalMatrix = mpiArgs->startRowInTotalMatrix; 
        
        for (g = 0; g < arguments->num_matrices; g++)
        {
            //Initialisieren der Ränder
            for (uint64_t i = 0; i < rows - 1; i++)
            {
                Matrix[g][i][0] = 1.0 - (h * (startRowInTotalMatrix + i - 1));
                Matrix[g][i][N] = h * (startRowInTotalMatrix + i - 1);
            }
            
            if (mpiArgs->rank == ROOT_RANK)
            {
                //Initialisieren der ersten Reihe
                for (uint64_t i = 0; i < columns; i++)
                {
                    Matrix[g][0][i] = 1.0 - (h * i);
                }
                Matrix[g][0][N] = 0;
            }
            else if (mpiArgs->rank == (mpiArgs->numberOfProcesses - 1))
            {
                //Initialisieren der letzten Reihe
                for (uint64_t i = 0; i < columns; i++)
                {
                    Matrix[g][rows - 1][i] = h * i;
                }
                Matrix[g][N][0] = 0;
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
calculateMpiJacobi (struct calculation_arguments const* arguments, struct calculation_results *results, struct options const* options, struct mpi_calc_arguments* mpiArgs)
{
    int i, j;                                   /* local variables for loops  */
    int m1, m2;                                 /* used as indices for old and new matrices       */
    double star;                                /* four times center value minus 4 neigh.b values */
    double residuum;                            /* residuum of current iteration                  */
    double maxresiduum;                         /* maximum residuum value of a slave in iteration */
    double const h = arguments->h;
    int columns = mpiArgs->matrixColumns;
    uint64_t startRowInTotalMatrix = mpiArgs->startRowInTotalMatrix;
    
    //Die Startreihe ist immer 1, denn der Root-Rank soll die erste Matrixzeile nicht berechnen
    //und die anderen Prozesse haben in der ersten Reihe die letzte Reihe des Vorgängers gespeichert
    int startRow = 1;
   
    //Letzte Reihe der Teilmatrix (über diese wird später nicht iteriert, weil sie entweder die erste Reihe des Nachfolgers ist
    //oder sie ist im letzten Rank die letzte Reihe der Gesamt-Matrix) 
    int lastRow = mpiArgs->matrixRows - 1;
    
    //Die Berechnung stimmen bei dem Root-Rank bzw letzten Rank nicht, weil der Root-Rank keinen Vorgänger und
    //der letzte keinen Nachfolger hat. Im Folgenden wird ein Senden bzw. Empfangen bei diesen Fällen jedoch verhindert,
    //sodass kein Fehler da durch aufritt.
    int nextTarget = mpiArgs->rank + 1;
    int previousTarget = mpiArgs->rank - 1;
    
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
        
        /***   Senden und Empfangen der Zeilen   ***/
        
        //Wenn nicht Root-Rank, dann erste eigene Reihe an Vorgänger senden
        if (mpiArgs->rank != ROOT_RANK)
        {
            MPI_Send(Matrix_In[1], mpiArgs->matrixColumns, MPI_DOUBLE, previousTarget, TAG_PREVIOUS_ROW, MPI_COMM_WORLD);
        }

        //Wenn nicht letzter Rank, dann letzte eigene Reihe an Nachfolger senden
        if (mpiArgs->rank != (mpiArgs->numberOfProcesses - 1))
        {
            MPI_Send(Matrix_In[lastRow - 1], mpiArgs->matrixColumns, MPI_DOUBLE, nextTarget, TAG_NEXT_ROW, MPI_COMM_WORLD);
        }

        //Wenn nicht Root-Rank, dann letzte Reihe vom Vorgänger Empfangen
        if (mpiArgs->rank != ROOT_RANK)
        {
            MPI_Recv(Matrix_In[0], mpiArgs->matrixColumns, MPI_DOUBLE,previousTarget, TAG_NEXT_ROW, MPI_COMM_WORLD, &mpiArgs->status);
        }
            
        //Wenn nicht letzter Rank, erst Reihe vom Nachfolger empfangen
        if (mpiArgs->rank != (mpiArgs->numberOfProcesses - 1))
        {
            MPI_Recv(Matrix_In[lastRow], mpiArgs->matrixColumns, MPI_DOUBLE,nextTarget, TAG_PREVIOUS_ROW, MPI_COMM_WORLD, &mpiArgs->status);
        }
        
        
        /* over all rows */
        for (i = startRow; i < lastRow; i++)
        {
            double fpisin_i = 0.0;
            
            if (options->inf_func == FUNC_FPISIN)
            {
                fpisin_i = fpisin * sin(pih * (double)(startRowInTotalMatrix + i - startRow));
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
        
        //Das MaxReisduum aller Prozesse ermitteln
        double maxresiduumOverAllProcesses;
        MPI_Allreduce(&maxresiduum, &maxresiduumOverAllProcesses, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        
        maxresiduum = maxresiduumOverAllProcesses;
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
calculateMpiGauss (struct calculation_arguments const* arguments, struct calculation_results *results, struct options const* options, struct mpi_calc_arguments* mpiArgs)
{
    //Initialisierung ist fast identisch mit Jacobi, außer dass nur eine Matrix benötigt wird.
    
    int i, j;                                   /* local variables for loops  */
    int m1, m2;                                 /* used as indices for old and new matrices */
    double star;                                /* four times center value minus 4 neigh.b values */
    double residuum;                            /* residuum of current iteration                  */
    double maxresiduum, prevMaxResiduum;        /* maximum residuum value of a slave in iteration */
    double const h = arguments->h;
    int columns = mpiArgs->matrixColumns;
    uint64_t startRowInTotalMatrix = mpiArgs->startRowInTotalMatrix;
    
    //Die Startreihe ist immer 1, denn der Root-Rank soll die erste Matrixzeile nicht berechnen
    //und die anderen Prozesse haben in der ersten Reihe die letzte Reihe des Vorgängers gespeichert
    int startRow = 1;
    
    //Letzte Reihe der Teilmatrix (über diese wird später nicht iteriert, weil sie entweder die erste Reihe des Nachfolgers ist
    //oder sie ist im letzten Rank die letzte Reihe der Gesamt-Matrix)
    int lastRow = mpiArgs->matrixRows - 1;
    
    //Die Berechnung stimmen bei dem Root-Rank bzw letzten Rank nicht, weil der Root-Rank keinen Vorgänger und
    //der letzte keinen Nachfolger hat. Im Folgenden wird ein Senden bzw. Empfangen bei diesen Fällen jedoch verhindert,
    //sodass kein Fehler da durch aufritt.
    int nextTarget = mpiArgs->rank + 1;
    int previousTarget = mpiArgs->rank - 1;
    
    int stopSend = 0;
    int stopReceived = 0; 
    double pih = 0.0;
    double fpisin = 0.0;
    int term_iteration = options->term_iteration;
    
    m1 = 0;
    m2 = 0;
    double** Matrix_In = arguments->Matrix[m1];
    double** Matrix_Out = arguments->Matrix[m2];
    
    if (options->inf_func == FUNC_FPISIN)
    {
        pih = PI * h;
        fpisin = 0.25 * TWO_PI_SQUARE * h * h;
    }
    
    MPI_Request stopSignal;
    if (options->termination == TERM_PREC && mpiArgs->rank == ROOT_RANK)
    {
	MPI_Irecv(&stopReceived, 1, MPI_INT, (mpiArgs->numberOfProcesses - 1), TAG_SEND_TERMINATION, MPI_COMM_WORLD, &stopSignal);
    }
    else if (options->termination == TERM_PREC)
    {
        MPI_Irecv(&stopReceived, 1, MPI_INT, previousTarget, TAG_SEND_TERMINATION, MPI_COMM_WORLD, &stopSignal);
    }
    
    while (term_iteration > 0)
    { 
        if(mpiArgs->rank != ROOT_RANK)
        {
            //Warten auf die Letzte Zeile vom Vorgänger (Startsignal für unsere Iteration). Diese Zeile ist bereits vom Vorgänger berechnet worden.
            //MPI_Recv(LetzteZeileVomVorgaengerEmpfangen);
            MPI_Recv(Matrix_In[0], mpiArgs->matrixColumns, MPI_DOUBLE, previousTarget, TAG_SEND_LOWER_ROW, MPI_COMM_WORLD, &mpiArgs->status);
            
        }
        
        maxresiduum = 0;
        
        //Die Berechnung an sich ist eigentlich ziemlich identsich mit jacobi, außer dass Matrix_In and Matrix_Out jetzt die auf die gleiche Matrix zeigen
        /* over all rows */
        for (i = startRow; i < lastRow; i++)
        {
            double fpisin_i = 0.0;
            
            if (options->inf_func == FUNC_FPISIN)
            {
                fpisin_i = fpisin * sin(pih * (double)(startRowInTotalMatrix + i - startRow));
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
            
            //Nach Berechnung der ersten Zeile, wenn nicht Root-Rank, an den Vorgänger senden. Diese ist für unseren Vorgänger 'unverändert', weil dieser einer Iteration weiter ist
            if(i == startRow && mpiArgs->rank != ROOT_RANK)
            {
                MPI_Send(Matrix_In[1], mpiArgs->matrixColumns, MPI_DOUBLE, previousTarget, TAG_SEND_UPPER_ROW, MPI_COMM_WORLD);
            }
        }
        
        results->stat_iteration++;
        
        if (mpiArgs->rank != (mpiArgs->numberOfProcesses - 1))
        {
            //Nach Berechnung die letzte Zeile an den Nachfolger senden, damit dieser die nächste Iteration starten kann (Startsignal für Nachfolger)
            MPI_Send(Matrix_In[lastRow - 1], mpiArgs->matrixColumns, MPI_DOUBLE, nextTarget, TAG_SEND_LOWER_ROW, MPI_COMM_WORLD);

            //Erste Zeile vom Nachfolger empfangen, die er nach seiner ersten Berechnung sofort schickt.
            //Dies ist eine Vorbereitung für unsere nächste Iteration, denn die Erste Zeile vom Nachfolger ist für uns quasi unverändert.
            MPI_Recv(Matrix_In[lastRow], mpiArgs->matrixColumns, MPI_DOUBLE, nextTarget, TAG_SEND_UPPER_ROW, MPI_COMM_WORLD, &mpiArgs->status);
      	}
        
	/* check for stopping calculation, depending on termination method */
        //Das MaxResiduum muss nur weitergereicht werden, wenn die Präzision als Abbruchbedienung gewählt wurde, sonst ist es unnötig
        if (options->termination == TERM_PREC)
        {
            //MaxResiduum ermitteln
            if (mpiArgs->rank != ROOT_RANK)
            {
                //MPI_Recv(MaxResiduumVomVorgaengerEmpfangen);
                MPI_Recv(&prevMaxResiduum, 1, MPI_DOUBLE, previousTarget, TAG_SEND_RESIDUUM, MPI_COMM_WORLD, &mpiArgs->status);
                
		if(prevMaxResiduum > maxresiduum)
                {
                    maxresiduum = prevMaxResiduum;
                }
            }
            
            if(stopReceived)
            {
		term_iteration = 0;
       	 	
		if (mpiArgs->rank != (mpiArgs->numberOfProcesses - 1))
        	{
            		MPI_Send(&stopReceived, 1, MPI_INT, nextTarget, TAG_SEND_TERMINATION, MPI_COMM_WORLD);
        	}
            }
            
            //MaxResiduum weiter an Nachfolger senden, wenn nicht letzter Rank
            if (mpiArgs->rank != (mpiArgs->numberOfProcesses - 1))
            {
                MPI_Send(&maxresiduum, 1, MPI_DOUBLE, nextTarget, TAG_SEND_RESIDUUM, MPI_COMM_WORLD);
            }
            
	    if(mpiArgs->rank == (mpiArgs->numberOfProcesses - 1))
	    {
                //Wenn letzter Rank und Präzision erreicht, dann Stopsignal an Root
                if (maxresiduum < options->term_precision && !stopSend)
                {
                    stopSend = 1;
                    MPI_Send(&stopSend, 1, MPI_INT, ROOT_RANK, TAG_SEND_TERMINATION, MPI_COMM_WORLD);
		}
            }
        }
        else if (options->termination == TERM_ITER)
        {
            term_iteration--;
        }
    }
    
    //Auf alle Prozesse warten...
    MPI_Barrier(MPI_COMM_WORLD);
    
    //... und danach das MaxReisduum aller Prozesse ermitteln
    double maxresiduumOverAllProcesses;
    MPI_Allreduce(&maxresiduum, &maxresiduumOverAllProcesses, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    results->stat_precision = maxresiduumOverAllProcesses;
    
    results->m = m1;
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
DisplayMatrixMpi (struct calculation_arguments* arguments, struct calculation_results* results, struct options* options, int rank, int size, int from, int to)
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
    MPI_Comm_size(MPI_COMM_WORLD, &mpiArgs.numberOfProcesses);

    //AskParams um einen Parameter erweitert, damit nur der Root-Rank die Info am Anfang anzeigt
	AskParams(&options, argc, argv, mpiArgs.rank);

	initVariables(&arguments, &results, &options);

    //Nur wenn Jacobi und mehr als 1 Prozess
    if (options.method == METH_JACOBI && mpiArgs.numberOfProcesses > 1)
    {
        allocateMatricesMpi(&arguments, &mpiArgs);
        initMatricesMpi(&arguments, &options, &mpiArgs);
        
        gettimeofday(&start_time, NULL);
        calculateMpiJacobi(&arguments, &results, &options, &mpiArgs);
        gettimeofday(&comp_time, NULL);
        
        //Nur Root-Rank soll die Statistiken anzeigen
        if (mpiArgs.rank == ROOT_RANK)
        {
            displayStatistics(&arguments, &results, &options);
        }
        
        DisplayMatrixMpi(&arguments, &results, &options, mpiArgs.rank, mpiArgs.numberOfProcesses, mpiArgs.startRowInTotalMatrix, mpiArgs.startRowInTotalMatrix + mpiArgs.matrixRows -3);
    }
    else if(options.method == METH_GAUSS_SEIDEL && mpiArgs.numberOfProcesses > 1)
    {
        allocateMatricesMpi(&arguments, &mpiArgs);
        initMatricesMpi(&arguments, &options, &mpiArgs);
    
        gettimeofday(&start_time, NULL);
        calculateMpiGauss(&arguments, &results, &options, &mpiArgs);
        gettimeofday(&comp_time, NULL);
    
        //Nur Root-Rank soll die Statistiken anzeigen
        if (mpiArgs.rank == ROOT_RANK)
        {
            displayStatistics(&arguments, &results, &options);
        }
    
        DisplayMatrixMpi(&arguments, &results, &options, mpiArgs.rank, mpiArgs.numberOfProcesses, mpiArgs.startRowInTotalMatrix, mpiArgs.startRowInTotalMatrix + mpiArgs.matrixRows -3);
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
