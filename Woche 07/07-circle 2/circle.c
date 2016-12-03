#include <stdio.h>
#include <stdlib.h>
#include <time.h>

typedef struct
{
    int rank;
    int num_procs;
    int* bufsize;
    MPI_Status status;
} CommParams;

int*
init (int N, CommParams* params)
{
    //Anzahl der Zahlen pro Rank
    int chunkSize = completeSize / params->num_procs;
    int extra = completeSize % params->num_procs;
    
    //Anzahl der Buffsizes
    params->bufsize = malloc(sizeof(int) * params->num_procs);
    
    //Ermitteln der Buffsizes der Ranks
    for (int i = 0; i < params->num_procs - 1; i++)
    {
        params->bufsize[i] = chunkSize;
        
        if (extra > 0)
        {
            params->bufsize[i] += 1;
            extra -= 1;
        }
    }
    
    int* buf = malloc(sizeof(int) * N);

    srand(time(NULL));

    for (int i = 0; i < N; i++)
    {
      buf[i] = rand() % 25; //do not modify % 25
    }

    return buf;
}

int*
circle (int* buf)
{
    // Todo
    return buf;
}

int
main (int argc, char** argv)
{
    char arg[256];
    int N;
    int rank;
    int* buf;
    CommParams params;


    if (argc < 2)
    {
        printf("Arguments error!\n");
        return EXIT_FAILURE;
    }

    sscanf(argv[1], "%s", arg);

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &params.rank);
    MPI_Comm_size(MPI_COMM_WORLD, &params.num_procs);
    
    // Array length
    N = atoi(arg);
    buf = init(N, &params);

    // Todo: myrank
    rank = 0;

    printf("\nBEFORE\n");

    for (int i = 0; i < N; i++)
    {
        printf ("rank %d: %d\n", rank, buf[i]);
    }

    circle(buf);

    printf("\nAFTER\n");

    for (int j = 0; j < N; j++)
    {
        printf ("rank %d: %d\n", rank, buf[j]);
    }

    return EXIT_SUCCESS;
}
