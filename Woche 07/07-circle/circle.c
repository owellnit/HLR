#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <ctype.h>
#include <mpi.h>



int*
init (int N)
{
    // Todo
    int* buf = malloc(sizeof(int) * N);
    
    srand(time(NULL));
    
    for (int i = 0; i < N; i++)
    {
        buf[i] = rand() % 25; //do not modify % 25
    }
    
    return buf;
}

//Erweiterung des Groesse-parameters, damit man das letzte Element leichter swappen kann

int*
circle (int* buf,int size)
{
    int i;
    int temp;
    temp = buf[size-1];
    
    for(i = sizeof(buf)-1; i > 0;--i)
    {
        buf[i] = buf[i-1];
    }
    buf[0] = temp;
    
    return buf;
}

int
main (int argc, char** argv)
{
    
    MPI_Status status;
    MPI_Request send_request, recv_request;
    
    char arg[256];
    int N;
    int rank, nprocs;
    int* buff;
    int* recvbuff;
    
    
    if (argc < 2)
    {
        printf("Arguments error!\n");
        return EXIT_FAILURE;
    }
    
    sscanf(argv[1], "%s", arg);
    
    // Array length
    //String to int atoi(String) --> int in N
    N = atoi(arg);
    
    //Check if argument non negativ and not letter
    if(N <= 0 || isdigit(N))
    {
        printf("No negativ arguments or argument not a number!\n");
        return EXIT_FAILURE;
    }
    
    //Determine Procsnumber
    if(N%2 == 0)
    {
        nprocs = N/2;
    }
    else
    {
        nprocs = N/2+1;
    }
    
    buff = init(N);
    int bufsize = sizeof(buff);
    recvbuff=(int *)malloc(sizeof(int)*bufsize);
    
    
    //MPI Initialize
    MPI_Init(&argc,&argv);
    //Welchen Prozess ermitteln
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //Anzahl der Prozesse
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    
    
    //Printing before ********************************************
    printf("\nBEFORE\n");
    
    for (int i = 0; i < N; i++)
    {
        printf ("rank %d: %d\n", rank, buff[i]);
    }
    
    //int first = buff[0];
    //Send communication n-1,n,n+1 ******************************
    
    
    
    for(int ind = 0; ind < nprocs-1; ++ind)
    {
        int first = buff[0];
        int last = buff[N-1];
        if(first != last)
        {
            recvbuff = circle(buff,N);
            if(rank == 0)
            {
                MPI_Isend(buff,bufsize,MPI_INT,ind+1,0,MPI_COMM_WORLD,&send_request);
                MPI_Irecv(recvbuff,bufsize,MPI_INT,ind-1,MPI_ANY_TAG,MPI_COMM_WORLD,&recv_request);
                
            }
            else if(rank == nprocs - 1)
            {
                MPI_Isend(buff,bufsize,MPI_INT,0,0,MPI_COMM_WORLD,&send_request);
                MPI_Irecv(recvbuff,bufsize,MPI_INT,ind-1,MPI_ANY_TAG,MPI_COMM_WORLD,&recv_request);
            }
            else
            {
                MPI_Isend(buff,bufsize,MPI_INT,ind+1,0,MPI_COMM_WORLD,&send_request);
                MPI_Irecv(recvbuff,bufsize,MPI_INT,ind-1,MPI_ANY_TAG,MPI_COMM_WORLD,&recv_request);
            }
        }
        else
        {
            printf("Erstes Element erreicht");
            break;
        }
    }
    
    //'Wartepunkt' für die Threads
    MPI_Barrier(MPI_COMM_WORLD); 
    //Nach Rotationen *********************************************
    printf("\nAFTER\n");
    
    for (int j = 0; j < N; j++)
    {
        printf ("rank %d: %d\n", rank, recvbuff[j]);
    }
   MPI_Finalize();    
    
    return EXIT_SUCCESS;
}

