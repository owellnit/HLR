#include <stdio.h>
#include <mpi.h>
#include <time.h>
#include <unistd.h>

const HOSTENAME_LENGTH = 40;
const TIMESTAMP_LENGTH = 40;
const BUFFER_SIZE = HOSTENAME_LENGTH + TIMESTAMP_LENGTH;

int main(int argc, char **argv)
{
  int numberOfProcesses, actualProcess;
  MPI_Status status;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &actualProcess);
  MPI_Comm_size(MPI_COMM_WORLD, &numberOfProcesses);

  if (actualProcess == 0)
  {
  	char buff[BUFF_SIZE];
  	for (int i = 1; i < numberOfProcesses; i++)
  	{
  		MPI_Recv(buff, BUFF_SIZE, MPI_CHAR, i, 0, MPI_COMM_WORLD, &status);
  		printf("%s\n", buff); 
  	}
  }
  else
  {
  	char hostname[HOSTENAME_LENGTH];
    char timestamp[TIMESTAMP_LENGTH];
    char buffer[BUFFER_SIZE];
  	struct timeval time;

  	gethostname(hostname, MAX_HOSTNAME_LEN);
    gettimeofday(&time, NULL);
    
    strftime(time, sizeof(time), "%F %T", localtime(&ts.tv_sec));
      
    sprintf(buff, "%s: %s.%09ld", hostname, time, ts.tv_nsec);
    
    MPI_Send(buff, BUFF_SIZE, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  printf("Rang %d beendet jetzt!\n", actualProcess);
  ierr = MPI_Finalize();
}
