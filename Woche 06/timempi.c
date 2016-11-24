#include <stdio.h>
#include <mpi.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>

const int HOSTNAME_LENGTH = 40;
const int TIMESTAMP_LENGTH = 40;
const int BUFFER_SIZE = 80;

int main(int argc, char **argv)
{
  int numberOfProcesses, actualProcess;
  MPI_Status status;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &actualProcess);
  MPI_Comm_size(MPI_COMM_WORLD, &numberOfProcesses);

  if (actualProcess == 0)
  {
  	char buffer[BUFFER_SIZE];
  	for (int i = 1; i < numberOfProcesses; i++)
  	{
  		MPI_Recv(buffer, BUFFER_SIZE, MPI_CHAR, i, 0, MPI_COMM_WORLD, &status);
  		printf("%s\n", buffer); 
  	}
  }
  else
  {
  	char hostname[HOSTNAME_LENGTH];
    	char timestamp[TIMESTAMP_LENGTH];
    	char buffer[BUFFER_SIZE];
  	struct timeval time;
	time_t nowtime;

  	gethostname(hostname, HOSTNAME_LENGTH);
    	gettimeofday(&time, NULL);
    
	nowtime = time.tv_sec;
    	strftime(timestamp, sizeof(timestamp), "%F %T", localtime(&nowtime));
      
    	sprintf(buffer, "%s: %s.%06ld", hostname, timestamp, time.tv_usec);
    
    	MPI_Send(buffer, BUFFER_SIZE, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  printf("Rang %d beendet jetzt!\n", actualProcess);
  MPI_Finalize();
}
