#include <stdio.h>
#include <mpi.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>

#define MAX_HOSTNAME_LEN 50
#define MAX_TIME_LEN 50
#define BUFF_SIZE MAX_HOSTNAME_LEN + MAX_TIME_LEN

int main(int argc, char **argv)
{
  int ierr, num_procs, my_id;
  int mymicro = 1000000;
  int absoluteMinMicro;
  MPI_Status status;

  ierr = MPI_Init(&argc, &argv);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

  if (my_id == 0)
  {
    char buff[BUFF_SIZE];
    for (int i = 1; i < num_procs; i++)
    {
      ierr = MPI_Recv(buff, BUFF_SIZE, MPI_CHAR, i, 0, MPI_COMM_WORLD, &status);
      printf("%s\n", buff); 
    }
  }
  else
  {
    char hostname[MAX_HOSTNAME_LEN];
    char time[MAX_TIME_LEN];
    char buff[BUFF_SIZE];
    struct timeval tv;

    gethostname(hostname, MAX_HOSTNAME_LEN);
    gettimeofday(&tv, NULL);

    strftime(time, sizeof(time), "%F %T", localtime(&tv.tv_sec));
    sprintf(buff, "%s: %s.%06ld",hostname, time, tv.tv_usec);

    mymicro = tv.tv_usec;
    ierr = MPI_Send(buff, BUFF_SIZE, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
  }

  MPI_Reduce(&mymicro, &absoluteMinMicro, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);

  if (my_id == 0)
  {
    printf("%d\n", absoluteMinMicro);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  printf("Rang %d beendet jetzt!\n", my_id);
  ierr = MPI_Finalize();
}
