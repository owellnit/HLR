#include <stdio.h>
#include <mpi.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include <limits.h>

const int HOSTNAME_LENGTH = 40;
const int TIMESTAMP_LENGTH = 40;
const int BUFFER_SIZE = 80;

int main(int argc, char **argv)
{
    //Anzahl der gesamt Prozesse und der aktuelle Prozess
    int numberOfProcesses, actualProcess;
    MPI_Status status;
    
    //MPI initialisieren
    MPI_Init(&argc, &argv);
    //Die Nummer des aktuellen Prozesses ermitteln
    MPI_Comm_rank(MPI_COMM_WORLD, &actualProcess);
    
    //Die Anzahl aller Prozesse ermitteln
    MPI_Comm_size(MPI_COMM_WORLD, &numberOfProcesses);
    
    //Wenn der aktuelle Prozess die Nummer 0 hat, in den Lesemodus
    if (actualProcess == 0)
    {
        //Nachrichten der anderen Threads der Reihenfolge nach auslesen
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
        
        //Hostname und Time ermitteln
        gethostname(hostname, HOSTNAME_LENGTH);
        gettimeofday(&time, NULL);
        
        //Time in String wandeln
        nowtime = time.tv_sec;
        strftime(timestamp, sizeof(timestamp), "%F %T", localtime(&nowtime));
        
        //Ausgabe mit Hostname und Timestamp zusammenbauen
        sprintf(buffer, "%s: %s.%06ld", hostname, timestamp, time.tv_usec);
        
        //Daten über MPI senden
        MPI_Send(buffer, BUFFER_SIZE, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
    }
    
    //'Wartepunkt' für die Threads
    MPI_Barrier(MPI_COMM_WORLD);
    
    //Melden, dass der THread fertig ist
    printf("Rang %d beendet jetzt!\n", actualProcess);
    MPI_Finalize();
}
