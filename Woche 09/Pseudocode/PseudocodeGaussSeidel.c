static
void
calculateMpiGauss (struct calculation_arguments const* arguments, struct calculation_results *results, struct options const* options, struct mpi_calc_arguments* mpiArgs)
{
    //Initialisierung ist fast identisch mit Jacobi, außer dass nur eine Matrix benötigt wird.
    
    int i, j;                                   /* local variables for loops  */
    int m1, m2;                                 /* used as indices for old and new matrices */
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
    
    m1 = 0;
    m2 = 0;
    double** Matrix_In = arguments->Matrix[m1];
    double** Matrix_Out = arguments->Matrix[m2];
    
    if (options->inf_func == FUNC_FPISIN)
    {
        pih = PI * h;
        fpisin = 0.25 * TWO_PI_SQUARE * h * h;
    }

    while (term_iteration > 0)
    {
        //Wie das genau mit dem Stop laufen soll, steht noch in den Sternen wie der Star Wars Abspann :D
        //Wenn Abbruch nach Präzision und Root-Rank, dann irgendwie gucken, obs ein Stopsignal gibt,...
        if (PräzisionErreichtUndRootRank)
        {
            StopsignalEmpfangen;
        }
        //...wenn ja, dann abbrechen und Kommunikation sauber beenden
        if(StopsignalEmpfangen)
        {
            SauberesBeendenDerKommunikation();
        }
        
        if (nichtRootRank)
        {
            //Warten auf die Letzte Zeile vom Vorgänger (Startsignal für unsere Iteration). Diese Zeile ist bereits vom Vorgänger berechnet worden.
            MPI_Recv(LetzteZeileVomVorgaengerEmpfangen);
            
            //Wenn nötig, das MaxResiduum des Vorgängers empfangen
            if (options->termination == TERM_PREC)
            {
                MPI_Recv(MaxResiduumVomVorgaengerEmpfangen);
            }
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
            if (startRowUndNichtRootRank)
            {
                MPI_Send(ErsteZeileAnVorgaengerSenden);
            }
        }
        
        results->stat_iteration++;
        
        if (NichtLetzterRank)
        {
            //Nach Berechnung die letzte Zeile an den Nachfolger senden, damit dieser die nächste Iteration starten kann (Startsignal für Nachfolger)
            MPI_Send(LetzteZeileAnNachfolgerSenden);
            
            //Erste Zeile vom Vorgänger empfangen, die er nach seiner ersten Berechnung sofort schickt.
            //Dies ist eine Vorbereitung für unsere nächste Iteration, denn die Erste Zeile vom Nachfolger ist für uns quasi unverändert.
            MPI_Recv(ErsteZeileVomNachfolgerEmpfangen);
        }
        
        /* check for stopping calculation, depending on termination method */
        //Das MaxResiduum muss nur weitergereicht werden, wenn die Präzision als Abbruchbedienung gewählt wurde, sonst ist es unnötig
        if (options->termination == TERM_PREC)
        {
            //MaxResiduum ermitteln
            if (nichtRoot)
            {
                if(MaxResiduumVorgaenger > maxresiduum)
                {
                    maxresiduum = MaxResiduumVorgaenger;
                }
            }
            
            //MaxResiduum weiter an Nachfolger senden, wenn nicht letzter Rank
            if (NichtLetzterRank)
            {
                MPI_Send(MaxResiduumAnNachfolgerSenden);
            }
            
            //Wenn letzter Rank und Präzision erreicht, dann Stopsignal an Root
            if (PräzisionErreichtUndLetzterRank)
            {
                MPI_Send(StopsignalSenden);
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
    
    results->m = 0;
}
