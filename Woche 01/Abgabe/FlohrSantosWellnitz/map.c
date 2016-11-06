#include <stdio.h>

// Definieren Sie ein enum cardd
typedef enum {
    N = 1 << 0, /* 0b0000000000000001 */
    W = 1 << 1, /* 0b0000000000000010 */
    E = 1 << 2, /* 0b0000000000000010 */
    S = 1 << 3  /* 0b0000000000001000 */
} cardd;

// Definieren Sie ein 3x3-Array namens map, das Werte vom Typ cardd enthält
cardd map[3][3];

//Überprüft die korrekte Eingabe und liefert das Ergebnis zurück
int checkDir(cardd dir)
{
    return (dir == N || dir == W ||dir == E ||dir == S ||
            dir == (N|W) || dir == (N|E) ||dir == (S|W) ||dir == (S|E));
}

// Die Funktion set_dir soll an Position x, y den Wert dir in das Array map eintragen
// Überprüfen Sie x und y um mögliche Arrayüberläufe zu verhindern
// Überprüfen Sie außerdem dir auf Gültigkeit
void set_dir (int x, int y, cardd dir)
{
    //Überlauf verhindern
    if((x < 3 && x >= 0) && (y < 3 && y >= 0) && checkDir(dir))
    {
        //Himmelsrichtung in Matrix setzen
        map[x][y] = dir;
    }
}

//Liefert die Himmelsrichtung als char-Array (String )zurück
char* getDirectionString(cardd dir)
{
    //Fallunterscheidung, um die korrekte Himmelsrichtung als Char-Array (String) auszugeben
    switch((int)dir)
    {
        case N:
            return "N";
        case W:
            return "W";
        case E:
            return "E";
        case S:
            return "S";
        case N|W:
            return "NW";
        case N|E:
            return "NE";
        case S|W:
            return "SW";
        case S|E:
            return "SE";
        default:
            return "0";
    }
}


// Die Funktion show_map soll das Array in Form einer 3x3-Matrix ausgeben
void show_map (void)
{
    //Durchlauf der X-Koordinaten
    for(int x = 0; x < 3; x++)
    {
        //Status, ob eine Reihe mit zwei Chars startet (0 = false)
        int startsWithTwoChars = 0;
        
        //Durchlauf der Y-Koordniaten
        for(int y = 0; y < 3; y++)
        {
            char* direction;
            cardd dir;
            
            dir = map[x][y];
            
            direction = getDirectionString(dir);
            
            //Status, ob eine Reihe mit zwei Chars überprüfen (1 = true)
            if(dir == (N|W) || dir == (S|W))
            {
                startsWithTwoChars = 1;
            }
            
            //Himmelsrichtung ausgeben, je nach Gegebenheit formatieren
            if(y == 2)
            {
                printf("%s", direction);
            }
            else if (startsWithTwoChars)
            {
                printf("%s__", direction);
            }
            else
            {
                printf("%s___", direction);
            }
        }
        printf("\n");
    }
    
    printf("\n");
    printf("\n");
    printf("\n");
}

int main (void)
{
    // In dieser Funktion darf nichts verändert werden!
    set_dir(0, 1, N);
    set_dir(1, 0, W);
    set_dir(1, 4, W);
    set_dir(1, 2, E);
    set_dir(2, 1, S);
    
    
    show_map();
    
    set_dir(0, 0, N|W);
    set_dir(0, 2, N|E);
    set_dir(0, 2, N|S);
    set_dir(2, 0, S|W);
    set_dir(2, 2, S|E);
    set_dir(2, 2, E|W);
    
    show_map();
    
    return 0;
}
