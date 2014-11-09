#include <stdlib.h>
#include <stdio.h>

/* Represents each tile in the board */
typedef struct {
    char bug; /* 1 if there's a bug on this tile, 0 otherwise */
    double emission; /* positive when there's a heat source on this tile,
                        negative when there's a "cold source", and zero
                        otherwise. */
} Tile;

int getInteger(char*, char*);
double getFloat(char*, char*);
void printUsageAndExit(void);

int main(int argc, char *argv[]) {
    int w, h, n, s, nh, nc, t, np, i, j, x, y;
    double c, ph, pc;
    Tile **board;
    
    /* Must receive 11 arguments */
    if (argc != 12)
        printUsageAndExit();
    
    /* Getting argument values */
    w = getInteger(argv[1], "W");
    h = getInteger(argv[2], "H");
    n = getInteger(argv[3], "N");
    s = getInteger(argv[4], "S");
    c = getFloat(argv[5], "C");
    ph = getFloat(argv[6], "PH");
    nh = getInteger(argv[7], "NH");
    pc = getFloat(argv[8], "PC");
    nc = getInteger(argv[9], "NC");
    t = getInteger(argv[10], "T");
    np = getInteger(argv[11], "NP");
    
    /* Constructing the board */
    board = malloc(h * sizeof(Tile *));
    for (i = 0; i < h; i++) {
        board[i] = malloc(w * sizeof(Tile));
        for (j = 0; j < w; j++) {
            board[i][j].bug = 0;
            board[i][j].emission = 0.0;
        }
    }
    
    /* Generating bugs */
    srand(s);
    for (i = 0; i < n; i++) {
        x = rand() % w;
        y = rand() % h;
        if (board[y][x].bug) i--;
        else board[y][x].bug = 1;
    }
    
    /* Printing */
    for (i = 0; i < h; i++) {
        for (j = 0; j < w; j++)
            printf("%d ", board[i][j].bug);
        printf("\n");
    }
    
    /* Freeing memory */
    for (i = 0; i < h; i++)
        free(board[i]);
    free(board);
    
    return 0;
}

int getInteger(char *arg, char *name) {
    int x = atoi(arg);
    if (x <= 0) {
        printf("The value of %s must be a positive integer.\n\n", name);
        printUsageAndExit();
    }
    return x;
}

double getFloat(char *arg, char *name) {
    double x = atof(arg);
    if (x <= 0.0) {
        printf("The value of %s must be a positive number.\n\n", name);
        printUsageAndExit();
    }
    return x;
}

void printUsageAndExit() {
    printf("Usage: hot-bugs W H N S C PH NH PC NC T NP\n=====\n\
W\twidth of the board\n\
H\theight of the board\n\
N\tnumber of bugs\n\
S\tseed for the random number generator\n\
C\tconstant of heat emission by the bugs\n\
PH\tprobability of heat source appearance\n\
NH\tnumber of iterations the heat sources will last\n\
PC\tprobability of \"cold source\" appearance\n\
NC\tnumber of iterations the \"cold sources\" will last\n\
T\ttotal iterations for the simulation\n\
NP\tnumber of processors to be used\n");
    exit(1);
}
