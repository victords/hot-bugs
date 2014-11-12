#include <stdlib.h>
#include <stdio.h>

/* Represents each tile in the board */
typedef struct Tile *tile;
struct Tile {
    char bug; /* 1 if there's a bug on this tile, 0 otherwise */
    double emission; /* positive when there's a heat source on this tile,
                        negative when there's a "cold source", and zero
                        otherwise. */
};

int getInteger(char*, char*);
double getFloat(char*, char*);
void printUsageAndExit(void);
tile *getNeighbors(tile**, int, int, int, int);

int main(int argc, char *argv[]) {
    int w, h, n, s, nh, nc, t, np, i, j, x, y;
    double c, ph, pc;
    tile **board;

    /* getNeighbors test variables */
    int a, b;
    tile* neigh;
    
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
    board = malloc(h * sizeof(tile *));
    for (i = 0; i < h; i++) {
        board[i] = malloc(w * sizeof(tile));
        for (j = 0; j < w; j++) {
            board[i][j] = malloc(sizeof(struct Tile));
            board[i][j]->bug = 0;
            board[i][j]->emission = 0.0;
        }
    }
    
    /* Generating bugs */
    srand(s);
    for (i = 0; i < n; i++) {
        x = rand() % w;
        y = rand() % h;
        if (board[y][x]->bug) i--;
        else board[y][x]->bug = 1;
    }
    
    /* Printing */
    for (i = 0; i < h; i++) {
        for (j = 0; j < w; j++)
            printf("%d ", board[i][j]->bug);
        printf("\n");
    }

    /* getNeighbors test */
    scanf("%d %d", &a, &b);
    neigh = getNeighbors(board, a, b, w, h);
    
    for (j = 0; j < 6; j++)
        if(neigh[j] == NULL)
            printf(" ");
        else
            printf("%d ", neigh[j]->bug);
    printf("\n");

    /* Freeing memory */
    for (i = 0; i < h; i++)
        free(board[i]);
    free(board);
    
    return 0;
}

tile *getNeighbors(tile **board, int x, int y, int w, int h) {
    tile *neighbors = malloc(6 * sizeof(tile));

    if(x % 2 == 0) {
        if (x == 0)
            neighbors[0] = NULL;
        else
            neighbors[0] = board[x-1][y];
        if (x == 0 || y == w - 1)
            neighbors[1] = NULL;
        else    
            neighbors[1] = board[x-1][y+1];
        if (y == 0)
            neighbors[2] = NULL;
        else
            neighbors[2] = board[x][y-1];
        if (y == w - 1)
            neighbors[3] = NULL;
        else
            neighbors[3] = board[x][y+1];
        if (x == h - 1)
            neighbors[4] = NULL;
        else
            neighbors[4] = board[x+1][y];
        if (x == h - 1 || y == w - 1)
            neighbors[5] = NULL;
        else
            neighbors[5] = board[x+1][y+1];
    }
    else {
        if (y == 0)
            neighbors[0] = NULL;
        else
            neighbors[0] = board[x-1][y-1];
        neighbors[1] = board[x-1][y];
        if (y == 0)
            neighbors[2] = NULL;
        else
            neighbors[2] = board[x][y-1];
        if (y == w - 1)
            neighbors[3] = NULL;
        else
            neighbors[3] = board[x][y+1];
        if (x == h - 1 || y == 0)
            neighbors[4] = NULL;
        else
            neighbors[4] = board[x+1][y-1];
        if (x == h -1)
            neighbors[5] = NULL;
        else
            neighbors[5] = board[x+1][y];
    }

    return neighbors;
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
