#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <float.h>
#include <pthread.h>

/* Represents coordinates in the board */
typedef struct {
    unsigned int x;
    unsigned int y;
} Point;

/* Represents each tile in the board */
typedef struct Tile *tile;
struct Tile {
    char bug; /* 1 if there's a bug on this tile, 0 otherwise */
    double emission; /* positive when there's a heat source on this tile,
                        negative when there's a "cold source", and zero
                        otherwise. */
    unsigned int seed; /* seed for the heat/cold source on this tile */
};

void *tileLoop(void*);
tile *getNeighbors(tile**, int, int, int, int);
int getInteger(char*, char*);
double getDouble(char*, char*, double, double);
void printUsageAndExit(void);
void printBoard(int, int);

tile **board;
int w, h, n, s, nh, nc, t, np;
double c, tmin, tmax, ph, pc;
pthread_barrier_t tileBarrier, bugBarrier;

int main(int argc, char *argv[]) {
    int i, j, x, y;
    pthread_t *tileThreads;
    Point *p;
    
    /* Must receive 11 arguments */
    if (argc != 14)
        printUsageAndExit();
    
    /* Getting argument values */
    w = getInteger(argv[1], "W");
    h = getInteger(argv[2], "H");
    n = getInteger(argv[3], "N");
    s = getInteger(argv[4], "S");
    c = getDouble(argv[5], "C", 0.0, DBL_MAX);
    tmin = getDouble(argv[6], "TMIN", -DBL_MAX, DBL_MAX);
    tmax = getDouble(argv[7], "TMAX", -DBL_MAX, DBL_MAX);
    ph = getDouble(argv[8], "PH", 0.0, 1.0);
    nh = getInteger(argv[9], "NH");
    pc = getDouble(argv[10], "PC", 0.0, 1.0);
    nc = getInteger(argv[11], "NC");
    t = getInteger(argv[12], "T");
    np = getInteger(argv[13], "NP");
    
    if (tmin >= tmax) {
        printf("The value of TMIN must be strictly less than TMAX.\n");
        printUsageAndExit();
    }
    
    printf("%d %d %d %d %lf %lf %d %lf %d %d %d\n", w, h, n, s, c, ph, nh, pc, nc, t, np);

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

    /* Generating seeds and threads for each tile */
    tileThreads = malloc(w * h * sizeof(pthread_t));
    pthread_barrier_init(&tileBarrier, NULL, w * h);
    for (i = 0; i < w * h; i++) {
        x = i % w;
        y = i / w;
        board[y][x]->seed = s = ((y + 1) * s + x) % RAND_MAX;
        p = malloc(sizeof(Point));
        p->x = x;
        p->y = y;
        pthread_create(&tileThreads[i], NULL, tileLoop, p);
    }

    for (i = 0; i < w * h; i++)
        pthread_join(tileThreads[i], NULL);

    printBoard(w, h);

    /* Freeing memory */
    for (i = 0; i < h; i++)
        free(board[i]);
    free(board);
    
    return 0;
}

/* Loop executed by each tile until the simulation ends */
void *tileLoop(void *args) {
    Point p = *((Point*)args);
    int count, lifetime;
    for (count = 0; count < t; count++) {
        /* Heat and cold sources */
        if (board[p.y][p.x]->emission == 0) {
            int r = rand_r(&board[p.y][p.x]->seed);
            char hs = ((float)r / RAND_MAX) <= ph;
            r = rand_r(&board[p.y][p.x]->seed);
            char cs = ((float)r / RAND_MAX) <= pc;
            if (hs && !cs) {
                board[p.y][p.x]->emission = c;
                lifetime = nh;
            }
            if (cs && !hs) {
                board[p.y][p.x]->emission = -c;
                lifetime = nc;
            }
        } else {
            lifetime--;
            if (lifetime == 0)
                board[p.y][p.x]->emission = 0;
        }
        pthread_barrier_wait(&tileBarrier);
        
        /* Calculating the temperature */
        
    }
    return NULL;
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

double getDouble(char *arg, char *name, double lower_limit, double upper_limit) {
    double x = atof(arg);
    if (x == 0.0) {
        printf("The value of %s must be a non-zero number.\n\n", name);
        printUsageAndExit();
    }
    if (x < lower_limit) {
        printf("The value of %s must be greater than %.1lf.\n\n", name, lower_limit);
        printUsageAndExit();
    }
    if (x > upper_limit) {
        printf("The value of %s must be less than %.1lf.\n\n", name, upper_limit);
        printUsageAndExit();
    }
    return x;
}

void printUsageAndExit() {
    printf("Usage: hot-bugs W H N S C TMIN TMAX PH NH PC NC T NP\n=====\n\
W\twidth of the board\n\
H\theight of the board\n\
N\tnumber of bugs\n\
S\tseed for the random number generator\n\
C\tconstant of heat emission by the bugs\n\
TMIN\ttemperature below which the bugs will look for heat\n\
TMAX\ttemperature above which the bugs will look for cold\n\
PH\tprobability of heat source appearance\n\
NH\tnumber of iterations the heat sources will last\n\
PC\tprobability of \"cold source\" appearance\n\
NC\tnumber of iterations the \"cold sources\" will last\n\
T\ttotal iterations for the simulation\n\
NP\tnumber of processors to be used\n");
    exit(1);
}

void printBoard(int w, int h) {
    int i, j;
    for (i = 0; i < h; i++) {
        for (j = 0; j < w; j++)
            printf("%hhd %10lf %10u|", board[i][j]->bug, board[i][j]->emission, board[i][j]->seed);
        printf("\n");
    }
}
