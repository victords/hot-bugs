#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <float.h>
#include <pthread.h>

/* Represents coordinates in the board */
typedef struct {
    unsigned int i;
    unsigned int j;
} Point;

/* Represents each tile in the board */
typedef struct Tile *tile;
struct Tile {
    char bug;           /* 1 if there's a bug on this tile, 0 otherwise */
    double temperature; /* temperature on this tile */
    double emission;    /* positive when there's a heat source on this tile,
                           negative when there's a "cold source", and zero
                           otherwise. */
    unsigned int seed;  /* seed for the heat/cold source on this tile */
};

void *tileLoop(void*);
void *bugLoop(void*);
tile *getNeighbors(tile**, int, int, int, int);
int getInteger(char*, char*);
double getDouble(char*, char*, double, double);
void printUsageAndExit(void);
void printBoard(int, int);
void printFinalState(Point*);

tile **board;
Point *bugs;
int w, h, n, s, nh, nc, t, np;
double c, tmin, tmax, ph, pc;
pthread_barrier_t tileBarrier, bugBarrier;

int main(int argc, char *argv[]) {
    int i, j, ii, jj;
    pthread_t *tileThreads, *bugThreads;
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
    
    if (n >= w * h) {
        printf("The number of bugs must be less than the size of the board (W * H).\n");
        printUsageAndExit();
    }
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
    bugs = malloc(n * sizeof(Point));
    bugThreads = malloc(n * sizeof(pthread_t));
    pthread_barrier_init(&bugBarrier, NULL, w * h);
    for (i = 0; i < n; i++) {
        ii = rand() % h;
        jj = rand() % w;
        if (board[ii][jj]->bug) i--;
        else {
            board[ii][jj]->bug = 1;
            bugs[i].i = ii;
            bugs[i].j = jj;
            pthread_create(&bugThreads[i], NULL, bugLoop, NULL);
        }
    }

    /* Generating seeds and threads for each tile */
    tileThreads = malloc(w * h * sizeof(pthread_t));
    pthread_barrier_init(&tileBarrier, NULL, w * h);
    for (i = 0; i < w * h; i++) {
        ii = i / w;
        jj = i % w;
        board[ii][jj]->seed = s = ((ii + 1) * s + jj) % RAND_MAX;
        p = malloc(sizeof(Point));
        p->i = ii;
        p->j = jj;
        pthread_create(&tileThreads[i], NULL, tileLoop, p);
    }

    for (i = 0; i < w * h; i++)
        pthread_join(tileThreads[i], NULL);
    for (i = 0; i < n; i++)
        pthread_join(bugThreads[i], NULL);

    printBoard(w, h);
    printFinalState(bugs);

    /* Freeing memory */
    for (i = 0; i < h; i++)
        free(board[i]);
    free(board);
    
    return 0;
}

/* Loop executed by each tile until the simulation ends */
void *tileLoop(void *args) {
    Point p = *((Point*)args);
    int count, lifetime, i;
    double temperature;
    for (count = 0; count < t; count++) {
        /* Heat and cold sources */
        if (board[p.i][p.j]->emission == 0) {
            int r = rand_r(&board[p.i][p.j]->seed);
            char hs = ((float)r / RAND_MAX) <= ph;
            r = rand_r(&board[p.i][p.j]->seed);
            char cs = ((float)r / RAND_MAX) <= pc;
            if (hs && !cs) {
                board[p.i][p.j]->emission = c;
                lifetime = nh;
            }
            if (cs && !hs) {
                board[p.i][p.j]->emission = -c;
                lifetime = nc;
            }
        } else {
            lifetime--;
            if (lifetime == 0)
                board[p.i][p.j]->emission = 0;
        }
        pthread_barrier_wait(&tileBarrier);
        
        /* Calculating the temperature */
    }
    return NULL;
}

void *bugLoop(void *args) {
    int count = 0;
    for (count = 0; count < t; count++) {
        /* fazer algo */
    }
    return NULL;
}

tile *getNeighbors(tile **board, int i, int j, int w, int h) {
    tile *neighbors = malloc(6 * sizeof(tile));

    if(i % 2 == 0) {
        if (i == 0)
            neighbors[0] = NULL;
        else
            neighbors[0] = board[i-1][j];
        if (i == 0 || j == w - 1)
            neighbors[1] = NULL;
        else    
            neighbors[1] = board[i-1][j+1];
        if (j == 0)
            neighbors[2] = NULL;
        else
            neighbors[2] = board[i][j-1];
        if (j == w - 1)
            neighbors[3] = NULL;
        else
            neighbors[3] = board[i][j+1];
        if (i == h - 1)
            neighbors[4] = NULL;
        else
            neighbors[4] = board[i+1][j];
        if (i == h - 1 || j == w - 1)
            neighbors[5] = NULL;
        else
            neighbors[5] = board[i+1][j+1];
    }
    else {
        if (j == 0)
            neighbors[0] = NULL;
        else
            neighbors[0] = board[i-1][j-1];
        neighbors[1] = board[i-1][j];
        if (j == 0)
            neighbors[2] = NULL;
        else
            neighbors[2] = board[i][j-1];
        if (j == w - 1)
            neighbors[3] = NULL;
        else
            neighbors[3] = board[i][j+1];
        if (i == h - 1 || j == 0)
            neighbors[4] = NULL;
        else
            neighbors[4] = board[i+1][j-1];
        if (i == h -1)
            neighbors[5] = NULL;
        else
            neighbors[5] = board[i+1][j];
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

double getDouble(char *arg, char *name, double lowerLimit, double upperLimit) {
    double x = atof(arg);
    if (x == 0.0) {
        printf("The value of %s must be a non-zero number.\n\n", name);
        printUsageAndExit();
    }
    if (x < lowerLimit) {
        printf("The value of %s must be greater than %.1lf.\n\n", name, lowerLimit);
        printUsageAndExit();
    }
    if (x > upperLimit) {
        printf("The value of %s must be less than %.1lf.\n\n", name, upperLimit);
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

void printFinalState(Point *bugs) {
	int i;
	FILE* saida;

	/* Open file */
	saida = fopen("hotbugs.txt", "w");
	if(saida == NULL) {
		fprintf(stderr, "Erro! Arquivo de saída não pode ser criado.\n");
		exit(EXIT_FAILURE);
	}

	/* Print inside the file */
	for(i = 0; i < n; i++) {
		fprintf(saida, "Joaninha nº %d\n", i);
		fprintf(saida, "Posição (i, j): (%u, %u)\n", bugs[i].i, bugs[i].j);
		fprintf(saida, "Temperatura: %lf\n", board[bugs[i].i][bugs[i].j]->temperature);
		if (i < n - 1)
			fprintf(saida, "========================================\n");
	}

	/* Close file */
	fclose(saida);
}
