#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <pthread.h>

/* Represents coordinates in the board */
typedef struct {
    unsigned int i;
    unsigned int j;
} Point;

/* Represents heat or cold sources */
typedef struct {
    unsigned int i;
    unsigned int j;
    double emission;
} Source;

/* Represents each tile in the board */
typedef struct {
    char element;       /* 0 if the tile is empty, 1 if there's a bug, 2 if there's a heat source or 3 if there's a cold source */
    double temperature; /* temperature on this tile */
    unsigned int seed;  /* seed for the heat/cold source on this tile */
} Tile;

void *loop(void*);
Tile **getNeighbors(Tile***, int, int, int, int);
int getInteger(char*, char*);
double getDouble(char*, char*, double, double);
void printUsageAndExit(void);
void printBoard(int, int);
void printFinalState(Point*);

Tile ***board;
Point *bugs;
Source *sources;
int w, h, n, s, nh, nc, t, np, sourceCount, tilesPerThread, tilesRemainder, bugsPerThread, bugsRemainder;
double c, tmin, tmax, ph, pc;
pthread_mutex_t mutex;
pthread_barrier_t barrier;

int main(int argc, char *argv[]) {
    int i, j, k, size, bugsStartPoint, tilesStartPoint;
    pthread_t *threads;
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
    
    size = w * h;
    if (n >= size) {
        printf("The number of bugs must be less than the size of the board (W * H).\n");
        printUsageAndExit();
    }
    if (tmin >= tmax) {
        printf("The value of TMIN must be strictly less than TMAX.\n");
        printUsageAndExit();
    }

    tilesPerThread = size / np;
    tilesRemainder = size - (np * tilesPerThread);
    bugsPerThread = n / np;
    bugsRemainder = size - (np * bugsPerThread);
    
    /* printf("%d %d %d %d %lf %lf %d %lf %d %d %d\n", w, h, n, s, c, ph, nh, pc, nc, t, np); */

    /* Constructing the board */
    board = malloc(h * sizeof(Tile **));
    for (i = 0; i < h; i++) {
        board[i] = malloc(w * sizeof(Tile *));
        for (j = 0; j < w; j++) {
            board[i][j] = malloc(sizeof(Tile));
            board[i][j]->element = 0;
        }
    }
    
    /* Generating bugs */
    srand(s);
    bugs = malloc(n * sizeof(Point));
    for (k = 0; k < n; k++) {
        i = rand() % h;
        j = rand() % w;
        if (board[i][j]->element) k--;
        else {
            board[i][j]->element = 1;
            bugs[k].i = i;
            bugs[k].j = j;
        }
    }

    /* Generating seeds and heat/cold sources array */
    sources = malloc(size * sizeof(Source));
    sourceCount = 0;
    for (k = 0; k < size; k++) {
        i = k / w;
        j = k % w;
        board[i][j]->seed = s = ((i + 1) * s + j) % RAND_MAX;
    }
    
    threads = malloc(np * sizeof(pthread_t));
    pthread_mutex_init(&mutex, NULL);
    pthread_barrier_init(&barrier, NULL, np);
    bugsStartPoint = tilesStartPoint = 0;
    for (k = 0; k < np; k++) {
        printf("iniciando thread %d %d\n", bugsStartPoint, tilesStartPoint);
        /* Using a Point just to store the two start points which must be passed as parameters to the thread function */
        p = malloc(sizeof(Point));
        p->i = bugsStartPoint;
        p->j = tilesStartPoint;
        pthread_create(&threads[k], NULL, loop, p);
        bugsStartPoint += bugsPerThread + (k < bugsRemainder ? 1 : 0);
        tilesStartPoint += tilesPerThread + (k < tilesRemainder ? 1 : 0);
    }

    for (k = 0; k < np; k++) {
        printf("tentando join...\n");
        pthread_join(threads[k], NULL);
    }
printf("joinou\n");
    printBoard(w, h);
    printFinalState(bugs);

    /* Freeing memory */
    for (k = 0; k < h; k++)
        free(board[k]);
    free(board);
    
    return 0;
}

void *loop(void *args) {
    printf("entrou\n");
    Point p = *((Point *)args);
    int bugsStartPoint = p.i,
        bugsId = bugsStartPoint / bugsPerThread,
        bugsEndPoint = bugsStartPoint + (bugsPerThread + (bugsId < bugsRemainder ? 1 : 0)),
        tilesStartPoint = p.j,
        tilesId = tilesStartPoint / tilesPerThread,
        tilesEndPoint = tilesStartPoint + (tilesPerThread + (tilesId < tilesRemainder ? 1 : 0)),
        count, lifetime, i, j, k;
    // double temperature;
    free(args);
    for (count = 0; count < t; count++) {
        printf("iteracao: %d\n", count);
        for (k = tilesStartPoint; k < tilesEndPoint; k++) {
            i = k / w;
            j = k % w;
            printf("(%d, %d)\n", i, j);
            /* Heat and cold sources */
            if (board[i][j]->element == 0) {
                int r = rand_r(&board[i][j]->seed);
                char hs = ((float)r / RAND_MAX) <= ph;
                r = rand_r(&board[i][j]->seed);
                char cs = ((float)r / RAND_MAX) <= pc;
                if (hs && !cs) {
                    board[i][j]->element = 2;
                    pthread_mutex_lock(&mutex);
                    sources[sourceCount].i = i;
                    sources[sourceCount].j = j;
                    sources[sourceCount++].emission = c;
                    pthread_mutex_unlock(&mutex);
                    lifetime = nh;
                }
                if (cs && !hs) {
                    board[i][j]->element = 3;
                    pthread_mutex_lock(&mutex);
                    sources[sourceCount].i = i;
                    sources[sourceCount].j = j;
                    sources[sourceCount++].emission = -c;
                    pthread_mutex_unlock(&mutex);
                    lifetime = nc;
                }
            } else {
                lifetime--;
                if (lifetime == 0)
                    board[i][j]->element = 0;
            }
            pthread_barrier_wait(&barrier);
            
            /* Calculating the temperature */
        }
    }
    return NULL;
}

Tile **getNeighbors(Tile ***board, int i, int j, int w, int h) {
    Tile **neighbors = malloc(6 * sizeof(Tile *));

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

double getDistance(int i1, int j1, int i2, int j2) {
	double dx, dy;
	if ((i2 - i1) % 2 == 0)
		dx = abs(j2 - j1);
	else if (j1 <= j2) {
		if (i1 % 2 == 0) dx = fabs(j2 - j1 - 0.5);
		else dx = j2 - j1 + 0.5;
	} else {
		if (i2 % 2 == 0) dx = fabs(j1 - j2 - 0.5);
		else dx = j1 - j2 + 0.5;
	}
	dy = abs(i2 - i1) * sqrt(3) * 0.5;
	return sqrt(dx * dx + dy * dy);
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
            printf("%hhd %10u|", board[i][j]->element, board[i][j]->seed);
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
