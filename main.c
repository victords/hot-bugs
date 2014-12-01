#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <pthread.h>

/* Used to be passed as parameter for the threads */
typedef struct {
    int bugsStart;
    int bugsEnd;
    int tilesStart;
    int tilesEnd;
} Limits;

/* Represents bugs */
typedef struct {
    int i;
    int j;
    void *nextTile; /* candidate to next tile where this bug will move */
                    /* type is 'void*' to avoid cross-reference */
} Bug;

/* Represents heat or cold sources */
typedef struct {
    int i;
    int j;
    double emission;
} Source;

/* Represents each tile in the board */
typedef struct {
    int i;
    int j;
    char element;       /* 0 if the tile is empty, 1 if there's a bug, 2 if there's a heat source or 3 if there's a cold source */
    double temperature; /* temperature on this tile */
    unsigned int seed;  /* seed for the heat/cold source on this tile */
    int lifetime;       /* remaining lifetime of the heat/cold source on this tile (if there's any) */
    Bug *candidates[6]; /* pointers for bugs that want to move to this tile */
    int candidateCount; /* number of used elements of the above array */
} Tile;

void *loop(void*);
Tile **getNeighbors(int, int);
double getDistance(int, int, int, int);
void moveBug(Bug*);
int getInteger(char*, char*);
double getDouble(char*, char*, double, double);
void printUsageAndExit(void);
void printBoard(int, int);
void printFinalState();

Tile ***board;
Bug **bugs;
Source **sources;
int w, h, n, s, nh, nc, t, np, tilesPerThread, tilesRemainder, bugsPerThread, bugsRemainder, sourceCount;
double c, tmin, tmax, ph, pc;
pthread_mutex_t mutex;
pthread_barrier_t barrier;

int main(int argc, char *argv[]) {
    int i, j, k, size, bugsStart, tilesStart;
    pthread_t *threads;
    Limits *p;
    
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
    bugsRemainder = n - (np * bugsPerThread);
    
    /* printf("%d %d %d %d %lf %lf %d %lf %d %d %d\n", w, h, n, s, c, ph, nh, pc, nc, t, np); */

    /* Constructing the board */
    board = malloc(h * sizeof(Tile **));
    for (i = 0; i < h; i++) {
        board[i] = malloc(w * sizeof(Tile *));
        for (j = 0; j < w; j++) {
            board[i][j] = malloc(sizeof(Tile));
            board[i][j]->i = i;
            board[i][j]->j = j;
            board[i][j]->element = 0;
            board[i][j]->temperature = 0;
            board[i][j]->candidateCount = 0;
        }
    }

    /* Generating bugs */
    srand(s);
    bugs = malloc(n * sizeof(Bug *));
    for (k = 0; k < n; k++) {
        i = rand() % h;
        j = rand() % w;
        if (board[i][j]->element) k--;
        else {
            board[i][j]->element = 1;
            bugs[k] = malloc(sizeof(Bug));
            bugs[k]->i = i;
            bugs[k]->j = j;
            bugs[k]->nextTile = NULL;
        }
    }
    
    /* Generating seeds */
    for (k = 0; k < size; k++) {
        i = k / w;
        j = k % w;
        board[i][j]->seed = s = ((i + 1) * s + j) % RAND_MAX;
    }
    
    /* Initializing heat/cold sources array */
    sources = malloc(size * sizeof(Source *));
    sourceCount = 0;

    /* Generating threads */
    threads = malloc(np * sizeof(pthread_t));
    pthread_mutex_init(&mutex, NULL);
    pthread_barrier_init(&barrier, NULL, np);
    bugsStart = tilesStart = 0;
    for (k = 0; k < np; k++) {
        p = malloc(sizeof(Limits));
        p->bugsStart = bugsStart;
        p->tilesStart = tilesStart;
        bugsStart += bugsPerThread + (k < bugsRemainder ? 1 : 0);
        tilesStart += tilesPerThread + (k < tilesRemainder ? 1 : 0);
        p->bugsEnd = bugsStart;
        p->tilesEnd = tilesStart;
        pthread_create(&threads[k], NULL, loop, p);
    }

    for (k = 0; k < np; k++)
        pthread_join(threads[k], NULL);

    // printBoard(w, h);
    printFinalState();

    /* Freeing memory */
    free(threads);

    for (k = 0; k < n; k++)
        free(bugs[k]);
    free(bugs);

    for (k = 0; k < sourceCount; k++)
        free(sources[k]);
    free(sources);

    for (i = 0; i < h; i++) {
        for (j = 0; j < w; j++)
            free(board[i][j]);
        free(board[i]);
    }
    free(board);
    
    return 0;
}

void *loop(void *args) {
    Limits p = *((Limits *)args);
    int count, i, j, k, l, m;
    double d, temperature, temp, temp2;
    char draw;
    Tile *next, **neighbors;
    Bug **candidates;
    free(args);
    //printf("starting thread (%d, %d)/(%d, %d)\n", p.tilesStart, p.tilesEnd, p.bugsStart, p.bugsEnd);
    for (count = 0; count < t; count++) {
        // printf("iteration %d, sourceCount %d\n", count, sourceCount);
        /* Heat and cold sources */
        for (k = p.tilesStart; k < p.tilesEnd; k++) {
            i = k / w;
            j = k % w;
            if (board[i][j]->element == 0) {
                int r = rand_r(&board[i][j]->seed);
                char hs = ((float)r / RAND_MAX) <= ph;
                r = rand_r(&board[i][j]->seed);
                char cs = ((float)r / RAND_MAX) <= pc;
                if (hs && !cs) {
                    board[i][j]->element = 2;
                    board[i][j]->lifetime = nh;
                    pthread_mutex_lock(&mutex);
                    sources[sourceCount] = malloc(sizeof(Source));
                    sources[sourceCount]->i = i;
                    sources[sourceCount]->j = j;
                    sources[sourceCount++]->emission = c;
                    pthread_mutex_unlock(&mutex);
                }
                if (cs && !hs) {
                    board[i][j]->element = 3;
                    board[i][j]->lifetime = nc;
                    pthread_mutex_lock(&mutex);
                    sources[sourceCount] = malloc(sizeof(Source));
                    sources[sourceCount]->i = i;
                    sources[sourceCount]->j = j;
                    sources[sourceCount++]->emission = -c;
                    pthread_mutex_unlock(&mutex);
                }
            } else if (board[i][j]->element > 1) {
                board[i][j]->lifetime--;
                if (board[i][j]->lifetime == 0) {
                    board[i][j]->element = 0;
                    pthread_mutex_lock(&mutex);
                    for (l = 0; l < sourceCount; l++) {
                        if (sources[l]->i == i && sources[l]->j == j) {
                            free(sources[l]);
                            sourceCount--;
                            for (; l < sourceCount; l++)
                                sources[l] = sources[l + 1];
                        }
                    }
                    pthread_mutex_unlock(&mutex);
                }
            }
        }
        pthread_barrier_wait(&barrier);
        
        /* Calculating the temperatures */
        for (k = p.tilesStart; k < p.tilesEnd; k++) {
            i = k / w;
            j = k % w;
            board[i][j]->temperature = 0;
            board[i][j]->candidateCount = 0;
            for (l = 0; l < n; l++) {
                if (i == bugs[l]->i && j == bugs[l]->j) continue;
                d = getDistance(i, j, bugs[l]->i, bugs[l]->j);
                board[i][j]->temperature += c / (d * d);
            }
            for (l = 0; l < sourceCount; l++) {
                if (i == sources[l]->i && j == sources[l]->j) continue;
                d = getDistance(i, j, sources[l]->i, sources[l]->j);
                board[i][j]->temperature += sources[l]->emission / (d * d);
            }
        }
        pthread_barrier_wait(&barrier);

        /* Computing bugs' candidates for movement */
        for (k = p.bugsStart; k < p.bugsEnd; k++) {
            temperature = board[bugs[k]->i][bugs[k]->j]->temperature;
            if (temperature > tmax || temperature < tmin) {
                neighbors = getNeighbors(bugs[k]->i, bugs[k]->j);
                for (l = 0; l < 6; l++) {
                    if (neighbors[l]) {
                        temp = neighbors[l]->temperature;
                        m = l;
                        break;
                    }
                }
                if (temperature > tmax) {
                    for (l = m; l < 6; l++)
                        if (neighbors[l] && neighbors[l]->temperature < temp) {
                            temp = neighbors[l]->temperature;
                            m = l;
                        }
                } else {
                    for (l = m; l < 6; l++)
                        if (neighbors[l] && neighbors[l]->temperature > temp) {
                            temp = neighbors[l]->temperature;
                            m = l;
                        }
                }
                if (neighbors[m]->element == 0) {
                    bugs[k]->nextTile = neighbors[m];
                    pthread_mutex_lock(&mutex);
                    neighbors[m]->candidates[neighbors[m]->candidateCount++] = bugs[k];
                    pthread_mutex_unlock(&mutex);
                }
            }
        }
        pthread_barrier_wait(&barrier);
        
        /* Computing movements */
        for (k = p.bugsStart; k < p.bugsEnd; k++) {
            temperature = board[bugs[k]->i][bugs[k]->j]->temperature;
            next = (Tile *)bugs[k]->nextTile;
            if (next) {
                if (next->candidateCount == 1)
                    moveBug(bugs[k]);
                else {
                    candidates = next->candidates;
                    temp = fabs(next->temperature - board[candidates[0]->i][candidates[0]->j]->temperature);
                    m = 0;
                    draw = 0;
                    for (l = 1; l < next->candidateCount; l++) {
                        temp2 = fabs(next->temperature - board[candidates[l]->i][candidates[l]->j]->temperature);
                        if (temp2 > temp) {
                            temp = temp2;
                            m = l;
                            draw = 0;
                        }
                        else if (temp2 == temp) draw = 1;
                    }
                    if (draw) m = -1;
                    for (l = 0; l < next->candidateCount; l++) {
                        if (l == m)
                            moveBug(candidates[l]);
                        else 
                            candidates[l]->nextTile = NULL;
                    }
                }
            }
        }
        pthread_barrier_wait(&barrier);
    }
    return NULL;
}

Tile **getNeighbors(int i, int j) {
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

void moveBug(Bug *bug) {
    Tile *next;
    pthread_mutex_lock(&mutex);
    next = (Tile *)bug->nextTile;
    if (next) {
        if (next->element == 0) {
            board[bug->i][bug->j]->element = 0;
            bug->i = next->i;
            bug->j = next->j;
            next->element = 1;
        }
        bug->nextTile = NULL;
    }
    pthread_mutex_unlock(&mutex);
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
            printf("%hhd|", board[i][j]->element);
        printf("\n");
    }
}

void printFinalState() {
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
		fprintf(saida, "Joaninha nº %5d | ", i);
		fprintf(saida, "Posição: (%4u, %4u) | ", bugs[i]->i, bugs[i]->j);
		fprintf(saida, "Temperatura: %10.3lf\n", board[bugs[i]->i][bugs[i]->j]->temperature);
	}

	/* Close file */
	fclose(saida);
}
