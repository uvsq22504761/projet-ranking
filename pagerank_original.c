#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define EPSILON 1e-6
#define MAX_ITER 2000
#define ALPHA 0.85

typedef struct Cell {
    int index;
    double val;
    struct Cell *next;
} Cell;

void ajout_arc(Cell **P, double prob, int i, int j) {
    Cell *new_c;
    new_c = (Cell *) malloc(sizeof(Cell));
    if (!new_c) {
        fprintf(stderr, "Malloc de cellule raté\n");
        return;
    }
    new_c->index = j;
    new_c->val = prob;
    new_c->next = P[i];
    P[i] = new_c;
}

// lecture de la matrice en format MatrixMarket, stockage sous forme de tableau de listes chainées
Cell **lecture_matrice_market(char *filename, int *N_out, double **f_out) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Erreur lecture fichier\n");
        exit(1);
    }

    // ignorer les lignes de commentaires commençant par %
    char line[1024];
    while (fgets(line, sizeof(line), file)) {
        if (line[0] != '%') break;
    }

    // lire N, N, M
    int N, M, tmp;
    if (sscanf(line, "%d %d %d", &N, &tmp, &M) != 3) {
        fprintf(stderr, "Erreur lecture en-tête\n");
        exit(1);
    }
    *N_out = N;
    fprintf(stderr, "Lecture de %d noeuds, %d arcs...\n", N, M);

    // initialiser P (tableau de listes chaînées)
    Cell **P = calloc(N, sizeof(Cell *));

    // compter le degré sortant de chaque noeud
    int *deg = calloc(N, sizeof(int));
    int *froms = malloc(M * sizeof(int));
    int *tos   = malloc(M * sizeof(int));

    for (int k = 0; k < M; k++) {
        if (k % 1000000 == 0) {
            fprintf(stderr, "\rArcs lus : %d / %d", k, M);
        }
        int i, j;
        if (fscanf(file, "%d %d", &i, &j) != 2) {
            fprintf(stderr, "Erreur lecture arc %d\n", k);
            exit(1);
        }
        froms[k] = i - 1; // passage en 0-indexé
        tos[k] = j - 1;
        deg[i - 1]++;
    }

    fprintf(stderr, "\rArcs lus : %d / %d\n", M, M);

    // construire f : f[i] = 1 si noeud sans arc sortant (dangling), 0 sinon
    *f_out = malloc(N * sizeof(double));
    double *f = *f_out;
    for (int i = 0; i < N; i++) {
        f[i] = (deg[i] == 0) ? 1.0 : 0.0;
    }

    // construire les listes chaînées avec probabilités uniformes
    fprintf(stderr, "Construction des listes chainees...\n");
    for (int k = 0; k < M; k++) {
        int i = froms[k];
        int j = tos[k];
        double prob = 1.0 / deg[i];
        ajout_arc(P, prob, i, j);
    }

    free(deg);
    free(froms);
    free(tos);
    fprintf(stderr, "Pret, debut des iterations.\n");
    fclose(file);
    return P;
}

void multiplier_axP(double *pi_resultat, double *x, Cell **P, int N, double alpha) {
    for (int j = 0; j < N; j++) {
        if (x[j] == 0.0) continue;
        Cell *c = P[j];
        while (c) {
            pi_resultat[c->index] += x[j] * c->val * alpha;
            c = c->next;
        }
    }
}

double multiplier_xf(double *x, double *f, int N) {
    double s = 0.0;
    for (int i = 0; i < N; i++) {
        s += (x[i] * f[i]);
    }
    return s;
}

void multiplier_e(double *e, double *x, double *f, int N, double alpha) {
    double val = (alpha * multiplier_xf(x, f, N) + (1.0 - alpha)) / N;
    for (int i = 0; i < N; i++) {
        e[i] = val;
    }
}

void multiplier(double *y, double *x, Cell **P, double *f, int N, double alpha) {
    for (int i = 0; i < N; i++) {
        y[i] = 0.0;
    }
    multiplier_axP(y, x, P, N, alpha);
    double *e;
    e = (double *) malloc(sizeof(double) * N);
    multiplier_e(e, x, f, N, alpha);
    for (int i = 0; i < N; i++) {
        y[i] = y[i] + e[i];
    }
    free(e);
}

double norme(double *pi_pair, double *pi_impair, int N) {
    double result = 0;
    for (int i = 0; i < N; i++) {
        result += fabs(pi_pair[i] - pi_impair[i]);
    }
    return result;
}

void iterer(Cell **P, double *f, int N, double eps, int max_iter, double alpha) {
    double *pi_pair;
    double *pi_impair;
    pi_pair = (double *) malloc(N * sizeof(double));
    pi_impair = (double *) malloc(N * sizeof(double));
    if (!pi_pair || !pi_impair) {
        fprintf(stderr, "Allocation de pi ratée\n");
        exit(1);
    }

    // init pi
    for (int i = 0; i < N; i++) {
        pi_pair[i] = 1.0 / N;
    }

    int iter = 0;

    do {
        multiplier(pi_impair, pi_pair, P, f, N, alpha);
        multiplier(pi_pair, pi_impair, P, f, N, alpha);

        fprintf(stderr, "Iteration %d, norme = %.2e\n", iter, norme(pi_pair, pi_impair, N));

        iter++;
        if (iter >= max_iter) {
            printf("(max itérations atteintes)\n");
            break;
        }
    } while (norme(pi_pair, pi_impair, N) > eps);

    // top 10 pour vérifier
    int top[10] = {0};
    for (int i = 1; i < N; i++) {
        for (int k = 0; k < 10; k++) {
            if (pi_pair[i] > pi_pair[top[k]]) {
                for (int l = 9; l > k; l--) top[l] = top[l-1];
                top[k] = i;
                break;
            }
        }
    }
    fprintf(stderr, "Top 10 noeuds :\n");
    double sum = 0;
    for (int i = 0; i < N; i++) sum += pi_pair[i];
    for (int k = 0; k < 10; k++)
        fprintf(stderr, "Node %d : %.10f\n", top[k], pi_pair[top[k]]);
    fprintf(stderr, "Somme = %.6f\n", sum);

    free(pi_pair);
    free(pi_impair);
}

int main(int argv, char** args)  {
    if (argv == 1 || argv > 2) {
        fprintf(stderr, "Trop ou pas assez d'arguments : écrire suivi du nom du fichier\n");
        exit(1);
    }
    int N;
    double *f = NULL;

    Cell **P = lecture_matrice_market(args[1], &N, &f);

    iterer(P, f, N, EPSILON, MAX_ITER, ALPHA);

    for (int j = 0; j < N; j++) {
        Cell *c = P[j];
        while (c) {
            Cell *tmp = c;
            c = c->next;
            free(tmp);
        }
    }
    free(P);

    return 0;
}