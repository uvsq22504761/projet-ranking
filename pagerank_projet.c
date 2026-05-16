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
    int *tos = malloc(M * sizeof(int));

    for (int k = 0; k < M; k++) {
        if (k % 1000000 == 0) { // affichage de progression
            fprintf(stderr, "\rArcs lus : %d / %d", k, M);
        }
        int i, j;
        if (fscanf(file, "%d %d", &i, &j) != 2) {
            fprintf(stderr, "Erreur lecture arc %d\n", k);
            exit(1);
        }
        froms[k] = i - 1;
        tos[k] = j - 1;
        deg[i - 1]++;
    }

    fprintf(stderr, "\rArcs lus : %d / %d\n", M, M);

    // construire f : f[i] = 1 si noeud sans arc sortant, 0 sinon
    *f_out = malloc(N * sizeof(double));
    double *f = *f_out;
    for (int i = 0; i < N; i++) {
        f[i] = (deg[i] == 0) ? 1.0 : 0.0;
    }

    // construire les listes chaînées avec probabilités uniformes
    fprintf(stderr, "Construction des listes chainees\n");
    for (int k = 0; k < M; k++) {
        int i = froms[k];
        int j = tos[k];
        double prob = 1.0 / deg[i];
        ajout_arc(P, prob, j, i); // arcs entrants
    }

    free(deg);
    free(froms);
    free(tos);
    fprintf(stderr, "Pret, debut des iterations.\n");
    fclose(file);
    return P;
}

double norme(double *x, double *y, int N) {
    double result = 0.0;
    for (int i = 0; i < N; i++) {
        result += fabs(x[i] - y[i]);
    }
    return result;
}

double *iterer(Cell **P, double *f, int N, double eps, int max_iter,
               double alpha, int *iterations) {

    double *pi     = (double *) malloc(N * sizeof(double));
    double *pi_old = (double *) malloc(N * sizeof(double));
    if (!pi || !pi_old) {
        fprintf(stderr, "Allocation de pi ratée\n");
        exit(1);
    }

    for (int i = 0; i < N; i++){
        pi[i] = 1.0 / N;
    }
    int iter = 0;
    double norm_val;

    do {
        // copie de l'état
        for (int i = 0; i < N; i++)
            pi_old[i] = pi[i];
        // ceci reste de la même façon
        double xf = 0.0;
        for (int i = 0; i < N; i++){
                xf += pi[i] * f[i];
        }
        double val = (alpha * xf + (1.0 - alpha)) / N;
        // parcours ascendant 
        for (int i = 0; i < N; i++) {
            // Changement
            // Un seul parcours de P[i] : s et Pii recuperees ensemble
            double s   = 0.0;
            double Pii = 0.0;
            Cell *c = P[i];
            while (c) {
                if (c->index == i) {
                    // arc de i vers lui-même 
                    Pii = c->val;
                } else {
                    // ça veut dire index diff 
                    s += pi[c->index] * c->val;
                }
                c = c->next;
            }
            s *= alpha;

            double Gii    = alpha * Pii + val;
            pi[i] = (s + val) / (1.0 - Gii);
        }

        // renormalisation — nécessaire car Gauss-Seidel ne préserve pas ||pi||=1
        double sum = 0.0;
        for (int i = 0; i < N; i++) sum += pi[i];
        for (int i = 0; i < N; i++) pi[i] /= sum;

        // Changement
        // norme calculée une seule fois et pour agiliser stocker swur norme val 
        // 
        norm_val = norme(pi, pi_old, N);
        fprintf(stderr, "Iteration %d, norme = %.2e\n", iter, norm_val);

        iter++;
        if (iter >= max_iter) {
            fprintf(stderr, "(max itérations atteintes)\n");
            break;
        }

    } while (norm_val > eps);   // réutilise norm_val déjà calculé

    *iterations = iter;
    free(pi_old);
    return pi;
}

// fonction pour enlever l'extension et le chemin d'un ficchier, gardant juste le nom de base
char *strip_file(char *filename) {
    char *base = filename;
    for (char *p = filename; *p; p++) { // dernier / (linux) ou \ (windows)
        if (*p == '/' || *p == '\\') {
            base = p + 1;
        }
    }
    char *end = base;
    char *dot = NULL;
    for (char *p = base; *p; p++) {
        if (*p == '.') {
            dot = p;
        }
        end = p;
    }
    size_t len = dot ? (size_t) (dot - base) : (size_t) (end - base + 1);
    char *result = malloc(len + 1);
    if (!result) {
        fprintf(stderr, "Malloc de nom de fichier stipped raté\n");
        exit(1);
    }
    for (size_t i = 0; i < len; i++) {
        result[i] = base[i];
    }
    result[len] = '\0';
    return result;
}

// sauvegarder dans un fichier txt pour le script python
void sauvegarder_txt(char *fichier_matrice, double *pi, int N, double alpha, double epsilon, int iterations) {
    char *fichier_matrice_stripped = strip_file(fichier_matrice);
    char filename[256];
    snprintf(filename, sizeof(filename), "results_%s_alpha_%.2f.txt", fichier_matrice_stripped, alpha);

    FILE *file = fopen(filename, "w");
    if (!file) {
        fprintf(stderr, "Fichier de résultat ne s'ouvre pas\n");
        return;
    }

    fprintf(file, "# Resultats de Pagerank Gauss-Seidel\n");
    fprintf(file, "# source : %s\n", fichier_matrice);
    fprintf(file, "# alpha : %.6f\n", alpha);
    fprintf(file, "# epsilon : %.2e\n", epsilon);
    fprintf(file, "# iterations : %d\n", iterations);
    fprintf(file, "# nombre de noeuds : %d\n", N);
    for (int i = 0; i < N; i++) {
        fprintf(file, "%d %.10f\n", i, pi[i]);
    }

    fclose(file);
    fprintf(stdout, "Les résultats sont sauvegardés dans le fichier %s\n", filename);
}

int main(int argv, char** args)  {
    // vérif nombre arguments
    if (argv < 2 || argv > 3) {
        fprintf(stderr, "Trop ou pas assez d'arguments : écrire suivi du nom du fichier et de la valeur d'alpha\n");
        exit(1);
    }

    // définir la valeur de alpha depuis l'input
    double alpha = ALPHA;
    if (argv == 3) {
        alpha = atof(args[2]);
        if (alpha <= 0.0 || alpha >= 1.0) {
            fprintf(stderr, "Pas la bonne valeur de alpha\n");
            exit(1);
        }
    }
    fprintf(stdout, "Alpha = %.6f\n", alpha);

    int N;
    double *f = NULL;

    // lecture du fichier contenant la matrice et construction de P
    Cell **P = lecture_matrice_market(args[1], &N, &f);

    // algorithme principal
    int iterations = 0;
    double *pi = iterer(P, f, N, EPSILON, MAX_ITER, alpha, &iterations);

    // résultats sauvegardés pour l'analyse
    sauvegarder_txt(args[1], pi, N, alpha, EPSILON, iterations);

    // libération mémoire
    for (int j = 0; j < N; j++) {
        Cell *c = P[j];
        while (c) {
            Cell *tmp = c;
            c = c->next;
            free(tmp);
        }
    }
    free(P);
    free(f);
    free(pi);

    return 0;
}