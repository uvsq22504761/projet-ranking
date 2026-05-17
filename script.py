import matplotlib.pyplot as plt
import subprocess
from pathlib import Path

repo_build = Path("build")
fichier_matrice = Path("wb-cs-stanford.mtx")
alphas = [0.7, 0.75, 0.8, 0.85, 0.9, 0.95]
epsilons = [1e-6, 1e-7, 1e-8, 1e-9, 1e-10]

subprocess.run(["make"])

for epsilon in epsilons:
    iterations = []
    for alpha in alphas:
        # on exécute pagerank avec les paramètres fichier et alpha, et les résultats sont stockés dans build/
        subprocess.run([str((repo_build / "pagerank").resolve()), str(fichier_matrice.resolve()), str(alpha), str(epsilon)], cwd=repo_build)
        fichier_resultat = repo_build / f"results_{fichier_matrice.stem}_alpha_{alpha:.2f}_eps_{epsilon:.0e}.txt"
        temps_tot = None; 
        with open(fichier_resultat) as f:
            for line in f:
                # on récupère le nombre d'itérations qui ira en ordonnée du graphique
                if line.startswith("# iterations"):
                    iterations.append(int(line.split(":")[1].strip()))
                if line.startswith("# temps"):
                    temps_tot = float(line.split(":")[1].strip())
                    

    # construction des graphiques, sauvegarde et affichage
    plt.plot(alphas, iterations)
    plt.xlabel("alpha")
    plt.ylabel("nombre d'itérations")
    plt.ylim(0, 180)
    plt.title(f"Pagerank Gauss-Seidel . epsilon = {epsilon:.0e}, temps total = {temps_tot:.4f}s")
    plt.grid(True)
    plt.savefig(repo_build / f"graphique_eps_{epsilon:.0e}.png") # aussi stocké dans build/
    plt.close()