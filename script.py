import matplotlib.pyplot as plt
import subprocess
from pathlib import Path

repo_build = Path("build")
fichier_matrice = Path("wb-cs-stanford.mtx")
alphas = [0.7, 0.75, 0.8, 0.85, 0.9, 0.95]
iterations = []

subprocess.run(["make"])

for alpha in alphas:
    # on exécute pagerank avec les paramètres fichier et alpha, et les résultats sont stockés dans build/
    subprocess.run([str((repo_build / "pagerank").resolve()), str(fichier_matrice.resolve()), str(alpha)], cwd=repo_build)
    fichier_resultat = repo_build / f"results_{fichier_matrice.stem}_alpha_{alpha:.2f}.txt"
    with open(fichier_resultat) as f:
        for line in f:
            # on récupère le nombre d'itérations qui ira en ordonnée du graphique
            if line.startswith("# iterations"):
                iterations.append(int(line.split(":")[1].strip()))
                break

# construction du graphique, sauvegarde et affichage
plt.plot(alphas, iterations)
plt.xlabel("alpha")
plt.ylabel("nombre d'itérations")
plt.title("Pagerank Gauss-Seidel en fonction de alpha")
plt.grid(True)
plt.savefig(repo_build / "graphique.png") # aussi stocké dans build/
plt.show()