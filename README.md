Le fichier `pagerank_original.c` correspond au code original de PageRank étudié en TD. Une fois compilé, pour l'exécuter, il faut ajouter en argument le nom du fichier contenant la matrice sur laquelle on veut travailler au format MatrixMarket.

Le fichier `pagerank_projet.c` correspond au réel code du projet avec la version de Gauss-Seidel ascendant.

Pour compiler le projet, taper `make` dans le terminal. Cela créera un répertoire build qui contient l'exécutable. Pour nettoyer le répertoire temporaire `build/`, taper `make clean`.
Pour exécuter le projet, taper le nom de l'exécutable suivi du nom du fichier au format MatrixMarket contenant la matrice sur laquelle travailler, suivi possiblement de la valeur alpha (par défaut elle vaut 0.85 si rien n'est entré).

Nous avons le fichier `script.py` qui permet d'exécuter plusieurs fois pagerank sur une même matrice avec des valeurs d'alpha différentes. Pour le lancer, taper `python3 script.py` dans le terminal. Le graphique construit sera affiché, mais également sauvegardé dans le répertoire temporaire `build/`.