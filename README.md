Le fichier `pagerank_original.c` correspond au code original de PageRank étudié en TD. Une fois compilé, pour l'exécuter, il faut ajouter en argument le nom du fichier contenant la matrice sur laquelle on veut travailler au format MatrixMarket.

Le fichier `pagerank_projet.c` correspond au réel code du projet avec la version de Gauss-Seidel ascendant.

Pour compiler le projet, taper `make` dans le terminal. Cela créera un répertoire build qui contient l'exécutable. Pour nettoyer le répertoire temporaire `build/`, taper `make clean`.
Pour exécuter le projet, taper le nom de l'exécutable suivi du nom du fichier au format MatrixMarket contenant la matrice sur laquelle travailler, suivi possiblement de la valeur alpha (par défaut elle vaut 0.85 si rien n'est entré), suivie possiblement de la valeur de epsilon (autrement sa valeur par défaut est 1e-6). Il faut donc entrer : `./build/pagerank NOM_FICHIER.mtx [alpha] [epsilon]`.

Nous avons le fichier `script.py` qui permet d'exécuter plusieurs fois pagerank (version Gauss-Seidel) sur une même matrice avec des valeurs d'alpha différentes. Il est implémenté de sorte à exécuter ceci pour plusieurs valeurs de epsilon fixées. On obtient donc autant de graphiques que de valeurs de epsilon sur lesquelles on veut travailler. Pour le lancer, taper `python3 script.py` dans le terminal. Les graphiques construits seront sauvegardés dans le répertoire temporaire `build/`.

Notre fichier `script_2.py` est similaire, sauf qu'il permet d'exécuter pagerank version initiale.


Pour changer la matrice avec laquelle on travaille : depuis le script python il faut changer dans le code le chemin du fichier, depuis l'exécutable du C c'est dans les arguments passés à la suite de l'exécutable.