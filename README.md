# reconstruire des tableaux 

L'objectif ici est de reconstruire un tableau préalablement traité par différentes méthodes de protection des données:

- la suppression de cases
- la perturbation aléatoire.

L'idée est d'estimer la probabilité qu'une case est en réalité sensible.

Avec certaines méthodes de perturbation aléatoire, il est possible de calculer une probabilité a posteriori à partir du mécanisme de perturbation.
Pour les tableaux auxquels sont appliquées des méthodes de blanchiment, cela nécessite de reconstruire l'ensemble des solutions possibles (pour des exemples de petite taille) 
et de calculer le nombre de solutions possibles dans lesquelles une case est confidentielle.  

