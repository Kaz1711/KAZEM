# Simulation de la végétation d’un corridor sous ligne électrique

Ce projet modélise l'évolution de la végétation dans un corridor sous ligne électrique à l’aide d’un modèle de transitions entre différents états de végétation. 

Le corridor est représenté par des parcelles pouvant être dans quatre états :
- sol nu (Barren)
- herbes (Grasses)
- buisson type 1 (Shrub_1)
- buisson type 2 (Shrub_2)


## Objectifs de gestion du corridor

Le modèle cherche à vérifier si les conditions suivantes peuvent être atteintes à l'équilibre :
- 20 % des parcelles doivent être végétalisées et parmi elles : 
  - 30 % doivent être des herbes
  - 70 % doivent être des buissons
    - le type de buisson le moins abondant doit représenter au moins 30 % des buissons.

## Modèle

L’évolution du système est décrite par une matrice de transition entre les états de végétation.
Chaque ligne représente l'état actuel d'une parcelle et chaque colonne la probabilité de transition vers un autre état à la génération suivante.
Les probabilités sont normalisées pour que la somme de chaque ligne soitégale à 1.

Deux types de simulations sont utilisés :
- Déterministe : les transitions suivent directement les probabilités.
- Stochastique : les transitions sont tirées aléatoirement à l’aide d’une distribution multinomiale.

## Structure du code

Fonctions principales :
- `check_transition_matrix!` : vérifie que les probabilités de transition somment à 1
- `check_function_arguments` : vérifie la cohérence des dimensions de la matrice
- `_sim_stochastic!` : simulation stochastique avec distribution multinomiale
- `_sim_determ!` : simulation déterministe
- `simulation` : fonction principale qui exécute le modèle
- `check_success` : évalue si les objectifs écologiques sont atteints

## Limites du modèle

- il ne modélise pas l’influence de la densité de végétation
- il ignore les perturbations environnementales

## Dépendances

Le projet utilise Julia
Les packages suivants sont nécessaires :
- CairoMakie
- Distributions
- Random
