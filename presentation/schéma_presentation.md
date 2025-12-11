# Schéma général de présentation

accroche : jeudi 4 décembre à 2 h 50 du matin a eu lieu un séisme de magnitude 2,5, le plus gros depuis la mise en exploitation du site, ce qui a conduit à la fermeture provisoire de la centrale de Rittershoffen

## Introduction

Transition : échec sur les corrélations donc on cherche une modélisation pour mieux comprendre le système + on observe des tendances que l'on veut retrouver

Description des objectifs : proposer un schéma de scénario d'injection appuyé qui s'appuie sur une modélisation simple mais capable de reproduire le comportement des séismes

Démarche : on part du rapport de l'INERIS

définition intuitive d'un séisme : avec claquement de doigt

## Modèle géomécanique

passage à la modélisation patin ressort : une loi d'élasticité mais aussi une loi de friction pour rendre compte de l'élasticité

modèle de friction loi de dieterich

Expliciter : c'est un modèle déjà utilisé pour étudier la sismicité naturelle que nous allons utiliser pour expliquer la sismicité induite, on doit donc ajouter une composante

## Modèle hydraulique

commencer par le cercle de Mohr

détailler les différents modèles et ce qu'ils apportent

## Implémentation

Architecture des codes

Etablissement du délai

Présentation rapide de RKF

## Discussion des résultats

étude paramétrique en fonction du retard

effet du moment d'injection et date du premier séisme après l'injection : ce n'est sûrement pas qu'une simple translation

discuter l'accélération des cycles et leur diminution d'amplitude

## Conclusion

proposer un scénario d'injection

identifier les paramètres pertinents de contrôle


Réponses de pierre :

l'asymétrie dans le second modèle s'explique soit par des asymétries introduites du fait de la discrétisation, soit de l'instabilité numériquee

ce second modèle permet d'expliciter deux échelles de temps, celle de la période de glissement asismique,  et celle du glissement sismique, le tout sur une seule et même faille

privilégier quelques fois des représentations en fonction du nombre d'itération pour repérer les zones de nucléations et de fractures.