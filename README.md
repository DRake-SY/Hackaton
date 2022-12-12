# Hackaton reproductibilité

Ce projet a pour but une analyse RNAseq reproductible par la création d'un workflow Snakemake. Ce dernier fonctionne uniquement avec les données spécifiées dans le Snakemake

Pour utiliser ce workflow, il vous faudra une machine virtuelle BioPipes de 16 coeurs et 64 Go ram disponible sur le cloud IFB


## Utilisation du workflow

### Installation des dépendances

Après avoir déployé la VM, il faudra utiliser la commande 

```
conda activate snakemake
```
pour activer l'environnment conda de la VM contenant l'outil Snakemake puis utiliser la commande 

```
conda install singularity
```
afin d'utiliser la commande singulatity du Snakefile

### Lancement du workflow

Après avoir installé les dépendances nécéssaires au bon fonctionnement 

```
snakemake -s Snakefile --cores 2 --use-singularity --force all
```
