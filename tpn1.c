#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#include <pthread.h>

#define SIZE_READ      100

/*
*On suppose que la taille maximale de la mémoire avec laquelle
*nous pouvons travailler est de 500Mo = 500 millions d'octets
*Ce qui nous amenera a appliquée une méthode de gestion
*optimisée de la mémoire, pour la réalisation de ce TP
*/
#define MAX_MEMORY_SIZE 500000000

/*
*Taille maximale que peut avoir un tableau de reads en tenant
*compte de MAX_MEMORY_SIZE
*/
#define MAX_READS_TABLE_SIZE MAX_MEMORY_SIZE / 100

/*
*Cette fonction a pour role de lire nombreDeGenomesALire lectures (GENOME)
*@Param tableau représente le tableau dans lequel les génomes lus sont stockés
*@Param nombreDeGenomesALire représente le nombre de génomes à lire
*@Param fichier représente le fichier à lire pour recuperer les lectures
*retourne le nombre de lecture lues
*/
int lectureDeGenomes(char** tableau, int nombreDeGenomesALire, FILE * fichier);

/* fonction utilisateur de comparaison fournie a qsort() */
static int compare (void const *a, void const *b);

/*
*Cette fonction a pour role de trier un tableau donné
*/
void trierTableauDeLectures(char** tableau, int taille);

/*
*Après avoir lus une quantité données de lectures à partir du fichier fasta.reads,
*et après les avoir triées dans un tableau,
*on stocke le contenu du tableau trié dans un nouveau fichier pour libérer la mémoire
*/
void ecrireTableauDeLecturesDansUnFichier(char** tableau, int taille, int* nombreDeFichiers);

/*
*on ajoute le contenu du tableaux resultats (resultas du tri-fusion) dans le fichierResultats
*on vide le tableau resultats
*/
void ecrireTableauResultatsDansUnFichier(char** tableau, int taille, FILE* fichierResultats);

/*
*Ici on fait le tri fusion des fichiers de lectures triées précédemment créés.
*On lui fourni le nombre de fichiers créés et le nombre maximal de lectures que peut contenir ce fichier
*/
void triFusion(int nombreDeFichiers, int nombreMaxDeLectures);

/*
*On verifie ici si le tri-fusion est terminé.
*Elle retourne 1 (faux) si l'un des tableaux contient encore des lectures  (indices[indiceDuTableau] != 0)
*Elle retourne 0 si tous les tableaux ne contiennent plus de lectures  (indices[indiceDuTableau] == 0)
*/
int triFusionTermine(int * indices, int nombreDeFichiers);

/*
*Un nombre nombreDeFichiers de tableau de tableaux de lectures etant créés dans triFusion(),
*On parcours ces differents tableaux contenus dans le atbleau lectures pour
*retourner l'indice du tableau (indice dans lectures) ayant la plus petite des lecture
*/
int indiceMinTableauLectures(char*** lectures, int nombreDeFichiers, int* indices);

/*
*On verifie pour chaque tableau de lectures si toutes ces lectures ont été parcourues
*si indices[indiceDuTableau] est egale à 1 alors toutes les lectures du tableau ont été parcourues
*Alors on remplace les lectures de ce tableau pas des chaines de caractères vides
*Puis on remplace le maximum possible par de nouvelles lectures à partir du fichier dont les précédentes avaient
*été extraites.
*/
void restaurationDeLectures(char*** lectures, int nombreDeFichiers, int nombreMaxDeLectures, int* indices, FILE** files);

int main (int argc, char *argv[]) {
    FILE * fichierDesLectures = NULL;
    char** tableauDeLectures;
    /*
    *Chaque fichier lu represente jusqu'a MAX_READS_TABLE_SIZE genomes lus dans le fichier de depart
    */
    int nombreDeFichiers = 0;

    /*
    *On ouvre le fichier des lectures en mode lecture
    */
    fichierDesLectures = fopen("reads.fasta", "r");
    
    int resultatDeLaLectureDeGenomes = 0;
    int taille = 0;
    do {
        tableauDeLectures = (char**) calloc(MAX_READS_TABLE_SIZE + 1, sizeof(char*));
        taille = lectureDeGenomes(tableauDeLectures, MAX_READS_TABLE_SIZE, fichierDesLectures);
        if (taille != 0) {
            printf("0\n");
            trierTableauDeLectures(tableauDeLectures, taille);
            printf("1\n");
            ecrireTableauDeLecturesDansUnFichier(tableauDeLectures, taille, &nombreDeFichiers);
            printf("2\n");
            nombreDeFichiers++; //on incremente le nombre de fichiers de sauvegardes de lectures triées
            free(tableauDeLectures); //on libere la memoire du tableau
        }
    } while(taille == MAX_READS_TABLE_SIZE);
    //nous avons ainsi créé n fichiers de lectures, tous triés
    //nous allons maintenant effectué le tri-fusion de ces fichiers
    //une partie de chacun de ces fichiers sera traité en mémoire, nous allons d'abord déterminé ce nombre
    printf("tris separés terminés (%d fichiers ont été créés)\n", nombreDeFichiers);
    int nombreMaxDeLectures = MAX_READS_TABLE_SIZE / (nombreDeFichiers + 1); //on crée autant de tableau que de fichiers
    //plus le tableau resultats où seront stockés temporairement une partie du resultat du tri-fusion
    printf("Lancement du tri fusion\n");
    triFusion(nombreDeFichiers, nombreMaxDeLectures);

    //
}
int lectureDeGenomes(char** tableau, int nombreDeGenomesALire, FILE * fichier) {
    char READ[SIZE_READ];
    int indice = 0;
    while(fgets(READ, SIZE_READ + 1, fichier) != NULL) {
        if (READ[0] == 'A' || READ[0] == 'T' || READ[0] == 'G' || READ[0] == 'C') {
            *(tableau + indice) = (char*) malloc(SIZE_READ + 1);
            strcpy(*(tableau + indice), READ);
            indice++;
            if (indice == nombreDeGenomesALire) {
                return indice;
            }
        }
    }
    return indice;
}

static int compare (void const *a, void const *b) {
   char const *const *pa = a;
   char const *const *pb = b;
   return strcmp (*pa, *pb);
}

void trierTableauDeLectures(char** tableau, int taille) {
    qsort (tableau, taille, sizeof *tableau, compare);
}

void ecrireTableauDeLecturesDansUnFichier(char** tableau, int taille, int* nombreDeFichiers) {
    FILE * tempFile = NULL;
    char nomFichier[20];
    sprintf(nomFichier, "%d.txt", *(nombreDeFichiers) + 1);
    tempFile = fopen(nomFichier, "w");
    for(int i = 0; i < taille; i++) {
        fprintf(tempFile, "%s\n", *(tableau + i));
    }
    fclose(tempFile);
}
//on ajoute le contenu du tableaux resultats (resultas du tri-fusion) dans le fichierResultats
//on vide le tableau resultats
void ecrireTableauResultatsDansUnFichier(char** tableau, int taille, FILE* fichierResultats) {
    for(int i = 0; i < taille; i++) {
        fprintf(fichierResultats, "%s\n", *(tableau + i));
        *(tableau + i) = "";
    }
}

void triFusion(int nombreDeFichiers, int nombreMaxDeLectures) {
    printf("debut du tri-fusion\n");
    char*** lectures;
    char** resultats;
    FILE** files;
    printf("tableaux créés, debut des allocations de memoire\n");
    files = (FILE**) calloc(nombreDeFichiers + 1, sizeof(char*));
    lectures = (char***) calloc(nombreDeFichiers + 1, sizeof(char**));
    int* indices = (int*) calloc(nombreDeFichiers + 1, sizeof(int));
    resultats = (char**) calloc(nombreMaxDeLectures + 1, sizeof(char*));
    printf("tableaux créés, fin des allocations de memoire\n");
    FILE* fichierResultats = fopen("resultats.txt", "w");
    printf("debut initialisation des tableaux\n");
    for(int i = 0; i < nombreDeFichiers; i++) {
        char nomFichier[20];
        sprintf(nomFichier, "%d.txt", i + 1);
        files[i] = (FILE*) malloc(sizeof(FILE*) + 1);
        files[i] = fopen(nomFichier, "r");
        lectures[i] = (char**) calloc(nombreMaxDeLectures + 1, sizeof(char*));
        indices[i] = lectureDeGenomes(lectures[i], nombreMaxDeLectures, files[i]);
    }
    printf("fin initialisation des tableaux\n");
    printf("debut initialisation du tableau des resultats\n");
    for (int i = 0; i < nombreMaxDeLectures; i++) {
        resultats[i] = (char*) malloc(sizeof(SIZE_READ + 1));
    }
    printf("fin initialisation du tableau des resultats\n");
    while(triFusionTermine(indices, nombreDeFichiers) == 1) {
        printf("debut etape tri-fusion\n");
        for(int i = 0; i < nombreMaxDeLectures; i++) {
            restaurationDeLectures(lectures, nombreDeFichiers, nombreMaxDeLectures, indices, files);
            int indiceMin = indiceMinTableauLectures(lectures, nombreDeFichiers, indices);
            *(resultats + i) = lectures[indiceMin][nombreMaxDeLectures - indices[indiceMin]];
            indices[indiceMin]--;
        }
        printf("fin etape tri-fusion\n");
        printf("debut ecriture des resultats\n");
        ecrireTableauResultatsDansUnFichier(resultats, nombreMaxDeLectures, fichierResultats);
        printf("fin ecriture des resultats\n");
    }
    fclose(fichierResultats); // a la fin du tri-fusion on ferme le fichierResultats, qui était ouvert en mode ecriture 
    //on libere la mémoire précédemment allouée
    printf("fin du tri-fusion\n");
    printf("liberation de la memoire allouée\n");
    free(resultats);
    for(int i = 0; i < nombreDeFichiers; i++) {
        free(lectures[i]);
    }
    free(lectures);
    free(files);
    free(indices);
}
int triFusionTermine(int * indices, int nombreDeFichiers) {
    int termine = 0;
    for(int i = 0; i < nombreDeFichiers; i++) {
        if(indices[i] != 0) {
            return 1; // 1 -> pas terminé
        }
    }
    return 0; // 0 -> terminé
}
int indiceMinTableauLectures(char*** lectures, int nombreDeFichiers, int* indices) {
    int indiceMin = -1; // on ne connait pas encore le tableau ayant la plus petite lecture (ordre alphabetique)
    for(int i = 0; i < nombreDeFichiers; i++) {
        //on recherche l'indice du premier tableau non vide
        //s'il est nonvide de lectures, on affecte l'indice de ce tableau (i) à indiceMin
        if (indiceMin == -1 && indices[i] != 0) {
            indiceMin = i;
        } else if (indiceMin != -1 && indices[i] != 0) { //si indiceMin est non vide, on verifie si la lecture du tableau d'indice indiceMin
        //est plus grande que la lecture du tableau d'indice i, alors indiceMin prend la valeur de l'indice i
            int cmp = strcmp(lectures[indiceMin][indices[indiceMin]], lectures[i][indices[i]]);
            if(cmp > 0) {
                indiceMin = i;
            }
        }
    }
    return indiceMin; // on retourne indiceMin apres avoir parcouru tout le tableau

}
void restaurationDeLectures(char*** lectures, int nombreDeFichiers, int nombreMaxDeLectures, int* indices, FILE** files) {
    for(int i = 1; i < nombreDeFichiers; i++) {
        //si indice (indices[i]) descend a 1 alors toutes les lectures chargéss ont ete parcourus, il faut lire de nouveau dans le fichier pour charger de nouvelles lecture
        if (indices[i] == 1) {
            // on remplace toutes les lectures du tableau par des chaines vides
            for(int j = 0; j < nombreMaxDeLectures; j++) {
                lectures[i][j] = "";
            }
            //on remplace les valeurs du tableau de lectures par le maximum de lectures possibles à
            //partir du fichier dont les précédentes lectures avaient été extraites
            indices[i] = lectureDeGenomes(lectures[i], nombreMaxDeLectures, files[i]);
        }
    }
}
