#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>

#define SIZE_READ      100

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
int triFusionTermine(int * indices, int * tailles, int nombreDeFichiers);

/*
*Un nombre nombreDeFichiers de tableau de tableaux de lectures etant créés dans triFusion(),
*On parcours ces differents tableaux contenus dans le atbleau lectures pour
*retourner l'indice du tableau (indice dans lectures) ayant la plus petite des lecture
*/
int indiceMinTableauLectures(char*** lectures, int nombreDeFichiers, int* indices, int* tailles);

/*
*On verifie pour chaque tableau de lectures si toutes ces lectures ont été parcourues
*si indices[indiceDuTableau] est egale à 1 alors toutes les lectures du tableau ont été parcourues
*Alors on remplace les lectures de ce tableau pas des chaines de caractères vides
*Puis on remplace le maximum possible par de nouvelles lectures à partir du fichier dont les précédentes avaient
*été extraites.
*/
void restaurationDeLectures(char*** lectures, int nombreDeFichiers, int nombreMaxDeLectures, int* indices, int* tailles, FILE** files);

/*
*On supprime du tableau toutes les valeurs qui n'existent qu'en une seule instance
*/
void filtrageTableau(char** lectures, int taille);

/*
*On supprime du fichier toutes les valeurs qui n'existent qu'en une seule instance
*Pour ce faire, on prend du fichier le maximum de valeurs possible
*Puis on les filtre avec la fonction filtrageTableau
*/
void filtrageFichier(FILE* fichier, int nombreMaxDeLectures);
/*
*retourne le temps en secondes
*/
double getTime(void); /**/

int main (int argc, char *argv[]) {
    // lecture des parametres
    if (argc != 2)
    {
      printf ("USAGE : %s  <TAILLE MAXIMALE AUTORISEE (RAM en Mo)>\n",argv[0]);
      exit(0);
    }
    /*
    *On suppose que la taille maximale de la mémoire avec laquelle
    *nous pouvons travailler est de MAX_MEMORY_SIZE Mo
    *Ce qui nous amenera a appliquée une méthode de gestion
    *optimisée de la mémoire, pour la réalisation de ce TP
    */
    int MAX_MEMORY_SIZE;

    /*
    *Taille maximale que peut avoir un tableau de reads en tenant
    *compte de MAX_MEMORY_SIZE
    */
    int MAX_READS_TABLE_SIZE;
    sscanf(argv[1], "%d",&MAX_MEMORY_SIZE);
    MAX_MEMORY_SIZE = MAX_MEMORY_SIZE * 1000000;
    MAX_READS_TABLE_SIZE = MAX_MEMORY_SIZE / 100;

    double tempsInitial = getTime();
    printf("debut du programme\t\t\t\t[time = 0sec]\n");
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
    printf("debut de la lecture du fichier reads.fasta\t[time = %1.2lf sec]\n", getTime()-tempsInitial);
    tableauDeLectures = (char**) calloc(MAX_READS_TABLE_SIZE + 1, sizeof(char*));
    do {
        printf("debut de la lecture et du tri n°%d\t\t[time = %1.2lf sec]\n", nombreDeFichiers+1, getTime()-tempsInitial);
        taille = lectureDeGenomes(tableauDeLectures, MAX_READS_TABLE_SIZE, fichierDesLectures);
        if (taille != 0) {
            trierTableauDeLectures(tableauDeLectures, taille);
            // dans le cas où une seule lecture du fichier fasta ne suffit pas pour recuperer toutes les lectures
            // on sauvegarde alors chaque tableau lecture trié dans un fichier sur le disque
            if (nombreDeFichiers == 0 && (taille < MAX_READS_TABLE_SIZE || taille == 0)) {
                nombreDeFichiers++; //on incremente le nombre de fichiers de sauvegardes de lectures triées
                printf("fin de la lecture et du tri\t\t\t[time = %1.2lf sec]  tout le contenu est chargé en mémoire\n", getTime()-tempsInitial);
                break;
            } else {
                ecrireTableauDeLecturesDansUnFichier(tableauDeLectures, taille, &nombreDeFichiers);
                nombreDeFichiers++; //on incremente le nombre de fichiers de sauvegardes de lectures triées
                //free(tableauDeLectures); //on libere la memoire du tableau
            }
        }
        printf("fin de la lecture et du tri n°%d\t\t\t[time = %1.2lf sec]\n", nombreDeFichiers, getTime()-tempsInitial);
    } while(taille == MAX_READS_TABLE_SIZE);
    fclose(fichierDesLectures); // on ferme le fichier fasta.reads
    if(nombreDeFichiers > 1) { //cas où il y'a eu necessité de créer des fichiers intermédiaires
        //nous avons ainsi créé n fichiers de lectures, tous triés
        //nous allons maintenant effectué le tri-fusion de ces fichiers
        //une partie de chacun de ces fichiers sera traité en mémoire, nous allons d'abord déterminé ce nombre
        int nombreMaxDeLectures = MAX_READS_TABLE_SIZE / (nombreDeFichiers + 1); //on crée autant de tableau que de fichiers
        //plus le tableau resultats où seront stockés temporairement une partie du resultat du tri-fusion
        free(tableauDeLectures); //on libere la memoire du tableau
        printf("fin de la lecture du fichier fasta.reads\t[time = %1.2lf sec]\n", getTime()-tempsInitial);
        printf("debut du tri fusion\t\t\t\t[time = %1.2lf sec]\n", getTime()-tempsInitial);
        triFusion(nombreDeFichiers, nombreMaxDeLectures);
        //les données ayant maintenant été triées et stockées dans le ,fichier resultats.txt
        //Nous allons maintenant les filtrer en supprimant celles qui n'ont qu'une seule occurence
        //Mais avant cela nous allons d'abord effacer les fichiers de resultats intermediaires créés sur le disque
        for (int i = 1; i <= nombreDeFichiers; i++) {
            char nomFichier[20];
            sprintf(nomFichier, "%d.txt", i);
            remove(nomFichier);
        }
        printf("fin du tri fusion\t\t\t\t[time = %1.2lf sec]\n", getTime()-tempsInitial);
        printf("debut du filtrage\t\t\t\t[time = %1.2lf sec]\n", getTime()-tempsInitial);
        FILE* fichierDesLectures = fopen("resultats.txt", "r");
        //FILE* fichierFinal = fopen("Sort_reads.txt", "w");
        filtrageFichier(fichierDesLectures, MAX_READS_TABLE_SIZE);
        remove("resultats.txt");
        //fclose(fichierFinal);
        fclose(fichierDesLectures);
        printf("fin du filtrage\t\t\t\t\t[time = %1.2lf sec]\n", getTime()-tempsInitial);
    } else { //cas où il n y'a pas eu necessité de créer des fichiers intermédiaires
        printf("fin de la lecture du fichier fasta.reads\t[time = %1.2lf sec]\n", getTime()-tempsInitial);
        printf("debut du filtrage du tableau\t\t\t[time = %1.2lf sec]\n", getTime()-tempsInitial);
        FILE* fichierFinal = fopen("Sort_reads.txt", "w");
        filtrageTableau(tableauDeLectures, taille);
        ecrireTableauResultatsDansUnFichier(tableauDeLectures, taille, fichierFinal);
        free(tableauDeLectures);
        fclose(fichierFinal);
        printf("fin du filtrage du tableau\t\t\t[time = %1.2lf sec]\n", getTime()-tempsInitial);
    }
    printf("\n\nfin du programme, consultez le resultat dans le fichier Sort_reads.txt\n");
    printf("-----------------------------------------\n");
    printf("|Temps total d'exécution: %1.2lf secondes|\n", getTime()-tempsInitial);
    printf("-----------------------------------------\n");
}
int lectureDeGenomes(char** tableau, int nombreDeGenomesALire, FILE * fichier) {
    char READ[SIZE_READ]; //chaine de caractère de même longueur qu'une lecture
    int taille = 0; //le nombre de lectures lues
    while(fgets(READ, SIZE_READ + 1, fichier) != NULL) { //on parcours le fichier
        if (READ[0] == 'A' || READ[0] == 'T' || READ[0] == 'G' || READ[0] == 'C') {
            //si une ligne commence par A, T, G ou C
            *(tableau + taille) = (char*) malloc(SIZE_READ + 1); //on alloue la memoire necessaire 
            strcpy(*(tableau + taille), READ); //ajoute la ligne au tableau
            taille++; // on incremente taille pour signifier qu'un nouvel element a ete lu et ajouté au tableau
            if (taille == nombreDeGenomesALire) { //si la taille maximale de lecture autorisé est atteinte
                return taille; //on retourne la taille
            }
        }
    }
    return taille; //on retourne la taille
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
    //on crée un fichier dont le nom est un chiffre inferieur au nombreDeFichiers à créer
    FILE * tempFile = NULL;
    char nomFichier[20];
    sprintf(nomFichier, "%d.txt", *(nombreDeFichiers) + 1);
    tempFile = fopen(nomFichier, "w");
    for(int i = 0; i < taille; i++) {
        if (strcmp(tableau[i], "") != 0) {
            fprintf(tempFile, "%s\n", *(tableau + i)); //on copie toute les valeurs du tableau dans le fichier
            *(tableau + i) = "";
        }
    }
    fclose(tempFile);
}
//on ajoute les lectures non-vides du tableaux dans le fichier
//on vide le tableau resultats
void ecrireTableauResultatsDansUnFichier(char** tableau, int taille, FILE* fichier) {
    for(int i = 0; i < taille; i++) {
        //printf("i=%d\n",i);
        if (strcmp(*(tableau + i), "") != 0) {
            fprintf(fichier, "%s\n", *(tableau + i));
            *(tableau + i) = "";
        }
    }
}

void triFusion(int nombreDeFichiers, int nombreMaxDeLectures) {
    char*** lectures; //tableau de tableaux de chaines de carateres
    //chaque tableau (lectures[i]) servira a stocké une partie du fichier portant le meme nom que son index
    char** resultats; //c'est ici que seront stockés les resultats temporaires du tri-fusion
    FILE** files; //les fichiers numerotés seront stocké dans ce tableau pour leurs utilisation ultérieures
    //allocations de memoire
    files = (FILE**) calloc(nombreDeFichiers + 1, sizeof(char*));
    lectures = (char***) calloc(nombreDeFichiers + 1, sizeof(char**));
    //le tableau d'entiers indices est un tableau d'entier qui a un indice i donné
    //contient le nombre d'éléments non encore sauvegardés du tableau lectures[i]
    //il nous permettra de savoir quel est le plus petit element non encore traité d'un tableau donné
    int* indices = malloc(nombreDeFichiers * sizeof(int));
    int* tailles = malloc(nombreDeFichiers * sizeof(int)); //les tailles des differents des tableau contenu dans lecture
    resultats = (char**) calloc(nombreMaxDeLectures + 1, sizeof(char*));
    FILE* fichierResultats = fopen("resultats.txt", "w"); //c'est dans ce fichiers que sera sauvegardé le
    //resultat du tri fusion
    //initialisation des tableaux précedemment créés
    for(int i = 0; i < nombreDeFichiers; i++) {
        char nomFichier[20];
        sprintf(nomFichier, "%d.txt", i + 1);
        files[i] = (FILE*) malloc(sizeof(FILE*) + 1);
        files[i] = fopen(nomFichier, "r");
        lectures[i] = (char**) calloc(nombreMaxDeLectures + 1, sizeof(char*));
        //indices[i] = lectureDeGenomes(lectures[i], nombreMaxDeLectures, files[i]);
        indices[i] = 0;
        tailles[i] = lectureDeGenomes(lectures[i], nombreMaxDeLectures, files[i]);
    }
    //on fait l'allocation memoire du tableau des resultats
    for (int i = 0; i < nombreMaxDeLectures; i++) {
        resultats[i] = (char*) malloc(SIZE_READ + 1);
    }
    //tri fusion
    //on parcours lectures
    //on recupere l'indice du tableau ayant le plus petit element non encore traité
    //on recuperes l'element concerné et on décremente le nombre d'élements non encore traité du tableau en question
    //on ecrit l'element en question dans le tableau resultats
    //quand le tableau resultats est vide, on le vide et on deverse son contenu dans le fichier resultat
    while(triFusionTermine(indices, tailles, nombreDeFichiers) == 1) {
        for(int i = 0; i < nombreMaxDeLectures; i++) {
            restaurationDeLectures(lectures, nombreDeFichiers, nombreMaxDeLectures, indices, tailles, files);
            int indiceMin = indiceMinTableauLectures(lectures, nombreDeFichiers, indices, tailles);
            if (indiceMin == -1) {
                break;
            }
            *(resultats + i) = lectures[indiceMin][indices[indiceMin]];
            lectures[indiceMin][indices[indiceMin]] = "";
            *(indices + indiceMin) = *(indices + indiceMin) + 1;
        }
        ecrireTableauResultatsDansUnFichier(resultats, nombreMaxDeLectures, fichierResultats);
    }
    fclose(fichierResultats); // a la fin du tri-fusion on ferme le fichierResultats, qui était ouvert en mode ecriture 
    //on libere la mémoire précédemment allouée
    free(resultats);
    for(int i = 0; i < nombreDeFichiers; i++) {
        free(lectures[i]);
    }
    free(lectures);
    free(files);
    free(indices);
}
int triFusionTermine(int* indices, int* tailles, int nombreDeFichiers) {
    int termine = 0;
    for(int i = 0; i < nombreDeFichiers; i++) {
        if(tailles[i] != 0) {
            return 1; // 1 -> pas terminé
        }
    }
    return 0; // 0 -> terminé
}
int indiceMinTableauLectures(char*** lectures, int nombreDeFichiers, int* indices, int* tailles) {
    int indiceMin = -1; // on ne connait pas encore le tableau ayant la plus petite lecture (ordre alphabetique)
    for(int i = 0; i < nombreDeFichiers; i++) {
        //on recherche l'indice du premier tableau non vide
        //s'il est nonvide de lectures, on affecte l'indice de ce tableau (i) à indiceMin
        if (tailles[i] != 0) {
            indiceMin = i;
            break;
        }        
    }
    for(int i = indiceMin; i<nombreDeFichiers; i++){
        if (indiceMin != -1 && tailles[i] != 0) { //si indiceMin est non vide, on verifie si la lecture du tableau d'indice indiceMin
        //est plus grande que la lecture du tableau d'indice i, alors indiceMin prend la valeur de l'indice i
            if(strcmp(lectures[i][*(indices+i)], lectures[indiceMin][*(indices+indiceMin)]) < 0) {
                indiceMin = i;
            }
        }
    }
    return indiceMin; // on retourne indiceMin apres avoir parcouru tout le tableau

}
void restaurationDeLectures(char*** lectures, int nombreDeFichiers, int nombreMaxDeLectures, int* indices, int* tailles, FILE** files) {
    //printf("debut de restauration\n");
    for(int i = 0; i < nombreDeFichiers; i++) {
        //si indice (indices[i]) descend a 1 alors toutes les lectures chargéss ont ete parcourus, il faut lire de nouveau dans le fichier pour charger de nouvelles lecture
        if (indices[i] == tailles[i]) {
            //on remplace les valeurs du tableau de lectures par le maximum de lectures possibles à
            //partir du fichier dont les précédentes lectures avaient été extraites
            tailles[i] = lectureDeGenomes(lectures[i], nombreMaxDeLectures, files[i]);
            if (tailles[i] != 0) {
                indices[i] = 0;
            }
        }
    }
}
void filtrageFichier(FILE* fichier, int nombreMaxDeLectures) {
    char** lectures; //c'est dans ce tableau que seront stockés les lectures tirées de fichier
    int taille = 0; //on initialise à 0 pour supposer que le fichier est vide
    FILE* fichierFinal = fopen("Sort_reads.txt", "w"); //le fichier de sortie (trié+filtré)
    char dernierElement[SIZE_READ] = ""; //puisque les lectures seront filtrées par bloc a partir du fichier,
    //Si on ne peut statuer sur la validité du dernier element (s'il n'a q'une occurence alors que son successeur se trouve dans le prochain bloc)
    //on le sauvegarde pour le comparer par la suite au premier element du dernier bloc
    //dans le cas ou successeur est lui-même deja valide dans son bloc, on n'oublie la valeur dernierElement
    //sinon on l'ecrit directement dans le fichier
    do {
        lectures= (char**) calloc(nombreMaxDeLectures + 1, sizeof(char*)); //allocation de mémoire
        taille = lectureDeGenomes(lectures, nombreMaxDeLectures, fichier); //les lectures sont extraites de fichier, 
        //stockées dans lecture et leur nombre est retourné dans taille
        if (taille != 0) {
            //on s'assure de verifier que les dernier elements de bloc de filtrage seront bien ajoutés s'ils sont
            // valides
            if (taille > 1) {
                if(strcmp(lectures[0], dernierElement) == 0 && strcmp(lectures[0], lectures[1]) != 0) {
                    fprintf(fichierFinal, "%s\n", dernierElement);
                }
                if (strcmp(lectures[taille-2], lectures[taille-1]) == 0)
                {
                    strcpy(dernierElement, "");
                } else if (strcmp(lectures[taille-2], lectures[taille-1]) != 0) {
                    strcpy(dernierElement, lectures[taille-1]);
                }
            }
            
            filtrageTableau(lectures, taille); //si non vide alors filtrage du tableau
            ecrireTableauResultatsDansUnFichier(lectures, taille, fichierFinal);
        }
        free(lectures);
    } while(taille != 0);
    fclose(fichierFinal);
}
void filtrageTableau(char** lectures, int taille) {
    int selected = 0;
    int verifie = 1;
    for(int i = 1; i < taille; i++) {
        if (strcmp(lectures[selected], lectures[i]) == 0) {
            lectures[i] = "";
            verifie = 0;
        } else {
            if (verifie == 1) {
                lectures[selected] = "";
            }
            selected = i;
            verifie = 1;
        }
    }
}
double getTime(void) /**/
{
  struct timeval t;
  gettimeofday(&t, NULL);
  return (1.0e-6*t.tv_usec + t.tv_sec);
}
