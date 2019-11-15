
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#include <pthread.h>

#define SIZE_READ      100

int main (int argc, char *argv[])
{

  FILE *fg, *fr;
  char *GENOME, *READ, c;
  int i,j,x,pos_erreur;
  int SIZE_GENOME;
  int COVERAGE;
  char int2nt[4] = {'A','C','G','T'};
  int NB_READS;

  // lecture des parametres
  if (argc != 3)
    {
      printf ("USAGE : %s  <genome size> <coverage>\n",argv[0]);
      exit(0);
    }
  sscanf (argv[1],"%d",&x);
  SIZE_GENOME = x * 1000000;
  sscanf (argv[2],"%d",&COVERAGE);
  NB_READS = (SIZE_GENOME / SIZE_READ) * COVERAGE;


  // Generation aleatoire du genome

  // reseravtion place en memoire pour le genome
  GENOME = (char *) malloc(SIZE_GENOME);
  READ = (char *) malloc(SIZE_READ);
  // ouverture du fichier qui contiendra le texte du genome
  fg = fopen("genome.fasta","w");
  // ecriture du commentaire demande par le format fasta
  fprintf (fg,">genome");

  for (i=0; i<SIZE_GENOME; i++)
    {
      // retour chariot tous les 80 caracteres
      if ((i%80)==0) fprintf (fg,"\n");
      // generation aleatoire d'un caractere A, C, G ou T
      c = int2nt[rand()%4];
      // ecriture de ce caratere dans le fichier
      fprintf (fg,"%c",c);
      // memorisation du caratere dans le tableau genome
      GENOME[i] = c;
    }
  // ecriture dernier retour chariot
  fprintf(fg,"\n");
  // fermeture du ficher genome.fasta
  fclose(fg);

  // 2 - generation jeu de lectures

  fr = fopen("reads.fasta","w");
  for (i=0; i<NB_READS; i++)
    {
      // on prend au hazard une position du genome
      x = rand()%(SIZE_GENOME-SIZE_READ);
      // on recopie dans READ la portion du genome a la position x
      strncpy(READ,&GENOME[x],SIZE_READ);
      // on ajoute du bruit (erreur de sequencage)
      // max = 1 erreur de substitution
      pos_erreur = -1;
      for (j=0; j<SIZE_READ; j++)
	{
	  if ((rand()%400) == 0)
	    {
	      c = int2nt[rand()%4];
	      if (c != READ[j])
		{
		  READ[j] = c;
		  pos_erreur = j;
		  break;
		}
	    }
	}
      // on ecrit la lecture au format fasta
      if (pos_erreur >= 0)
	// on indique la position de la lecture et la position de l'erreur dans la lecture
	fprintf(fr,">read_%d %d %d\n%s\n",i,x,pos_erreur,READ);
      else
	// on indique seulement la position de la lecture
	fprintf(fr,">read_%d %d\n%s\n",i,x,READ);

    }
  // fermeture du fichier reads.fasta
  fclose(fr);
  // liberationde la memoire
  free(GENOME);
  free(READ);

}
