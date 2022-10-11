#include "yxq.h"
#include "tools.h"
/*Read in plain formated prediction 
  and print it out in FASTA format.
  arguments: rawpred [secstr.fasta] 
*/
int main(int argc, char* argv[])
{
   if(argc<2)
   {
      fprintf(stderr,"Usage: charpred rawpred [secstr.fasta]\n");
      return 0;
   }
   char *seq,*pred,*anno;
   real *rawpred;
   int  len,i;
   FILE *fp1, *fp2=NULL;
   fp1=fopen(argv[1],"r");
   if(argc>2) fp2=fopen(argv[2],"r"); 
   while(len=readplain(fp1,&rawpred)/3)
   {
      if(fp2) readfasta(fp2,&anno,&seq);
      pred=(char*)malloc((len+1)*sizeof(char));
      for(i=0; i<len; i++)
      {
         pred[i] = getstr(maxind(rawpred+i*3,3));
      }
      pred[i]='\0';
      if(fp2) writefasta(stdout,anno,pred);
      else writefasta(stdout,NULL,pred);
   
     /* free up space */
      free(pred);
      if(fp2) {free(anno); free(seq);}
      free(rawpred);
   }
   fclose(fp1);
   if(fp2) fclose(fp2);
   return 0;
}

