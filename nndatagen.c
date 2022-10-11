#include "yxq.h"
#include "tools.h"
void printUsage()
{
   printf("Usage: nndatagen [-str fnmstr] [-pssm fnmpssm] [-nn1o fnmnn1output] ");
   printf("-o output -w -linear[-sigmoid]\n");
}
int main(int argc, char** argv) 
{
   if(argc < 2 || !strcmp(argv[1],"-")) 
   {
      printUsage();
      return 0;
   }

   FILE *fp1=NULL, *fp2=NULL, *fp3=NULL, *fpout=NULL;
   int  i,j,k,l,w=0;
   bool blinear  = false, 
        bsigmoid = false,
        btrain   = false,
        bsecnn   = false;
   for(i=1; i<argc; i++)
   {
      if(!strcmp(argv[i],"-str")) 
      {
         fp1=fopen(argv[i+1],"r");
         btrain=true;
      }
      else if(!strcmp(argv[i],"-pssm"))
      {
         fp2=fopen(argv[i+1],"r");
      }
      else if(!strcmp(argv[i],"-nn1o"))
      {
         fp3=fopen(argv[i+1],"r");
         bsecnn=true;
      }
      else if(!strcmp(argv[i],"-o"))
         fpout=fopen(argv[i+1],"w");
      else if(!strcmp(argv[i],"-w"))
         w=atoi(argv[i+1]);
      else if(!strcmp(argv[i],"-linear"))
         blinear=true;
      else if(!strcmp(argv[i],"-sigmoid"))
         bsigmoid=true;
   }
   if(!w || !fpout || (!blinear && !bsigmoid && !bsecnn) 
      || (blinear && bsigmoid) || (!fp2 && !bsecnn))
   {
      printUsage();
      return 1;
   }
   if(bsecnn) blinear=bsigmoid=false;

   char *anno,*str, ch;
   real *pssm, x[21]; /* pssm stores two kinds of data */
   int   len, ncols=0;
   if(bsecnn) {len=readplain(fp3,&pssm)/3; ncols=3;}
   else if(fp2) {len=readpssm(fp2,&pssm); ncols=20;}
   while(len)
   {
      if(btrain) 
      {
         readfasta(fp1,&anno,&str);
         for(j=0; j<len; j++)
         {
            str[j]=getstr(transtr(str[j]));
         }
      }
      for(i=0; i<len; i++)
      {
         for(j=0; j<ncols; j++)
         {
            if(blinear)
               pssm[i*ncols+j] = linear(pssm[i*ncols+j],7,-7);
            if(bsigmoid)
               pssm[i*ncols+j] = sigmoid(pssm[i*ncols+j]);
         }
      }
      for(i=0; i<len; i++)
      {
         for(j=0; j<w; j++)
         {
            k=i-(w-1)/2+j;
            for(l=0; l<ncols+1; l++) x[l]=0.0;
            if(k<0||k>=len)
            {
               x[ncols]=1.0;
            }
            else
            {
               for(l=0; l<ncols; l++) x[l] = pssm[k*ncols+l];
            }
            for(l=0; l<ncols+1; l++) fprintf(fpout,"%.3f,",x[l]);
         }
         if(btrain) 
         {
            ch=str[i];
            fprintf(fpout,"%1d,%1d,%1d\n",ch=='H',ch=='E',ch=='C');
         }
         else fprintf(fpout,",,,\n");
      }
      fprintf(fpout,"//\n");

     /* free up space */
      if(btrain) {free(anno),free(str);}
      free(pssm);
      if(bsecnn) len=readplain(fp3,&pssm)/3;
      else if(fp2) len=readpssm(fp2,&pssm);
   }
   if(btrain) fclose(fp1);
   if(fp2) fclose(fp2);
   if(bsecnn) fclose(fp3);
   fclose(fpout);
   return 0;
}
