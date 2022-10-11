#include "yxq.h"
#include "tools.h"

void printUsage()
{
      printf("Usage: dbndatagen [-str fnmstr] -pssm fnmpssm ");
      printf("-o fnmout -laa -lss -dmax -linear[-sigmoid] [-R]\n");
      fflush(stdout);
}
int ipow(int x, int y)
{
   int i,v=1;
   for(i=0; i<y; i++) v=v*x;
   return v;
}
/* 1 based for matlab programs */
int main(int argc, char** argv)
{
   if(argc<2 || !strcmp(argv[1],"-"))
   {
      printUsage();
      return 0;
   }
   
   FILE *fp1, *fp2, *fpout;
   int   laa, lss, dmax, i,j,k,l,m;
   bool  breverse = false,
         blinear  = false,
         bsigmoid = false,
         btrain   = false;
   for(i=1,j=0; i<argc; i++)
   {
      if(!strcmp(argv[i],"-str"))
      {
         fp1 = fopen(argv[i+1],"r");
         btrain = true;
      }
      else if(!strcmp(argv[i],"-pssm"))
      {
         fp2 = fopen(argv[i+1],"r");
         j++;
      }
      else if(!strcmp(argv[i],"-o"))
      {
         fpout = fopen(argv[i+1],"w");
         j++;
      }
      else if(!strcmp(argv[i],"-laa"))
      {
         laa = atoi(argv[i+1]);
         j++;
      }
      else if(!strcmp(argv[i],"-lss"))
      {
         lss = atoi(argv[i+1]);
         j++;
      }
      else if(!strcmp(argv[i],"-dmax"))
      {
         dmax = atoi(argv[i+1]);
         j++;
      }
      else if(!strcmp(argv[i],"-sigmoid"))
      {
         bsigmoid = true;
         j++;
      }
      else if(!strcmp(argv[i],"-linear"))
      {
         blinear = true;
         j++;
      }
      else if(!strcmp(argv[i],"-R"))
         breverse = true;
   }
   if(j<6)
   {
      printUsage();
      return 0;
   }
   if(blinear && bsigmoid)
   {
      printf("[Error] Can't set both linear and sigmoid transformation\n");
      printUsage();
      return 1;
   }

   int        nd=lss,
              nnds=(laa>0?1:0)+lss+4,
              maxlenseg = -1,
              len,
              nsegs,
              nseqs=0;
   sSegment **pseg;
   real      *pssm,
              val[21], /* R's value stored */
              t;
   char      *anno,
             *str,
              ch;
   while(len=readpssm(fp2,&pssm))
   { 
      if(btrain)
      {
         readfasta(fp1,&anno,&str);
         for(i=0; i<len; i++)
         {
            str[i]=getstr(transtr(str[i]));
         }
      }
      if(breverse)
      {
         for(i=0; i<len/2; i++)
         {
            if(btrain)
            {
               ch=str[i]; 
               str[i]=str[len-i-1]; 
               str[len-i-1]=ch;
            }
            for(j=0; j<20; j++)
            { 
               t=pssm[i*20+j];
               pssm[i*20+j]=pssm[(len-i-1)*20+j];
               pssm[(len-i-1)*20+j]=t;
            }
         }
      }
      for(i=0; i<len; i++)
      {
         for(j=0; j<20; j++)
         {
            if(blinear)
               pssm[i*20+j] = linear(pssm[i*20+j],7,-7);
            if(bsigmoid)
               pssm[i*20+j] = sigmoid(pssm[i*20+j]); 
         }
      }
      if(btrain)
      {
         nsegs=split(str,&pseg);
         for(i=0; i<nsegs; i++)
         {
             if(pseg[i]->end - pseg[i]->start + 1 > maxlenseg)
                maxlenseg = pseg[i]->end - pseg[i]->start + 1;
             for(j=pseg[i]->start; j<=pseg[i]->end; j++)
             {
                k=pseg[i]->end-j+1;
                fprintf(fpout,"%d,%d,%d",transtr(str[j])+1,k>=dmax?dmax:k,k==1?2:1);
                for(k=0; k<lss; k++)
                {
                   l=j-lss+k;
                   if(l<0) m=4;
                   else m=transtr(str[l])+1;
//		   m=m+dd[k]*ipow(4,lss-k-1);
                   fprintf(fpout,",%d",m);
                }
//		if(lss>0) fprintf(fpout,",%d",m+1);
                for(k=0; k<laa; k++)
                {
                   l=j-laa+k;
                   if(l<0)
                   {
                      for(m=0; m<20; m++) val[m]=0;
                      val[20]=1;
                   }
                   else
                   {
                      for(m=0; m<20; m++) val[m]=pssm[l*20+m];
                      val[20]=0;
                   }
                   for(m=0; m<20+1; m++) fprintf(fpout,",%.3f",val[m]);
                }
                for(k=0; k<20; k++) fprintf(fpout,",%.3f",pssm[j*20+k]);
                fprintf(fpout,"\n");
             }
         }
         for(i=0; i<nsegs; i++) free(pseg[i]);
         free(pseg);
      }
      else
      {
         for(i=0; i<len; i++)
         {
            for(j=0; j<2+lss; j++) fprintf(fpout,",");
            for(j=0; j<laa; j++)
            {
               k=i-laa+j;
               if(k<0)
               {
                  for(l=0; l<20; l++) val[l]=0;
                  val[20]=1;
               }
               else
               {
                  for(l=0; l<20; l++) val[l]=pssm[k*20+l];
                  val[20]=0;
               }
               for(l=0; l<20+1; l++) fprintf(fpout,",%.3f",val[l]);
            }
            for(j=0; j<20; j++) fprintf(fpout,",%.3f",pssm[i*20+j]);
            fprintf(fpout,"\n");
         }
      }
      fprintf(fpout,"//\n");    
      nseqs++;
     /* free space */
      free(pssm);
      if(btrain) {free(anno), free(str);}
   }
   fclose(fp2);
   fclose(fpout);
   if(btrain) fclose(fp1);
   return 0;
} 
