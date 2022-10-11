#include "yxq.h"
#include "tools.h"
/* Usage: combine file1 file2 [file3...] */
int main(int argc, char *argv[])
{
   if(argc<3)
   {
      printf("Usage: combine file1 file2 [file3 ...]\n");
      return 0;
   }
   int nfiles=argc-1,i,j,len;
   real *data, *dataT;
   FILE *fp[nfiles];

   for(i=0; i<nfiles; i++) fp[i]=fopen(argv[i+1],"r");
   while(len=readplain(fp[0],&data)/3)
   {
      for(i=1; i<nfiles; i++)
      {
        readplain(fp[i],&dataT); 
        for(j=0; j<len; j++) 
        {
           normalize(dataT+j*3,3);
           data[j*3]  += dataT[j*3];
           data[j*3+1]+= dataT[j*3+1];
           data[j*3+2]+= dataT[j*3+2];
        }
        free(dataT);
      }
      for(i=0; i<len; i++) 
         printf("%.3f,%.3f,%.3f\n",data[i*3]/nfiles,
                 data[i*3+1]/nfiles,data[i*3+2]/nfiles);
      printf("//\n");
      free(data);
   }
   for(i=0; i<nfiles; i++) fclose(fp[i]); 
   return 0;
}
