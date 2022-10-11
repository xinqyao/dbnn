#include "tools.h"

/* DSSP definition */
int transtr(char ch) {
   int n;
   switch(ch) {
      case 'H':
      case 'G':
         n=0;
         break;
      case 'E':
      case 'B':
         n=1;
         break;
      case 'C':
      case 'I':
      case 'S':
      case 'T':
      case '_':
      case '?':
      case '-':
         n=2;
         break;
      default:
         n=2;
         printf("[Warning] Unknown letter for secondary structure\n");
   }
   return n;
}
int tranaa(char ch) 
{
   int n;
   switch(ch) 
   {
      case 'A':
      case 'X':
         n=0;
         break;
      case 'R':
         n=1;
         break;
      case 'N':
         n=2;
         break;
      case 'B':
      case 'D':
         n=3;      
         break;
      case 'C':
         n=4;
         break;
      case 'Q':
      case 'Z':
         n=5;
         break;
      case 'E':
         n=6;
         break;
      case 'G':
         n=7;
         break;
      case 'H':
         n=8;
         break;
      case 'I':
         n=9;
         break;
      case 'L':
         n=10;
         break;
      case 'K':
         n=11;
         break;
      case 'M':
         n=12;
         break;
      case 'F':
         n=13;
         break;
      case 'P':
         n=14;
         break;
      case 'S':
         n=15;
         break;
      case 'T':
         n=16;
         break;
      case 'W':
         n=17;
         break;
      case 'Y':
         n=18;
         break;
      case 'V':
         n=19;
         break;
      default:
         n=0;
         printf("[Warning] Unknown letter for amino acid\n");
   }
   return n;
}

char getaa(int n) 
{
   char ch;
   switch(n) 
   {
      case 0:
         ch='A';
         break;
      case 1:
         ch='R';
         break;
      case 2:
         ch='N';
         break;
      case 3:
         ch='D';
         break;
      case 4:
         ch='C';
         break;
      case 5:
         ch='Q';
         break;
      case 6:
         ch='E';
         break;
      case 7:
         ch='G';
         break;
      case 8:
         ch='H';
         break;
      case 9:
         ch='I';
         break;
      case 10:
         ch='L';
         break;
      case 11:
         ch='K';
         break;
      case 12:
         ch='M';
         break;
      case 13:
         ch='F';
         break;
      case 14:
         ch='P';
         break;
      case 15:
         ch='S';
         break;
      case 16:
         ch='T';
         break;
      case 17:
         ch='W';
         break;
      case 18:
         ch='Y';
         break;
      case 19:
         ch='V';
         break;
      default:
         ch='A';
         printf("[Warning] Bad index for amino acid (0~19)\n");
   }
   return ch;
}

char getstr(int n)
{
   char ch;
   switch(n) 
   {
      case 0:
         ch='H';
         break;
      case 1:
         ch='E';
         break;
      case 2:
         ch='C';
         break;
      default:
         ch='C';
         printf("[Warning] Bad index for secondary structure (0~2)\n");
   }
   return ch;
}

/*split a sequence into segments */
/* return the number of segments */
/* 0 based */
int split(char *seq, sSegment ***pSeg)
{
   int i=0,j,k=0;
   int len=strlen(seq);

   *pSeg=NULL;
   while(i<len)
   {
      if(is_space(seq[i]) || seq[i] == '\n') 
      {
         i++;
         continue;
      }
      j=i+1;
      if(k==0) *pSeg=(sSegment**)malloc(sizeof(sSegment*));
      else *pSeg=(sSegment**)realloc(*pSeg,(k+1)*sizeof(sSegment*));
      while(seq[j] == seq[i]) j++;
      (*pSeg)[k] = (sSegment*)malloc(sizeof(sSegment));
      (*pSeg)[k]->start = i;
      (*pSeg)[k]->end   = j-1;
      (*pSeg)[k]->type  = seq[i];
      k++;
      i=j;
   }
   return k;
}

/* read one field of data of plain format*/
real readdata(FILE *fp)
{
   int i=0;
   char c,buff[256];
   real dat=0;
   c=fgetc(fp);
   while(!feof(fp) && c!='\n' && c!=',')
   {
      if(~is_space(c)) buff[i++]=c;
      c=fgetc(fp);
   }
   buff[i]='\0';
   if(!strcmp(buff,"//")) dat=INF+1;  /* a new sequence */
   else if(i==0 && c==',') dat=-INF-1; /* missing data */
   else if(i==0 && !feof(fp)) dat=readdata(fp); /* continue to read */
   else if(i>0) 
   {
      dat=atof(buff);
      if(dat>=INF)
      {
         printf("[Warning] Data overflow\n");
         dat=0.0;
      }
      else if(dat<=-INF)
      {
         printf("[Warning] Data underflow\n");
         dat=0.0;
      }
   }
   return dat; /* if feof, dat=0 */
}

/* read a block of data of plain format */
/* if successful returns number of fields, 
   otherwise return false */
int readplain(FILE *fp, real **data)
{
   int i=0;
   real dat;
   
   *data=NULL; 
   dat=readdata(fp); 
   while(!feof(fp))
   {
      if(dat>=INF)
      {
         return i;
      }
      else
      {
         if(i==0) *data=(real*)malloc(sizeof(real));
         else *data=(real*)realloc(*data,(i+1)*sizeof(real));
         (*data)[i]=dat;
         i++;
         if(i>=MAXSIZEARRAY) 
         {
            printf("[Warning] Too many data for one sequence\n");
            return false;
         }
      }
      dat=readdata(fp);
   }
   return false;
}

/* read PSSM data once a sequence */
/* if successful returns the length of sequence,
   otherwise returns false */ 
/* this routine is not so robust */
int readpssm(FILE *fp,real **pssm)
{
   char buff[256];
   int i,j,j2;

   *pssm=NULL;
   fgets(buff,5,fp);
   while(!feof(fp) && strcmp(buff,"Last")) 
      fgets(buff,5,fp);
   if (feof(fp)) return false;

   /* if find a line started with "Last", it assumes 
      that the following lines are automatically arranged
      in the default PSSM format, with two consecutive
      newlines indicating the end of PSSM data
   */

   j=0; 
   while(!feof(fp) && fgetc(fp)!='\n');   /* ignore the "Last" line */
   while(!feof(fp) && fgetc(fp)!='\n');   /* ignore the Caption line */
   if(feof(fp)) return false;
   while(fpeek(fp)!='\n')    
   {
      fscanf(fp,"%d%s",&j2,buff);
      if(j2!=j+1)
      {
         printf("[Warning] Bad PSSM file\n");
         return false; 
      }
      if(j==0) *pssm=(real*)malloc(20*sizeof(real));
      else *pssm=(real*)realloc(*pssm,(j+1)*20*sizeof(real));
#ifdef __DOUBLE_ACCURACY
      for(i=0; i<20; i++) fscanf(fp,"%lf",*pssm+j*20+i);
#else
      for(i=0; i<20; i++) fscanf(fp,"%f", *pssm+j*20+i);
#endif
      j=j2;
      if(j>=MAXSEQLEN) 
      {
         printf("[Warning] Too long sequence\n");
         return false;
      }
      while(!feof(fp) && fgetc(fp)!='\n');
      if(feof(fp)) return false;
   }
   return j;
}

/* read a sequence in FASTA format 
   if success it returns the length of sequence,
   otherwise returns false */
/* not so good... */
int readfasta(FILE *fp,char **anno,char **seq)
{
   char ch;
   int  k;
   
   *anno=NULL, *seq=NULL;
   while(!feof(fp) && fgetc(fp)!='>');
   if(feof(fp)) return false;
   *anno=(char*)malloc(MAXSEQLEN*sizeof(char));
   *seq =(char*)malloc(MAXSEQLEN*sizeof(char));
   k=0; 
   while(((*anno)[k]=fgetc(fp))!='\n')
   {
      k++;
      if(k>=MAXSEQLEN) {printf("[Warning] Too long annotation\n");k--;break;}
   }
   (*anno)[k]='\0';
   k=0;
   while(!feof(fp) && fpeek(fp)!='>')
   {
      ch=fgetc(fp);
      if((ch>='A' && ch<='Z'||ch=='-'||ch=='?'||ch=='_'||ch=='.'))
      {
         (*seq)[k++]=ch; 
      }
      else if(ch>='a' && ch<='z')
      {
         (*seq)[k++]=ch-'a'+'A';
      }
      if(k>=MAXSEQLEN) {printf("[Warning] Too long sequence\n");return false;}
   }
   (*seq)[k]='\0';
   return k;
}
void writefasta(FILE *fp,char *anno, char *seq)
{
   int j=0;
   if(anno)
      fprintf(fp,">%s\n",anno);
   else
      fprintf(fp,">\n");
   while(seq[j])
   {
      fprintf(fp,"%c",seq[j]);
      j++;
      if(j%70 == 0) fprintf(fp,"\n");
   }
   if(j%70) fprintf(fp,"\n");
}
real sigmoid(real x)
{
   return (1.0/(1.0+exp(-x)));
}
real linear(real x, real max, real min)
{
   real k=1.0/(max-min);
   real c=-min/(max-min);
   real y=k*x+c;
   return y;
}

/* return the index of maximum value (0 based) */
int maxind(real *vector, int size)
{
   real maxv;
   int i,j=2;
   maxv = -INF;
   for(i=0; i<size; i++)
   {
      if(vector[i]>maxv)
      {
         maxv=vector[i];
         j=i;
      }
   }
   return j;
}
/* 0 based */
int minind(real *vector, int size)
{
   real minv;
   int i,j=2;
   minv = INF;
   for(i=0; i<size; i++)
   {
      if(vector[i]<minv)
      {
         minv = vector[i];
         j=i;
      }
   }
   return j;
}

bool is_space(char ch)
{
   if(ch==' ' || ch=='\t' || ch=='\b') return true;
   return false;
}
char fpeek(FILE *fp)
{
   char ch;
   ch = getc(fp);
   ungetc(ch,fp);
   return ch;
}
void normalize(real *data, int num)
{
   int i;
   real sum=0;
   for(i=0; i<num; i++) sum+=data[i];
   for(i=0; i<num; i++) data[i]/=sum;
}
