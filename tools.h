#ifndef __TOOLS_H
#define __TOOLS_H
#include "yxq.h"
bool is_space(char);
char fpeek(FILE *);
int readpssm(FILE*,real**);
real readdata(FILE*); /* used by readplain */
int readplain(FILE*,real**);
typedef struct {
   int start;      /* 0 based */
   int end;        /* 0 based */
   char type;
} sSegment;
int transtr(char);
int tranaa(char);
char getaa(int);
char getstr(int);
int  readfasta(FILE *,char **,char **);
void writefasta(FILE *,char *,char *);
real sigmoid(real);
real linear(real,real,real);
int minind(real*,int);
int maxind(real*,int);
int split(char*,sSegment***);
void normalize(real *, int);
#endif
