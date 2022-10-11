#ifndef __YXQ_H
#define __YXQ_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define PI       3.14159265358979323846
#define true     1
#define false    0
/*#define __DOUBLE_ACCURACY*/
#ifdef __DOUBLE_ACCURACY
typedef double   real;
#else
typedef float    real;
#endif
typedef int      bool;
#define INF          HUGE_VAL
#define MAXSEQLEN    2500
#define MAXSIZEARRAY 800000   /* limit for DBN input 
			         data of one sequence */
#define MAXSIZE      27600000 /* limit for NN input data */
#endif
