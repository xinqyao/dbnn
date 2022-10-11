/************************* Neural Network for ***********************
 **************** Protein Secondary Structure Prediction ************ 
 *
 * Usage: nn N H M input weight [prediction]
 *     N: the number of units of the input layer
 *     H: the number of units of the hidden layer
 *     M: the number of units of the output layer
 *
 *     Input file should be in FLAT format: each line contains
 *     the profiles (within the window) and the secstr encoded 
 *     in 0-based for each residue site; sequences are separated
 *     by "//".
 *
 *     Weight file stores the weights learned from the training
 *     set; if no training data provided, the weight file will 
 *     be regard as input. The weights are saved in natural order
 *     as  arranged in the network.
 *
 *     Prediction file is also in FLAT format: each line contains
 *     the values of the output layer for each residue site; 
 *     sequences are searated by "//"; do not provide the prediction
 *     file if you only want to train the network.
 *
 * I give credit to K. Kutza, who wrote the Backpropagation Network
 * Simulator, the basis of this program.
 * 
 * Xinqiu Yao
 * Jan 13, 2006
 */

/******************************************************************************
                            D E C L A R A T I O N S
 ******************************************************************************/
#include "yxq.h"
#include "tools.h"

#define NOT           !
#define AND           &&
#define OR            ||

#define MIN_REAL      -INF
#define MAX_REAL      +INF
#define MIN(x,y)      ((x)<(y) ? (x) : (y))
#define MAX(x,y)      ((x)>(y) ? (x) : (y))

#define LO            0.1
#define HI            0.9
#define BIAS          1

#define sqr(x)        ((x)*(x))

typedef struct {                     /* A LAYER OF A NET:                     */
        int           Units;         /* - number of units in this layer       */
        real*         Output;        /* - output of ith unit                  */
        real*         Error;         /* - error term of ith unit              */
        real**        Weight;        /* - connection weights to ith unit      */
        real**        WeightSave;    /* - saved weights for stopped training  */
        real**        dWeight;       /* - last weight deltas for momentum     */
} LAYER;

typedef struct {                     /* A NET:                                */
        LAYER**       Layer;         /* - layers of this net                  */
        LAYER*        InputLayer;    /* - input layer                         */
        LAYER*        OutputLayer;   /* - output layer                        */
        real          Alpha;         /* - momentum factor                     */
        real          Eta;           /* - learning rate                       */
        real          Gain;          /* - gain of sigmoid function            */
        real          Error;         /* - total net error                     */
} NET;

/******************************************************************************
        R A N D O M S   D R A W N   F R O M   D I S T R I B U T I O N S
 ******************************************************************************/

void InitializeRandoms()
{
  srand(4711);
}

int RandomEqualINT(int Low, int High)
{
  return rand() % (High-Low+1) + Low;
}

real RandomEqualREAL(real Low, real High)
{
  return ((real) rand() / RAND_MAX) * (High-Low) + Low;
}

/******************************************************************************
               A P P L I C A T I O N - S P E C I F I C   C O D E
 ******************************************************************************/
#define MAXNTRAIN     20
#define NUM_LAYERS    3

int                   Units[NUM_LAYERS],
                      N,H,M,
                      TRAIN_LWB, TRAIN_UPB, TRAIN_NUM,
                      TEST_LWB, TEST_UPB, TEST_NUM,
                      EVAL_LWB, EVAL_UPB, EVAL_NUM,
                      NUM, LEN;
char                  Fin[256], Fout1[256], Fout2[256];
bool                  DEBUG=false;
real                 *Profile,
                      TestError;

/* parse arguments, count the number of sequences, 
   and set the training or evaluation boundaries */
int InitializeApplication(int npara, char ** para)
{
  int i,j;
  FILE *fp;
  bool bTrain=true;
  real *profileT;
 
  if(npara >= 8)
  {
    if(!strcmp(para[7],"-debug")) DEBUG = true;
  }
  N = atoi(para[1]);
  H = atoi(para[2]);
  M = atoi(para[3]);
  strcpy(Fin, para[4]);
  strcpy(Fout1,para[5]);
  if(npara>=7) strcpy(Fout2,para[6]);
  else strcpy(Fout2,"preds");
  Units[0] = N, Units[1] = H, Units[2] = M;
  NUM = 0;
 
  fp=fopen(Fin,"r");
  while(LEN=readplain(fp,&profileT)) 
  {
     /* if exists any missing data, then regards as prediction*/
     if(NUM==0) for(i=0; i<LEN; i++) if(profileT[i]<=-INF) {bTrain=false; break;}
     free(profileT);
     NUM+=LEN/(N+M);
  }
  fclose(fp);

  if(NUM==0)
  {
     printf("[Error] No data avaiable, please check you input file\n");
     return false;
  }
  if(bTrain) 
  {
     /* read the whole data in */
     if(NUM*(N+M)>=MAXSIZE)
     {
        printf("[Error] Too many data for training\n");
        return false;
     }
     Profile=(real*)malloc(NUM*(N+M)*sizeof(real));
     fp=fopen(Fin,"r");
     j=0;
     while(LEN=readplain(fp,&profileT))
     {
        for(i=0; i<LEN; i++) Profile[j+i]=profileT[i];
        j+=LEN;
        free(profileT);
     }
     fclose(fp);
     TRAIN_LWB = 0;
     TRAIN_UPB = NUM/7*6;
     TRAIN_NUM = TRAIN_UPB-TRAIN_LWB+1;
     TEST_LWB  = TRAIN_UPB+1;
     TEST_UPB  = NUM-1;
     TEST_NUM  = TEST_UPB-TEST_LWB+1;
     EVAL_NUM  = 0;
     EVAL_LWB  = 0;
     EVAL_UPB  = -1;
  }
  else 
  {
     EVAL_LWB  = 0;
     EVAL_UPB  = NUM-1;
     EVAL_NUM  = NUM; 
     TRAIN_NUM = 0;
     TRAIN_LWB = 0;
     TRAIN_UPB = -1;
     TEST_NUM  = 0;
     TEST_LWB  = 0;
     TEST_UPB  = -1;
  }
  if(DEBUG)
  {
    printf("\nN=%-6dH=%-6dM=%-6d\ninputs=%s\n\nweights=%s\npreds=%s\n\nTRAIN_LWB=%-8dTRAIN_UPB=%-8dTRAIN_NUM=%-8d\nTEST_LWB=%-8dTEST_UPB=%-8dTEST_NUM=%-8d\nEVAL_LWB=%-8dEVAL_UPB=%-8dEVAL_NUM=%-8d\nNUM=%-8d\n\n",N,H,M,Fin,Fout1,Fout2,TRAIN_LWB,TRAIN_UPB,TRAIN_NUM,TEST_LWB,TEST_UPB,TEST_NUM,EVAL_LWB,EVAL_UPB,EVAL_NUM,NUM);
      return false;
  }
  return true;
}

void WriteOut(NET * Net, FILE *fp)
{
  int i,j,k,l;
  if(EVAL_NUM > 0) 
  {
     for(j=0; j<LEN; j++)
     {
        for(k=0; k<M-1; k++) fprintf(fp,"%.3f,",Profile[j*(M+N)+N+k]);
        fprintf(fp,"%.3f\n",Profile[j*(M+N)+N+M-1]);
     }
     fprintf(fp,"//\n");
  }
  if(NUM-EVAL_NUM>0) 
  {
     fp=fopen(Fout1,"w");
     for (l=1; l<NUM_LAYERS; l++) {
        for (i=1; i<=Net->Layer[l]->Units; i++) {
           for (j=0; j<=Net->Layer[l-1]->Units; j++) {
              fprintf(fp,"%.3f\n",Net->Layer[l]->Weight[i][j]);
           }
        }
     }
     fclose(fp);
  }
}
  
/******************************************************************************
                          I N I T I A L I Z A T I O N
 ******************************************************************************/

void GenerateNetwork(NET* Net)
{
  int l,i;

  Net->Layer = (LAYER**) calloc(NUM_LAYERS, sizeof(LAYER*));

  for (l=0; l<NUM_LAYERS; l++) {
    Net->Layer[l] = (LAYER*) malloc(sizeof(LAYER));

    Net->Layer[l]->Units      = Units[l];
    Net->Layer[l]->Output     = (real*)  calloc(Units[l]+1, sizeof(real));
    Net->Layer[l]->Error      = (real*)  calloc(Units[l]+1, sizeof(real));
    Net->Layer[l]->Weight     = (real**) calloc(Units[l]+1, sizeof(real*));
    Net->Layer[l]->WeightSave = (real**) calloc(Units[l]+1, sizeof(real*));
    Net->Layer[l]->dWeight    = (real**) calloc(Units[l]+1, sizeof(real*));
    Net->Layer[l]->Output[0]  = BIAS;

    if (l != 0) {
      for (i=1; i<=Units[l]; i++) {
        Net->Layer[l]->Weight[i]     = (real*) calloc(Units[l-1]+1,
sizeof(real));
        Net->Layer[l]->WeightSave[i] = (real*) calloc(Units[l-1]+1,
sizeof(real));
        Net->Layer[l]->dWeight[i]    = (real*) calloc(Units[l-1]+1,
sizeof(real));
      }
    }
  }
  Net->InputLayer  = Net->Layer[0];
  Net->OutputLayer = Net->Layer[NUM_LAYERS - 1];
  Net->Alpha       = 0.9;
  Net->Eta         = 0.005;
  Net->Gain        = 1;
}

void RandomWeights(NET* Net)
{
  int l,i,j;

  for (l=1; l<NUM_LAYERS; l++) {
    for (i=1; i<=Net->Layer[l]->Units; i++) {
      for (j=0; j<=Net->Layer[l-1]->Units; j++) {
        Net->Layer[l]->Weight[i][j] = RandomEqualREAL(-0.5, 0.5);
      }
    }
  }
}

void SetInput(NET* Net, real* Input)
{
  int i;

  for (i=1; i<=Net->InputLayer->Units; i++) {
    Net->InputLayer->Output[i] = Input[i-1];
  }
}

void GetOutput(NET* Net, real* Output)
{
  int i;

  for (i=1; i<=Net->OutputLayer->Units; i++) {
    Output[i-1] = Net->OutputLayer->Output[i];
  }
}

/******************************************************************************
            S U P P O R T   F O R   S T O P P E D   T R A I N I N G
 ******************************************************************************/
void LoadWeights(NET* Net)
{
   int l,i,j;
   real w;
   FILE *fp=fopen(Fout1,"r");
   for(l=1; l<NUM_LAYERS; l++) {
      for(i=1;i<=Net->Layer[l]->Units; i++) {
         for(j=0; j<=Net->Layer[l-1]->Units; j++) {
#ifdef __DOUBLE_ACCURACY
            fscanf(fp,"%lf",&w);
#else
            fscanf(fp,"%f",&w);
#endif
            Net->Layer[l]->Weight[i][j] = w;
         }
      }
   }
   fclose(fp);
}

void SaveWeights(NET* Net)
{
  int l,i,j;

  for (l=1; l<NUM_LAYERS; l++) {
    for (i=1; i<=Net->Layer[l]->Units; i++) {
      for (j=0; j<=Net->Layer[l-1]->Units; j++) {
        Net->Layer[l]->WeightSave[i][j] = Net->Layer[l]->Weight[i][j];
      }
    }
  }
}

void RestoreWeights(NET* Net)
{
  int l,i,j;
  for (l=1; l<NUM_LAYERS; l++) {
    for (i=1; i<=Net->Layer[l]->Units; i++) {
      for (j=0; j<=Net->Layer[l-1]->Units; j++) {
        Net->Layer[l]->Weight[i][j] = Net->Layer[l]->WeightSave[i][j];
      }
    }
  }
}

/******************************************************************************
                     P R O P A G A T I N G   S I G N A L S
 ******************************************************************************/

void PropagateLayer(NET* Net, LAYER* Lower, LAYER* Upper)
{
  int  i,j;
  real Sum;

  for (i=1; i<=Upper->Units; i++) {
    Sum = 0;
    for (j=0; j<=Lower->Units; j++) {
      Sum += Upper->Weight[i][j] * Lower->Output[j];
    }
    Upper->Output[i] = 1 / (1 + exp(-Net->Gain * Sum));
  }
}

void PropagateNet(NET* Net)
{
  int l;

  for (l=0; l<NUM_LAYERS-1; l++) {
    PropagateLayer(Net, Net->Layer[l], Net->Layer[l+1]);
  }
}

/******************************************************************************
                  B A C K P R O P A G A T I N G   E R R O R S
 ******************************************************************************/

void ComputeOutputError(NET* Net, real* Target)
{
  int  i;
  real Out, Err;

  Net->Error = 0;
  for (i=1; i<=Net->OutputLayer->Units; i++) {
    Out = Net->OutputLayer->Output[i];
    Err = Target[i-1]-Out;
    Net->OutputLayer->Error[i] = Net->Gain * Out * (1-Out) * Err;
    Net->Error += 0.5 * sqr(Err);
  }
}

void BackpropagateLayer(NET* Net, LAYER* Upper, LAYER* Lower)
{
  int  i,j;
  real Out, Err;

  for (i=1; i<=Lower->Units; i++) {
    Out = Lower->Output[i];
    Err = 0;
    for (j=1; j<=Upper->Units; j++) {
      Err += Upper->Weight[j][i] * Upper->Error[j];
    }
    Lower->Error[i] = Net->Gain * Out * (1-Out) * Err;
  }
}

void BackpropagateNet(NET* Net)
{
  int l;

  for (l=NUM_LAYERS-1; l>1; l--) {
    BackpropagateLayer(Net, Net->Layer[l], Net->Layer[l-1]);
  }
}

void AdjustWeights(NET* Net)
{
  int  l,i,j;
  real Out, Err, dWeight;

  for (l=1; l<NUM_LAYERS; l++) {
    for (i=1; i<=Net->Layer[l]->Units; i++) {
      for (j=0; j<=Net->Layer[l-1]->Units; j++) {
        Out = Net->Layer[l-1]->Output[j];
        Err = Net->Layer[l]->Error[i];
        dWeight = Net->Layer[l]->dWeight[i][j];
        Net->Layer[l]->Weight[i][j] += Net->Eta * Err * Out + Net->Alpha *
dWeight;
        Net->Layer[l]->dWeight[i][j] = Net->Eta * Err * Out;
      }
    }
  }
}

/******************************************************************************
                      S I M U L A T I N G   T H E   N E T
 ******************************************************************************/

void SimulateNet(NET* Net, real* Input, real* Output, real* Target, bool
Training)
{
  SetInput(Net, Input);
  PropagateNet(Net);
  GetOutput(Net, Output);
  ComputeOutputError(Net, Target);
  if (Training) {
    BackpropagateNet(Net);
    AdjustWeights(Net);
  }
}

void TrainNet(NET* Net, int Epochs)
{
  real Output[M];
  int i, site;
  for(i=0; i<Epochs*TRAIN_NUM; i++) 
  {
     site=RandomEqualINT(TRAIN_LWB,TRAIN_UPB);
     SimulateNet(Net, Profile+site*(N+M), Output, Profile+site*(N+M)+N, true);
  }
}

void TestNet(NET* Net)
{
  real Output[M];
  int i;
  
  TestError=0;
  for(i=TEST_LWB; i<=TEST_UPB; i++)
  {
     SimulateNet(Net, Profile+i*(N+M), Output, Profile+i*(N+M)+N, false);
     TestError+=Net->Error;
  }
}

void EvaluateNet(NET* Net,real *Input)
{
    int i;
    for(i=0; i<LEN; i++)
       SimulateNet(Net, Input+i*(N+M), Input+i*(N+M)+N, Input+i*(N+M)+N, false);
}

/******************************************************************************
                                    M A I N
 ******************************************************************************/

int main(int argc, char** argv)
{
  if(argc < 6) 
  {
     printf("Usage: nn N H M input weight [prediction]\n");
     printf("\tN\tnumber of units in input layer\n");
     printf("\tH\tnumber of units in hidden layer\n");
     printf("\tM\tnumber of units in output layer\n");
     return 0;
  }
  
  NET  Net;
  bool bStop;	/* stop training */
  real MinTestError;
  int  i,NumTrain=0;
  FILE *fp, *fp2; 
 
  if(!InitializeApplication(argc,argv)) return 1;
  
  InitializeRandoms();
  GenerateNetwork(&Net);
  RandomWeights(&Net);
  
  bStop = false;
  if(NUM-EVAL_NUM==0)  /* No data for training */
  {
     fp=fopen(Fin,"r");
     fp2=fopen(Fout2,"w");
     LoadWeights(&Net);
     while(LEN=readplain(fp,&Profile)/(N+M))
     {
        EvaluateNet(&Net,Profile);
        WriteOut(&Net,fp2);
        free(Profile);
     }
     fclose(fp);
     fclose(fp2);
  }
  else
  {
    MinTestError = MAX_REAL;
    do {
      TrainNet(&Net,8);
      TestNet(&Net);
      NumTrain++;
      printf("TestError=%10.3f\tMinTestError=%10.3f\n",TestError,MinTestError);
      if (TestError < MinTestError) {
        MinTestError = TestError;
        printf(" - saving Weights ...\n");
        SaveWeights(&Net);
      }
      else if (TestError > 1.0* MinTestError) {
        bStop = true;
        printf(" - stopping Training and restoring Weights ...\n");
        RestoreWeights(&Net);
      }
      if( NOT bStop AND NumTrain>=MAXNTRAIN) {
        printf(" - stopping Training and restoring Weights ...\n");
        bStop=true;
        RestoreWeights(&Net);
      }
    } while (NOT bStop);
    WriteOut(&Net,NULL);
  }
  return 0;
}
