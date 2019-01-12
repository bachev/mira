/*********************************************************
  under.c
  --------------------------------------------------------
  generated at Tue Feb 22 14:20:40 2000
  by snns2c ( Bernward Kett 1995 ) 
*********************************************************/

#include <math.h>

#define Act_Logistic(sum, bias)  ( (sum+bias<10000.0) ? ( 1.0/(1.0 + exp(-sum-bias) ) ) : 0.0 )
#define NULL (void *)0

typedef struct UT {
          float act;         /* Activation       */
          float Bias;        /* Bias of the Unit */
          int   NoOfSources; /* Number of predecessor units */
   struct UT   **sources; /* predecessor units */
          float *weights; /* weights from predecessor units */
        } UnitType, *pUnit;

  /* Forward Declaration for all unit types */
  static UnitType Units[22];
  /* Sources definition section */
  static pUnit Sources[] =  {
Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 6, Units + 7, Units + 8, Units + 9, Units + 10, 
Units + 11, Units + 12, Units + 13, 
Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 6, Units + 7, Units + 8, Units + 9, Units + 10, 
Units + 11, Units + 12, Units + 13, 
Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 6, Units + 7, Units + 8, Units + 9, Units + 10, 
Units + 11, Units + 12, Units + 13, 
Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 6, Units + 7, Units + 8, Units + 9, Units + 10, 
Units + 11, Units + 12, Units + 13, 
Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 6, Units + 7, Units + 8, Units + 9, Units + 10, 
Units + 11, Units + 12, Units + 13, 
Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 6, Units + 7, Units + 8, Units + 9, Units + 10, 
Units + 11, Units + 12, Units + 13, 
Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 6, Units + 7, Units + 8, Units + 9, Units + 10, 
Units + 11, Units + 12, Units + 13, 
Units + 14, Units + 15, Units + 16, Units + 17, Units + 18, Units + 19, Units + 20, 

  };

  /* Weigths definition section */
  static float Weights[] =  {
0.553260, -1.862410, -0.419620, -1.422950, -1.062190, 0.021970, -0.587490, 1.960180, -1.057200, 2.469110, 
-1.020990, -1.980180, 0.205530, 
-1.571930, 2.284010, -0.423110, 1.452060, 1.396540, 0.134870, 1.147750, -0.689870, -1.250520, -0.525170, 
-0.561510, 0.284070, 0.308750, 
2.502770, -2.599970, -0.599820, -3.754130, -0.685920, -1.750990, -1.216380, 1.042120, -0.621250, 0.114690, 
0.135980, -0.473180, 0.088790, 
0.474500, 0.119810, -0.734830, -0.773710, -0.232180, 0.792410, 0.616230, 0.481070, -0.143960, -0.276020, 
-0.072460, -0.024980, 0.468570, 
-1.551220, 1.893560, 0.632360, 1.763980, 1.350610, 1.449080, 0.044910, -0.723790, 0.377790, -0.103080, 
-0.082100, 0.311310, -0.187220, 
0.358510, -0.225730, -1.043730, -0.108010, -0.462990, -0.750880, -1.351200, 0.834550, -0.261140, -0.031780, 
0.195460, -1.134300, 0.049470, 
0.697110, -0.726470, 0.241840, -1.310910, -1.031490, 1.039500, -0.450670, 0.940720, -1.617490, -0.517930, 
-0.647590, -0.334740, 0.169390, 
-3.681990, 3.893300, -5.278480, -0.080400, 3.726740, -0.805490, -1.940030, 

  };

  /* unit definition section (see also UnitType) */
  static UnitType Units[22] = 
  {
    { 0.0, 0.0, 0, NULL , NULL },
    { /* unit 1 (Old: 1) */
      0.0, -0.489050, 0,
       &Sources[0] , 
       &Weights[0] , 
      },
    { /* unit 2 (Old: 2) */
      0.0, 0.314600, 0,
       &Sources[0] , 
       &Weights[0] , 
      },
    { /* unit 3 (Old: 3) */
      0.0, 0.777830, 0,
       &Sources[0] , 
       &Weights[0] , 
      },
    { /* unit 4 (Old: 4) */
      0.0, 0.155320, 0,
       &Sources[0] , 
       &Weights[0] , 
      },
    { /* unit 5 (Old: 5) */
      0.0, 0.837780, 0,
       &Sources[0] , 
       &Weights[0] , 
      },
    { /* unit 6 (Old: 6) */
      0.0, -0.712880, 0,
       &Sources[0] , 
       &Weights[0] , 
      },
    { /* unit 7 (Old: 7) */
      0.0, 0.002400, 0,
       &Sources[0] , 
       &Weights[0] , 
      },
    { /* unit 8 (Old: 8) */
      0.0, 0.042370, 0,
       &Sources[0] , 
       &Weights[0] , 
      },
    { /* unit 9 (Old: 9) */
      0.0, -0.381620, 0,
       &Sources[0] , 
       &Weights[0] , 
      },
    { /* unit 10 (Old: 10) */
      0.0, -0.669190, 0,
       &Sources[0] , 
       &Weights[0] , 
      },
    { /* unit 11 (Old: 11) */
      0.0, 0.074200, 0,
       &Sources[0] , 
       &Weights[0] , 
      },
    { /* unit 12 (Old: 12) */
      0.0, -0.155280, 0,
       &Sources[0] , 
       &Weights[0] , 
      },
    { /* unit 13 (Old: 13) */
      0.0, -0.067950, 0,
       &Sources[0] , 
       &Weights[0] , 
      },
    { /* unit 14 (Old: 14) */
      0.0, 1.118810, 13,
       &Sources[0] , 
       &Weights[0] , 
      },
    { /* unit 15 (Old: 15) */
      0.0, -0.335500, 13,
       &Sources[13] , 
       &Weights[13] , 
      },
    { /* unit 16 (Old: 16) */
      0.0, 1.441040, 13,
       &Sources[26] , 
       &Weights[26] , 
      },
    { /* unit 17 (Old: 17) */
      0.0, -0.264910, 13,
       &Sources[39] , 
       &Weights[39] , 
      },
    { /* unit 18 (Old: 18) */
      0.0, -1.040220, 13,
       &Sources[52] , 
       &Weights[52] , 
      },
    { /* unit 19 (Old: 19) */
      0.0, 0.488900, 13,
       &Sources[65] , 
       &Weights[65] , 
      },
    { /* unit 20 (Old: 20) */
      0.0, 0.774590, 13,
       &Sources[78] , 
       &Weights[78] , 
      },
    { /* unit 21 (Old: 21) */
      0.0, 0.673250, 7,
       &Sources[91] , 
       &Weights[91] , 
      }

  };



int under(float *in, float *out, int init)
{
  int member, source;
  float sum;
  enum{OK, Error, Not_Valid};
  pUnit unit;


  /* layer definition section (names & member units) */

  static pUnit Input[13] = {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 12, Units + 13}; /* members */

  static pUnit Hidden1[7] = {Units + 14, Units + 15, Units + 16, Units + 17, Units + 18, Units + 19, Units + 20}; /* members */

  static pUnit Output1[1] = {Units + 21}; /* members */

  static int Output[1] = {21};

  for(member = 0; member < 13; member++) {
    Input[member]->act = in[member];
  }

  for (member = 0; member < 7; member++) {
    unit = Hidden1[member];
    sum = 0.0;
    for (source = 0; source < unit->NoOfSources; source++) {
      sum += unit->sources[source]->act
             * unit->weights[source];
    }
    unit->act = Act_Logistic(sum, unit->Bias);
  };

  for (member = 0; member < 1; member++) {
    unit = Output1[member];
    sum = 0.0;
    for (source = 0; source < unit->NoOfSources; source++) {
      sum += unit->sources[source]->act
             * unit->weights[source];
    }
    unit->act = Act_Logistic(sum, unit->Bias);
  };

  for(member = 0; member < 1; member++) {
    out[member] = Units[Output[member]].act;
  }

  return(OK);
}
