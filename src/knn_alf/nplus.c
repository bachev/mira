/*********************************************************
  nplus.c
  --------------------------------------------------------
  generated at Tue Feb 29 16:27:52 2000
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
  static UnitType Units[13];
  /* Sources definition section */
  static pUnit Sources[] =  {
Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 6, Units + 7, 
Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 6, Units + 7, 
Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 6, Units + 7, 
Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 6, Units + 7, 
Units + 8, Units + 9, Units + 10, Units + 11, 

  };

  /* Weigths definition section */
  static float Weights[] =  {
0.383750, -3.785780, 2.561380, -0.080500, 1.967470, 1.961440, -0.237310, 
-0.069180, -2.196420, 1.659500, 0.093840, 1.988610, 1.169960, -0.494460, 
0.005600, -3.731170, 1.729160, 0.602080, 2.120740, 1.771150, 0.536360, 
0.041810, 3.964270, -2.570720, -0.242730, -2.557900, -0.515340, -0.791130, 
4.330270, 2.795560, 3.835630, -5.290350, 

  };

  /* unit definition section (see also UnitType) */
  static UnitType Units[13] = 
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
      0.0, 0.191920, 7,
       &Sources[0] , 
       &Weights[0] , 
      },
    { /* unit 9 (Old: 9) */
      0.0, -0.027840, 7,
       &Sources[7] , 
       &Weights[7] , 
      },
    { /* unit 10 (Old: 10) */
      0.0, 1.090890, 7,
       &Sources[14] , 
       &Weights[14] , 
      },
    { /* unit 11 (Old: 11) */
      0.0, -0.679520, 7,
       &Sources[21] , 
       &Weights[21] , 
      },
    { /* unit 12 (Old: 12) */
      0.0, -2.739330, 4,
       &Sources[28] , 
       &Weights[28] , 
      }

  };



int nplus(float *in, float *out, int init)
{
  int member, source;
  float sum;
  enum{OK, Error, Not_Valid};
  pUnit unit;


  /* layer definition section (names & member units) */

  static pUnit Input[7] = {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 6, Units + 7}; /* members */

  static pUnit Hidden1[4] = {Units + 8, Units + 9, Units + 10, Units + 11}; /* members */

  static pUnit Output1[1] = {Units + 12}; /* members */

  static int Output[1] = {12};

  for(member = 0; member < 7; member++) {
    Input[member]->act = in[member];
  }

  for (member = 0; member < 4; member++) {
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
