/*********************************************************
  ncall.c
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
  static UnitType Units[16];
  /* Sources definition section */
  static pUnit Sources[] =  {
Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 6, Units + 7, Units + 8, Units + 9, 
Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 6, Units + 7, Units + 8, Units + 9, 
Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 6, Units + 7, Units + 8, Units + 9, 
Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 6, Units + 7, Units + 8, Units + 9, 
Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 6, Units + 7, Units + 8, Units + 9, 
Units + 10, Units + 11, Units + 12, Units + 13, Units + 14, 

  };

  /* Weigths definition section */
  static float Weights[] =  {
0.075160, 0.289060, 1.797810, 2.813980, -4.205520, 1.215330, -0.564600, -0.650220, -5.272390, 
-0.476600, -0.110200, -1.722730, 0.000330, -1.553920, -3.613230, -0.213780, 0.125970, 5.025960, 
0.378910, 0.193040, 0.897440, 0.460310, 1.169750, 0.594330, -3.628080, 0.546550, 1.518100, 
0.200730, -0.338270, -1.073700, -1.008920, 3.622840, 0.364810, 3.394980, 0.468640, 1.995030, 
0.086310, 0.127610, 0.348910, 0.460080, 0.401440, -0.783350, 1.333270, -0.157220, -0.957470, 
-6.503880, 4.887620, 3.858670, 4.254280, -1.907580, 

  };

  /* unit definition section (see also UnitType) */
  static UnitType Units[16] = 
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
      0.0, 3.062330, 9,
       &Sources[0] , 
       &Weights[0] , 
      },
    { /* unit 11 (Old: 11) */
      0.0, -1.751130, 9,
       &Sources[9] , 
       &Weights[9] , 
      },
    { /* unit 12 (Old: 12) */
      0.0, 0.822930, 9,
       &Sources[18] , 
       &Weights[18] , 
      },
    { /* unit 13 (Old: 13) */
      0.0, -2.552420, 9,
       &Sources[27] , 
       &Weights[27] , 
      },
    { /* unit 14 (Old: 14) */
      0.0, -0.124220, 9,
       &Sources[36] , 
       &Weights[36] , 
      },
    { /* unit 15 (Old: 15) */
      0.0, -3.718900, 5,
       &Sources[45] , 
       &Weights[45] , 
      }

  };



int ncall(float *in, float *out, int init)
{
  int member, source;
  float sum;
  enum{OK, Error, Not_Valid};
  pUnit unit;


  /* layer definition section (names & member units) */

  static pUnit Input[9] = {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 6, Units + 7, Units + 8, Units + 9}; /* members */

  static pUnit Hidden1[5] = {Units + 10, Units + 11, Units + 12, Units + 13, Units + 14}; /* members */

  static pUnit Output1[1] = {Units + 15}; /* members */

  static int Output[1] = {15};

  for(member = 0; member < 9; member++) {
    Input[member]->act = in[member];
  }

  for (member = 0; member < 5; member++) {
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
