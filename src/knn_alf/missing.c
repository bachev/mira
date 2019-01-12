/*********************************************************
  missing.c
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
  static UnitType Units[24];
  /* Sources definition section */
  static pUnit Sources[] =  {
Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 6, Units + 7, Units + 8, Units + 9, Units + 10, 
Units + 11, Units + 12, Units + 13, Units + 14, 
Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 6, Units + 7, Units + 8, Units + 9, Units + 10, 
Units + 11, Units + 12, Units + 13, Units + 14, 
Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 6, Units + 7, Units + 8, Units + 9, Units + 10, 
Units + 11, Units + 12, Units + 13, Units + 14, 
Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 6, Units + 7, Units + 8, Units + 9, Units + 10, 
Units + 11, Units + 12, Units + 13, Units + 14, 
Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 6, Units + 7, Units + 8, Units + 9, Units + 10, 
Units + 11, Units + 12, Units + 13, Units + 14, 
Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 6, Units + 7, Units + 8, Units + 9, Units + 10, 
Units + 11, Units + 12, Units + 13, Units + 14, 
Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 6, Units + 7, Units + 8, Units + 9, Units + 10, 
Units + 11, Units + 12, Units + 13, Units + 14, 
Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 6, Units + 7, Units + 8, Units + 9, Units + 10, 
Units + 11, Units + 12, Units + 13, Units + 14, 
Units + 15, Units + 16, Units + 17, Units + 18, Units + 19, Units + 20, Units + 21, Units + 22, 

  };

  /* Weigths definition section */
  static float Weights[] =  {
-0.775630, -0.307030, 0.901090, 1.059520, 0.464540, 2.556430, 0.778640, -0.515940, 1.113650, -0.365000, 
0.442630, 0.402400, 0.250000, 1.656440, 
1.751540, -0.941550, 0.289640, -1.163610, 0.530390, -2.488400, -1.327890, 0.230870, -0.332220, -0.560510, 
-0.820210, -0.222230, -0.122980, -1.876830, 
-2.313840, 0.778810, -0.242500, 1.313870, 0.570160, 1.063080, 0.568190, 0.466490, 1.388440, 1.172740, 
1.038820, 1.453360, 1.017650, 0.589880, 
0.737940, 0.445220, 0.118930, 0.260750, -0.384090, -2.409260, -0.409700, 0.109810, -0.228920, -0.618760, 
-0.571430, -0.720400, 0.655780, -0.991180, 
1.335600, -0.946470, 0.647550, -1.144130, -0.461660, -0.119990, -0.872790, -0.436990, 0.763900, -0.094560, 
0.763040, -1.064020, 0.694070, -0.948490, 
-1.370480, 0.912410, -0.487610, -0.094780, 0.573910, -0.296660, -0.105720, 0.799720, 0.210590, -0.473000, 
0.020200, -0.805560, -0.731080, 0.846070, 
1.297920, -0.843400, -0.680370, -1.605840, -0.299790, 0.423510, -0.685570, -0.474020, 1.371400, -1.193960, 
0.417590, 0.694080, -1.058760, -1.100260, 
-1.032550, 0.179860, -0.212590, 1.056550, 0.238420, 0.712240, 0.993500, 0.220890, 1.039400, -0.637820, 
0.498250, 0.604910, 0.665120, 1.932210, 
3.253960, -3.963880, 3.573960, -2.533980, -2.308270, 1.057980, -2.778150, 2.272780, 

  };

  /* unit definition section (see also UnitType) */
  static UnitType Units[24] = 
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
      0.0, 0.948220, 0,
       &Sources[0] , 
       &Weights[0] , 
      },
    { /* unit 15 (Old: 15) */
      0.0, -0.346610, 14,
       &Sources[0] , 
       &Weights[0] , 
      },
    { /* unit 16 (Old: 16) */
      0.0, 1.026970, 14,
       &Sources[14] , 
       &Weights[14] , 
      },
    { /* unit 17 (Old: 17) */
      0.0, -1.325400, 14,
       &Sources[28] , 
       &Weights[28] , 
      },
    { /* unit 18 (Old: 18) */
      0.0, -0.536750, 14,
       &Sources[42] , 
       &Weights[42] , 
      },
    { /* unit 19 (Old: 19) */
      0.0, 1.015300, 14,
       &Sources[56] , 
       &Weights[56] , 
      },
    { /* unit 20 (Old: 20) */
      0.0, -0.407300, 14,
       &Sources[70] , 
       &Weights[70] , 
      },
    { /* unit 21 (Old: 21) */
      0.0, 0.036180, 14,
       &Sources[84] , 
       &Weights[84] , 
      },
    { /* unit 22 (Old: 22) */
      0.0, 0.475610, 14,
       &Sources[98] , 
       &Weights[98] , 
      },
    { /* unit 23 (Old: 23) */
      0.0, -1.036240, 8,
       &Sources[112] , 
       &Weights[112] , 
      }

  };



int missing(float *in, float *out, int init)
{
  int member, source;
  float sum;
  enum{OK, Error, Not_Valid};
  pUnit unit;


  /* layer definition section (names & member units) */

  static pUnit Input[14] = {Units + 1, Units + 2, Units + 3, Units + 4, Units + 5, Units + 6, Units + 7, Units + 8, Units + 9, Units + 10, Units + 11, Units + 12, Units + 13, Units + 14}; /* members */

  static pUnit Hidden1[8] = {Units + 15, Units + 16, Units + 17, Units + 18, Units + 19, Units + 20, Units + 21, Units + 22}; /* members */

  static pUnit Output1[1] = {Units + 23}; /* members */

  static int Output[1] = {23};

  for(member = 0; member < 14; member++) {
    Input[member]->act = in[member];
  }

  for (member = 0; member < 8; member++) {
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
