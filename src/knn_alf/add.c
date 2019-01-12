/*********************************************************
  add.c
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
-0.618930, -0.189080, -0.794380, 0.009860, -0.169890, -0.007730, 0.062230, 0.337980, -0.890290, 0.856050, 
-0.732770, -0.951400, -0.644090, 
3.038520, 0.513190, -2.829950, -0.453670, 1.527360, -0.368030, 1.822900, -1.213060, -0.609130, 0.671600, 
-0.957500, -0.702360, 0.099150, 
-2.142120, -0.086490, 0.832910, -1.062650, -0.253880, 0.212890, -1.218620, 0.655910, 0.805110, 0.439000, 
0.718680, 0.457260, 0.790070, 
2.667970, -0.083000, -2.680620, -0.825920, 0.236360, -0.073970, 1.453210, -0.132850, -0.036460, 0.325540, 
-0.032890, -0.116150, 0.058830, 
2.670230, -0.156380, -1.282520, 0.064340, 1.516800, 0.015360, 0.424450, -0.046210, -0.695470, -0.508280, 
-0.527550, -0.400790, -0.867030, 
-2.120440, 0.939510, 1.186680, 0.673390, -0.654040, 0.309290, -1.979610, 1.123930, -0.305430, -0.262850, 
0.419640, -0.821440, -0.161700, 
-2.115660, 0.517020, 3.336870, -0.838710, -1.784470, 1.051590, -1.677760, 0.821500, -1.404430, -0.704730, 
-0.058540, -0.544560, 0.435740, 
0.127140, 4.886230, -2.728880, 3.674980, 3.079150, -2.975730, -5.054220, 

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
      0.0, 0.925570, 13,
       &Sources[0] , 
       &Weights[0] , 
      },
    { /* unit 15 (Old: 15) */
      0.0, 0.502870, 13,
       &Sources[13] , 
       &Weights[13] , 
      },
    { /* unit 16 (Old: 16) */
      0.0, 0.429900, 13,
       &Sources[26] , 
       &Weights[26] , 
      },
    { /* unit 17 (Old: 17) */
      0.0, 0.044410, 13,
       &Sources[39] , 
       &Weights[39] , 
      },
    { /* unit 18 (Old: 18) */
      0.0, -0.241020, 13,
       &Sources[52] , 
       &Weights[52] , 
      },
    { /* unit 19 (Old: 19) */
      0.0, 0.345880, 13,
       &Sources[65] , 
       &Weights[65] , 
      },
    { /* unit 20 (Old: 20) */
      0.0, 0.731500, 13,
       &Sources[78] , 
       &Weights[78] , 
      },
    { /* unit 21 (Old: 21) */
      0.0, -0.195490, 7,
       &Sources[91] , 
       &Weights[91] , 
      }

  };



int add(float *in, float *out, int init)
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
