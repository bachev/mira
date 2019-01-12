/*********************************************************
  add.c
  --------------------------------------------------------
  generated at Tue Feb 22 14:20:39 2000
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
-0.672380, -0.105020, -0.664100, 0.027140, -0.186750, -0.003710, 0.025020, 0.536790, -0.861780, 0.964860, 
-0.687690, -0.895390, -0.675340, 
4.974400, 0.945640, -2.006370, -0.005400, 1.085460, -0.200780, 3.171800, -2.240540, -0.498900, -0.803070, 
0.160160, -2.488470, -0.579550, 
-2.233090, -0.500910, -0.163610, -1.228330, 0.595270, 0.152740, -1.206640, 1.373060, 0.281470, 0.462890, 
0.048750, 1.315130, 0.747070, 
3.221040, -0.246640, -2.066770, -0.695040, -0.214200, -0.034200, 1.909330, -0.213160, -0.175820, 0.252840, 
0.250910, -1.456660, -0.532660, 
2.286180, -0.096540, -0.110770, 0.126850, 0.945610, 0.126490, 0.230370, 0.514270, -0.585760, -1.667200, 
0.100310, -1.311160, -0.885550, 
-1.270310, 0.314090, -0.631960, 0.540150, -0.131940, -0.197690, -1.637230, 2.394210, -0.514470, 1.214290, 
0.282420, -0.621710, -0.729490, 
-2.098290, 0.609410, 2.450490, -0.883100, -0.878350, 1.470870, -1.695120, -0.293530, -0.751880, -1.275710, 
-0.965290, 1.689720, 2.344350, 
0.513810, 7.012270, -3.013900, 4.056260, 2.679840, -2.671080, -4.449240, 

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
      0.0, 0.960020, 13,
       &Sources[0] , 
       &Weights[0] , 
      },
    { /* unit 15 (Old: 15) */
      0.0, 1.399300, 13,
       &Sources[13] , 
       &Weights[13] , 
      },
    { /* unit 16 (Old: 16) */
      0.0, 0.098240, 13,
       &Sources[26] , 
       &Weights[26] , 
      },
    { /* unit 17 (Old: 17) */
      0.0, 0.306210, 13,
       &Sources[39] , 
       &Weights[39] , 
      },
    { /* unit 18 (Old: 18) */
      0.0, -0.115990, 13,
       &Sources[52] , 
       &Weights[52] , 
      },
    { /* unit 19 (Old: 19) */
      0.0, 0.079310, 13,
       &Sources[65] , 
       &Weights[65] , 
      },
    { /* unit 20 (Old: 20) */
      0.0, 0.642690, 13,
       &Sources[78] , 
       &Weights[78] , 
      },
    { /* unit 21 (Old: 21) */
      0.0, -0.171910, 7,
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
