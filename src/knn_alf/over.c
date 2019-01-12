/*********************************************************
  over.c
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
-0.918370, 0.370610, 1.419830, 0.291790, -0.305740, 0.457350, -1.196380, -0.175770, -2.659160, 0.733330, 
0.310130, -3.870330, -1.625070, 
0.915850, 0.127310, -0.645950, -0.631970, 1.637260, 0.471990, 0.687560, 0.017850, -1.299840, 0.266150, 
-0.362100, -0.473790, -0.285090, 
-3.364300, 0.211250, 0.091130, -1.333520, 1.250740, 2.097380, -0.773040, 0.390560, 3.853620, 1.748080, 
1.149280, 0.170110, 0.911580, 
6.011590, -1.431990, -5.931900, -1.080570, -0.754630, 1.945610, 2.263600, 3.397350, 0.691530, 1.166110, 
-0.297710, 1.999240, -1.013310, 
4.052380, -2.356590, -1.414310, -1.499220, 4.025910, 2.307680, 0.387890, 1.207810, 0.986170, -1.624860, 
-0.144260, 1.310330, -0.974570, 
-0.749680, 0.732430, 0.560220, 0.714560, -0.222730, -1.934050, -1.542080, 4.318600, -3.202300, 0.713290, 
-0.229660, -0.556860, -1.880610, 
-4.263560, -0.677420, 3.185130, 0.011430, -2.323510, 0.208840, -1.990960, 2.176030, 2.153380, -2.045750, 
-2.845090, -1.417130, 2.885740, 
-5.088340, 0.941440, -4.849660, 7.412120, 6.326800, -3.988370, -6.462060, 

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
      0.0, 2.274920, 13,
       &Sources[0] , 
       &Weights[0] , 
      },
    { /* unit 15 (Old: 15) */
      0.0, 0.058120, 13,
       &Sources[13] , 
       &Weights[13] , 
      },
    { /* unit 16 (Old: 16) */
      0.0, 0.127750, 13,
       &Sources[26] , 
       &Weights[26] , 
      },
    { /* unit 17 (Old: 17) */
      0.0, -0.708540, 13,
       &Sources[39] , 
       &Weights[39] , 
      },
    { /* unit 18 (Old: 18) */
      0.0, -1.144540, 13,
       &Sources[52] , 
       &Weights[52] , 
      },
    { /* unit 19 (Old: 19) */
      0.0, 1.236020, 13,
       &Sources[65] , 
       &Weights[65] , 
      },
    { /* unit 20 (Old: 20) */
      0.0, 0.280990, 13,
       &Sources[78] , 
       &Weights[78] , 
      },
    { /* unit 21 (Old: 21) */
      0.0, 0.303350, 7,
       &Sources[91] , 
       &Weights[91] , 
      }

  };



int over(float *in, float *out, int init)
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
