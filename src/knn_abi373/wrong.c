/*********************************************************
  wrong.c
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
0.821270, 0.347730, 0.951210, 0.551740, -4.569540, -3.122010, -2.151940, 1.162300, -0.200580, 
0.101760, 1.291510, -0.537760, 0.579040, -1.891680, -2.147570, -1.018570, 0.626710, 0.417660, 
0.314540, -0.419920, 0.827890, 0.502720, 0.689150, -0.366510, -0.380200, 0.439680, -0.586430, 
-0.508890, 0.122640, 0.612790, 0.650820, 2.615490, -2.414350, 4.559760, -1.095140, 2.415640, 
0.945880, 0.178610, 0.147690, 0.138760, 2.212590, 0.391190, 1.831230, -0.247900, -1.589140, 
-5.603230, -2.497710, -0.042000, 4.901480, 1.988140, 

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
      0.0, 1.216290, 9,
       &Sources[0] , 
       &Weights[0] , 
      },
    { /* unit 11 (Old: 11) */
      0.0, 0.125010, 9,
       &Sources[9] , 
       &Weights[9] , 
      },
    { /* unit 12 (Old: 12) */
      0.0, 0.658340, 9,
       &Sources[18] , 
       &Weights[18] , 
      },
    { /* unit 13 (Old: 13) */
      0.0, -0.461420, 9,
       &Sources[27] , 
       &Weights[27] , 
      },
    { /* unit 14 (Old: 14) */
      0.0, 0.005310, 9,
       &Sources[36] , 
       &Weights[36] , 
      },
    { /* unit 15 (Old: 15) */
      0.0, -1.649980, 5,
       &Sources[45] , 
       &Weights[45] , 
      }

  };



int wrong(float *in, float *out, int init)
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
