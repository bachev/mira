/*********************************************************
  over.c
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
-2.741110, 1.612550, 0.072220, -2.176600, -1.328590, -0.173110, -0.838530, -1.464620, -2.374660, 0.752570, 
1.035680, -4.081110, -4.421090, 
1.279820, 0.714940, 0.188530, -0.114850, 0.831470, 1.272000, 0.642760, -0.609680, 0.154040, -0.424740, 
-1.289250, -0.169390, 1.875400, 
0.635110, -0.223790, -1.144410, -0.142940, 0.624480, -0.047670, 0.097170, 0.496940, -0.239030, 0.096900, 
0.711610, 1.794680, -0.902940, 
3.426980, -1.518050, -3.700240, 1.552130, 0.070020, 3.895670, 1.627850, 0.568500, -0.625040, 0.650830, 
1.830140, 4.092140, -0.703900, 
2.080610, -1.479060, -1.120990, 1.407860, 1.297010, -0.358070, -0.229900, 1.225030, 0.218510, 0.106970, 
0.149620, 2.196170, -1.392770, 
-1.665130, 1.297160, -0.468820, -0.695840, 0.958910, 0.120550, -1.660320, 2.721870, 0.967870, 0.750830, 
2.654700, -4.732220, 0.049560, 
-2.213550, 0.843870, 3.702150, -2.971100, -0.736910, -0.177500, -1.545490, 1.149680, -3.489310, -0.644770, 
-0.902890, -4.929760, -0.187020, 
-7.344120, 2.232830, 2.247220, 7.348330, 4.070420, -6.094320, -8.280490, 

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
      0.0, 1.321120, 13,
       &Sources[0] , 
       &Weights[0] , 
      },
    { /* unit 15 (Old: 15) */
      0.0, 0.483940, 13,
       &Sources[13] , 
       &Weights[13] , 
      },
    { /* unit 16 (Old: 16) */
      0.0, 0.761950, 13,
       &Sources[26] , 
       &Weights[26] , 
      },
    { /* unit 17 (Old: 17) */
      0.0, -1.351390, 13,
       &Sources[39] , 
       &Weights[39] , 
      },
    { /* unit 18 (Old: 18) */
      0.0, -1.215800, 13,
       &Sources[52] , 
       &Weights[52] , 
      },
    { /* unit 19 (Old: 19) */
      0.0, 0.959950, 13,
       &Sources[65] , 
       &Weights[65] , 
      },
    { /* unit 20 (Old: 20) */
      0.0, 2.362970, 13,
       &Sources[78] , 
       &Weights[78] , 
      },
    { /* unit 21 (Old: 21) */
      0.0, 0.337930, 7,
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
