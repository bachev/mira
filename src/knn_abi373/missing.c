/*********************************************************
  missing.c
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
-0.879130, -0.178880, 0.202790, 0.053380, 0.399030, 0.820630, 0.827890, -0.949970, -0.189170, -1.172030, 
-0.474610, 0.573670, -0.085110, 1.781750, 
0.706790, -0.645120, 0.895300, 0.418910, 0.455110, -0.674400, -0.908570, 0.849720, -0.218110, -0.064390, 
0.445130, 0.340180, -0.196660, -0.358650, 
-2.539240, 1.638450, -0.254560, -0.333820, 0.991600, 0.337970, 0.944570, 1.358890, 1.423090, -0.048890, 
0.260500, 2.955450, 0.195520, 2.134110, 
1.812730, -0.545800, 0.843810, 1.472980, -0.751050, -0.279610, -0.991680, 0.264490, 0.394910, -0.563050, 
0.045960, -2.766530, 0.073620, -3.544490, 
0.324910, -0.676920, 0.696690, -0.615960, -0.343840, 0.738860, -0.481600, -0.087610, 0.062830, -0.405430, 
0.838180, -0.469810, 1.337330, 0.372410, 
-1.301730, 0.745030, -0.742470, -0.291550, 0.506000, -0.538930, -0.175510, 0.728360, 0.505400, -0.231330, 
0.072880, -0.748380, -0.417540, 0.661810, 
0.952460, -0.523330, -0.728740, -0.747490, -0.313570, -0.382590, -0.581780, 0.106290, 0.536330, -0.972230, 
0.570760, 0.482180, 0.144510, -0.747070, 
-1.725080, 0.740740, -0.450920, 0.089550, 0.448510, 0.226260, 1.352380, -0.158780, 0.383830, -0.325020, 
0.471630, 1.679200, 0.635480, 2.876520, 
2.355600, -1.511380, 4.810840, -5.275760, 0.310820, 0.394520, -1.042470, 3.578630, 

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
      0.0, -0.631700, 14,
       &Sources[0] , 
       &Weights[0] , 
      },
    { /* unit 16 (Old: 16) */
      0.0, 0.888890, 14,
       &Sources[14] , 
       &Weights[14] , 
      },
    { /* unit 17 (Old: 17) */
      0.0, -1.164350, 14,
       &Sources[28] , 
       &Weights[28] , 
      },
    { /* unit 18 (Old: 18) */
      0.0, -0.155440, 14,
       &Sources[42] , 
       &Weights[42] , 
      },
    { /* unit 19 (Old: 19) */
      0.0, 0.676200, 14,
       &Sources[56] , 
       &Weights[56] , 
      },
    { /* unit 20 (Old: 20) */
      0.0, -0.496060, 14,
       &Sources[70] , 
       &Weights[70] , 
      },
    { /* unit 21 (Old: 21) */
      0.0, 0.000710, 14,
       &Sources[84] , 
       &Weights[84] , 
      },
    { /* unit 22 (Old: 22) */
      0.0, 0.266120, 14,
       &Sources[98] , 
       &Weights[98] , 
      },
    { /* unit 23 (Old: 23) */
      0.0, -1.084810, 8,
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
