/*********************************************************
  missing.h
  --------------------------------------------------------
  generated at Tue Feb 29 16:27:52 2000
  by snns2c ( Bernward Kett 1995 ) 
*********************************************************/

extern int missing(float *in, float *out, int init);

static struct {
  int NoOfInput;    /* Number of Input Units  */
  int NoOfOutput;   /* Number of Output Units */
  int(* propFunc)(float *, float*, int);
} missingREC = {14,1,missing};
