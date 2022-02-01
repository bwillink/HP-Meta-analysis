/* file env_gradient.c */
#include <R.h>
#include <math.h>
static double parms[9];

#define b parms[0]
#define rho parms[1]
#define K parms[2]
#define B parms[3]
#define BE parms[4]
#define d parms[5]
#define dE parms[6]
#define ENV parms[7]
#define v parms[8]



/* initializer */
void initmod(void (* odeparms)(int *, double *))
{
int N=9;
odeparms(&N, parms);
}

/* Derivatives and 2 output variables */
void derivs (int *neq, double *t, double *y, double *ydot,
double *yout, int *ip)
{
if (ip[0] <1) error("nout should be at least 1");

double S = y[0];
double I = y[1];

ydot[0] = b*(S + rho*I)*(1 - (S+I)/K) - B*(1 + BE*ENV)*S*I - d*(1 + dE*ENV)*S;
ydot[1] =  B*(1 + BE*ENV)*S*I - (d + v)*(1 + dE*ENV)*I;

yout[0] = 0;

}

/* END file env_gradient.c */ 
