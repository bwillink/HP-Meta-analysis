/* file food_gradient.c */
#include <R.h>
#include <math.h>
static double parms[10];

#define r parms[0]
#define K parms[1]
#define fM parms[2]
#define hr parms[3]
#define ht parms[4]
#define e parms[5]
#define rho parms[6]
#define BM parms[7]
#define d parms[8]
#define v parms[9]


/* initializer */
void initmod(void (* odeparms)(int *, double *))
{
int N=10;
odeparms(&N, parms);
}

/* Derivatives */
void derivs (int *neq, double *t, double *y, double *ydot,
double *yout, int *ip)
{


double S = y[0];
double I = y[1];
double R = y[2];

ydot[0] = e*fM*R/(hr + R)*(S + rho*I) - BM*R/(ht + R)*S*I - d*S;
ydot[1] =  BM*R/(ht + R)*S*I - (d + v)*I;
ydot[2] = r*R*(1 - R/K) - fM*R/(hr + R)*(S+I);


}

/* END file food_gradient.c */ 
