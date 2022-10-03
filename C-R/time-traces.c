// simulate time behaviour for given parameterisations

#include <stdio.h>
#include <math.h>

#define MAXT 144

#include <stdlib.h>
#define RND drand48()

// structure for model parameters
typedef struct
{
  double a, b, k, mu, lambda, num, nuc, diff;
  int noise;
  double initp;
} Params;

// expression signal function, takes environmental demand and current supply and returns desired amount of gene expression
double m(double demand, double supply)
{
  if(supply > demand) return 0;
  return (demand-supply);
}

int main(void)
{
  double xm, dxm, xc, dxc, p, dp;
  double xmstar;
  FILE *fp;
  Params P;
  Params Q[20];
  int nexpt;
  int i;
  double t, dt = 0.01;
  int alpha;
  double cost;
  double nu;
  int noise;
  int costfn = 1;
  
  nexpt = 0;

  // this is just a set of specific experiments (parameterisations) to simulate
  // nexpt, the experiment counter, is incremented each time
  
  // static env, p = 1, high nuc low diff
  Q[nexpt].initp = 1; Q[nexpt].lambda = 1; Q[nexpt].num = 0.1; Q[nexpt].nuc = 1; Q[nexpt].diff = 0.1;
  Q[nexpt].noise = 0; Q[nexpt].a = 1; Q[nexpt].b = 0; Q[nexpt].k = 0; nexpt++;

  // static env, p = 1, low nuc high diff
  Q[nexpt].initp = 1; Q[nexpt].lambda = 1; Q[nexpt].num = 0.1; Q[nexpt].nuc = 0; Q[nexpt].diff = 10;
  Q[nexpt].noise = 0; Q[nexpt].a = 1; Q[nexpt].b = 0; Q[nexpt].k = 0; nexpt++;

  // static env, p = 0.9, low nuc high diff
  Q[nexpt].initp = 0.9; Q[nexpt].lambda = 1; Q[nexpt].num = 0.1; Q[nexpt].nuc = 0; Q[nexpt].diff = 10;
  Q[nexpt].noise = 0; Q[nexpt].a = 1; Q[nexpt].b = 0; Q[nexpt].k = 0; nexpt++;

  // dyn env, low diff, k = 1
  Q[nexpt].initp = 0.9; Q[nexpt].lambda = 1; Q[nexpt].num = 0.1; Q[nexpt].nuc = 0; Q[nexpt].diff = 0.1;
  Q[nexpt].noise = 0; Q[nexpt].a = 1; Q[nexpt].b = 1; Q[nexpt].k = 1; nexpt++;

  // dyn env, low diff, k = 4
  Q[nexpt].initp = 0.9; Q[nexpt].lambda = 1; Q[nexpt].num = 0.1; Q[nexpt].nuc = 0; Q[nexpt].diff = 0.1;
  Q[nexpt].noise = 0; Q[nexpt].a = 1; Q[nexpt].b = 1; Q[nexpt].k = 4; nexpt++;

  // dyn env, high diff, k = 1
  Q[nexpt].initp = 0.9; Q[nexpt].lambda = 1; Q[nexpt].num = 0.1; Q[nexpt].nuc = 0; Q[nexpt].diff = 1;
  Q[nexpt].noise = 0; Q[nexpt].a = 1; Q[nexpt].b = 1; Q[nexpt].k = 1; nexpt++;

  // dyn env, high diff, k = 16
  Q[nexpt].initp = 0.9; Q[nexpt].lambda = 1; Q[nexpt].num = 0.1; Q[nexpt].nuc = 0; Q[nexpt].diff = 1;
  Q[nexpt].noise = 0; Q[nexpt].a = 1; Q[nexpt].b = 1; Q[nexpt].k = 16; nexpt++;

  // white noise, low nuc
  Q[nexpt].initp = 0.9; Q[nexpt].lambda = 1; Q[nexpt].num = 0.1; Q[nexpt].nuc = 0.1; Q[nexpt].diff = 1;
  Q[nexpt].noise = 1; Q[nexpt].a = 1; Q[nexpt].b = 0; Q[nexpt].k = 0; nexpt++;

  // white noise, high nuc
  Q[nexpt].initp = 0.9; Q[nexpt].lambda = 1; Q[nexpt].num = 0.1; Q[nexpt].nuc = 1; Q[nexpt].diff = 1;
  Q[nexpt].noise = 1;  Q[nexpt].a = 1; Q[nexpt].b = 0;   Q[nexpt].k = 0; nexpt++;

  // red noise, low nuc
  Q[nexpt].initp = 0.9; Q[nexpt].lambda = 1; Q[nexpt].num = 0.1; Q[nexpt].nuc = 0.1; Q[nexpt].diff = 1;
  Q[nexpt].noise = 2; Q[nexpt].a = 1; Q[nexpt].b = 0; Q[nexpt].k = 0; nexpt++;

  // red noise, high nuc
  Q[nexpt].initp = 0.9; Q[nexpt].lambda = 1; Q[nexpt].num = 0.1; Q[nexpt].nuc = 1; Q[nexpt].diff = 1;
  Q[nexpt].noise = 2;  Q[nexpt].a = 1; Q[nexpt].b = 0;   Q[nexpt].k = 0; nexpt++;

  // red noise, low nuc, small step
  Q[nexpt].initp = 0.9; Q[nexpt].lambda = 1; Q[nexpt].num = 0.1; Q[nexpt].nuc = 0.1; Q[nexpt].diff = 1;
  Q[nexpt].noise = 3; Q[nexpt].a = 1; Q[nexpt].b = 0; Q[nexpt].k = 0; nexpt++;

  // red noise, low nuc, tiny step
  Q[nexpt].initp = 0.9; Q[nexpt].lambda = 1; Q[nexpt].num = 0.1; Q[nexpt].nuc = 0.1; Q[nexpt].diff = 1;
  Q[nexpt].noise = 4; Q[nexpt].a = 1; Q[nexpt].b = 0; Q[nexpt].k = 0; nexpt++;

  // open file for output
  fp = fopen("traces-gen.csv", "w");
  fprintf(fp, "expt,a,k,mu,diff,num,nuc,alpha,t,xmstar,xm,xc,p,cost\n");

  for(i = 0; i < nexpt; i++)
    {
      P.num = Q[i].num; P.nuc = Q[i].nuc; P.diff = Q[i].diff; P.lambda = Q[i].lambda;
      P.a = Q[i].a; P.b = Q[i].b; P.k = Q[i].k;
      P.mu = 0; P.noise = Q[i].noise; P.initp = Q[i].initp;
      
      printf("%f %f %f\n", P.a, P.k, P.mu);
      // loop through encoding strategy
      for(alpha = 0; alpha <= 1; alpha++)
	{
	  // initial conditions
	  xm = xc = 0; p = P.initp; cost = 0;
	  // euler time loop
	  for(t = 0; t < MAXT*2; t += dt)
	    {
	      // environmental demand
	      if(P.noise == 0) {  xmstar = P.a*(1+P.b*sin(t * P.k*2.*3.14159/MAXT)); }
	      if(P.noise == 1) { xmstar = RND*2; }
	      if(P.noise == 2) { xmstar += (RND-0.5)*0.1; if(xmstar < 0) xmstar = 0; if(xmstar > 2) xmstar = 2; }
	      if(P.noise == 3) { xmstar += (RND-0.5)*0.01; if(xmstar < 0) xmstar = 0; if(xmstar > 2) xmstar = 2; }
	      if(P.noise == 4) { xmstar += (RND-0.5)*0.001; if(xmstar < 0) xmstar = 0; if(xmstar > 2) xmstar = 2; }
	      // equations of motion depending on encoding strategy
	      if(alpha == 0)
		{
		  dxm = dt* ( P.lambda*m(xmstar, xm)*p - P.num*xm );
		  dxc = dt* ( 0 );
		  dp  = dt* ( -P.mu*p );
		}
	      else
		{
		  dxm = dt* ( P.diff*xc - P.num*xm );
		  dxc = dt* ( P.lambda*m(xmstar, xm) - P.diff*xc - P.nuc*xc );
		  dp  = dt* ( -P.mu*p );
		}
	      // update state and compute cost
	      xm += dxm; xc += dxc; p += dp;
	      	      // output to file
	      if(t > MAXT)
		{
	      if(costfn == 0) cost += (xm < xmstar ? xmstar-xm : 0);
	      else if(costfn == 1) cost += (xm < xmstar ? xmstar-xm : xm-xmstar);
	      else if(costfn == 2) cost += (xmstar-xm)*(xmstar-xm);

		  fprintf(fp, "%i,%f,%f,%f,%f,%f,%f,%i,%f,%f,%f,%f,%f,%f\n", i, P.a, P.k, P.mu, P.diff, P.num, P.nuc, alpha, t-MAXT, xmstar, xm, xc, p, cost);
		}
	    }
	}
    }
  
  fclose(fp);

  
  return 0;
}
