// scan parameter space of model, returning costs of both encoding strategies after a simulated time window
// structures control different experiments (which parameters will be looped through, and how)

#include <stdio.h>
#include <math.h>

#include <stdlib.h>
#define RND drand48()

#define MAXT 144

// structure for model parameters
typedef struct
{
  double a, b, k, mu, lambda, num, nuc, diff;
  double initp;
  int noise;
} Params;

// structure for controlling loops through parameters in a given experiment
typedef struct
{
  // numbers of different values to simulate
  int na, nb, nk, nnum, nnuc, ndiff;
  // specific values to simulate
  double a[100], b[100], k[100], num[100], nuc[100], diff[100];
} Loops;

// structure for counters through the above loops
typedef struct
{
  int ac, bc, kc, numc, nucc, diffc;
} Counters;

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
  Loops E;
  Counters C;
  double t, dt = 0.01;
  int alpha;
  double cost, totalcost;
  int nsample;
  char fstr[100];
  int expt;
  int i;
  int costfn;
  int sample;

  // loop through experiments
  for(expt = 0; expt <= 11; expt++)
    {
      P.mu = 0; P.a = 1; P.b = 1; P.k = 1;
      P.lambda = 1; P.num = 1; P.nuc = 1; P.diff = 1;

      // this sets up the parameters for an experiment. the first is commented
      // originally we had costfn = 0 throughout and initp = 0.75 as the damage condition.
      // now we switch to costfn = 1, and everything favours NU a bit more. try initp = 0.9?
      switch(expt)
	{
	case 0:
	  // set 0 (static env)
	  // initialise the number of values for a to zero
	  // populate the list of values of a to simulate, incrementing counter each time
	  // here we just have one a value so it's pretty trivial...
	  E.na    = 0; for(i = 0; i < 1;  i++) { E.a[i] = 1;     E.na++; }
	  E.nb    = 0; for(i = 0; i < 1;  i++) { E.b[i] = 0;     E.nb++; }
	  E.nk    = 0; for(i = 0; i < 1;  i++) { E.k[i] = 0;     E.nk++; }
	  // but here we have a set of different values for nuc...
	  E.nnuc  = 0; for(i = 0; i < 20;  i++) { E.nuc[i] = 0.05*i;     E.nnuc++; }
	  E.nnum  = 0; for(i = 0; i < 20;  i++) { E.num[i] = 0.05*i;     E.nnum++; }
	  // and here we have a set of explicitly given diff values...
	  E.diff[0] = 0.01; E.diff[1] = 0.1; E.diff[2] = 1; E.diff[3] = 10; E.ndiff = 4;
	  P.initp = 1; costfn = 1; P.noise = 0;
	  break;
  
	case 1:
	  // set 1 (as set zero with damaged oDNA)
	  E.na    = 0; for(i = 0; i < 1;  i++) { E.a[i] = 1;     E.na++; }
	  E.nb    = 0; for(i = 0; i < 1;  i++) { E.b[i] = 0;     E.nb++; }
	  E.nk    = 0; for(i = 0; i < 1;  i++) { E.k[i] = 0;     E.nk++; }
	  E.nnuc  = 0; for(i = 0; i < 20;  i++) { E.nuc[i] = 0.05*i;     E.nnuc++; }
	  E.nnum  = 0; for(i = 0; i < 20;  i++) { E.num[i] = 0.05*i;     E.nnum++; }
	  E.diff[0] = 0.01; E.diff[1] = 0.1; E.diff[2] = 1; E.diff[3] = 10; E.ndiff = 4;
	  P.initp = 0.9; costfn = 1; P.noise = 0;
	  break;

	case 2:
	  // set 2 (changing env)
	  E.na    = 0; for(i = 0; i < 1;  i++) { E.a[i] = 1;     E.na++; }
	  E.nb    = 0; for(i = 0; i <= 10; i++) { E.b[i] = 0.1*i; E.nb++; }
	  E.nk    = 0; for(i = 0; i <= 20; i++) { E.k[i] = i;     E.nk++; }
	  E.nnuc  = 0; for(i = 0; i < 3;  i++) { E.nuc[i] = 0.1*i;     E.nnuc++; }
	  E.nnum  = 0; for(i = 0; i < 3;  i++) { E.num[i] = 0.1*i;     E.nnum++; }
	  E.diff[0] = 0.01; E.diff[1] = 0.1; E.diff[2] = 1; E.diff[3] = 10; E.ndiff = 4;
	  P.initp = 1; costfn = 1; P.noise = 0;
	  break;

	case 3:
	  // set 3 (as set 2 with damaged oDNA)
	  E.na    = 0; for(i = 0; i < 1;  i++) { E.a[i] = 1;     E.na++; }
	  E.nb    = 0; for(i = 0; i <= 10; i++) { E.b[i] = 0.1*i; E.nb++; }
	  E.nk    = 0; for(i = 0; i <= 20; i++) { E.k[i] = i;     E.nk++; }
	  E.nnuc  = 0; for(i = 0; i < 3;  i++) { E.nuc[i] = 0.1*i;     E.nnuc++; }
	  E.nnum  = 0; for(i = 0; i < 3;  i++) { E.num[i] = 0.1*i;     E.nnum++; }
	  E.diff[0] = 0.01; E.diff[1] = 0.1; E.diff[2] = 1; E.diff[3] = 10; E.ndiff = 4;
	  P.initp = 0.9; costfn = 1; P.noise = 0;

	case 4:
	  // set 4 (as set 3 with different cost function)
	  E.na    = 0; for(i = 0; i < 1;  i++) { E.a[i] = 1;     E.na++; }
	  E.nb    = 0; for(i = 0; i <= 10; i++) { E.b[i] = 0.1*i; E.nb++; }
	  E.nk    = 0; for(i = 0; i <= 20; i++) { E.k[i] = i;     E.nk++; }
	  E.nnuc  = 0; for(i = 0; i < 3;  i++) { E.nuc[i] = 0.1*i;     E.nnuc++; }
	  E.nnum  = 0; for(i = 0; i < 3;  i++) { E.num[i] = 0.1*i;     E.nnum++; }
	  E.diff[0] = 0.01; E.diff[1] = 0.1; E.diff[2] = 1; E.diff[3] = 10; E.ndiff = 4;
	  P.initp = 0.9; costfn = 0; P.noise = 0;
	  break;

	case 5:
	  // as set 3 -- different lambda
	  E.na    = 0; for(i = 0; i < 1;  i++) { E.a[i] = 1;     E.na++; }
	  E.nb    = 0; for(i = 0; i <= 10; i++) { E.b[i] = 0.1*i; E.nb++; }
	  E.nk    = 0; for(i = 0; i <= 20; i++) { E.k[i] = i;     E.nk++; }
	  E.nnuc  = 0; for(i = 0; i < 3;  i++) { E.nuc[i] = 0.1*i;     E.nnuc++; }
	  E.nnum  = 0; for(i = 0; i < 3;  i++) { E.num[i] = 0.1*i;     E.nnum++; }
	  E.diff[0] = 0.01; E.diff[1] = 0.1; E.diff[2] = 1; E.diff[3] = 10; E.ndiff = 4;
	  P.initp = 0.9; costfn = 1; P.noise = 0;
	  P.lambda = 10;
	  break;

        case 6:
	  // as set 3 -- different lambda
	  E.na    = 0; for(i = 0; i < 1;  i++) { E.a[i] = 1;     E.na++; }
	  E.nb    = 0; for(i = 0; i <= 10; i++) { E.b[i] = 0.1*i; E.nb++; }
	  E.nk    = 0; for(i = 0; i <= 20; i++) { E.k[i] = i;     E.nk++; }
	  E.nnuc  = 0; for(i = 0; i < 3;  i++) { E.nuc[i] = 0.1*i;     E.nnuc++; }
	  E.nnum  = 0; for(i = 0; i < 3;  i++) { E.num[i] = 0.1*i;     E.nnum++; }
	  E.diff[0] = 0.01; E.diff[1] = 0.1; E.diff[2] = 1; E.diff[3] = 10; E.ndiff = 4;
	  P.initp = 0.9; costfn = 1; P.noise = 0;
	  P.lambda = 0.1;
	  break;

 	case 7:
	  // white noise
	  E.na    = 0; for(i = 0; i < 1;  i++) { E.a[i] = 1;     E.na++; }
	  E.nb    = 0; for(i = 0; i < 1;  i++) { E.b[i] = 0;     E.nb++; }
	  E.nk    = 0; for(i = 0; i < 1;  i++) { E.k[i] = 0;     E.nk++; }
	  E.nnuc  = 0; for(i = 0; i < 10;  i++) { E.nuc[i] = 0.1*i;     E.nnuc++; }
	  E.nnum  = 0; for(i = 0; i < 10;  i++) { E.num[i] = 0.1*i;     E.nnum++; }
	  E.diff[0] = 0.01; E.diff[1] = 0.1; E.diff[2] = 1; E.diff[3] = 10; E.ndiff = 4;
	  P.initp = 0.9; costfn = 1; P.noise = 1;
	  break;

   	case 8:
	  // red noise
	  E.na    = 0; for(i = 0; i < 1;  i++) { E.a[i] = 1;     E.na++; }
	  E.nb    = 0; for(i = 0; i < 1;  i++) { E.b[i] = 0;     E.nb++; }
	  E.nk    = 0; for(i = 0; i < 1;  i++) { E.k[i] = 0;     E.nk++; }
	  E.nnuc  = 0; for(i = 0; i < 10;  i++) { E.nuc[i] = 0.1*i;     E.nnuc++; }
	  E.nnum  = 0; for(i = 0; i < 10;  i++) { E.num[i] = 0.1*i;     E.nnum++; }
	  E.diff[0] = 0.01; E.diff[1] = 0.1; E.diff[2] = 1; E.diff[3] = 10; E.ndiff = 4;
	  P.initp = 0.9; costfn = 1; P.noise = 2;
	  break;

  	case 9:
	  // set 9 (as set 3 with different again (quadratic) cost function)
	  E.na    = 0; for(i = 0; i < 1;  i++) { E.a[i] = 1;     E.na++; }
	  E.nb    = 0; for(i = 0; i <= 10; i++) { E.b[i] = 0.1*i; E.nb++; }
	  E.nk    = 0; for(i = 0; i <= 20; i++) { E.k[i] = i;     E.nk++; }
	  E.nnuc  = 0; for(i = 0; i < 3;  i++) { E.nuc[i] = 0.1*i;     E.nnuc++; }
	  E.nnum  = 0; for(i = 0; i < 3;  i++) { E.num[i] = 0.1*i;     E.nnum++; }
	  E.diff[0] = 0.01; E.diff[1] = 0.1; E.diff[2] = 1; E.diff[3] = 10; E.ndiff = 4;
	  P.initp = 0.9; costfn = 2; P.noise = 0;
	  break;

     	case 10:
	  // red noise small step
	  E.na    = 0; for(i = 0; i < 1;  i++) { E.a[i] = 1;     E.na++; }
	  E.nb    = 0; for(i = 0; i < 1;  i++) { E.b[i] = 0;     E.nb++; }
	  E.nk    = 0; for(i = 0; i < 1;  i++) { E.k[i] = 0;     E.nk++; }
	  E.nnuc  = 0; for(i = 0; i < 10;  i++) { E.nuc[i] = 0.1*i;     E.nnuc++; }
	  E.nnum  = 0; for(i = 0; i < 10;  i++) { E.num[i] = 0.1*i;     E.nnum++; }
	  E.diff[0] = 0.01; E.diff[1] = 0.1; E.diff[2] = 1; E.diff[3] = 10; E.ndiff = 4;
	  P.initp = 0.9; costfn = 1; P.noise = 3;
	  break;


    	case 11:
	  // red noise tiny step
	  E.na    = 0; for(i = 0; i < 1;  i++) { E.a[i] = 1;     E.na++; }
	  E.nb    = 0; for(i = 0; i < 1;  i++) { E.b[i] = 0;     E.nb++; }
	  E.nk    = 0; for(i = 0; i < 1;  i++) { E.k[i] = 0;     E.nk++; }
	  E.nnuc  = 0; for(i = 0; i < 10;  i++) { E.nuc[i] = 0.1*i;     E.nnuc++; }
	  E.nnum  = 0; for(i = 0; i < 10;  i++) { E.num[i] = 0.1*i;     E.nnum++; }
	  E.diff[0] = 0.01; E.diff[1] = 0.1; E.diff[2] = 1; E.diff[3] = 10; E.ndiff = 4;
	  P.initp = 0.9; costfn = 1; P.noise = 4;
	  break;


	}

      // open file for output
      sprintf(fstr, "costsscans-%i.csv", expt);
      fp = fopen(fstr, "w");
      fprintf(fp, "a,b,k,mu,diff,num,nuc,alpha,cost\n");
  
      // loop through total demand
      for(C.ac = 0; C.ac < E.na; C.ac++)
	{
	  P.a = E.a[C.ac];
	  // loop through environmental amplitude
	  for(C.bc = 0; C.bc < E.nb; C.bc++)
	    {
	      P.b = E.b[C.bc];
	      // loop through environmental frequency
	      for(C.kc = 0; C.kc < E.nk; C.kc++)
		{
		  P.k = E.k[C.kc];
		  // loop through transport
		  for(C.diffc = 0; C.diffc < E.ndiff; C.diffc++)
		    {
		      P.diff = E.diff[C.diffc];
		      // loop through degradation
		      for(C.numc = 0; C.numc <= E.nnum; C.numc++)
			{
			  P.num = E.num[C.numc];
			  for(C.nucc = 0; C.nucc < E.nnuc; C.nucc++)
			    {
			      P.nuc = E.nuc[C.nucc];
			      // loop through encoding strategy
			      for(alpha = 0; alpha <= 1; alpha++)
				{
				  // if we're using random environments, take 1000 samples, otherwise 1 (deterministic)
				  totalcost = 0;
				  if(P.noise == 0) nsample = 1;
				  else nsample = 1000;
				  
				  for(sample = 0; sample < nsample; sample++)
				    {
				      // initial conditions
				      xm = xc = 0; p = P.initp; cost = 0;
				      xmstar = 1;
				      // euler time loop
				      for(t = 0; t < MAXT*2; t += dt)
					{
					  // environmental demand -- different models
					  if(P.noise == 0) {  xmstar = P.a*(1+P.b*sin(t * P.k*2.*3.14159/MAXT)); }
					  if(P.noise == 1) { xmstar = RND*2; }
					  if(P.noise == 2) { xmstar += (RND-0.5)*0.1; if(xmstar < 0) xmstar = 0; if(xmstar > 2) xmstar = 2; }
					  if(P.noise == 3) { xmstar += (RND-0.5)*0.01; if(xmstar < 0) xmstar = 0; if(xmstar > 2) xmstar = 2; }
					  if(P.noise == 4) { xmstar += (RND-0.5)*0.001; if(xmstar < 0) xmstar = 0; if(xmstar > 2) xmstar = 2; };
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
					  if(t > MAXT)
					    {
					      if(costfn == 0) cost += (xm < xmstar ? xmstar-xm : 0);
					      else if(costfn == 1) cost += (xm < xmstar ? xmstar-xm : xm-xmstar);
      					      else if(costfn == 2) cost += (xmstar-xm)*(xmstar-xm);
					    }
					  // output to file
					}
				      totalcost += cost;
				    }
				  // output results to file and screen 
				  fprintf(fp, "%f,%f,%f,%f,%f,%f,%f,%i,%f\n", P.a, P.b, P.k, P.mu, P.diff, P.num, P.nuc, alpha, totalcost/nsample);
				  printf("Expt %i : %f %f %f %f %f %f %f %i\n", expt, P.a, P.b, P.k, P.mu, P.diff, P.num, P.nuc, alpha);

				}
			    } 
			} 
		    }
		}
	    }
	}
      
      fclose(fp);
    }
  
  return 0;
}
