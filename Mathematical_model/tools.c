/***************************************************************************/
/*                                                                         */
/*  Program: tools.C                                                       */
/*  Version: 2.0                                                           */
/*  By: Mette Olufsen, Math-Tech                                           */
/*  Date: 14. Jan. 1997                                                    */
/*                                                                         */
/*  This module includes auxiliary functions such as error-handling,       */
/*                      a version of the Numerical Recipes in C routine    */
/*  for solving a nonlinear set of equations using Newton Raphson's met-   */
/*  hod, and finally one for solving a one-dimensional nonlinear equation. */
/*  Since these functions are taken directly from Numerical Recipes they   */
/*  will not be commented further.                                         */
/*                                                                         */
/*  The program tools uses a number of vector and matrix handling functions*/
/*  that are available in the library nrutil, therefore I have included    */
/*  nrutil.h. Apart from that only standard c-libraries stdio.h, stdlib.h, */
/*  and math.h are used.                                                   */
/*                                                                         */
/***************************************************************************/

// $Id: tools.C,v 1.5 2010-10-20 15:06:06 mette Exp $

#include "tools.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <algorithm>

using namespace std;

// Standard parameters.
const double EPS = 1.0e-20;

// The general error handling of this program prints out an error message
// and terminates the program.
void error (const char *s, const char *u)
{
  fprintf (stdout, "ERROR (%s): %s\n", s, u);
  //exit (0);
}

// The remaining functions are from the Numerical Recipes in C library.
// They are all used with Newton Raphson's method when solving a set of
// non-linear equations, and when solving a one-dimensional nonlinear equation,
// also using Newton Raphson's method.

void lubksb (double *a[], int n, int indx[], double b[])
{
  int ii = 0;

  for (int i=0; i<n; i++)
  {
    int ip = indx[i];
    double sum = b[ip-1];
    b[ip-1] = b[i];
    if (ii != 0) {
      for (int j=ii-1; j<i; j++)
      {
        sum = sum - a[i][j]*b[j];
      }
    }
    else if (sum) ii=i+1;
    b[i] = sum;
  }
  for (int i=n-1; i>=0; i--)
  {
    double sum = b[i];
    for (int j=i+1; j<n; j++)
    {
      sum = sum - a[i][j]*b[j];
    }
    b[i] = sum / a[i][i];
  }
}


void ludcmp (double *a[], int n, int indx[], double *d)
{
  int imax = -1;
  double vv[n];
  *d   = 1.0;
    
    //Debugging:MJC
//      fprintf(stdout,"A VALUES\n");
//      for (int i2=0; i2<n; i2++)
//      {
//          for (int j2=0; j2<n; j2++)
//          {
//              fprintf(stdout,"%lf ",a[i2][j2]);
//          }
//          fprintf(stdout,"\n");
//      }
//    fprintf(stdout,"\n\n");
    
  for (int i=0; i<n; i++)
  {
    double big = 0.0;
    for (int j=0; j<n; j++)
    {
	big = max(big, fabs(a[i][j]));
    }
    if (big == 0.0) {
      const char *error_text = "Singular matrix in routine LUDCMP";
      fprintf (stdout, "Numerical Recipes run-time error...\n");
      fprintf (stdout, "%s\n",error_text);
      fprintf (stdout, "...now exiting to system...\n");
        //Debugging:MJC
        fprintf(stdout,"A VALUES\n");
        for (int i=0; i<n; i++)
        {
            for (int j=0; j<n; j++)
            {
                fprintf(stdout,"%lf ",a[i][j]);
            }
            fprintf(stdout,"\n");
        }
      fprintf(stdout,"\n\n");
      exit (1);
    }
    vv[i] = 1.0/big;
  }

  for (int j=0; j<n; j++)
  {
    for (int i=0; i<j; i++)
    {
      double sum = a[i][j];
      for (int k=0; k<i; k++)
      {
        sum = sum - a[i][k]*a[k][j];
      }
      a[i][j] = sum;
    }

    double big = 0.0;
    for (int i=j; i<n; i++)
    {
      double sum = a[i][j];
      for (int k=0; k<j; k++)
      {
	sum = sum - a[i][k]*a[k][j];
      }
      a[i][j] = sum;

      double dum = vv[i]*fabs(sum);
      if (dum >= big) {
	big  = dum;
	imax = i;
      }
    }

    if (j != imax) {
      for (int k=0; k<n; k++)
      {
	swap(a[imax][k], a[j][k]);
      }
      *d       = -(*d);
      vv[imax] = vv[j];
    }
    indx[j] = imax+1;
    if (a[j][j] == 0.0) a[j][j] = EPS;
    if (j+1 != n) {
      double dum = 1.0 / (a[j][j]);
      for (int i=j+1; i<n; i++)
      {
        a[i][j] = a[i][j]*dum;
      }
    }
  }
}


// Takes one step with Newton Raphson's method. Assumes to find zeros
// for a multi-dimensional problem.
int zero (double* x, int n, double tolx, double tolf,
          double fvec[], double *fjac[])
{

  int indx[n];
  double p[n];
  double errf = 0.0;
  for (int i=0; i<n; i++)
  {
    errf = errf + fabs(fvec[i]);
  }
  if (errf <= tolf)
  {
    return (1);
  }

  for (int i=0; i<n; i++) {
    p[i] = - fvec[i];
      // MJC: For large gradient, take smaller steps
//      if (p[i]>1.0 || p[i]<-1.0) {
//          fprintf(stdout,"step: %lf",p[i]);
//          p[i] = p[i]/fabs(p[i]);
//      }
  }
//    fprintf(stdout,"\n");
  double d; // unused return value
  ludcmp (fjac, n, indx, &d);
  lubksb (fjac, n, indx, p);
  double errx = 0.0;
  for (int i=0; i<n; i++)
  {
    errx = errx + fabs(p[i]);
    x[i] = x[i] + p[i];
  }
  if (errx <= tolx)
  {
    return (1);
  }
  return (2);
}


// Takes one step with Newton Raphson's method. Assumes to find zeros for
// a one-dimensional problem.
bool zero_1d (double *x, double f, double df, double tolx)
{
  double dx  = f/df;
  *x  = *x-dx;
  //if (fabs(dx) < tolx) return(1); else   // Original statement.
  return (fabs(dx) < tolx && fabs(f) < tolx);
}


void linesearch (bool armijo_flag, int n, double *xb, double *xold, double fvec[], double *fjac[], double *pln, double SSE_old, double slope)
{
//    // ADD NEWTON PARAMETERS w/ LINE SEARCH HERE
//        double tolx = 1e-8;
//        double tol_res = 1e-5;
//        double tol_min = 1e-8;
//        double stpmaxMAX = 100.0;
//        double alam, alamin, a, b, den, disc, f2, rhs1, rhs2, stpmax, sum, temp, test, tmplam;
//        double alam2 = 0.0;
////        double fvec[n], pln[n], grad[n], fold[n];
//        double ALF = 1.0e-4;
//        double SSE_new;
//    
//    test=0.0;
//    for (int i=0; i<n; i++) {
//        temp=abs(pln[i])/fmax(abs(xold[i]),1.0);
//        if (temp > test) test=temp;
//    }
//    
//    alamin=tolx/test;
//    alam=1.0;
//    f2 = 0.0;
//
//    
//        // Calculate SSE
//        SSE_new = 0.0;
//        for (int i=0; i<n; i++)
//        {
//            SSE_new+=fvec[i]*fvec[i];
//            if (isnan(fvec[i])) fprintf(stdout,"i: %d      fvec %lf\n",i,fvec[i]);
//        }
//        SSE_new*=0.5;
//        
//        if (alam<alamin) {
//            for (int i=0; i<n; i++) xb[i] = xold[i];
//            armijo_flag = true;
//        }
//        else if (SSE_new <= SSE_old +ALF*alam*slope)
//        {
//            armijo_flag=true;
//        }
//        else
//        {
//            if (alam == 1.0)
//            {
//                tmplam = -slope/(2.0*(SSE_new-SSE_old-slope));
//            }
//            else
//            {
//                rhs1 = SSE_new-SSE_old-alam*slope;
//                rhs2 = f2-SSE_old-alam2*slope;
//                a    = (rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
//                b    = (-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
//                if (a == 0.0) tmplam = -slope/(2.0*b);
//                else {
//                    disc=b*b-3.0*a*slope;
//                    if (disc < 0.0) tmplam=0.9*alam;
//                    else if (b <= 0.0) tmplam=(-b+sqrt(disc))/(3.0*a);
//                    else tmplam=-slope/(b+sqrt(disc));
//                    }
//                    if (tmplam>0.9*alam)
//                        tmplam=0.9*alam;
//            }
//            if (isnan(tmplam))
//             {
//                 printf(" -  tmplam not a number\n");
////                     for (int i=0; i<12; i++) fprintf(stdout,"    %lf    ",fvec[i]);
////                     abort();
//             }
//            alam2=alam;
//            f2 = SSE_new;
//            alam=fmax(tmplam,0.1*alam);
//        }
//    fprintf(stdout,"Armijo: %d",armijo_flag);
    
    }



