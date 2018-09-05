
#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>
#include "CodeMLoptimizer.h"
#include "CompilerHints.h"

#ifdef USE_LAPACK
#include "blas.h"
#endif


inline static double min2(double a, double b) {return (a < b) ? a : b;}
inline static double max2(double a, double b) {return (a > b) ? a : b;}
inline static double square(double a) {return a*a;}
inline static void   zero(double x[], int n) {memset(x, 0, n*sizeof(double));}
inline static void   xtoy(const double x[], double y[], int n) {memcpy(y, x, n*sizeof(double));}

inline static void identityMatrix(double* x, int n) 
{
	memset(x, 0, n*n*sizeof(double));
	for(int i = 0; i < n; i++) x[i*n + i] = 1.0;
}

static inline double norm(const double x[], int n) 
{
#ifdef USE_LAPACK
	return dnrm2_(&n, x, &I1);
#else
	double t = 0;

	for(int i = 0; i < n; ++i) t += square(x[i]);

	return sqrt(t);
#endif
}

static inline double innerp(const double x[], const double y[], int n) 
{
#ifdef USE_LAPACK
	return ddot_(&n, x, &I1, y, &I1);
#else
	double t = 0.;

	for(int i=0; i < n; ++i) t += x[i] * y[i];

	return t;
#endif
}

static inline double distance(const double* RESTRICT x, const double* RESTRICT y, int n) 
{
	double t = 0;

	for(int i = 0; i < n; ++i) t += square(x[i] - y[i]);

	return sqrt(t);
}


double Ming2::minimizeFunction(std::vector<double>& aVars)
{
	int np = static_cast<int>(aVars.size());
	std::vector<double> space(np*(np*2+9+2));
	std::vector<int> ispace(2*np);

	bool sy = std::cout.sync_with_stdio(true);

	//mAlwaysCenter = false;
	double lnL = 0;
	int sts = ming2(mTrace ? stdout : NULL, &lnL, &aVars[0], &mLowerBound[0], &mUpperBound[0], &space[0], &ispace[0], mRelativeError, np);
	if(sts < 0 && mVerbose > 0) std::cout << "Check ming2 convergence" << std::endl;
	std::cout.sync_with_stdio(sy);

	return -lnL;
}


int Ming2::ming2(FILE *fout, double *f,	double x[], const double xl[], const double xu[], double space[], int ispace[], double rel_error, int n)
{
    /* n-variate minimization with bounds using the BFGS algorithm
         g0[n] g[n] p[n] x0[n] y[n] s[n] z[n] H[n*n] C[n*n] tv[2*n]
         xmark[n],ix[n]
       Size of space should be (check carefully?)
          #define spaceming2(n) ((n)*((n)*2+9+2)*sizeof(double))
       nfree: # free variables
       xmark[i]=0 for inside space; -1 for lower boundary; 1 for upper boundary.
       x[] has initial values at input and returns the estimates in return.
       ix[i] specifies the i-th free parameter

    */
    int i, j, i1, i2, it, fail = 0, nfree;
    int Ngoodtimes = 2, goodtimes = 0;
    double small = 1.e-30, sizep0 = 0;	/* small value for checking |w|=0 */
    double w, v, am, h, maxstep = 8;
	int iround;

	// Sanity check.
    if(n == 0) return 0;

	// Prepare the temporary variables
    double* g0 = space;
    double* g  = g0 + n;
    double* p  = g + n;
    double* x0 = p + n;
    double* y  = x0 + n;
    double* s  = y + n;
    double* z  = s + n;
    double* H  = z + n;
    double* C  = H + n * n;
	double* tv = C + n * n;
    //int* xmark = (int *) (tv + 2 * n);
	int* xmark = ispace;
    int* ix    = ispace + n;

	// Initialize the two integer workspaces
    //for(i = 0; i < n; i++)
    //{
    //    xmark[i] = 0;
    //    ix[i] = i;
    //}
	memset(ispace, 0, 2*n*sizeof(int));

    for(i = 0, nfree = 0; i < n; i++)
    {
        if(x[i] <= xl[i])
        {
            x[i] = xl[i];
            xmark[i] = -1;
            //continue;
        }
        else if(x[i] >= xu[i])
        {
            x[i] = xu[i];
            xmark[i] = 1;
           // continue;
        }
		else
			ix[nfree++] = i;
    }

    if(mNoisy > 2 && nfree < n && n < 50)
    {
        printf("\n");
        for(j=0; j < n; ++j) printf(" %9.6f", x[j]);
        printf("\n\n");
        for(j=0; j < n; ++j) printf(" %9.6f", xl[j]);
        printf("\n\n");
        for(j=0; j < n; ++j) printf(" %9.6f", xu[j]);
        printf("\n");

        if(nfree < n && mNoisy >= 3)
        {
            printf("warning: ming2, %d params at boundary.", n - nfree);
        }
    }

    //double f0 = *f = (*fun) (x, n);	++mNumFunCall;
	double f0 = *f = -mModel->computeLikelihood(x, n, mTraceFun);
    xtoy(x, x0, n);
    double sizep = 99.;

    if(mNoisy > 2)
    {
        printf("\nInitial: fx = %12.6f\nx =", f0);
        for(i=0; i < n; ++i) printf(" %8.6f", x[i]);
        printf("\n");
    }

    gradientB(n, x0, f0, g0, tv, xmark, sizep);

    identityMatrix(H, nfree);

    for(iround = 0; iround < mMaxIterations; ++iround)
    {
		// Check if the optimization can be stopped in advance due to LRT non satisfied
		if(mStopIfBigger && *f < mThreshold) throw FastCodeMLEarlyStopLRT();

        if(fout)
        {
            fprintf(fout, "\n%3d %7.4f %13.6f  x: ", iround, sizep0, f0);
            for(j=0; j < n; ++j) fprintf(fout, "%8.6f  ", x0[j]);
			fprintf(fout, "\n");
            fflush(fout);
        }

        for(i = 0, zero(p, n); i < nfree; i++)
            for(j=0; j < nfree; ++j)
				p[ix[i]] -= H[i * nfree + j] * g0[ix[j]];

        sizep0 = sizep;
        sizep = norm(p, n);	/* check this */

        for (i = 0, am = maxstep; i < n; i++)  	/* max step length */
        {
            if (p[i] > 0 && (xu[i] - x0[i]) / p[i] < am)
            {
                am = (xu[i] - x0[i]) / p[i];
            }
            else if (p[i] < 0 && (xl[i] - x0[i]) / p[i] < am)
            {
                am = (xl[i] - x0[i]) / p[i];
            }
        }

        if(iround == 0)
        {
			// First round
            h = fabs(2 * f0 * .01 / innerp(g0, p, n));	/* check this?? */
            h = min2(h, am / 2000.0);

        }
        else
        {
            h = norm(s, nfree) / sizep;
            h = max2(h, am / 500.0);
        }

        h = max2(h, 1e-5);
        h = min2(h, am / 5.0);
        *f = f0;
        double alpha = LineSearch2(f, x0, p, h, am, min2(1e-3, rel_error), tv, iround, n);	/* n or nfree? */

        if(alpha <= 0)
        {
            if(fail)
            {
                if(mAlwaysCenter)
                {
                    iround = mMaxIterations;
                    break;
                }
                else
                {
                    mAlwaysCenter = true;
                    identityMatrix(H, n);
                    fail = 1;
                }
            }
            else
            {
                if(mNoisy > 2)
                {
                    printf(".. ");
                }

                identityMatrix(H, nfree);
                fail = 1;
            }
        }
        else
        {
            fail = 0;
            for(i=0; i < n; ++i) x[i] = x0[i] + alpha * p[i];
            w = min2(2., rel_error * 1000.);

            if(rel_error < 1e-4 && rel_error > 1e-6)
            {
                w = 0.01;
            }

            if (iround == 0 || sizep < sizep0 || (sizep < .001 && sizep0 < .001))
            {
                goodtimes++;
            }
            else
            {
                goodtimes = 0;
            }

            if ((n == 1 || goodtimes >= Ngoodtimes) && sizep < (rel_error > 1e-5 ? 1 : .001) && H_end(x0, x, f0, *f, rel_error, rel_error, n))
            {
                break;
            }
        }

		gradientB(n, x, *f, g, tv, xmark, sizep);

        /* modify the working set */
        for (i = 0; i < n; i++)  	/* add constraints, reduce H */
        {
            if (xmark[i])
            {
                continue;
            }

            if (fabs(x[i] - xl[i]) < 1e-6 && -g[i] < 0)
            {
                xmark[i] = -1;
            }
            else if (fabs(x[i] - xu[i]) < 1e-6 && -g[i] > 0)
            {
                xmark[i] = 1;
            }

            if (xmark[i] == 0)
            {
                continue;
            }

            xtoy(H, C, nfree * nfree);

            for (it = 0; it < nfree; it++)
                if (ix[it] == i)
                {
                    break;
                }

            for (i1 = it; i1 < nfree - 1; i1++)
            {
                ix[i1] = ix[i1 + 1];
            }

            for (i1 = 0, nfree--; i1 < nfree; i1++)
                for(i2=0; i2 < nfree; ++i2)
					H[i1 *nfree + i2] = C[(i1 + (i1 >= it)) * (nfree + 1) + i2 + (i2 >= it)];
        }

        for (i = 0, it = 0, w = 0; i < n; ++i)  	/* delete a constraint, enlarge H */
        {
            if (xmark[i] == -1 && -g[i] > w)
            {
                it = i;
                w = -g[i];
            }
            else if (xmark[i] == 1 && -g[i] < -w)
            {
                it = i;
                w = g[i];
            }
        }

        if (w > 10 * sizep / nfree)  	/* *** */
        {
            xtoy(H, C, nfree * nfree);

			//for(i1=0; i1 < nfree; ++i1)
			//	for(i2=0; i2 < nfree; ++i2)
			//		H[i1 * (nfree + 1) + i2] = C[i1 * nfree + i2];
			for(i1=0; i1 < nfree; ++i1) memcpy(&H[i1 * (nfree + 1)], &C[i1 * nfree], nfree*sizeof(double));


			for(i1=0; i1 < (nfree+1); ++i1) H[i1 * (nfree + 1) + nfree] = H[nfree * (nfree + 1) + i1] = 0;

            H[(nfree + 1) * (nfree + 1) - 1] = 1;
            xmark[it] = 0;
            ix[nfree++] = it;
        }

        if (mNoisy > 2)
        {
            printf(" | %d/%d", n - nfree, n);
        }

        for (i = 0, f0 = *f; i < nfree; i++)
        {
            y[i] = g[ix[i]] - g0[ix[i]];
            s[i] = x[ix[i]] - x0[ix[i]];
        }

		//for(i=0; i < n; ++i)
  //      {
  //          g0[i] = g[i];
  //          x0[i] = x[i];
  //      }
		memcpy(g0, g, n*sizeof(double)); memcpy(x0, x, n*sizeof(double));
        /* renewal of H varies with different algorithms   */
		/* BFGS */

        for (i = 0, w = v = 0.; i < nfree; i++)
        {
            for (j = 0, z[i] = 0.; j < nfree; j++)
            {
                z[i] += H[i * nfree + j] * y[j];
            }

            w += y[i] * z[i];
            v += y[i] * s[i];
        }

        if(fabs(v) < small)
        {
            identityMatrix(H, nfree);
            fail = 1;
            continue;
        }

		for(i=0; i < nfree; ++i)
			for(j=0; j < nfree; ++j)
				H[i*nfree + j] += ((1 + w / v) * s[i] * s[j] - z[i] * s[j] - s[i] * z[j]) / v;

    } // end of for(iround, mMaxIterations)

    /* try to remove this after updating LineSearch2() */
    //*f = (*fun) (x, n);	++mNumFunCall;
	*f = -mModel->computeLikelihood(x, n, mTraceFun);

    if(mNoisy > 2)
    {
        printf("\n");
    }

    if(iround == mMaxIterations)
    {
        if (fout)
        {
            fprintf(fout, "\ncheck convergence! (max number of iterations reached)\n");
        }

        return -1;
    }

    if(nfree == n)
    {
        xtoy(H, space, n * n);	/* H has variance matrix, or inverse of Hessian */
        return 1;
    }

    return 0;
}

bool Ming2::H_end(const double x0[], const double x1[], double f0, double f1, double e1, double e2, int n) const
{
	// Himmelblau termination rule. Return true for stop, false otherwise.
    double r = norm(x0, n);

    if(r < e2) r = 1;

    r *= e1;

    if(distance(x1, x0, n) >= r) return false;

    r = fabs(f0);

    if(r < e2) r = 1;

    r *= e1;

	return (fabs(f1 - f0) < r);
}

double Ming2::fun_LineSearch(double t, const double x0[], const double p[], double x[], int n)
{
    for(int i=0; i < n; ++i) x[i] = x0[i] + t * p[i];
    //return ((*fun) (x, n));
	return -mModel->computeLikelihood(x, n, mTraceFun);
}



double Ming2::LineSearch2(double *f, const double x0[], const double p[], double step, double limit, double e, double space[], int iround, int n)
{
    /* linear search using quadratic interpolation
       from x0[] in the direction of p[],
                    x = x0 + a*p        a ~(0,limit)
       returns (a).    *f: f(x0) for input and f(x) for output

       x0[n] x[n] p[n] space[n]

       adapted from Wolfe M. A.  1978.  Numerical methods for unconstrained
       optimization: An introduction.  Van Nostrand Reinhold Company, New York.
       pp. 62-73.
       step is used to find the bracket and is increased or reduced as necessary,
       and is not terribly important.
    */
    int ii = 0, maxround = 10, status, i, nsymb = 0;
    double *x = space, factor = 4, small = 1e-10, smallgapa = 0.2;
    double a0, a1, a2, a3, a4 = -1, a5, a6, f0, f1, f2, f3, f4 = -1, f5, f6;

    /* look for bracket (a1, a2, a3) with function values (f1, f2, f3)
       step length step given, and only in the direction a>=0
    */

    if (mNoisy > 2)
        printf("\n%3d h-m-p %7.4f %6.4f %8.4f ", iround + 1, step, limit, norm(p, n));

    if (step <= 0 || limit < small || step >= limit)
    {
        if (mNoisy > 2)
            printf("\nh-m-p:%20.8e%20.8e%20.8e %12.6f\n", step, limit, norm(p, n), *f);

        return 0;
    }

    a0 = a1 = 0;
    f1 = f0 = *f;
    a2 = a0 + step;
    f2 = fun_LineSearch(a2, x0, p, x, n);

    if (f2 > f1)  		/* reduce step length so the algorithm is decreasing */
    {
        for (;;)
        {
            step /= factor;

            if (step < small)
            {
                return 0;
            }

            a3 = a2;
            f3 = f2;
            a2 = a0 + step;
            f2 = fun_LineSearch(a2, x0, p, x, n);

            if (f2 <= f1)
            {
                break;
            }

            if (mNoisy > 2)   //CMV added correct #if
            {
                putchar('-');
                nsymb++;
            }
        }
    }
    else  			/* step length is too small? */
    {
        for (;;)
        {
            step *= factor;

            if (step > limit)
            {
                step = limit;
            }

            a3 = a0 + step;
            f3 = fun_LineSearch(a3, x0, p, x, n);

            if (f3 >= f2)
            {
                break;
            }

            if (mNoisy > 2)
            {
                putchar('+');
                nsymb++;
            }

            a1 = a2;
            f1 = f2;
            a2 = a3;
            f2 = f3;

            if(step >= limit)
            {
                if(mNoisy > 2) //CMV added correct #if
				{
                    for (; nsymb < 5; nsymb++)
                    {
                        printf(" ");
                    }

                    printf(" %12.6f%3c %6.4f", *f = f3, 'm', a3);
                    //printf(" %12.6f%3c %6.4f %5d", *f = f3, 'm', a3, mNumFunCall);
				}

                *f = f3;
                return a3;
            }
        }
    }

    /* iteration by quadratic interpolation, fig 2.2.9-10 (pp 71-71) */
    for(ii = 0; ii < maxround; ii++)
    {
        /* a4 is the minimum from the parabola over (a1,a2,a3)  */
        a4 = (a2 - a3) * f1 + (a3 - a1) * f2 + (a1 - a2) * f3;

        if (fabs(a4) > 1e-100)
            a4 = ((a2 * a2 - a3 * a3) * f1 + (a3 * a3 - a1 * a1) * f2 + (a1 * a1 - a2 * a2) * f3) / (2 * a4);

        if (a4 > a3 || a4 < a1)  	/* out of range */
        {
            a4 = (a1 + a2) / 2;
            status = 'N';
        }
        else
        {
            if ((a4 <= a2 && a2 - a4 > smallgapa * (a2 - a1))
            || (a4 > a2 && a4 - a2 > smallgapa * (a3 - a2)))
            {
                status = 'Y';
            }
            else
            {
                status = 'C';
            }
        }

        f4 = fun_LineSearch(a4, x0, p, x, n);

        if (mNoisy > 2) //CMV added correct #if
        {
            putchar(status);
        }

        if (fabs(f2 - f4) < e * (1 + fabs(f2)))
        {
            if (mNoisy > 2) //CMV added correct #if
                for (nsymb += ii + 1; nsymb < 5; nsymb++)
                {
                    putchar(' ');
                }
            break;
        }

        /* possible multiple local optima during line search */
        if (mNoisy > 2 && ((a4 < a2 && f4 > f1) || (a4 > a2 && f4 > f3)))   //CMV added correct #if
        {
            printf("\n\na %12.6f %12.6f %12.6f %12.6f", a1, a2, a3, a4);
            printf("\nf %12.6f %12.6f %12.6f %12.6f\n", f1, f2, f3, f4);

            for (a5 = a1; a5 <= a3; a5 += (a3 - a1) / 20)
            {
                printf("\t%.6e ", a5);

                if (n < 5)
                {
                    for(i=0; i < n; ++i) printf("\t%.6f", x0[i] + a5 * p[i]);
                }

                printf("\t%.6f\n", fun_LineSearch(a5, x0, p, x, n));
            }

            puts("Linesearch2 a4: multiple optima?");
        }

        if (a4 <= a2)  		/* fig 2.2.10 */
        {
            if (a2 - a4 > smallgapa * (a2 - a1))
            {
                if (f4 <= f2)
                {
                    a3 = a2;
                    a2 = a4;
                    f3 = f2;
                    f2 = f4;
                }
                else
                {
                    a1 = a4;
                    f1 = f4;
                }
            }
            else
            {
                if (f4 > f2)
                {
                    a5 = (a2 + a3) / 2;
                    f5 = fun_LineSearch(a5, x0, p, x, n);

                    if (f5 > f2)
                    {
                        a1 = a4;
                        a3 = a5;
                        f1 = f4;
                        f3 = f5;
                    }
                    else
                    {
                        a1 = a2;
                        a2 = a5;
                        f1 = f2;
                        f2 = f5;
                    }
                }
                else
                {
                    a5 = (a1 + a4) / 2;
                    f5 = fun_LineSearch(a5, x0, p, x, n);

                    if (f5 >= f4)
                    {
                        a3 = a2;
                        a2 = a4;
                        a1 = a5;
                        f3 = f2;
                        f2 = f4;
                        f1 = f5;
                    }
                    else
                    {
                        a6 = (a1 + a5) / 2;
                        f6 = fun_LineSearch(a6, x0, p, x, n);

                        if (f6 > f5)
                        {
                            a1 = a6;
                            a2 = a5;
                            a3 = a4;
                            f1 = f6;
                            f2 = f5;
                            f3 = f4;
                        }
                        else
                        {
                            a2 = a6;
                            a3 = a5;
                            f2 = f6;
                            f3 = f5;
                        }
                    }
                }
            }
        }
        else  		/* fig 2.2.9 */
        {
            if (a4 - a2 > smallgapa * (a3 - a2))
            {
                if (f2 >= f4)
                {
                    a1 = a2;
                    a2 = a4;
                    f1 = f2;
                    f2 = f4;
                }
                else
                {
                    a3 = a4;
                    f3 = f4;
                }
            }
            else
            {
                if (f4 > f2)
                {
                    a5 = (a1 + a2) / 2;
                    f5 = fun_LineSearch(a5, x0, p, x, n);

                    if (f5 > f2)
                    {
                        a1 = a5;
                        a3 = a4;
                        f1 = f5;
                        f3 = f4;
                    }
                    else
                    {
                        a3 = a2;
                        a2 = a5;
                        f3 = f2;
                        f2 = f5;
                    }
                }
                else
                {
                    a5 = (a3 + a4) / 2;
                    f5 = fun_LineSearch(a5, x0, p, x, n);

                    if (f5 >= f4)
                    {
                        a1 = a2;
                        a2 = a4;
                        a3 = a5;
                        f1 = f2;
                        f2 = f4;
                        f3 = f5;
                    }
                    else
                    {
                        a6 = (a3 + a5) / 2;
                        f6 = fun_LineSearch(a6, x0, p, x, n);

                        if (f6 > f5)
                        {
                            a1 = a4;
                            a2 = a5;
                            a3 = a6;
                            f1 = f4;
                            f2 = f5;
                            f3 = f6;
                        }
                        else
                        {
                            a1 = a5;
                            a2 = a6;
                            f1 = f5;
                            f2 = f6;
                        }
                    }
                }
            }
        }
    }

    if (f2 > f0 && f4 > f0)
    {
        a4 = 0;
    }

    if (f2 <= f4)
    {
        *f = f2;
        a4 = a2;
    }
    else
    {
        *f = f4;
    }

    if (mNoisy > 2)
    {
        printf(" %12.6f%3d %6.4f", *f, ii, a4);
        //printf(" %12.6f%3d %6.4f %5d", *f, ii, a4, mNumFunCall);
    }

    return a4;
}


void Ming2::gradientB(int n, const double x[], double f0, double g[], double space[], const int xmark[], double sizep) const
{
    /* f0=fun(x) is always provided.
       xmark=0: central; 1: upper; -1: down
    */
    int i, j;
    double *x0 = space, *x1 = space + n, eh;	/* eh0=1e-6 || 1e-7 */

    for(i=0; i < n; ++i)
    {
        eh = mDeltaForGradient * (fabs(x[i]) + 1);

        if(xmark[i] == 0 && (mAlwaysCenter || sizep < 1))  	/* central */
        {
            for(j=0; j < n; ++j) x0[j] = x1[j] = x[j];
            eh = pow(eh, .67);
            x0[i] -= eh;
            x1[i] += eh;
            //g[i] = ((*fun) (x1, n) - (*fun) (x0, n)) / (eh * 2.0);
			g[i] = (-mModel->computeLikelihood(x1, n, mTraceFun) + mModel->computeLikelihood(x0, n, mTraceFun)) / (eh * 2.0);
        }
        else  		/* forward or backward */
        {
            //for(j=0; j < n; ++j)  x1[j] = x[j];
			memcpy(x1, x, n*sizeof(double));

            if(xmark[i])
            {
                eh *= -xmark[i];
            }

            x1[i] += eh;
            g[i] = (-mModel->computeLikelihood(x1, n, mTraceFun) - f0) / eh;
        }
    }
}

