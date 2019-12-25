#include "LOSSolver.h"



void LOS::multyMatrixVector(vectorD x, vectorD res)
{
	for (int i = 0; i < n; ++i) {
		int gi = ig[i], gi_1 = ig[i + 1];
		res[i] = di[i] * x[i];
		for (int j = gi; j < gi_1; ++j) {
			int column = jg[j];
			res[i] += gl[j] * x[column];
			res[column] += gu[j] * x[i];
		}
	}
}

double LOS::scal(vectorD a, vectorD b)
{
	double s = 0.0;
	for (int i = 0; i < n; i++)
		s += a[i] * b[i];
	return s;
}

double LOS::norm(vectorD a)
{
	return sqrt(scal(a, a));
}


LOS::LOS(int locn, int* lig, int* ljg, vectorD ld, vectorD ggu, vectorD ggl, vectorD lb)
{
	FILE* f;
	fopen_s(&f, "options.txt", "r");
	fscanf_s(f, "%d%lf", &maxiter, &eps);
	n = locn;
	ig = lig;
	jg = ljg;
	di = ld;
	gl = ggl; gu = ggu;
	b = lb;
	mv = new double[n];
	z = new double[n];
	r = new double[n];
	p = new double[n];
	x = new double[n];
	diag = new double[n];
}

LOS::~LOS()
{
	delete[] mv;
	delete[] z;
	delete[] r;
	delete[] p;
	delete[] x;
	delete[] diag;
}

void LOS::solve()
{
	int count = 0;
	for (int i = 0; i < n; ++i)
	{
		x[i] = 1;
	}
	multyMatrixVector(x, mv);
	for (int i = 0; i < n; ++i)
	{
		r[i] = b[i] - mv[i];
		z[i] = r[i];
	}
	multyMatrixVector(z, p);
	double sr = scal(r, r);
	while (sr > eps&& count <= maxiter)
	{
		double pp = scal(p, p);
		double ak = scal(p, r) / pp;
		for (int i = 0; i < n; ++i)
		{
			x[i] = x[i] + ak * z[i];
			r[i] = r[i] - ak * p[i];
		}
		multyMatrixVector(r, mv);
		double bk = -scal(p, mv) / pp;
		for (int i = 0; i < n; ++i)
		{
			z[i] = r[i] + bk * z[i];
			p[i] = mv[i] + bk * p[i];
		}
		sr = sqrt(scal(r, r));
		++count;
		//printf_s("%d: %.2lf", count, sr);
	}
	for (int i = 0; i < n; ++i)
	{
		b[i] = x[i];
	}
}
