#pragma once
#include <cstdio>
#include <cmath>
#define myfopen(file, path, mode) fopen_s(file, (path).c_str(), mode)
typedef double* vectorD;

class LOS {



private:
	double eps;
	int *ig, *jg;
	vectorD di, gu, gl;
	vectorD b;
	vectorD x;
	vectorD mv;
	vectorD z;
	vectorD r;
	vectorD p;
	vectorD diag;
	int n, maxiter;
public:
	void multyMatrixVector(vectorD x, vectorD res);
	double scal(vectorD x, vectorD y);
	double norm(vectorD a);
	void solve();
	LOS(int locn, int* lig, int* ljg, vectorD ld, vectorD ggu, vectorD ggl, vectorD lb);
	~LOS();
};