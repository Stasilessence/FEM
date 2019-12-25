#pragma once
#include <string>
#include <vector>
#include "LOSSolver.h"
using namespace std;
struct numberNode
{
	int x[9];
	int num_sqr;
};
struct point
{
	double x;
	double y;
	point();
};
struct firstBoundCondition
{
	int el;
	double value;
};
struct secondBoundCondition
{
	int nodes[3];
	int numFunction;
};
struct thirdBoundCondition
{
	int nodes[3];
	int numFunction;
	double b;
};
class FEM
{
private:
	static const int countBasicFuncs = 9;
	void initMatrix();
	double basicFunc(double x, double xpStart, double xpMid, double xpEnd, int num);
	void initGM();
	int globalIndex(int elem, int num);
	int getNumberElement(double x, double y);
	double lambda(int elem, double x, double y);
	double gamma(int elem, double x, double y);
	double function(double x, double y);
	double funcsKr2(int number, double x, double y);
	double funcsKr3(int number, double x, double y);
	void calcFLocal(int elem);
	void calcGM(int elem);
	double bisquar(int def_elem, int number, double x, double y);
	void addElem(int i, int j, double a);
	void ukr1(int number);
	void ukr2(int number);
	void ukr3(int number);
	numberNode* nvtr;
	point* xy;
	firstBoundCondition* nvk1;
	secondBoundCondition* nvk2;
	thirdBoundCondition* nvk3;
	double GLocal[9][9];
	double MLocal[9][9];
	double Gxy[3][3];
	double M[3][3];
	double Mxy[3][3];
	double fLocal[9];
	double bLocal[9];
	int* ig, * jg;
	double* ggl, * ggu, * di;
	int size_gg;
	int nk1;
	int nk2;
	int nk3;
	double* b;
	int nuz, nel;
public:
	FEM(string direct);
	void solveSystem();
	double solveInPoint(double x, double y);
	void output(string dir);
	~FEM();
};


