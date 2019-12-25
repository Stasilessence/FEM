#include "fragmentationGrid.h"
#include "FEM.h"

#define TEST 2


int main()
{
	string dir = "test10";
	fragmentationGrid A;
	A.process(dir);
	FEM f(dir);
	f.solveSystem();
	f.output(dir);
}