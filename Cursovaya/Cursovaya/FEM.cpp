#include "FEM.h"
point::point()
{
	x = 0;
	y = 0;
}
FEM::FEM(string direct)
{
	string dir = direct + "/";
	FILE* f;
	myfopen(&f, dir + "xy.txt", "r");
	fscanf_s(f, "%d", &nuz);
	xy = new point[nuz];
	for (int i = 0; i < nuz; ++i)
	{
		fscanf_s(f, "%lf%lf", &xy[i].x, &xy[i].y);
	}
	fclose(f); 
	myfopen(&f, dir + "nvtr.txt", "r");
	fscanf_s(f, "%d", &nel);
	nvtr = new numberNode[nel];
	for (int i = 0; i < nel; ++i)
	{
		for (int j = 0; j < 9; ++j) {
			fscanf_s(f, "%d", &nvtr[i].x[j]);
		}
	}
	fclose(f);
	myfopen(&f, dir + "nvk.txt", "r");
	fscanf_s(f, "%d", &nk1);
	nvk1 = new firstBoundCondition[nk1];
	for (int i = 0; i < nk1; ++i)
		fscanf_s(f, "%d%lf", &nvk1[i].el, &nvk1[i].value);
	fclose(f);
	myfopen(&f, dir + "nvk1.txt", "r");
	fscanf_s(f, "%d", &nk2);
	nvk2 = new secondBoundCondition[nk2];
	for (int i = 0; i < nk2; ++i)
	{
		fscanf_s(f, "%d%d%d%d", &nvk2[i].nodes[0], &nvk2[i].nodes[1], &nvk2[i].nodes[2], &nvk2[i].numFunction);
	}
	fclose(f);
	myfopen(&f, dir + "nvk2.txt", "r");
	fscanf_s(f, "%d", &nk3);
	nvk3 = new thirdBoundCondition[nk3];
	for (int i = 0; i < nk3; ++i)
		fscanf_s(f, "%d%d%d%d%lf", &nvk3[i].nodes[0], &nvk3[i].nodes[1], &nvk3[i].nodes[2], &nvk3[i].numFunction, &nvk3[i].b);
	fclose(f);
	b = new double[nuz];
	for (int i = 0; i < nuz; ++i)
		b[i] = 0;
	initMatrix();
	initGM();
}

void FEM::solveSystem()
{
	double g = gamma(0, 0, 0);
	for (int i = 0; i < nel; ++i)
	{
		calcGM(i);
		for (int k = 0; k < 9; ++k)
		{
			for (int j = 0; j < 9; ++j)
			{
				addElem(nvtr[i].x[k], nvtr[i].x[j], GLocal[k][j] + g * MLocal[k][j]);
			}
		}
		calcFLocal(i);
		for (int k = 0; k < 9; ++k)
		{
			b[nvtr[i].x[k]] += bLocal[k];
		}
	}
	for (int i = 0; i < nk2; ++i)
	{
		ukr2(i);
	}
	for (int i = 0; i < nk3; ++i)
	{
		ukr3(i);
	}
	for (int i = 0; i < nk1; ++i)
	{
		ukr1(i);
	}
	LOS solver(nuz, ig, jg, di, ggu, ggl, b);
	solver.solve();
}

void FEM::output(string dir)
{
	FILE* f;
	myfopen(&f, dir + "/out.txt", "w");
	for (int i = 0; i < nuz; ++i) {
		fprintf(f, "%lf %lf %.12le\n", xy[i].x, xy[i].y, b[i]);
	}
	fclose(f);
}

FEM::~FEM()
{
	delete[]xy;
	delete[]nvtr;
	delete[]nvk1;
	delete[]nvk2;
	delete[]nvk3;
	delete[]b;
	delete[]ig;
	delete[]jg;
	delete[]ggl;
	delete[]ggu;
	delete[]di;
}

void FEM::initMatrix()
{
	di = new double[nuz];
	for (int i = 0; i < nuz; ++i)
	{
		di[i] = 0;
	}
	ig = new int[nuz + 1];
	vector<vector<int>> list(nuz);
	int listsize = 0;
	for (int elem = 0; elem < nel; elem++)
	{
		for (int i = 0; i < countBasicFuncs; ++i) {
			int k = nvtr[elem].x[i];
			for (int j = i + 1; j < countBasicFuncs; j++) {
				int ind1 = k;
				int ind2 = nvtr[elem].x[j];
				if (ind2 < ind1) {
					ind1 = ind2;
					ind2 = k;
				}
				int iaddr = list[ind2].size();
				if (iaddr == 0) {
					listsize++;
					list[ind2].push_back(ind1);
				}
				else
				{
					int l = 0;
					for (l = 0; l < iaddr && list[ind2][l] < ind1; l++);
					if (l == iaddr)
					{
						listsize++;
						list[ind2].push_back(ind1);
					}
					else
					{
						if (list[ind2][l] > ind1) {
							listsize++;
							list[ind2].emplace(list[ind2].begin() + l, ind1);
						}
					}
				}
			}
		}
	}
	jg = new int[listsize];
	ig[0] = 0;
	size_gg = listsize;
	int size_sum = 0;
	for (int i = 0; i < nuz; ++i)
	{
		int size = list[i].size();
		ig[i + 1] = ig[i] + list[i].size();
		for (int j = 0; j < size; j++)
			jg[size_sum + j] = list[i][j];
		size_sum += size;
	}
	ggl = new double[listsize];
	ggu = new double[listsize];
	for (int i = 0; i < listsize; ++i)
	{
		ggl[i] = 0; 
		ggu[i] = 0; 
	}
	list.~vector();
}

double FEM::funcsKr2(int numFunc, double x, double y)
{
	switch (numFunc) 
	{
	case 1: {return x / (y * y); }
	case 2: {return -x / (y * y); }
	}
}

double FEM::funcsKr3(int numFunc, double x, double y)
{
	switch (numFunc) 
	{
	case 1: {return x / y - x / (y * y); }
	case 2: {return x / y + x / (y * y); }
	}
}

void FEM::ukr1(int number)
{
	int elem = nvk1[number].el;
	di[elem] = 1e+16;
	b[elem] = nvk1[number].value * 1e+16;
}

void FEM::ukr2(int number)
{
	int* uzel = nvk2[number].nodes;
	double h = xy[uzel[2]].y - xy[uzel[0]].y;
	if (h == 0)
	{
		h = xy[uzel[2]].x - xy[uzel[0]].x;
	}
	double buff[3];
	for (int i = 0; i < 3; ++i) {
		buff[i] = 0;
		for (int j = 0; j < 3; j++)
			buff[i] += funcsKr2(nvk2[number].numFunction, xy[uzel[j]].x, xy[uzel[j]].y) * Mxy[i][j];
		b[uzel[i]] += buff[i] * h / 30.;
	}
}

void FEM::ukr3(int number)
{
	int* uzel = nvk3[number].nodes;
	double h = xy[uzel[2]].y - xy[uzel[0]].y;
	double beta = nvk3[number].b;
	if (h == 0) h = xy[uzel[2]].x - xy[uzel[0]].x;
	double buff[3];
	for (int i = 0; i < 3; ++i) {
		buff[i] = 0;
		for (int j = 0; j < 3; j++)
		{
			buff[i] += funcsKr3(nvk3[number].numFunction, xy[uzel[j]].x, xy[uzel[j]].y) * Mxy[i][j];
		}
		b[uzel[i]] += buff[i] * h * beta / 30.;
	}
	for (int k = 0; k < 3; k++)
	{
		for (int j = 0; j < 3; j++)
		{
			addElem(nvk3[number].nodes[k], nvk3[number].nodes[j], Mxy[k][j] * beta * h / 30.0);
		}
	}
}

void FEM::initGM()
{
	Gxy[0][0] =  11.; Gxy[0][1] = -12.; Gxy[0][2] =  1.;
	Gxy[1][0] = -12.; Gxy[1][1] =  16.; Gxy[1][2] = -4.;
	Gxy[2][0] =   1.; Gxy[2][1] =  -4.; Gxy[2][2] =  3.;

	Mxy[0][0] =  7.; Mxy[0][1] =  4.; Mxy[0][2] = -1.;
	Mxy[1][0] =  4.; Mxy[1][1] = 16.; Mxy[1][2] =   0;
	Mxy[2][0] = -1.; Mxy[2][1] =  0.; Mxy[2][2] =  1.;

	M[0][0] =  4.; M[0][1] =  2.; M[0][2] = -1.;
	M[1][0] =  2.; M[1][1] = 16.; M[1][2] =  2.;
	M[2][0] = -1.; M[2][1] =  2.; M[2][2] =  4.;	
}

double FEM::basicFunc(double x, double xpStart, double xpMid, double xpEnd, int num)
{
	double hx1 = xpMid - xpStart;
	double hx2 = xpEnd - xpStart;
	double hx3 = xpEnd - xpMid;
	switch (num)
	{
	case 0: {return((xpMid - x) * (xpEnd - x)) / (hx1 * hx2); }
	case 1: {return((x - xpStart) * (xpEnd - x)) / (hx1 * hx3); }
	case 2: {return((x - xpStart) * (x - xpMid)) / (hx2 * hx3); }
	}
}

int FEM::globalIndex(int elem, int num)
{
	return nvtr[elem].x[num];
}

double FEM::solveInPoint(double x, double y)
{
	double u = 0;
	int numElem = getNumberElement(x, y);
	for (int i = 0; i < countBasicFuncs; ++i) {
		u += b[globalIndex(numElem, i)] * bisquar(numElem, i, x, y);
	}
	return u;
}

int FEM::getNumberElement(double x, double y)
{
	for (int i = 0; i < nel; ++i)
	{
		int x0 = globalIndex(i, 0), x8 = globalIndex(i, 8);
		if (xy[x0].x <= x && xy[x0].y <= y && xy[x8].x >= x && xy[x8].y >= y)
			return i;
	}
	return -1;
}

double FEM::lambda(int elem, double x, double y)
{
	return 1;
	//return x + y;
	//return x * x;
}

double FEM::gamma(int elem, double x, double y)
{
	return 1;
}

double FEM::function(double x, double y)
{
	//TEST1	
	//return x + y;
	//return x + y - 2;
	//TEST2
	//return x * x + y * y - 4;
	//return x * x + y * y - 4 * x - 2 * y;
	//TEST3
	//return -y * y * y - 6 * x * x * y - 9 * y * y * x + x * y * y * y;
	//TEST4
	//return -x;
	//TEST5
	return x / y - 2 * x / (y * y * y);
}

void FEM::calcFLocal(int elem)
{	
	int node;
	for (int i = 0; i < countBasicFuncs; ++i)
	{
		node = globalIndex(elem , i);
		fLocal[i] = function(xy[node].x, xy[node].y);
	}
	for (int i = 0; i < countBasicFuncs; ++i)
	{
		bLocal[i] = 0;
		for (int j = 0; j < countBasicFuncs; j++)
		{
			bLocal[i] += MLocal[i][j] * fLocal[j];
		}
	}
}

void FEM::calcGM(int elem)
{
	int xy0 = globalIndex(elem, 0);
	int xy2 = globalIndex(elem, 2);
	int xy6 = globalIndex(elem, 6);
	int xy8 = globalIndex(elem, 8);
	double hx = xy[xy8].x - xy[xy0].x;
	double hy = xy[xy8].y - xy[xy0].y;
	double lk0k0 = lambda(elem, xy[xy0].x, xy[xy0].y);
	double lk1k0 = lambda(elem, xy[xy2].x, xy[xy2].y);
	double lk0k1 = lambda(elem, xy[xy6].x, xy[xy6].y);
	double lk1k1 = lambda(elem, xy[xy8].x, xy[xy8].y);
	int nui, mui, nuj, muj;
	for (int i = 0; i < countBasicFuncs; ++i) {
		nui = (i) % 3;
		mui = (i) / 3;
		for (int j = 0; j < countBasicFuncs; ++j)
		{
			nuj = (j) % 3;
			muj = (j) / 3;
			GLocal[i][j] = ((hy / hx) * (lk0k0 * Gxy[nui][nuj] * M[mui][muj] + lk1k0 * Gxy[2 - nui][2 - nuj] * M[mui][muj]
				+ lk0k1 * Gxy[nui][nuj] * M[2 - mui][2 - muj] + lk1k1 * Gxy[2 - nui][2 - nuj] * M[2 - mui][2 - muj])
				+ (hx / hy) * (lk0k0 * Gxy[mui][muj] * M[nui][nuj] + lk1k0 * Gxy[2 - mui][2 - muj] * M[nui][nuj]
				+ lk0k1 * Gxy[mui][muj] * M[2 - nui][2 - nuj] + lk1k1 * Gxy[2 - mui][2 - muj] * M[2 - nui][2 - nuj])) / 360.;
			MLocal[i][j] = hx * hy * (M[mui][muj] * M[nui][nuj]) / 900.;
		}
	}
}

double FEM::bisquar(int elem, int num, double x, double y)
{
	int xySt = globalIndex(elem, 0);
	int xyMid = globalIndex(elem, 4);
	int xyEnd = globalIndex(elem, 8);
	double xpEnd = xy[xyEnd].x, xpSt = xy[xySt].x, xpMid = xy[xyMid].x;
	double ypEnd = xy[xyEnd].y, ypSt = xy[xySt].y, ypMid = xy[xyMid].y;
	switch (num)
	{
	case 0: {return basicFunc(x, xpSt, xpMid, xpEnd, 0) * basicFunc(y, ypSt, ypMid, ypEnd, 0); }
	case 1: {return basicFunc(x, xpSt, xpMid, xpEnd, 1) * basicFunc(y, ypSt, ypMid, ypEnd, 0); }
	case 2: {return basicFunc(x, xpSt, xpMid, xpEnd, 2) * basicFunc(y, ypSt, ypMid, ypEnd, 0); }
	case 3: {return basicFunc(x, xpSt, xpMid, xpEnd, 0) * basicFunc(y, ypSt, ypMid, ypEnd, 1); }
	case 4: {return basicFunc(x, xpSt, xpMid, xpEnd, 1) * basicFunc(y, ypSt, ypMid, ypEnd, 1); }
	case 5: {return basicFunc(x, xpSt, xpMid, xpEnd, 2) * basicFunc(y, ypSt, ypMid, ypEnd, 1); }
	case 6: {return basicFunc(x, xpSt, xpMid, xpEnd, 0) * basicFunc(y, ypSt, ypMid, ypEnd, 2); }
	case 7: {return basicFunc(x, xpSt, xpMid, xpEnd, 1) * basicFunc(y, ypSt, ypMid, ypEnd, 2); }
	case 8: {return basicFunc(x, xpSt, xpMid, xpEnd, 2) * basicFunc(y, ypSt, ypMid, ypEnd, 2); }
	};
	return 0;
}

void FEM::addElem(int i, int j, double a)
{
	if (i == j)
	{
		di[i] += a;
	}
	else
	{
		if (i < j)
		{
			int ind;
			for (ind = ig[j]; ind < ig[j + 1]; ++ind)
			{
				if (jg[ind] == i) break;
			}
			ggu[ind] = ggu[ind] + a;
		}
		else
		{
			int ind;
			for (ind = ig[i]; ind < ig[i + 1]; ++ind)
			{
				if (jg[ind] == j)
					break;
			}
			ggl[ind] = ggl[ind] + a;
		}
	}
}
