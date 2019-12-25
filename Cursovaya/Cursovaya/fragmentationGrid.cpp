#include "fragmentationGrid.h"

void fragmentationGrid::read(string dir)
{
	FILE* f;
	fopen_s(&f, (dir + "/x.txt").c_str(), "r");
	fscanf_s(f, "%d", &howx);
	xw = new double[howx];
	ihx = new int[howx - 1];
	for (int i = 0; i < howx; ++i) {
		fscanf_s(f, "%lf", &xw[i]);
	}
	fclose(f);
	fopen_s(&f, (dir + "/y.txt").c_str(), "r");
	fscanf_s(f, "%d", &howy);
	yw = new double[howy];
	ihy = new int[howy - 1];
	for (int i = 0; i < howy; ++i) {
		fscanf_s(f, "%lf", &yw[i]);
	}
	fclose(f);
	fopen_s(&f, (dir + "/hxhy.txt").c_str(), "r");
	for (int i = 0; i < howx - 1; ++i) {
		fscanf_s(f, "%d", &ihx[i]);
		ihx[i] *= 2;
		hx = (xw[i + 1] - xw[i]) / ihx[i];
		for (int j = 0; j < ihx[i]; ++j)
		{
			x.push_back(xw[i] + j * hx);
		}
	}
	x.push_back(xw[howx - 1]);
	for (int i = 0; i < howy - 1; ++i) {
		fscanf_s(f, "%d", &ihy[i]);
		ihy[i] *= 2;
		hy = (yw[i + 1] - yw[i]) / ihy[i];
		for (int j = 0; j < ihy[i]; j++)
		{
			y.push_back(yw[i] + j * hy);
		}
	}
	y.push_back(yw[howy - 1]);
	xsize = x.size();
	ysize = y.size();
	fclose(f);
}

double fragmentationGrid::NodeValue(double x, double y)
{
	//return x + y;
	//return x * x + y * y;
	//return x * y * y * y;
	return x / y;
}

void fragmentationGrid::process(string dir)
{
	read(dir);
	FILE* f;

	fopen_s(&f, (dir + "/xy.txt").c_str(), "w");
	fprintf(f, "%d\n", ysize * xsize);
	for (int i = 0; i < ysize; i++)
		for (int j = 0; j < xsize; j++)
			fprintf(f, " %0.15lg %0.15lg\n", x[j], y[i]);
	fclose(f);
	fopen_s(&f, (dir + "/nvtr.txt").c_str(), "w");
	fprintf(f, "%d\n", (xsize / 2) * (ysize / 2));

	for (int i = 0; i < (xsize / 2) * (ysize / 2); i++) {
		int k = K(i);
		fprintf(f, "%d %d %d ", k, k + 1, k + 2);
		fprintf(f, "%d %d %d ", k + 2 * xsize / 2, k + 2 * xsize / 2 + 1, k + 2 * xsize / 2 + 2);
		fprintf(f, "%d %d %d\n", k + 2 * (2 * xsize / 2), k + 2 * (2 * xsize / 2) + 1, k + 2 * (2 * xsize / 2) + 2);
	}
	fclose(f);
	fopen_s(&f, (dir + "/nvk.txt").c_str(), "w");
	fprintf(f, "%d\n", xsize * 2 + ysize * 2 - 4);

	for (int i = 0; i < xsize; i++)
	{
		fprintf(f, "%d %le\n", i, NodeValue(x[i], y[0]));
	}
	for (int i = 0; i < xsize; i++)
	{
		fprintf(f, "%d %le\n", i + (2 * xsize / 2) * (ysize - 1), NodeValue(x[i], y[ysize - 1]));
	}
	for (int i = 1; i < ysize - 1; i++)
	{
		fprintf(f, "%d %le\n", i * (2 * xsize / 2), NodeValue(x[0], y[i]));
	}
	for (int i = 1; i < ysize - 1; i++)
	{
		fprintf(f, "%d %le\n", (i + 1) * (2 * xsize / 2) - 1, NodeValue(x[xsize - 1], y[i]));
	}
	fclose(f);
	delete[]xw;
	delete[]yw;
	delete[]ihx;
	delete[]ihy;
	x.~vector();
	y.~vector();
}

int fragmentationGrid::K(int i)
{
	return 2 * ((i) / (xsize / 2)) * (xsize) + 2 * (i % (xsize / 2));
}
