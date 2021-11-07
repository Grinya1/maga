#include <iostream>
#include <fstream>
#include <stdio.h>
#include <omp.h>
#include<conio.h>
#include<math.h>
//#include <vcl.h>
#include<string.h>

#define _CRT_SECURE_NO_WARNINGS
#define N 12
#define M 5
#define Mconc 10
#define Mreac 8

using namespace std;

/*
множитель перед omega с чертой
*/
double konst(double tmp, int ii);

/*
тепловой эффект реакции -- последнее слагаемое в уравнении тепла
*/

/*double summa(double yy, double tmp)
{
	for (int m = 1; m < M; m++)
		delta_teta[m] = qp[m] / (T_op * c_p);
	double sum = 0;
	sum += delta_teta[1] * konst(tmp, 1) * yy;
	sum += delta_teta[2] * konst(tmp, 2) * yy;
	sum += delta_teta[3] * konst(tmp, 3);
	sum += delta_teta[4] * konst(tmp, 4) * yy;
	return sum;
	}
double sum_mu(double yy, double temp)
{
	double sum = 0.;
	sum += delt[1] * konst(temp, 1) * yy;
	sum += delt[2] * konst(temp, 2) * yy;
	sum += delt[3] * konst(temp, 3);
	sum += delt[4] * konst(temp, 4) * yy;
	return sum;
	}*/

//double teta1, teta2, koefa, koefb, koefc; //расчет комплексов на поверхности зерна Боденштейн


double Sk(double qc);

double Rc(double qc);

double S(double qc);

double w1(double qc, double teta1, double teta2, double yy1, double z1, double z2, double tempK);

double w2(double qc, double teta1, double teta2, double yy1, double z1, double z2, double tempK);

double w3(double qc, double teta1, double teta2, double yy1, double z1, double z2, double tempK);

double w4(double qc, double teta1, double teta2, double yy1, double z1, double z2, double tempK);

double w5(double qc, double teta1, double teta2, double yy1, double z1, double z2, double tempK);

double w6(double qc, double teta1, double teta2, double yy1, double z1, double z2, double tempK);

double w7(double qc, double teta1, double teta2, double yy1, double z1, double z2, double tempK);

double fqc(double qc, double teta1, double teta2, double yy1, double z1, double z2, double tempK);

double fz1(double qc, double teta1, double teta2, double yy1, double z1, double z2, double tempK);

double fz2(double qc, double teta1, double teta2, double yy1, double z1, double z2, double tempK);

double fy1(double qc, double teta1, double teta2, double yy1, double z1, double z2, double tempK);

double fy2(double qc, double teta1, double teta2, double yy1, double z1, double z2, double tempK);

double fy3(double qc, double teta1, double teta2, double yy1, double z1, double z2, double tempK);

double fy4(double qc, double teta1, double teta2, double yy1, double z1, double z2, double tempK);

double ft1(double qc, double teta1, double teta2, double yy1, double z1, double z2, double tempK);

double ft2(double qc, double teta1, double teta2, double yy1, double z1, double z2, double tempK);

void print_to_file(const char *fname, double array[], int size);

struct Results {
	double teta[N];
	double y[Mconc][N];
};

Results calculation();
