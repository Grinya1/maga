#include <iostream>
//#include <dir.h>
#include <fstream>
#include <stdio.h>
#include <omp.h>
#include<conio.h>
#include<math.h>
//#include <vcl.h>
#include<string.h>

#define _CRT_SECURE_NO_WARNINGS
#define N 120
#define M 5
#define Mconc 10
#define Mreac 8

using namespace std;

/*
множитель перед omega с чертой
*/
double konst(double tmp, int ii);

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
