#include "headers.h"

double h = 1.e-1;
double k_op[] = {0, 0.2423, 0.8772, 0.6165, 5.2358, 4.3495, 1.5637, 5.3238}; // 1/c при t_op
double E_op[] = {0, 99100.0, 90800.0, 52800.0, 50400.0, 27300.0, 164500.0, 42700.0}; // J/mole при t_op
double qp[] = {0, 83700.0, 394000.0, 67600.0, 303000.0, 311000.0}; // J/mole тепловые эффекты скоростей стадий
double x[] = { 0,0.02,0,0,0 };
double delt[] = { 0,-1.0,0,1.0,0 };//последняя строка матрицы стехиометрических коэффициентов

double Rz = 0.005; // grain radius, m

//double Dz = 0.000005; // effective diffusion factor, m2/sec

double T_op = 793.0; // Basic temperature, Calvin 520C

//double T0 = 400.0; // начальная температура

double c_p = 33.47; // J/(mole*Cal)

double tau_k = 4.8; // contact time, sec

//double tau = 4.8/1000000.; // contact time, sec

//double TMAX = 10000*tau; //время, до которого ведется расчет

double qc0 = 0.05; // the initial and current content of coke

// on the catalyst, g/g

double gam0 = 6.0; // (c_k/c_p) 500-700 t/m3

double ek = 0.3;

double beta = 0.0115;

double tile; // parameter Tile

//Здесь Язовцева начинается

double tau = 1.e-10;
double c0 = 17.0;
double TMAX = 100000 * tau;
double Rc0 = 0.00005;
double Dz = 0.000005;
double T0 = 400.0;//C
double T0_inside = 300.0;//C

double z1, z2; // фазовые переменные

double Mc = 0.000012;
double roC = 1.6;

double Sk0 = 3 * qc0 / (Rc0 * roC);
double phi0 = gam0 * Sk0 / c0;

const double t10 = 0.12;
const double t20 = 0.0;
const double z10 = 0.015;
const double z20 = 0.0;
const double yy10 = 1.0;	//o2 
const double yy20 = 0.0;	//co2  
const double yy30 = 0.0;	//co  
const double yy40 = 0.0;	//h2o 

double b_j[Mreac];
double delta_teta[M];
double teta[N], teta_new[N];//температура T/T_op
double y[Mconc][N], y_new[Mconc][N];//концентрации
double mu[N]; // стефановский поток
double ist_y[M]; // источниковый член в уравнении концентраций
double f[Mconc]; // источниковый член в уравнении концентраций

/*
множитель перед omega с чертой
*/
double konst(double tmp, int ii) {
	//double sk = pow((qc / qc0), 2.0 / 3.0);
	for (int m = 1; m < Mreac; m++)
		b_j[m] = E_op[m] / (8.314 * T_op);
	double aa = k_op[ii] * exp(b_j[ii] - b_j[ii] / tmp);
	//aa = aa * tau_k * sk * tile;
	aa = aa * tau_k * phi0;
	return aa;
}

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


double Sk(double qc) {
	return Sk0 * pow(qc / qc0, 2.0 / 3.0);
}

double Rc(double qc) {
	return Rc0 * pow(qc / qc0, 1.0 / 3.0);
}

double S(double qc) {
	return pow(qc / qc0, 2.0 / 3.0);
}

double w1(double qc, double teta1, double teta2, double yy1, double z1, double z2, double tempK) {
	return konst(tempK, 1) * (1 - teta2 - teta1) * (1 - teta2 - teta1) * yy1;
}

double w2(double qc, double teta1, double teta2, double yy1, double z1, double z2, double tempK) {
	return konst(tempK, 2) * teta2 * yy1;
}

double w3(double qc, double teta1, double teta2, double yy1, double z1, double z2, double tempK) {
	return konst(tempK, 3) * teta2;
}

double w4(double qc, double teta1, double teta2, double yy1, double z1, double z2, double tempK) {
	return konst(tempK, 4) * teta1 * yy1;
}

double w5(double qc, double teta1, double teta2, double yy1, double z1, double z2, double tempK) {
	return konst(tempK, 5) * teta2 * teta2;
}

double w6(double qc, double teta1, double teta2, double yy1, double z1, double z2, double tempK) {
	return konst(tempK, 6) * roC / Rc(qc) * (teta1 / 6 - z1);
}

double w7(double qc, double teta1, double teta2, double yy1, double z1, double z2, double tempK) {
	return konst(tempK, 7) * roC / Rc(qc) * (4 * teta2 / 3 - z2);
}

double fqc(double qc, double teta1, double teta2, double yy1, double z1, double z2, double tempK) {
	return -Mc * Sk(qc) / phi0 *
		(w2(qc, teta1, teta2, yy1, z1, z2, tempK) + w3(qc, teta1, teta2, yy1, z1, z2, tempK) + w5(qc, teta1, teta2, yy1, z1, z2, tempK));
}

double fz1(double qc, double teta1, double teta2, double yy1, double z1, double z2, double tempK) {
	return Sk(qc) / qc / phi0 *
		(w6(qc, teta1, teta2, yy1, z1, z2, tempK) + z1 * Mc *
		(w2(qc, teta1, teta2, yy1, z1, z2, tempK) + w3(qc, teta1, teta2, yy1, z1, z2, tempK) + w5(qc, teta1, teta2, yy1, z1, z2, tempK)));
}

double fz2(double qc, double teta1, double teta2, double yy1, double z1, double z2, double tempK) {
	return Sk(qc) / qc / phi0 *
		(w7(qc, teta1, teta2, yy1, z1, z2, tempK) + z2 * Mc *
		(w2(qc, teta1, teta2, yy1, z1, z2, tempK) + w3(qc, teta1, teta2, yy1, z1, z2, tempK) + w5(qc, teta1, teta2, yy1, z1, z2, tempK)));
}

double fy1(double qc, double teta1, double teta2, double yy1, double z1, double z2, double tempK) {
	return -S(qc) * (w1(qc, teta1, teta2, yy1, z1, z2, tempK) + w2(qc, teta1, teta2, yy1, z1, z2, tempK) + w4(qc, teta1, teta2, yy1, z1, z2, tempK));
}

double fy2(double qc, double teta1, double teta2, double yy1, double z1, double z2, double tempK) {
	return S(qc) * (w3(qc, teta1, teta2, yy1, z1, z2, tempK));
}

double fy3(double qc, double teta1, double teta2, double yy1, double z1, double z2, double tempK) {
	return S(qc) * (w2(qc, teta1, teta2, yy1, z1, z2, tempK) + w5(qc, teta1, teta2, yy1, z1, z2, tempK));
}

double fy4(double qc, double teta1, double teta2, double yy1, double z1, double z2, double tempK) {
	return S(qc) * (w4(qc, teta1, teta2, yy1, z1, z2, tempK));
}

double ft1(double qc, double teta1, double teta2, double yy1, double z1, double z2, double tempK) {
	return -(S(qc) * w4(qc, teta1, teta2, yy1, z1, z2, tempK) + Sk(qc) * w6(qc, teta1, teta2, yy1, z1, z2, tempK) / phi0);
}

double ft2(double qc, double teta1, double teta2, double yy1, double z1, double z2, double tempK) {
	return S(qc) * (2 * w1(qc, teta1, teta2, yy1, z1, z2, tempK) - w3(qc, teta1, teta2, yy1, z1, z2, tempK)
		+ w4(qc, teta1, teta2, yy1, z1, z2, tempK) - 2 * w5(qc, teta1, teta2, yy1, z1, z2, tempK))
		- Sk(qc) * w7(qc, teta1, teta2, yy1, z1, z2, tempK) / phi0;
}
