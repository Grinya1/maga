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

double beta_1 = 0.0115;

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

void print_to_file(const char *fname, double array [], int size) {
	// вывод массива в текстовый файл
	FILE* file;

	if ((file = fopen(fname, "w")) == NULL) {
		exit(1);
	}
	for (int i = 0; i < size; i++) {
		fprintf(file, "%5.2d \t %11.6f \n", i, array[i]);
	}
}

Results calculation() {
	double teta[N], teta_new[N];//температура T/T_op
	double y[Mconc][N], y_new[Mconc][N];//концентраци
	double ck = c_p * gam0;
	double teta_0, e, beta_0;

	teta_0 = (T0 + 273.0) / T_op;
	tile = Rz * Rz / (Dz * tau_k);
	e = 1./(tile * ek);
	beta_0 = beta_1 * Rz / Dz;

	//начальные условия
	for (int i = 0; i < N; i++)	{
		//teta[i] = teta_0;
		mu[i] = 0.;
	}
	teta[1] = teta_0;
	teta[0] = teta_0;
	double TT=(T0_inside + 273) / T_op;
	for (int i = 2; i < N; i++)	{
		teta[i] = TT;
	}
	for (int m = 0; m < M; m++) { 
		for (int i = 0; i < N; i++) {
			y[m][i] = 0.;
			y[1][i] = 1.0;
		}
	} //y[1][1] = 1.0;
	for (int i = 0; i < N; i++) {
		y[5][i] = t10; 
		y[6][i] = t20;
		y[7][i] = z10; 
		y[8][i] = z20; 
		y[9][i] = qc0;
	}

	double t = 0.;

	while (t < TMAX) {
		t = t + tau;
		/*cout << w1(y[9][0], y[5][0], y[6][0], y[1][0], y[7][0], y[8][0], teta[0]) << endl;
		cout << w2(y[9][0], y[5][0], y[6][0], y[1][0], y[7][0], y[8][0], teta[0]) << endl;
		cout << w3(y[9][0], y[5][0], y[6][0], y[1][0], y[7][0], y[8][0], teta[0]) << endl;
		cout << w4(y[9][0], y[5][0], y[6][0], y[1][0], y[7][0], y[8][0], teta[0]) << endl;
		cout << w5(y[9][0], y[5][0], y[6][0], y[1][0], y[7][0], y[8][0], teta[0]) << endl;
		cout << w6(y[9][0], y[5][0], y[6][0], y[1][0], y[7][0], y[8][0], teta[0]) << endl;
		cout << w7(y[9][0], y[5][0], y[6][0], y[1][0], y[7][0], y[8][0], teta[0]) << endl;*/
		//printf("%5.2f \n", t);
		//граничные условия
		mu[0] = mu[1];
	   // mu[N - 1] = mu[N-2];
		for (int m = 0; m < M; m++) {
			y[m][0] = y[m][1];
			//y[m][N - 1] = y[m][N-2]+h*beta_0*(y[m][N-1]-x[m]);
			y[m][N - 1] = y[m][N - 2];
			
		}
		//cout << y[1][0] << "  " << y[1][N - 2] << endl;
		teta[1] = teta_0;
		//teta[0] = teta[1]-h*teta_0;
		teta[0] = teta[1];
		//teta[N - 1] = teta[N - 2]+beta_0*h*(teta_0-teta[N-1]);
		teta[N - 1] = teta[N - 2];

		// считаю стефановский поток
		mu[1] = 0.;
		for (int i = 2; i <= N - 1; i++) {
			//mu[i] = (1 / (i * h*i*h)) * ((i-1)*h*(i-1)*h*mu[i-1]+(h*tile*c0/gam0)*(-3*konst(teta[i], 1)*y[1][i] + konst(teta[i], 3)*y[1][i] + 2*konst(teta[i], 4)*y[2][i])* (i - 1) * h * (i - 1) * h);
			mu[i] = (1 / (i * h * i * h)) * ((i - 1) * h * (i - 1) * h * mu[i - 1] + (h * tile * c0 / gam0) *
				(-w1(y[9][i - 1], y[5][i - 1], y[6][i - 1], y[1][i - 1], y[7][i - 1], y[8][i - 1], teta[i - 1])
				+ w3(y[9][i - 1], y[5][i - 1], y[6][i - 1], y[1][i - 1], y[7][i - 1], y[8][i - 1], teta[i - 1])
				+ w5(y[9][i - 1], y[5][i - 1], y[6][i - 1], y[1][i - 1], y[7][i - 1], y[8][i - 1], teta[i - 1]))
				* (i - 1) * h * (i - 1) * h);
		}

		//считаю концентрации
	/*	for (int m = 0; m < M; m++) {
			for (int i = 1; i < N-1; i++) {
				ist_y[1] = -(konst(teta[i], 1) + konst(teta[i], 2) + konst(teta[i], 4)) * y[1][i];
				ist_y[2] = konst(teta[i], 4) * y[1][i];
				ist_y[3] = konst(teta[i], 2);
				ist_y[4] = konst(teta[i], 3) * y[1][i];
				y_new[m][i] = y[m][i]+(tau*e/((i+0.5)*h* (i + 0.5) * h))*((0.5*((i+1)*h*(i+1)*h+ (i) * h * (i) * h)*(y[m][i+1]-y[m][i])/h- 0.5 * ((i) * h * (i) * h + (i-1)* h * (i-1)* h) * (y[m][i] - y[m][i-1]) / h)/h-
					Rz*((i+1)*h*(i+1)*h*mu[i+1]*y[m][i+1] - (i - 1) * h * (i - 1) * h * mu[i - 1] * y[m][i - 1])/(2*h))+
					e*ist_y[m]*tau;
			}
		}*/

		//считаю концентрации
      for (int m = 0; m < M; m++) {
		for (int i = 1; i < N-1; i++) {
			/*f[1] = -(konst(teta[i], 1) + konst(teta[i], 2) + konst(teta[i], 3)) * y[1][i];
			f[2] = konst(teta[i], 1) * y[1][i] - konst(teta[i], 4) * y[2][i];
			f[3] = konst(teta[i], 2) * y[1][i];
			f[4] = konst(teta[i], 3) * y[1][i] + konst(teta[i], 4) * y[2][i];*/

			f[1] = fy1(y[9][i], y[5][i], y[6][i], y[1][i], y[7][i], y[8][i], teta[i]);
			f[2] = fy2(y[9][i], y[5][i], y[6][i], y[1][i], y[7][i], y[8][i], teta[i]);
			f[3] = fy3(y[9][i], y[5][i], y[6][i], y[1][i], y[7][i], y[8][i], teta[i]);
			f[4] = fy4(y[9][i], y[5][i], y[6][i], y[1][i], y[7][i], y[8][i], teta[i]);
			
			y_new[m][i] = y[m][i] + (tau * e) * (1.0 / ((i - 0.5) * h * (i - 0.5) * h)) * (1.0 / h) * (h * i * h * i * (y[m][i + 1] - y[m][i]) / h - h * (i - 1) * h * (i - 1) * (y[m][i] - y[m][i - 1]) / h)
				- (tau * e) * (1.0 / ((i - 0.5) * h * (i - 0.5) * h)) * (1.0 / (h)) * (h * i * h * i * mu[i] * 0.5 * (y[m][i + 1] + y[m][i]) - h * (i - 1) * h * (i - 1) * mu[i - 1] * (y[m][i] + y[m][i - 1]))
				+ f[m] * tau / ek;
		}
	  }

		//считаю температуру
	/*	for (int i = 1; i < N - 1; i++) {
			teta_new[i] = teta[i] + (tau / (gam*(i + 0.5) * h * (i + 0.5) * h)) * (0.5 * ((i + 1) * h * (i + 1) * h + (i)* h * (i)* h) * (teta[i + 1] - teta[i]) / h - 0.5 * ((i)* h * (i)* h + (i - 1) * h * (i - 1) * h) * (teta[i] - teta[i - 1]) / h) / h+
				tau*summa(y[1][i], teta[i])/gam;
		}*/

		//считаю температуру

		for (int i = 1; i < N - 1; i++) {			
			f[0] = c0 / ck / T_op * (qp[1] * w1(y[9][i], y[5][i], y[6][i], y[1][i], y[7][i], y[8][i], teta[i])
			+ qp[2] * w2(y[9][i], y[5][i], y[6][i], y[1][i], y[7][i], y[8][i], teta[i])
			+ qp[3] * w3(y[9][i], y[5][i], y[6][i], y[1][i], y[7][i], y[8][i], teta[i])
			+ qp[4] * w4(y[9][i], y[5][i], y[6][i], y[1][i], y[7][i], y[8][i], teta[i])
			+ qp[5] * w5(y[9][i], y[5][i], y[6][i], y[1][i], y[7][i], y[8][i], teta[i]));

			teta_new[i] = teta[i] + (tau / tile) * (1.0 / ((i - 0.5) * h * (i - 0.5) * h)) * (1.0 / h) * (h * i * h * i * (teta[i + 1] - teta[i]) / h - h * (i - 1) * h * (i - 1) * (teta[i] - teta[i - 1]) / h) + tau * f[0];
			//cout << f[0] << endl;
		}
	
		for (int i = 1; i < N-1; i++) teta[i] = teta_new[i];

		for (int m = 0; m < M; m++) {
			for (int i = 1; i < N-1; i++) {
				y[m][i] = y_new[m][i];
			}
		}

		for (int m = 5; m < 10; m++) {
			for (int i = 1; i < N - 1; i++) {
				f[5] = ft1(y[9][i], y[5][i], y[6][i], y[1][i], y[7][i], y[8][i], teta[i]);
				f[6] = ft2(y[9][i], y[5][i], y[6][i], y[1][i], y[7][i], y[8][i], teta[i]);
				f[7] = fz1(y[9][i], y[5][i], y[6][i], y[1][i], y[7][i], y[8][i], teta[i]);
				f[8] = fz2(y[9][i], y[5][i], y[6][i], y[1][i], y[7][i], y[8][i], teta[i]);
				f[9] = fqc(y[9][i], y[5][i], y[6][i], y[1][i], y[7][i], y[8][i], teta[i]);
				y_new[m][i] = y[m][i] + tau * f[m];
			}
		}

		//пока для вывода

		//граничные условия

		for (int m = 0; m < M; m++) {
			y[m][0] = y[m][1];
			//y[m][N - 1] = y[m][N - 2] + h * beta_0 * ( x[m]+y[m][N - 2]);
			y[m][N - 1] = y[m][N - 2];
		}

		//teta[0] = teta[1]-h*teta_0;
		teta[0] = teta[1];
		//teta[N - 1] = teta[N - 2] + beta_0 * h * (teta_0 - teta[N - 2]);
		teta[N - 1] = teta[N - 2];

		for (int m = 5; m < 10; m++) {
			for (int i = 0; i < N-1; i++) {
				y[m][i] = y_new[m][i];
			}
		}
		
		//cout << f[9] << endl;
	}

	Results results;
	for (int i=0; i<N; i++){
		results.teta[i] = teta[i] * T_op - 273.0;
	}
	for (int i=0; i<Mconc; i++){
		for (int j=0; j<N; j++){
			results.y[i][j] = y[i][j];
		}
	}
	return results;
}
