#include "headers.h"

double h = 1.e-1;
double k_op[] = {0, 0.2423, 0.8772, 0.6165, 5.2358, 4.3495, 1.5637, 5.3238}; // 1/c при t_op
double E_op[] = {0, 99100.0, 90800.0, 52800.0, 50400.0, 27300.0, 164500.0, 42700.0}; // J/mole при t_op
double qp[] = {0, 83700.0, 394000.0, 67600.0, 303000.0, 311000.0}; // J/mole тепловые эффекты скоростей стадий
double x[] = { 0,0.02,0,0,0 };
double delt[] = { 0,-1.0,0,1.0,0 };//последняя строка матрицы стехиометрических коэффициентов
double Rz = 0.005; // grain radius, m
double T_op = 793.0; // Basic temperature, Calvin 520C
double c_p = 33.47; // J/(mole*Cal)
double tau_k = 4.8; // contact time, sec
double qc0 = 0.05; // the initial and current content of coke
double gam0 = 6.0; // (c_k/c_p) 500-700 t/m3
double ek = 0.3;
double beta_1 = 0.0115;
double tile; // parameter Tile
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
	for (int m = 1; m < Mreac; m++)
		b_j[m] = E_op[m] / (8.314 * T_op);
	double aa = k_op[ii] * exp(b_j[ii] - b_j[ii] / tmp);
	aa = aa * tau_k * phi0;
	return aa;
}

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

Results calculation(int numproc, int rank) {
	double teta[N], teta_new[N];//температура T/T_op
	double y[Mconc][N], y_new[Mconc][N];//концентраци
	double ck = c_p * gam0;
	double teta_0, e, beta_0;

	teta_0 = (T0 + 273.0) / T_op;
	tile = Rz * Rz / (Dz * tau_k);
	e = 1./(tile * ek);
	beta_0 = beta_1 * Rz / Dz;
	double tau_e=tau*e;
	double hh=h*h;
	double h_tile_c0_gam0=h * tile * c0 / gam0;
	double h_1=1.0 / h;
		//double tau_ek=tau/ek;
	//начальные условия
	
	for (int i = 0; i < N; i++)	{	
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
	} 
	for (int i = 0; i < N; i++) {
		y[5][i] = t10; 
		y[6][i] = t20;
		y[7][i] = z10; 
		y[8][i] = z20; 
		y[9][i] = qc0;
	}
	for (int o=0; o<100000; o++){
		
		//граничные условия
		mu[0] = mu[1];
	  
		for (int m = 0; m < M; m++) {
			y[m][0] = y[m][1];
			y[m][N - 1] = y[m][N - 2];
		}
		
		teta[1] = teta_0;
		teta[0] = teta[1];
		teta[N - 1] = teta[N - 2];
		
        // считаю стефановский поток
	
        for (int i = 2; i <= N - 1; i++) {
            mu[i] = (1 / (i * hh * i)) * ((i-1) * hh * (i-1)* mu[i - 1] + h_tile_c0_gam0 *
            (-w1(y[9][i - 1], y[5][i - 1], y[6][i - 1], y[1][i - 1], y[7][i - 1], y[8][i - 1], teta[i - 1])
            + w3(y[9][i - 1], y[5][i - 1], y[6][i - 1], y[1][i - 1], y[7][i - 1], y[8][i - 1], teta[i - 1])
            + w5(y[9][i - 1], y[5][i - 1], y[6][i - 1], y[1][i - 1], y[7][i - 1], y[8][i - 1], teta[i - 1]))
            * (i-1) * hh * (i-1));
		}
	
		//считаю концентрации
		
        for (int m = rank; m < rank+1; m++) { 
            for (int i = 1; i < N-1; i++) {
                f[1] = fy1(y[9][i], y[5][i], y[6][i], y[1][i], y[7][i], y[8][i], teta[i]);
                f[2] = fy2(y[9][i], y[5][i], y[6][i], y[1][i], y[7][i], y[8][i], teta[i]);
                f[3] = fy3(y[9][i], y[5][i], y[6][i], y[1][i], y[7][i], y[8][i], teta[i]);
                f[4] = fy4(y[9][i], y[5][i], y[6][i], y[1][i], y[7][i], y[8][i], teta[i]);
                y_new[m][i] = y[m][i] + tau_e * (1.0 / ((i-0.5) * hh * (i-0.5))) * h_1* (i * hh * i * (y[m][i + 1] - y[m][i]) / h -(i-1) * hh * (i-1)* (y[m][i] - y[m][i - 1]) / h)
                    - tau_e * (1.0 / ((i-0.5) * hh * (i-0.5))) * h_1 * (i * hh * i* mu[i] * 0.5 * (y[m][i + 1] + y[m][i]) - (i-1) * hh * (i-1) * mu[i - 1] * (y[m][i] + y[m][i - 1]))
                    + f[m] *tau/ek;
            }
        }

        //считаю температуру
        double c0_ck_T_op=c0 / ck / T_op;
        double tau_tile=tau / tile;
        for (int i = 1; i < N - 1; i++) {
            
            f[0] = c0_ck_T_op * (qp[1] * w1(y[9][i], y[5][i], y[6][i], y[1][i], y[7][i], y[8][i], teta[i])
            + qp[2] * w2(y[9][i], y[5][i], y[6][i], y[1][i], y[7][i], y[8][i], teta[i])
            + qp[3] * w3(y[9][i], y[5][i], y[6][i], y[1][i], y[7][i], y[8][i], teta[i])
            + qp[4] * w4(y[9][i], y[5][i], y[6][i], y[1][i], y[7][i], y[8][i], teta[i])
            + qp[5] * w5(y[9][i], y[5][i], y[6][i], y[1][i], y[7][i], y[8][i], teta[i]));

            teta_new[i] = teta[i] + tau_tile * (1.0 / ((i - 0.5) * h * (i - 0.5) * h)) * (1.0 / h) * (h * i * h * i * (teta[i + 1] - teta[i]) / h - h * (i - 1) * h * (i - 1) * (teta[i] - teta[i - 1]) / h) + tau * f[0];
        }

        for (int i = 1; i < N-1; i++){ 
            teta[i] = teta_new[i];
        }
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
        //граничные условия

        for (int m = 0; m < M; m++) {
            y[m][0] = y[m][1];
            y[m][N - 1] = y[m][N - 2];
        }

        teta[0] = teta[1];
        teta[N - 1] = teta[N - 2];

        for (int m = 5; m < 10; m++) {
            for (int i = 0; i < N-1; i++) {
                y[m][i] = y_new[m][i];
            }
        }	
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