#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include "wolf_research.h"
#include "TNT/tnt_array1d.h"
#include "Steeptest_Descent.h"
#include "Newton.h"
#include "Gradient.h"
#include "BFGS.h"
#include "AugmentedLagrange.h"
using namespace std;



double Vf_test5(Array1D<double>& X, Array1D<double>& niu){

	double C = 4.0;
	double a = X[0];
	double b = X[1];
	double niu1 = niu[0];
	double niu2 = niu[1];
	double niu3 = niu[2];
	double niu4 = niu[3];

	double func = (a - 2)*(a - 2)*(a - 2)*(a - 2) + (a - 2 * b)*(a - 2 * b);//目标函数

	double tempF = niu1 + C*(a*a + b*b - 4);//第一个不等式约束，并加以处理
	double niuFunc1 = 0.5*C*((tempF > 0) ? (tempF*tempF) : 0.0) - niu1*niu1;//增广项

	double niuFunc2 = 0.5*C*(a > 0 ? a*a : 0.0) - niu2*niu2;//后两个不等式约束的增广项
	double niuFunc3 = 0.5*C*(b > 0 ? b*b : 0.0) - niu3*niu3;

	double niuFunc4 = niu4*(a*a - b) + 0.5*C*(a*a - b)*(a*a - b);//等式约束的增广项

	//	return X*X + x1*x1;
	return func + niuFunc1 + niuFunc2 + niuFunc3 + niuFunc4;

}

Array1D<double> Vg_test5(Array1D<double>& X, Array1D<double>& niu){
	double C = 4.0;
	double a = X[0];
	double b = X[1];
	double niu1 = niu[0];
	double niu2 = niu[1];
	double niu3 = niu[2];
	double niu4 = niu[4];

//	double func = (a - 2)*(a - 2)*(a - 2)*(a - 2) + (a - 2 * b)*(a - 2 * b);
	double funcA = 4.0*(a - 2)*(a - 2)*(a - 2) + 2.0*(a - 2 * b);//目标函数对a的求导
	double funcB = 4.0*(2.0*b - a);//目标函数对b的求导

	double tempF = niu1 + C*(a*a + b*b - 4);
	double tempFA = C*2.0*a;
	double tempFB = C*2.0*b;
	double niuFunc1A = 0.5*C*((tempF > 0) ? (2.0*tempF*tempFA) : 0.0);//niuFunc1对a，b的求导
	double niuFunc1B = 0.5*C*((tempF > 0) ? (2.0*tempF*tempFB) : 0.0);

	double niuFunc2A = 0.5*C*(a > 0 ? 2.0*a : 0.0);//niuFunc2对a，b的求导
	double niuFunc2B = 0.0;

	double niuFunc3 = 0.5*C*(b > 0 ? b*b : 0.0) - niu3*niu3;
	double niuFunc3A = 0.0;
	double niuFunc3B = 0.5*C*(b > 0 ? 2.0*b : 0.0);//niuFunc3对a，b的求导

	double niuFunc4 = niu4*(a*a - b) + 0.5*C*(a*a - b)*(a*a - b);
	double niuFunc4A = niu4*2.0*a + C*(a*a - b)*2.0*a;
	double niuFunc4B = niu4 - C*(a*a - b);//niuFunc4对a，b的求导

	Array1D<double> g(X.dim(), 0.0);
	g[0] = niuFunc1A + niuFunc2A + niuFunc3A + niuFunc4A;
	g[1] = niuFunc1B + niuFunc2B + niuFunc3B + niuFunc4B;
	return g;
}



Array1D<double> niuUpdate(Array1D<double>& X, Array1D<double>& niuPrevious){
	double C = 4.0;
	double a = X[0];
	double b = X[1];
	double niu1 = niuPrevious[0];
	double niu2 = niuPrevious[1];
	double niu3 = niuPrevious[2];
	double niu4 = niuPrevious[3];


	double tempF = niu1 + C*(a*a + b*b - 4);
	double niuFunc1 = (tempF > 0) ? tempF : 0.0 ;

	double niuFunc2 = a > 0 ? niu2+C*a : 0.0;
	double niuFunc3 = b > 0 ? niu3+C*b : 0.0;

	double niuFunc4 = niu4+ C*(a*a - b);

	
	Array1D<double> niuNew(niuPrevious.dim(), 0.0);

	niuNew[0] = niuFunc1;
	niuNew[1] = niuFunc2;
	niuNew[2] = niuFunc3;
	niuNew[3] = niuFunc4;
	return niuNew;
}


double Vf_test4(Array1D<double>& X){

	double a = X[0];
	double b = X[1];
	//	return X*X + x1*x1;
	return 10.0*(a - 1.0)*(a - 1.0) + (b + 1.0)*(b + 1.0)*(b + 1.0)*(b + 1.0);

}

double restraint(Array1D<double>& X, Array1D<double>& niu){
	double x1 = X[0];
	double x2 = X[1];
	double x3 = X[2];
	double d = X[3];//设计变量

	double C = 4.0;//增广拉格朗日的辅助变量
	double niu1 = niu[0];
	double niu2 = niu[1];
	double niu3 = niu[2];
	double niu4 = niu[3];
	double niu5 = niu[4];
	double niu6 = niu[5];
	double niu7 = niu[6];

	double restraintF1 = niu1 + C*(d - 0.3);//第一个不等式约束，并加以处理
	double niuFunc1 = 0.5*(((restraintF1 > 0) ? (restraintF1*restraintF1) : 0.0) - niu1*niu1) / C;//增广项

	double restraintF2 = niu2 + C*(1.8 - d);//第一个不等式约束，并加以处理
	double niuFunc2 = 0.5*(((restraintF2 > 0) ? (restraintF2*restraintF2) : 0.0) - niu2*niu2) / C;//增广项

	double restraintF3 = niu3 + C*(x1 - d - 0.3);//第一个不等式约束，并加以处理
	double niuFunc3 = 0.5*(((restraintF3 > 0) ? (restraintF3*restraintF3) : 0.0) - niu3*niu3) / C;//增广项

	double restraintF4 = niu4 + C*(x2 - d - 0.3);//第一个不等式约束，并加以处理
	double niuFunc4 = 0.5*(((restraintF4 > 0) ? (restraintF4*restraintF4) : 0.0) - niu4*niu4) / C;//增广项

	double restraintF5 = niu5 + C*(x3 - 0.3);//第一个不等式约束，并加以处理
	double niuFunc5 = 0.5*(((restraintF5 > 0) ? (restraintF5*restraintF5) : 0.0) - niu5*niu5) / C;//增广项

	double restraintF6 = niu6 + C*(asin(x3 / sqrt(x1*x1 + x2*x2 + x3*x3)) - M_PI / 9.0);//第一个不等式约束，并加以处理
	double niuFunc6 = 0.5*(((restraintF6 > 0) ? (restraintF6*restraintF6) : 0.0) - niu6*niu6) / C;//增广项

	double restraintF7 = niu7 + C*(M_PI / 4.0 - asin(x3 / sqrt(x1*x1 + x2*x2 + x3*x3)));//第一个不等式约束，并加以处理
	double niuFunc7 = 0.5*(((restraintF7 > 0) ? (restraintF7*restraintF7) : 0.0) - niu7*niu7) / C;//增广项

	return niuFunc1 + niuFunc2 + niuFunc3 + niuFunc4 + niuFunc5 + niuFunc6 + niuFunc7;
}

double relative_density(Array1D<double>& X, Array1D<double>& niu){
	double x1 = X[0];
	double x2 = X[1];
	double x3 = X[2];
	double d = X[3];


	double func_relative_density = M_PI*d*d*sqrt(x1*x1 + x2*x2 + x3*x3) / x1*x2*x3;//目标函数



	return func_relative_density + restraint(X,niu);

}
Array1D<double> delta_relative_density(Array1D<double>& X){
	double x1 = X[0];
	double x2 = X[1];
	double x3 = X[2];
	double d = X[3];

	Array1D<double> g(X.dim(), 0.0);
	g[0] = M_PI*d*d*(1.0 / sqrt(x1*x1 + x2*x2 + x3*x3) - sqrt(x1*x1 + x2*x2 + x3*x3) / x1*x1) / x2*x3;
	g[1] = M_PI*d*d*(1.0 / sqrt(x1*x1 + x2*x2 + x3*x3) - sqrt(x1*x1 + x2*x2 + x3*x3) / x2*x2) / x1*x3;
	g[2] = M_PI*d*d*(1.0 / sqrt(x1*x1 + x2*x2 + x3*x3) - sqrt(x1*x1 + x2*x2 + x3*x3) / x3*x3) / x1*x2;
	g[3] = M_PI*2.0*d*sqrt(x1*x1 + x2*x2 + x3*x3) / x1*x2*x3;
	return g;

}

double modulus(Array1D<double>& X, Array1D<double>& niu){
	double x1 = X[0];
	double x2 = X[1];
	double x3 = X[2];
	double d = X[3];
	double temp_length = sqrt(x1*x1 + x2*x2 + x3*x3);

	double A1 = 8.0*(x1*x1 + x2*x2)*temp_length / M_PI / (x3*x3 + d*d) + 8.0*pow(x1*x1 + x2*x2 + x3*x3, 1.5) / 3.0*M_PI*d*d*d*d;
	double temp_A2 = 4.0*(x1*x1 + x2*x2) / M_PI / (x3*x3 + d*d) + 4.0*(x1*x1 + x2*x2 + x3*x3) / (3.0*M_PI*d*d*d*d);
	double A2 = sqrt(x1*x1 + x2*x2) * temp_A2;
	double A3 = pow(x1*x1 + x2*x2 + x3*x3, 1.5) / (2.0*M_PI * (x3*x3 + d*d)) - temp_length / (M_PI*d*d);

	double A0 = (A2 + sqrt(A2*A2 - 4.0*A1*A3)) / A1;

	double E1 = 4.0*A0*(x1*x1 + x2*x2 + x3*x3)*sqrt(x1*x1 + x2*x2) / (3.0*M_PI*d*d*d*d);
	double E2 = (4.0*A0*sqrt(x1*x1 + x2*x2) - temp_length) / (M_PI*d*d);
	double E = x1*x2*(E1 - E2) / x3;

	return E + restraint(X,niu);
}
Array1D<double> delta_modulus(Array1D<double>& X){
	double x1 = X[0];
	double x2 = X[1];
	double x3 = X[2];
	double d = X[3];

	double temp_length = sqrt(x1*x1 + x2*x2 + x3*x3);
	double A1 = 8.0*(x1*x1 + x2*x2)*temp_length / M_PI / (x3*x3 + d*d) + 8.0*pow(x1*x1 + x2*x2 + x3*x3, 1.5) / (3.0*M_PI*d*d*d*d);
	double temp_A2 = 4.0*(x1*x1 + x2*x2) / M_PI / (x3*x3 + d*d) + 4.0*(x1*x1 + x2*x2 + x3*x3) / (3.0*M_PI*d*d*d*d);
	double A2 = sqrt(x1*x1 + x2*x2) * temp_A2;
	double A3 = pow(x1*x1 + x2*x2 + x3*x3, 1.5) / (2.0*M_PI * (x3*x3 + d*d)) - temp_length / (M_PI*d*d);

	double A0 = (A2 + sqrt(A2*A2 - 4.0*A1*A3)) / A1;
	//计算第一个
	double delta_A1_1 = 8.0*(2.0*x1*temp_length + x1*(x1*x1 + x2*x2) / temp_length) / (M_PI * (x3*x3 + d*d)) + 8.0*x1*sqrt(x1*x1 + x2*x2 + x3*x3) / (M_PI*d*d*d*d);
	double delta_A1_2 = 8.0*(2.0*x2*temp_length + x2*(x1*x1 + x2*x2) / temp_length) / (M_PI * (x3*x3 + d*d)) + 8.0*x2*sqrt(x1*x1 + x2*x2 + x3*x3) / (M_PI*d*d*d*d);
	double delta_A1_3 = 8.0*(x1*x1 + x2*x2)*(x3*(x3*x3 + d*d) / temp_length - 2.0*x3*temp_length) / (M_PI * (x3*x3 + d*d)) + 8.0*x3*sqrt(x1*x1 + x2*x2 + x3*x3) / (M_PI*d*d*d*d);
	double delta_A1_d = -16.0*(x1*x1 + x2*x2)*temp_length*d / (M_PI*(x3*x3 + d*d)*(x3*x3 + d*d)) - 32.0*(x1*x1 + x2*x2)*temp_length / (3.0*M_PI*d*d*d*d*d);


	double temp_delta_A2 = 4.0*(x1*x1 + x2*x2) / M_PI / (x3*x3 + d*d) + 4.0*(x1*x1 + x2*x2 + x3*x3) / (3.0*M_PI*d*d*d*d);
	double delta_A2_1 = x1*temp_delta_A2 + sqrt(x1*x1 + x2*x2)*(8.0*x1 / M_PI / (x3*x3 + d*d) + 8.0*x1/(3.0*M_PI*d*d*d*d));
	double delta_A2_2 = x2*temp_delta_A2 + sqrt(x1*x1 + x2*x2)*(8.0*x2 / M_PI / (x3*x3 + d*d) + 8.0*x2 / (3.0*M_PI*d*d*d*d));
	double delta_A2_3 = sqrt(x1*x1 + x2*x2)*(8.0*x3 / (3.0*M_PI*d*d*d*d) - 8.0*x3 *(x1*x1 + x2*x2) / (x3*x3 + d*d)*(x3*x3 + d*d));
	double delta_A2_d = -2.0*sqrt(x1*x1 + x2*x2)*4.0*(x1*x1 + x2*x2) / (M_PI*(x3*x3 + d*d)*(x3*x3 + d*d)) - 4.0*sqrt(x1*x1 + x2*x2)* 4.0*(x1*x1 + x2*x2 + x3*x3) / (3.0*M_PI*d*d*d*d*d);

	double delta_A3_1 = 1.5*x1*temp_length / (M_PI * (x3*x3 + d*d)) - x1 / (M_PI*d*d*temp_length);
	double delta_A3_2 = 1.5*x2*temp_length / (M_PI * (x3*x3 + d*d)) - x2 / (M_PI*d*d*temp_length);
	double delta_A3_3 = 1.5*x3*temp_length*(x3*x3 + d*d) - 2.0*x3*pow(x1*x1 + x2*x2 + x3*x3, 1.5) / (M_PI * (x3*x3 + d*d) * (x3*x3 + d*d)) - x3 / (M_PI*d*d*temp_length);
	double delta_A3_d = -pow(x1*x1 + x2*x2 + x3*x3, 1.5) / (M_PI * (x3*x3 + d*d) * (x3*x3 + d*d)) + 2.0*temp_length / (M_PI*d*d*d);



	double delta_A0_1 = (delta_A2_1 + (A2*delta_A2_1 - 2.0*delta_A1_1*A3 - 2.0*A1*delta_A3_1) / (sqrt(A2*A2 - 4.0*A1*A3))) / A1 - delta_A1_1*(A2 + sqrt(A2*A2 - 4.0*A1*A3)) / A1 / A1;
	double delta_A0_2 = (delta_A2_2 + (A2*delta_A2_2 - 2.0*delta_A1_2*A3 - 2.0*A1*delta_A3_2) / (sqrt(A2*A2 - 4.0*A1*A3))) / A1 - delta_A1_2*(A2 + sqrt(A2*A2 - 4.0*A1*A3)) / A1 / A1;
	double delta_A0_3 = (delta_A2_3 + (A2*delta_A2_3 - 2.0*delta_A1_3*A3 - 2.0*A1*delta_A3_3) / (sqrt(A2*A2 - 4.0*A1*A3))) / A1 - delta_A1_3*(A2 + sqrt(A2*A2 - 4.0*A1*A3)) / A1 / A1;
	double delta_A0_d = (delta_A2_d + (A2*delta_A2_d - 2.0*delta_A1_d*A3 - 2.0*A1*delta_A3_d) / (sqrt(A2*A2 - 4.0*A1*A3))) / A1 - delta_A1_d*(A2 + sqrt(A2*A2 - 4.0*A1*A3)) / A1 / A1;

	//计算第二个
	double delta_E1_1 = 4.0*(delta_A0_1*(x1*x1 + x2*x2 + x3*x3)*sqrt(x1*x1 + x2*x2) + 2.0*x1*A0*sqrt(x1*x1 + x2*x2) + x1*A0*(x1*x1 + x2*x2 + x3*x3) / sqrt(x1*x1 + x2*x2)) / (3.0*M_PI*d*d*d*d);
	double delta_E1_2 = 4.0*(delta_A0_2*(x1*x1 + x2*x2 + x3*x3)*sqrt(x1*x1 + x2*x2) + 2.0*x2*A0*sqrt(x1*x1 + x2*x2) + x2*A0*(x1*x1 + x2*x2 + x3*x3) / sqrt(x1*x1 + x2*x2)) / (3.0*M_PI*d*d*d*d);
	double delta_E1_3 = 4.0*(delta_A0_3*(x1*x1 + x2*x2 + x3*x3)*sqrt(x1*x1 + x2*x2) + 2.0*x3*A0*sqrt(x1*x1 + x2*x2) + x3*A0*(x1*x1 + x2*x2 + x3*x3) / sqrt(x1*x1 + x2*x2)) / (3.0*M_PI*d*d*d*d);

	double delta_E2_1 = 4.0*(delta_A0_1*sqrt(x1*x1 + x2*x2) + 2.0*A0*x1*sqrt(x1*x1 + x2*x2)) / (M_PI*d*d) - x1 / (M_PI*d*d*temp_length);
	double delta_E2_2 = 4.0*(delta_A0_2*sqrt(x1*x1 + x2*x2) + 2.0*A0*x2*sqrt(x1*x1 + x2*x2)) / (M_PI*d*d) - x2 / (M_PI*d*d*temp_length);
	double delta_E2_3 = 4.0*(delta_A0_3*sqrt(x1*x1 + x2*x2) + 2.0*A0*x3*sqrt(x1*x1 + x2*x2)) / (M_PI*d*d) - x3 / (M_PI*d*d*temp_length);

	double E1 = 4.0*A0*(x1*x1 + x2*x2 + x3*x3)*sqrt(x1*x1 + x2*x2) / (3.0*M_PI*d*d*d*d);
	double delta_E1_d = - 4.0*E1 / d;

	double E2 = (4.0*A0*sqrt(x1*x1 + x2*x2) - temp_length) / (M_PI*d*d);
	double delta_E2_d = - 2.0*E2 / d;

	double E = x1*x2*(E1 - E2) / x3;

	double delta_E_1 = x2*(E1 - E2) / x3 + x1*x2*(delta_E1_1 - delta_E2_1) / x3;
	double delta_E_2 = x1*(E1 - E2) / x3 + x1*x2*(delta_E1_2 - delta_E2_2) / x3;
	double delta_E_3 = x1*x2*(delta_E1_3 - delta_E2_3) / x3 - x1*x2*(E1 - E2) / (x3*x3);
	double delta_E_d = x1*x2*(delta_E1_d - delta_E2_d) / x3;



	Array1D<double> g(X.dim(), 0.0);
	g[0] = delta_E_1;
	g[1] = delta_E_2;
	g[2] = delta_E_3;
	g[3] = delta_E_d;
	return g;

}

double plastic_strength(Array1D<double>& X, Array1D<double>& niu){
	double x1 = X[0];
	double x2 = X[1];
	double x3 = X[2];
	double d = X[3];

	double temp_length = sqrt(x1*x1 + x2*x2 + x3*x3);
	double A1 = 8.0*(x1*x1 + x2*x2)*temp_length / M_PI / (x3*x3 + d*d) + 8.0*pow(x1*x1 + x2*x2 + x3*x3, 1.5) / (3.0*M_PI*d*d*d*d);
	double temp_A2 = 4.0*(x1*x1 + x2*x2) / M_PI / (x3*x3 + d*d) + 4.0*(x1*x1 + x2*x2 + x3*x3) / (3.0*M_PI*d*d*d*d);
	double A2 = sqrt(x1*x1 + x2*x2) * temp_A2;
	double A3 = pow(x1*x1 + x2*x2 + x3*x3, 1.5) / (2.0*M_PI * (x3*x3 + d*d)) - temp_length / (M_PI*d*d);

	double A0 = (A2 + sqrt(A2*A2 - 4.0*A1*A3)) / A1;

	return 4.0*x1*x2*temp_length*A0 / (M_PI*d*d*d) + restraint(X, niu);

}

Array1D<double> delta_plastic_strength(Array1D<double>& X){
	double x1 = X[0];
	double x2 = X[1];
	double x3 = X[2];
	double d = X[3];

	double temp_length = sqrt(x1*x1 + x2*x2 + x3*x3);
	double A1 = 8.0*(x1*x1 + x2*x2)*temp_length / M_PI / (x3*x3 + d*d) + 8.0*pow(x1*x1 + x2*x2 + x3*x3, 1.5) / (3.0*M_PI*d*d*d*d);
	double temp_A2 = 4.0*(x1*x1 + x2*x2) / M_PI / (x3*x3 + d*d) + 4.0*(x1*x1 + x2*x2 + x3*x3) / (3.0*M_PI*d*d*d*d);
	double A2 = sqrt(x1*x1 + x2*x2) * temp_A2;
	double A3 = pow(x1*x1 + x2*x2 + x3*x3, 1.5) / (2.0*M_PI * (x3*x3 + d*d)) - temp_length / (M_PI*d*d);

	double A0 = (A2 + sqrt(A2*A2 - 4.0*A1*A3)) / A1;

	double delta_A1_1 = 8.0*(2.0*x1*temp_length + x1*(x1*x1 + x2*x2) / temp_length) / (M_PI * (x3*x3 + d*d)) + 8.0*x1*sqrt(x1*x1 + x2*x2 + x3*x3) / (M_PI*d*d*d*d);
	double delta_A1_2 = 8.0*(2.0*x2*temp_length + x2*(x1*x1 + x2*x2) / temp_length) / (M_PI * (x3*x3 + d*d)) + 8.0*x2*sqrt(x1*x1 + x2*x2 + x3*x3) / (M_PI*d*d*d*d);
	double delta_A1_3 = 8.0*(x1*x1 + x2*x2)*(x3*(x3*x3 + d*d) / temp_length - 2.0*x3*temp_length) / (M_PI * (x3*x3 + d*d)) + 8.0*x3*sqrt(x1*x1 + x2*x2 + x3*x3) / (M_PI*d*d*d*d);
	double delta_A1_d = - 16.0*(x1*x1 + x2*x2)*temp_length*d / (M_PI*(x3*x3 + d*d)*(x3*x3 + d*d)) - 32.0*(x1*x1 + x2*x2)*temp_length / (3.0*M_PI*d*d*d*d*d);


	double temp_delta_A2 = 4.0*(x1*x1 + x2*x2) / M_PI / (x3*x3 + d*d) + 4.0*(x1*x1 + x2*x2 + x3*x3) / (3.0*M_PI*d*d*d*d);
	double delta_A2_1 = x1*temp_delta_A2 + sqrt(x1*x1 + x2*x2)*(8.0*x1 / M_PI / (x3*x3 + d*d) + 8.0*x1 / (3.0*M_PI*d*d*d*d));
	double delta_A2_2 = x2*temp_delta_A2 + sqrt(x1*x1 + x2*x2)*(8.0*x2 / M_PI / (x3*x3 + d*d) + 8.0*x2 / (3.0*M_PI*d*d*d*d));
	double delta_A2_3 = sqrt(x1*x1 + x2*x2)*(8.0*x3 / (3.0*M_PI*d*d*d*d) - 8.0*x3 *(x1*x1 + x2*x2) / (x3*x3 + d*d)*(x3*x3 + d*d));
	double delta_A2_d = - 2.0*sqrt(x1*x1 + x2*x2)*4.0*(x1*x1 + x2*x2) / (M_PI*(x3*x3 + d*d)*(x3*x3 + d*d)) - 4.0*sqrt(x1*x1 + x2*x2)* 4.0*(x1*x1 + x2*x2 + x3*x3) / (3.0*M_PI*d*d*d*d*d);

	double delta_A3_1 = 1.5*x1*temp_length / (M_PI * (x3*x3 + d*d)) - x1 / (M_PI*d*d*temp_length);
	double delta_A3_2 = 1.5*x2*temp_length / (M_PI * (x3*x3 + d*d)) - x2 / (M_PI*d*d*temp_length);
	double delta_A3_3 = 1.5*x3*temp_length*(x3*x3 + d*d) - 2.0*x3*pow(x1*x1 + x2*x2 + x3*x3, 1.5) / (M_PI * (x3*x3 + d*d) * (x3*x3 + d*d)) - x3 / (M_PI*d*d*temp_length);
	double delta_A3_d = - pow(x1*x1 + x2*x2 + x3*x3, 1.5) / (M_PI * (x3*x3 + d*d) * (x3*x3 + d*d)) + 2.0*temp_length / (M_PI*d*d*d);



	double delta_A0_1 = (delta_A2_1 + (A2*delta_A2_1 - 2.0*delta_A1_1*A3 - 2.0*A1*delta_A3_1) / (sqrt(A2*A2 - 4.0*A1*A3))) / A1 - delta_A1_1*(A2 + sqrt(A2*A2 - 4.0*A1*A3)) / A1 / A1;
	double delta_A0_2 = (delta_A2_2 + (A2*delta_A2_2 - 2.0*delta_A1_2*A3 - 2.0*A1*delta_A3_2) / (sqrt(A2*A2 - 4.0*A1*A3))) / A1 - delta_A1_2*(A2 + sqrt(A2*A2 - 4.0*A1*A3)) / A1 / A1;
	double delta_A0_3 = (delta_A2_3 + (A2*delta_A2_3 - 2.0*delta_A1_3*A3 - 2.0*A1*delta_A3_3) / (sqrt(A2*A2 - 4.0*A1*A3))) / A1 - delta_A1_3*(A2 + sqrt(A2*A2 - 4.0*A1*A3)) / A1 / A1;
	double delta_A0_d = (delta_A2_d + (A2*delta_A2_d - 2.0*delta_A1_d*A3 - 2.0*A1*delta_A3_d) / (sqrt(A2*A2 - 4.0*A1*A3))) / A1 - delta_A1_d*(A2 + sqrt(A2*A2 - 4.0*A1*A3)) / A1 / A1;


	Array1D<double> g(X.dim(), 0.0);
	g[0] = 4.0*x2*temp_length*A0 / (M_PI*d*d*d) + 4.0*x1*x2*x1 / (M_PI*d*d*d*temp_length) + 4.0*x1*x2*temp_length*delta_A0_1 / (M_PI*d*d*d);
	g[1] = 4.0*x1*temp_length*A0 / (M_PI*d*d*d) + 4.0*x1*x2*x2 / (M_PI*d*d*d*temp_length) + 4.0*x1*x2*temp_length*delta_A0_2 / (M_PI*d*d*d);
	g[2] = 4.0*x1*x2*x3 / (M_PI*d*d*d*temp_length) + 4.0*x1*x2*temp_length*delta_A0_3 / (M_PI*d*d*d);
	g[3] = 4.0*x1*x2*temp_length*delta_A0_d / (M_PI*d*d*d) - 12.0*x1*x2*temp_length*A0 / (M_PI*d*d*d*d);
	return g;

}



Array1D<double> Vg_test4(Array1D<double>& X){
	double a = X[0];
	double b = X[1];
	Array1D<double> g(X.dim(), 0.0);
	g[0] = 20.0*(a - 1.0);
	g[1] = 4.0*(b + 1.0)*(b + 1.0)*(b + 1.0);
	return g;
}
Array2D<double> H_test4(Array1D<double>& X){
	double a = X[0];
	double b = X[1];
	int d = X.dim();
	Array2D<double> h(d, d, 0.0);

	h[0][0] = 20.0;
	h[1][1] = 12.0*(b + 1.0)*(b + 1.0);
	h[0][1] = h[1][0] = 0.0;
	return h;
}

int main2(){
	Array1D<double> x(2, 0.0);
	x[0] = 0.0;
	x[1] = 0.0;
	Array1D<double> niu(4, 0.0);

	Gradient gra(Vf_test4,Vg_test4);
	cout << "********GRADIENT*******\n";
	double d = gra.Vsolve(x);
	cout << d << endl;
	cout << "*********BFGS*****\n";
	BFGS bfgs(Vf_test4, Vg_test4);
	cout << bfgs.Vsolve(x);


	Newton newt(Vf_test4, Vg_test4,H_test4);
	cout << "********NEWTON*******\n";
	cout << newt.Vsolve(x);

	Steepest_Decend sp(Vf_test4, Vg_test4);
	cout << "***************STEEP****** " << sp.Vsolve(x);

	AugmentedLagrange lagrange(Vf_test5, Vg_test5, niuUpdate);

	cout<<"*************lagrange "<<lagrange.Vsolve(x, niu);

	system("pause");
	return 0;

}
int main1(){
	Array1D<double> x(2, 0.0);
	x[0] = 1.0;
	x[1] = 1.0;
	Array1D<double> niu(4, 0.0);

	Gradient gra(Vf_test4, Vg_test4);
	cout << "********梯度法计算结果*******\n";
	double d = gra.Vsolve(x);
	cout << "\t" << d << endl;
	cout << "*********BFGS法计算结果*****\n";
	BFGS bfgs(Vf_test4, Vg_test4);
	cout << "\t" << bfgs.Vsolve(x) << endl;


	Newton newt(Vf_test4, Vg_test4, H_test4);
	cout << "********牛顿法计算结果*******\n";
	cout << "\t" << newt.Vsolve(x) << endl;

	Steepest_Decend sp(Vf_test4, Vg_test4);
	cout << "*********最速下降法计算结果****** \n";
	cout << "\t" << sp.Vsolve(x) << endl;

	AugmentedLagrange lagrange(Vf_test5, Vg_test5, niuUpdate);
	cout << "*********增广拉格朗日函数法计算结果****** \n";
	cout << "\t" << lagrange.Vsolve(x, niu) << endl;
	cout << -M_PI << endl;


	system("pause");
	return 0;

}

