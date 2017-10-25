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



Array1D<double> niuUpdate(Array1D<double>& X, Array1D<double>& niu){
	
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
	//约束函数对 X 的梯度

	double restraintF2 = niu2 + C*(1.8 - d);//第一个不等式约束，并加以处理

	double restraintF3 = niu3 + C*(x1 - d - 0.3);//第一个不等式约束，并加以处

	double restraintF4 = niu4 + C*(x2 - d - 0.3);//第一个不等式约束，并加以处理

	double restraintF5 = niu5 + C*(x3 - 0.3);//第一个不等式约束，并加以处理

	double restraintF6 = niu6 + C*(asin(x3 / sqrt(x1*x1 + x2*x2 + x3*x3)) - M_PI / 9.0);//第一个不等式约束，并加以处理

	double restraintF7 = niu7 + C*(M_PI / 4.0 - asin(x3 / sqrt(x1*x1 + x2*x2 + x3*x3)));//第一个不等式约束，并加以处理
	



	double niuFunc1 = (restraintF1 > 0) ? (restraintF1) : 0.0;

	double niuFunc2 = (restraintF2 > 0) ? (restraintF2) : 0.0;
	double niuFunc3 = (restraintF3 > 0) ? (restraintF3) : 0.0;

	double niuFunc4 = (restraintF4 > 0) ? (restraintF4) : 0.0;
	double niuFunc5 = (restraintF5 > 0) ? (restraintF5) : 0.0;
	double niuFunc6 = (restraintF6 > 0) ? (restraintF6) : 0.0;
	double niuFunc7 = (restraintF7 > 0) ? (restraintF7) : 0.0;
	
	Array1D<double> niuNew(niu.dim(), 0.0);

	niuNew[0] = niuFunc1;
	niuNew[1] = niuFunc2;
	niuNew[2] = niuFunc3;
	niuNew[3] = niuFunc4;
	niuNew[4] = niuFunc5;
	niuNew[5] = niuFunc6;
	niuNew[6] = niuFunc7;
	cout << "*************更新niu*************\n";
	cout << "niuNew[0] = niuFunc0 : " << niuFunc1 << endl;
	cout << "niuNew[1] = niuFunc1 : " << niuFunc2 << endl;
	cout << "niuNew[2] = niuFunc2 : " << niuFunc3 << endl;
	cout << "niuNew[3] = niuFunc3 : " << niuFunc4 << endl;
	cout << "niuNew[4] = niuFunc4 : " << niuFunc5 << endl;
	cout << "niuNew[5] = niuFunc5 : " << niuFunc6 << endl;
	cout << "niuNew[6] = niuFunc6 : " << niuFunc7 << endl;
	return niuNew;
}


Array1D<double> delta_restraint(Array1D<double>& X, Array1D<double>& niu){
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
	//约束函数对 X 的梯度

	double restraintF2 = niu2 + C*(1.8 - d);//第一个不等式约束，并加以处理

	double restraintF3 = niu3 + C*(x1 - d - 0.3);//第一个不等式约束，并加以处理

	double restraintF4 = niu4 + C*(x2 - d - 0.3);//第一个不等式约束，并加以处理

	double restraintF5 = niu5 + C*(x3 - 0.3);//第一个不等式约束，并加以处理

	double restraintF6 = niu6 + C*(asin(x3 / sqrt(x1*x1 + x2*x2 + x3*x3)) - M_PI / 9.0);//第一个不等式约束，并加以处理

	double restraintF7 = niu7 + C*(M_PI / 4.0 - asin(x3 / sqrt(x1*x1 + x2*x2 + x3*x3)));//第一个不等式约束，并加以处理

	double delta_restraintF6_x1 = -x1*x3 / (sqrt(x1*x1 + x2*x2)*(x1*x1 + x2*x2 + x3*x3));//计算第六个反三角正弦的倒数
	double delta_restraintF6_x2 = -x2*x3 / (sqrt(x1*x1 + x2*x2)*(x1*x1 + x2*x2 + x3*x3));
	double delta_restraintF6_x3 = 1 / sqrt(x1*x1 + x2*x2) - x3*x3 / (sqrt(x1*x1 + x2*x2)*(x1*x1 + x2*x2 + x3*x3));

	double g_restraint_1 = 0.5*((restraintF3 > 0) ? 2.0*C*restraintF3 : 0.0) / C + 0.5*((restraintF6 > 0) ? 2.0*C*delta_restraintF6_x1 : 0.0) / C + 0.5*((restraintF7 > 0) ? -C*delta_restraintF6_x1 : 0.0) / C;
	double g_restraint_2 = 0.5*((restraintF4 > 0) ? 2.0*C*restraintF4 : 0.0) / C + 0.5*((restraintF6 > 0) ? 2.0*C*delta_restraintF6_x2 : 0.0) / C + 0.5*((restraintF7 > 0) ? -C*delta_restraintF6_x2 : 0.0) / C;
	double g_restraint_3 = 0.5*((restraintF5 > 0) ? 2.0*C*restraintF5 : 0.0) / C + 0.5*((restraintF6 > 0) ? 2.0*C*delta_restraintF6_x3 : 0.0) / C + 0.5*((restraintF7 > 0) ? -C*delta_restraintF6_x3 : 0.0) / C;
	double g_restraint_d = 0.5*((restraintF1 > 0) ? 2.0*C*restraintF1 : 0.0) / C + 0.5*(restraintF2 > 0) ? (-2.0*C*restraintF2) : 0.0 / C +
		0.5*((restraintF3 > 0) ? -2.0*C*restraintF3 : 0.0) / C + 0.5*((restraintF4 > 0) ? -2.0*C*restraintF4 : 0.0) / C;

	Array1D<double> g(X.dim(), 0.0);
	g[0] = g_restraint_1;
	g[1] = g_restraint_2;
	g[2] = g_restraint_3;
	g[3] = g_restraint_d;
//	cout << "**********计算约束的梯度***********\n";
	/**
		cout << "g[0] = g_restraint_1 : " << g_restraint_1 << endl;
		cout << "g[1] = g_restraint_2 : " << g_restraint_2 << endl;
		cout << "g[2] = g_restraint_3 : " << g_restraint_3 << endl;
		cout << "g[3] = g_restraint_d : " << g_restraint_d << endl;
	
	*/

	return g;
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
	//约束函数对 X 的梯度

	double restraintF2 = niu2 + C*(1.8 - d);//第一个不等式约束，并加以处理
	double niuFunc2 = 0.5*(((restraintF2 > 0) ? (restraintF2*restraintF2) : 0.0) - niu2*niu2) / C;//增广项


	double restraintF3 = niu3 + C*(x1 - d - 0.3);//第一个不等式约束，并加以处理
	double niuFunc3 = 0.5*(((restraintF3 > 0) ? (restraintF3*restraintF3) : 0.0) - niu3*niu3) / C;//增广项


	double restraintF4 = niu4 + C*(x2 - d - 0.3);//第一个不等式约束，并加以处理
	double niuFunc4 = 0.5*(((restraintF4 > 0) ? (restraintF4*restraintF4) : 0.0) - niu4*niu4) / C;//增广项

	double restraintF5 = niu5 + C*(x3 - 0.6);//第一个不等式约束，并加以处理
	double niuFunc5 = 0.5*(((restraintF5 > 0) ? (restraintF5*restraintF5) : 0.0) - niu5*niu5) / C;//增广项


	double restraintF6 = niu6 + C*(asin(x3 / sqrt(x1*x1 + x2*x2 + x3*x3)) - M_PI / 9.0);//第一个不等式约束，并加以处理
	double niuFunc6 = 0.5*(((restraintF6 > 0) ? (restraintF6*restraintF6) : 0.0) - niu6*niu6) / C;//增广项

	double restraintF7 = niu7 + C*(M_PI / 4.0 - asin(x3 / sqrt(x1*x1 + x2*x2 + x3*x3)));//第一个不等式约束，并加以处理
	double niuFunc7 = 0.5*(((restraintF7 > 0) ? (restraintF7*restraintF7) : 0.0) - niu7*niu7) / C;//增广项
//	cout << "*********计算约束************\n";
	/************************************************************************/
	/* 
		cout << "niuFunc1 = " << niuFunc1 << endl;
		cout << "niuFunc2 = " << niuFunc2 << endl;
		cout << "niuFunc3 = " << niuFunc3 << endl;
		cout << "niuFunc4 = " << niuFunc4 << endl;
		cout << "niuFunc5 = " << niuFunc5 << endl;
		cout << "niuFunc6 = " << niuFunc6 << endl;
		cout << "niuFunc7 = " << niuFunc7 << endl;
	*/
	/************************************************************************/




	return niuFunc1 + niuFunc2 + niuFunc3 + niuFunc4 + niuFunc5 + niuFunc6 + niuFunc7;
}



double relative_density(Array1D<double>& X, Array1D<double>& niu){
	double x1 = X[0];
	double x2 = X[1];
	double x3 = X[2];
	double d = X[3];

	double func_relative_density = M_PI*d*d*sqrt(x1*x1 + x2*x2 + x3*x3) / (x1*x2*x3);//目标函数
//	cout << "*******计算想对密度*******\n";
//	cout <<"func_relative_density + restraint(X,niu) = "<< func_relative_density + restraint(X, niu) << endl;
	cout << "$$$$$$$$$$$$$$$$$$$ restraint(X,niu) = " << restraint(X, niu) << endl;
	return func_relative_density + restraint(X,niu);
}
Array1D<double> delta_relative_density(Array1D<double>& X, Array1D<double>& niu){
	double x1 = X[0];
	double x2 = X[1];
	double x3 = X[2];
	double d = X[3];

	Array1D<double> g(X.dim(), 0.0);
	Array1D<double> delta_res(X.dim(), 0.0);
	delta_res = delta_restraint(X, niu);
	g[0] = M_PI*d*d*(1.0 / sqrt(x1*x1 + x2*x2 + x3*x3) - sqrt(x1*x1 + x2*x2 + x3*x3) / x1/x1) / (x2*x3);
	g[1] = M_PI*d*d*(1.0 / sqrt(x1*x1 + x2*x2 + x3*x3) - sqrt(x1*x1 + x2*x2 + x3*x3) / x2/x2) / (x1*x3);
	g[2] = M_PI*d*d*(1.0 / sqrt(x1*x1 + x2*x2 + x3*x3) - sqrt(x1*x1 + x2*x2 + x3*x3) / x3/x3) / (x1*x2);
	g[3] = M_PI*2.0*d*sqrt(x1*x1 + x2*x2 + x3*x3) / (x1*x2*x3);
//	cout << "*******计算相对密度梯度*******\n";
//	cout << "g[0] + delta_restraint(X,niu) = " << g[0] << " + " << delta_res[0]<< endl;
	/*
	cout << "g[1] + delta_restraint(X,niu) = " << g[1] << " + " << delta_restraint(X, niu) << g[1] + delta_restraint(X, niu) << endl;
	cout << "g[2] + delta_restraint(X,niu) = " << g[2] << " + " << delta_restraint(X, niu) << g[2] + delta_restraint(X, niu) << endl;
	cout << "g[3] + delta_restraint(X,niu) = " << g[3] << " + " << delta_restraint(X, niu) << g[3] + delta_restraint(X, niu) << endl;
                                                                     
	/************************************************************************/
//	cout << g + delta_restraint(X, niu);
	g = g + delta_res;
//	cout << g;
	return g ;
}

double modulus(Array1D<double>& X, Array1D<double>& niu){
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

	double E1 = 4.0*A0*(x1*x1 + x2*x2 + x3*x3)*sqrt(x1*x1 + x2*x2) / (3.0*M_PI*d*d*d*d);
	double E2 = (4.0*A0*sqrt(x1*x1 + x2*x2) - temp_length) / (M_PI*d*d);
	double E = x1*x2*(E1 - E2) / x3;
	cout << "*********计算初始刚度***********\n";
//	cout << "E + restraint(X,niu) = " << E + restraint(X, niu) << endl;
	return E + restraint(X,niu);
}
Array1D<double> delta_modulus(Array1D<double>& X, Array1D<double>& niu){
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
	//计算A1对x1,x2,x3,d的偏导
	double delta_A1_1 = 8.0*(2.0*x1*temp_length + x1*(x1*x1 + x2*x2) / temp_length) / (M_PI * (x3*x3 + d*d)) + 8.0*x1*temp_length / (M_PI*d*d*d*d);
	double delta_A1_2 = 8.0*(2.0*x2*temp_length + x2*(x1*x1 + x2*x2) / temp_length) / (M_PI * (x3*x3 + d*d)) + 8.0*x2*temp_length / (M_PI*d*d*d*d);
	double delta_A1_3 = 8.0*(x1*x1 + x2*x2)*(x3*(x3*x3 + d*d) / temp_length - 2.0*x3*temp_length) / (M_PI * (x3*x3 + d*d)) + 8.0*x3*temp_length / (M_PI*d*d*d*d);
//	double delta_A1_d = - 16.0*(x1*x1 + x2*x2)*temp_length*d / (M_PI*(x3*x3 + d*d)*(x3*x3 + d*d)) - 32.0*(x1*x1 + x2*x2)*temp_length / (3.0*M_PI*d*d*d*d*d);
	double delta_A1_d = -16.0*(x1*x1 + x2*x2)*temp_length*d / (M_PI*(x3*x3 + d*d)*(x3*x3 + d*d)) - 32.0* pow(x1*x1 + x2*x2 + x3*x3, 1.5) / (3.0*M_PI*d*d*d*d*d);
	//计算A2对x1,x2,x3,d的偏导
	double temp_delta_A2 = 4.0*(x1*x1 + x2*x2) / M_PI / (x3*x3 + d*d) + 4.0*(x1*x1 + x2*x2 + x3*x3) / (3.0*M_PI*d*d*d*d);
	//改动delta_A2_1， / sqrt(x1*x1 + x2*x2)
	double delta_A2_1 = x1*temp_delta_A2 / sqrt(x1*x1 + x2*x2) + sqrt(x1*x1 + x2*x2)*(8.0*x1 / M_PI / (x3*x3 + d*d) + 8.0*x1 / (3.0*M_PI*d*d*d*d));
	double delta_A2_2 = x2*temp_delta_A2 / sqrt(x1*x1 + x2*x2) + sqrt(x1*x1 + x2*x2)*(8.0*x2 / M_PI / (x3*x3 + d*d) + 8.0*x2 / (3.0*M_PI*d*d*d*d));

	double delta_A2_3 = sqrt(x1*x1 + x2*x2)*(8.0*x3 / (3.0*M_PI*d*d*d*d) - 8.0*x3 *(x1*x1 + x2*x2) / (x3*x3 + d*d)*(x3*x3 + d*d));
	double delta_A2_d = - 2.0*d*sqrt(x1*x1 + x2*x2)*4.0*(x1*x1 + x2*x2) / (M_PI*(x3*x3 + d*d)*(x3*x3 + d*d)) - 4.0*d*sqrt(x1*x1 + x2*x2)* 4.0*(x1*x1 + x2*x2 + x3*x3) / (3.0*M_PI*d*d*d*d*d);
	//计算A3对x1,x2,x3,d的偏导
	double delta_A3_1 = 1.5*x1*temp_length / (M_PI * (x3*x3 + d*d)) - x1 / (M_PI*d*d*temp_length);
	double delta_A3_2 = 1.5*x2*temp_length / (M_PI * (x3*x3 + d*d)) - x2 / (M_PI*d*d*temp_length);
	//少了一对（）
	double delta_A3_3 = (1.5*x3*temp_length*(x3*x3 + d*d) - 2.0*x3*pow(x1*x1 + x2*x2 + x3*x3, 1.5)) / (M_PI * (x3*x3 + d*d) * (x3*x3 + d*d)) - x3 / (M_PI*d*d*temp_length);
	//少了个d
	double delta_A3_d = - d*pow(x1*x1 + x2*x2 + x3*x3, 1.5) / (M_PI * (x3*x3 + d*d) * (x3*x3 + d*d)) + 2.0*temp_length / (M_PI*d*d*d);



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
	cout << "********计算初始刚度梯度*********\n";
	/************************************************************************/
	/* 	cout << "delta_E_1 + delta_restraint(X,niu) = " << g[0] << " + " << delta_restraint(X, niu) << g[0] + delta_restraint(X, niu) << endl;
	cout << "delta_E_2 + delta_restraint(X,niu) = " << g[1] << " + " << delta_restraint(X, niu) << g[1] + delta_restraint(X, niu) << endl;
	cout << "delta_E_3 + delta_restraint(X,niu) = " << g[2] << " + " << delta_restraint(X, niu) << g[2] + delta_restraint(X, niu) << endl;
	cout << "delta_E_d + delta_restraint(X,niu) = " << g[3] << " + " << delta_restraint(X, niu) << g[3] + delta_restraint(X, niu) << endl;
                                                                     */
	/************************************************************************/

	Array1D<double> g(X.dim(), 0.0);
	g[0] = delta_E_1;
	g[1] = delta_E_2;
	g[2] = delta_E_3;
	g[3] = delta_E_d;
//	cout << g + delta_restraint(X, niu);
	return g + delta_restraint(X, niu);

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

	cout << "******计算塑性强度*******\n";
	cout << "A0 = " << A0 << endl;
	cout << "4.0*x1*x2*temp_length*A0 / (M_PI*d*d*d) + restraint(X, niu) = " << 4.0*x1*x2*temp_length*A0 / (M_PI*d*d*d) + restraint(X, niu) << endl;

	return 4.0*x1*x2*temp_length*A0 / (M_PI*d*d*d) + restraint(X, niu);

}

Array1D<double> delta_plastic_strength(Array1D<double>& X, Array1D<double>& niu){
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

	//计算A1对x1,x2,x3,d的偏导
	double delta_A1_1 = 8.0*(2.0*x1*temp_length + x1*(x1*x1 + x2*x2) / temp_length) / (M_PI * (x3*x3 + d*d)) + 8.0*x1*temp_length / (M_PI*d*d*d*d);
	double delta_A1_2 = 8.0*(2.0*x2*temp_length + x2*(x1*x1 + x2*x2) / temp_length) / (M_PI * (x3*x3 + d*d)) + 8.0*x2*temp_length / (M_PI*d*d*d*d);
	double delta_A1_3 = 8.0*(x1*x1 + x2*x2)*(x3*(x3*x3 + d*d) / temp_length - 2.0*x3*temp_length) / (M_PI * (x3*x3 + d*d)) + 8.0*x3*temp_length / (M_PI*d*d*d*d);
	//	double delta_A1_d = - 16.0*(x1*x1 + x2*x2)*temp_length*d / (M_PI*(x3*x3 + d*d)*(x3*x3 + d*d)) - 32.0*(x1*x1 + x2*x2)*temp_length / (3.0*M_PI*d*d*d*d*d);
	double delta_A1_d = -16.0*(x1*x1 + x2*x2)*temp_length*d / (M_PI*(x3*x3 + d*d)*(x3*x3 + d*d)) - 32.0* pow(x1*x1 + x2*x2 + x3*x3, 1.5) / (3.0*M_PI*d*d*d*d*d);
	//计算A2对x1,x2,x3,d的偏导
	double temp_delta_A2 = 4.0*(x1*x1 + x2*x2) / M_PI / (x3*x3 + d*d) + 4.0*(x1*x1 + x2*x2 + x3*x3) / (3.0*M_PI*d*d*d*d);
	//改动delta_A2_1， / sqrt(x1*x1 + x2*x2)
	double delta_A2_1 = x1*temp_delta_A2 / sqrt(x1*x1 + x2*x2) + sqrt(x1*x1 + x2*x2)*(8.0*x1 / M_PI / (x3*x3 + d*d) + 8.0*x1 / (3.0*M_PI*d*d*d*d));
	double delta_A2_2 = x2*temp_delta_A2 / sqrt(x1*x1 + x2*x2) + sqrt(x1*x1 + x2*x2)*(8.0*x2 / M_PI / (x3*x3 + d*d) + 8.0*x2 / (3.0*M_PI*d*d*d*d));

	double delta_A2_3 = sqrt(x1*x1 + x2*x2)*(8.0*x3 / (3.0*M_PI*d*d*d*d) - 8.0*x3 *(x1*x1 + x2*x2) / (x3*x3 + d*d)*(x3*x3 + d*d));
	double delta_A2_d = -2.0*d*sqrt(x1*x1 + x2*x2)*4.0*(x1*x1 + x2*x2) / (M_PI*(x3*x3 + d*d)*(x3*x3 + d*d)) - 4.0*d*sqrt(x1*x1 + x2*x2)* 4.0*(x1*x1 + x2*x2 + x3*x3) / (3.0*M_PI*d*d*d*d*d);
	//计算A3对x1,x2,x3,d的偏导
	double delta_A3_1 = 1.5*x1*temp_length / (M_PI * (x3*x3 + d*d)) - x1 / (M_PI*d*d*temp_length);
	double delta_A3_2 = 1.5*x2*temp_length / (M_PI * (x3*x3 + d*d)) - x2 / (M_PI*d*d*temp_length);
	//少了一对（）
	double delta_A3_3 = (1.5*x3*temp_length*(x3*x3 + d*d) - 2.0*x3*pow(x1*x1 + x2*x2 + x3*x3, 1.5)) / (M_PI * (x3*x3 + d*d) * (x3*x3 + d*d)) - x3 / (M_PI*d*d*temp_length);
	//少了个d
	double delta_A3_d = -d*pow(x1*x1 + x2*x2 + x3*x3, 1.5) / (M_PI * (x3*x3 + d*d) * (x3*x3 + d*d)) + 2.0*temp_length / (M_PI*d*d*d);



	double delta_A0_1 = (delta_A2_1 + (A2*delta_A2_1 - 2.0*delta_A1_1*A3 - 2.0*A1*delta_A3_1) / (sqrt(A2*A2 - 4.0*A1*A3))) / A1 - delta_A1_1*(A2 + sqrt(A2*A2 - 4.0*A1*A3)) / A1 / A1;
	double delta_A0_2 = (delta_A2_2 + (A2*delta_A2_2 - 2.0*delta_A1_2*A3 - 2.0*A1*delta_A3_2) / (sqrt(A2*A2 - 4.0*A1*A3))) / A1 - delta_A1_2*(A2 + sqrt(A2*A2 - 4.0*A1*A3)) / A1 / A1;
	double delta_A0_3 = (delta_A2_3 + (A2*delta_A2_3 - 2.0*delta_A1_3*A3 - 2.0*A1*delta_A3_3) / (sqrt(A2*A2 - 4.0*A1*A3))) / A1 - delta_A1_3*(A2 + sqrt(A2*A2 - 4.0*A1*A3)) / A1 / A1;
	double delta_A0_d = (delta_A2_d + (A2*delta_A2_d - 2.0*delta_A1_d*A3 - 2.0*A1*delta_A3_d) / (sqrt(A2*A2 - 4.0*A1*A3))) / A1 - delta_A1_d*(A2 + sqrt(A2*A2 - 4.0*A1*A3)) / A1 / A1;



	Array1D<double> g(X.dim(), 0.0);
	g[0] = 4.0*x2*temp_length*A0 / (M_PI*d*d*d) + 4.0*x1*x2*x1 / (M_PI*d*d*d*temp_length) + 4.0*x1*x2*temp_length*delta_A0_1 / (M_PI*d*d*d);
	g[1] = 4.0*x1*temp_length*A0 / (M_PI*d*d*d) + 4.0*x1*x2*x2 / (M_PI*d*d*d*temp_length) + 4.0*x1*x2*temp_length*delta_A0_2 / (M_PI*d*d*d);
	g[2] = 4.0*x1*x2*x3 / (M_PI*d*d*d*temp_length) + 4.0*x1*x2*temp_length*delta_A0_3 / (M_PI*d*d*d);
	g[3] = 4.0*x1*x2*temp_length*delta_A0_d / (M_PI*d*d*d) - 12.0*x1*x2*temp_length*A0 / (M_PI*d*d*d*d);
	/************************************************************************/
	/* 	cout << "delta_E_1 + delta_restraint(X,niu) = " << g[0] << " + " << delta_restraint(X, niu) << g[0] + delta_restraint(X, niu) << endl;
	cout << "delta_E_2 + delta_restraint(X,niu) = " << g[1] << " + " << delta_restraint(X, niu) << g[1] + delta_restraint(X, niu) << endl;
	cout << "delta_E_3 + delta_restraint(X,niu) = " << g[2] << " + " << delta_restraint(X, niu) << g[2] + delta_restraint(X, niu) << endl;
	cout << "delta_E_d + delta_restraint(X,niu) = " << g[3] << " + " << delta_restraint(X, niu) << g[3] + delta_restraint(X, niu) << endl;
                                                                     */
	/************************************************************************/
	cout << "********计算塑性强度梯度*******\n";
	cout << g + delta_restraint(X, niu);

	return g + delta_restraint(X, niu);

}

double evaluate_Func(Array1D<double>& x){

	Array1D<double> niu(7, 0.0);
	AugmentedLagrange lagrange1(relative_density,delta_relative_density, niuUpdate);
	AugmentedLagrange lagrange2(modulus, delta_modulus,	niuUpdate);
	AugmentedLagrange lagrange3(plastic_strength,delta_plastic_strength, niuUpdate);


	double dealPoint = 
		sqrt((relative_density(x, niu) -lagrange1.Vsolve(x, niu))*(relative_density(x, niu) -lagrange1.Vsolve(x, niu))
		+(modulus(x, niu) - lagrange2.Vsolve(x, niu))*(modulus(x, niu) - lagrange2.Vsolve(x, niu)) 
		+(plastic_strength(x, niu) - lagrange3.Vsolve(x,niu))*(plastic_strength(x, niu) - lagrange3.Vsolve(x, niu)));
	
	return dealPoint;
}
Array1D<double> delta_evaluate_Func(Array1D<double>& x){
	using namespace TNT;
	Array1D<double> niu1(7, 0.0);

	AugmentedLagrange lagrange1(relative_density,delta_relative_density, niuUpdate);
	cout << lagrange1.Vsolve(x, niu1)<<endl;
	AugmentedLagrange lagrange2(modulus, delta_modulus,niuUpdate);
	cout << lagrange2.Vsolve(x, niu1) << endl;
	AugmentedLagrange lagrange3(plastic_strength,delta_plastic_strength, niuUpdate);
	cout << lagrange3.Vsolve(x, niu1) << endl;

	Array1D<double> g(x.dim(), 0.0);
	g = delta_relative_density(x, niu1)*(relative_density(x,niu1) - lagrange1.Vsolve(x, niu1)) 
		+ delta_modulus(x, niu1)*(modulus(x, niu1) - lagrange2.Vsolve(x, niu1)) 
		+ delta_plastic_strength(x,niu1)*(plastic_strength(x, niu1) - lagrange3.Vsolve(x, niu1));
	return g / evaluate_Func(x);
}

int main(){
	Array1D<double> x(4, 0.0);
	x[0] = 0.6;
	x[1] = 0.6;
	x[2] = 0.6;
	x[3] = 0.3;

//	Gradient gra(evaluate_Func,delta_evaluate_Func);
//	cout<<gra.Vsolve(x);

	Array1D<double> niu1(7, 0.0);
	AugmentedLagrange lagrange1(modulus, delta_modulus, niuUpdate);
	cout << lagrange1.Vsolve(x, niu1) << endl;

	system("pause");
	return 0;

}

