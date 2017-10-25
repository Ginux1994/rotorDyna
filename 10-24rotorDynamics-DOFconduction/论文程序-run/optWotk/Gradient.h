#ifndef GRADIENT_H
#define GTADIENT_H

#pragma once
#include "TNT/tnt_array1d.h"
using namespace TNT;

class Gradient
{
public:
	Gradient();
	Gradient(double(*f_test_)(Array1D<double>&), Array1D<double>(*g_test_)(Array1D<double>&));
	Gradient(double(*f_test_)(Array1D<double>&, Array1D<double>&), Array1D<double>(*g_test_)(Array1D<double>&, Array1D<double>&));

	double Vsolve(Array1D<double>& Vx_current);
	double LangrangeSolve(Array1D<double>& Vx_current, Array1D<double>& niu);
protected:
public:
	/************************************************************************/
	/* WOLF一维搜索                                                                     */
	/************************************************************************/
	wolf_research wolf;
	/************************************************************************/
	/* 函数参数                                                                     */
	/************************************************************************/
	double(*f_test)(double&);
	double(*g_test)(double&);

	double(*Vf_test)(Array1D<double>&);
	Array1D<double>(*Vg_test)(Array1D<double>&);

	//Lagrange augment function
	double(*Lf_test)(Array1D<double>&, Array1D<double>&);
	Array1D<double>(*Lg_test)(Array1D<double>&, Array1D<double>&);
	//输出参数
	Array1D<double> Vx_optimal,Vx_current, Vg_current, Vd_current;
	double x_optimal;
	double f_optimal;
};

Gradient::Gradient(){

}

Gradient::Gradient(double(*f_test_)(Array1D<double>&), Array1D<double>(*g_test_)(Array1D<double>&)){
	Vf_test = f_test_;
	Vg_test = g_test_;


	wolf = wolf_research(f_test_, g_test_);
}
Gradient::Gradient(double(*f_test_)(Array1D<double>&, Array1D<double>&), Array1D<double>(*g_test_)(Array1D<double>&, Array1D<double>&)){
	Lf_test = f_test_;
	Lg_test = g_test_;


	wolf = wolf_research(f_test_, g_test_);
}


double Gradient::Vsolve(Array1D<double>& Vx_initial){

	int dim = Vx_initial.dim();

	Vx_optimal = Array1D<double>(dim, 0.0);
	f_optimal = 0.0;

	Array1D<double> Vx_current = Vx_initial;
	Array1D<double> Vg_current = Vg_test(Vx_current);
	

	Array1D<double> Vd_current = -Vg_current;
	
	double tolerance = 1e-5;


	wolf.Vsolve(Vx_current, Vd_current);
	int k = 0;
	Array1D<double> Vg_previous(dim, 0.0);
	Array1D<double> Vd_previous(dim, 0.0);
	double beta_current;
	while ((wolf.Vx_next - Vx_current).norm() > tolerance)
	{
		k++;
		Vg_previous = Vg_current;
		Vd_previous = Vd_current;

		Vx_current = wolf.Vx_next;
		Vg_current = Vg_test(Vx_current);
		
		if (k%dim==0)
		{
			Vd_current = -Vg_current;
		}
		else{
			beta_current = (Vg_current*Vg_current) / (Vd_previous*(Vg_current - Vg_previous));
			for (int i = 0; i < dim;i++)
			{
				Vd_current[i] = Vd_previous[i] * beta_current - Vg_current[i];
			}

			//Vd_current = -Vg_current + Vd_previous*beta_current;
		}

		wolf.Vsolve(Vx_current, Vd_current);
		if (Vg_test(wolf.Vx_next).norm()<tolerance) break;//梯度为0是退出
		cout << "++++++++++++第 " << k << " 次迭代++++++++";
	}

	Vx_optimal = wolf.Vx_next;
	f_optimal = wolf.f_next;


	return f_optimal;



}

double Gradient::LangrangeSolve(Array1D<double>& Vx_initial, Array1D<double>& Vniu_initial){
	using namespace TNT;
	int dim = Vx_initial.dim1();

	Vx_optimal = Array1D<double>(4, 0.0);
	Vx_current = Array1D<double>(4, 0.0);
	Vg_current = Array1D<double>(4, 0.0);
	Vd_current = Array1D<double>(4, 0.0);
	f_optimal = 0.0;

	Vx_current = Vx_initial;
	Vg_current = Lg_test(Vx_current, Vniu_initial);

	cout << "Vg_current = " << Vg_current;

	Vd_current = -Vg_current;

	double tolerance = 1e-5;


	wolf.LagrangeSolve(Vx_current, Vd_current, Vniu_initial);
	int k = 0;
	Array1D<double> Vg_previous(dim, 0.0);
	Array1D<double> Vd_previous(dim, 0.0);
	double beta_current;
	while ((wolf.Vx_next - Vx_current).norm() < tolerance)
	{
		k++;
		Vg_previous = Vg_current;
		Vd_previous = Vd_current;

		Vx_current = wolf.Vx_next;
		Vg_current = Lg_test(Vx_current, Vniu_initial);

		if (k%dim == 0)
		{
			Vd_current = -Vg_current;
		}
		else{
			beta_current = (Vg_current*Vg_current) / (Vd_previous*(Vg_current - Vg_previous));
			for (int i = 0; i < dim; i++)
			{
				Vd_current[i] = Vd_previous[i] * beta_current - Vg_current[i];
			}
//			Vd_current = -Vg_current + Vd_previous*beta_current;
		}

		wolf.LagrangeSolve(Vx_current, Vd_current, Vniu_initial);
		if (Lg_test(wolf.Vx_next, Vniu_initial).norm() < tolerance) break;//梯度为0是退出
	}

	Vx_optimal = wolf.Vx_next;
	f_optimal = wolf.f_next;


	return f_optimal;



}






















#endif