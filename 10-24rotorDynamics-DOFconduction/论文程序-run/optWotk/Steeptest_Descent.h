#ifndef STEEPTEST_DESCENT_H
#define STEEPTEST_DESCENT_H

#include "TNT/tnt_array1d.h"
#include "wolf_research.h"
using namespace std;
using namespace TNT;
class Steepest_Decend
{
public:

	Steepest_Decend(double(*f_test_)(double&), double(*g_test_)(double&));
	Steepest_Decend(double(*f_test_)(Array1D<double>&), Array1D<double>(*g_test_)(Array1D<double>&));

	void solve(double x_current);
	double Vsolve(Array1D<double>& Vx_current);

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

	Array1D<double> Vx_optimal;
	double x_optimal;
	double f_optimal;
	
};

Steepest_Decend::Steepest_Decend(double(*f_test_)(Array1D<double>&), Array1D<double>(*g_test_)(Array1D<double>&)){
	Vf_test = f_test_;
	Vg_test = g_test_;


	wolf = wolf_research(f_test_, g_test_);
}
Steepest_Decend::Steepest_Decend(double(*f_test_)(double&), double(*g_test_)(double&)){
	f_test = f_test_;
	g_test = g_test_;


	wolf = wolf_research(f_test, g_test);
}


void Steepest_Decend::solve(double x_initial){
	double x_current = x_initial;
	double g_current = g_test(x_current);
	double d_current = -g_current;
	double tolerance = 1e-6;


	wolf.solve(x_current, d_current);
	int k = 0;
	while (abs(wolf.x_next-x_current)<tolerance)
	{
		k++;
		x_current = wolf.x_next;
		g_current = g_test(x_current);

		d_current = -g_current;
		wolf.solve(x_current, d_current);

	}
	x_optimal = wolf.x_next;
	f_optimal = wolf.f_next;
}


double Steepest_Decend::Vsolve(Array1D<double>& Vx_initial){

	Vx_optimal = Array1D<double>(Vx_initial.dim(), 0.0);
	f_optimal = 0.0;

	Array1D<double> Vx_current = Vx_initial;
	Array1D<double> Vg_current = Vg_test(Vx_current);

	Array1D<double> Vd_current = -Vg_current;
	double tolerance = 1e-4;


	wolf.Vsolve(Vx_current, Vd_current);
	int k = 0;
	while ((wolf.Vx_next - Vx_current).norm() > tolerance)
	{
		k++;
		//temp = (wolf.Vx_next - Vx_current).norm();
		

		Vx_current = wolf.Vx_next;
		Vg_current = Vg_test(Vx_current);
		if (Vg_current.norm()<tolerance) break;//梯度为0是退出
		Vd_current = -Vg_current;

		wolf.Vsolve(Vx_current, Vd_current);
		if (Vg_test(wolf.Vx_next).norm()<tolerance) break;//梯度为0是退出
	}

	Vx_optimal = wolf.Vx_next;
	f_optimal = wolf.f_next;


	return f_optimal;



}





























#endif