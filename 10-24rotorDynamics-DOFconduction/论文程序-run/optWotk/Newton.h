#ifndef NEWTON_H
#define NEWTON_H

#include "TNT/tnt_array1d.h"
#include "TNT/tnt_array2d.h"
#include "TNT/tnt_array2d_utils.h"
#include "JAMA/jama_lu.h"
#include "wolf_research.h"
using namespace std;
using namespace TNT;
using namespace JAMA;
class Newton
{
public:


	Newton(double(*f_test_) (Array1D<double>&), Array1D<double> (*g_test_)(Array1D<double>&), Array2D<double>(*H_test_)(Array1D<double>&));


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
	//输出参数
	Array1D<double> Vx_optimal;
	double x_optimal;
	double f_optimal;
	//附加嗨森矩阵
	Array2D<double> (*H_test)(Array1D<double>&);
	Array2D<double> H_current;

};

Newton::Newton(double(*f_test_)(Array1D<double>&), Array1D<double>(*g_test_)(Array1D<double>&), Array2D<double>(*H_test_)(Array1D<double>&)){
	Vf_test = f_test_;
	Vg_test = g_test_;
	H_test = H_test_;

	wolf = wolf_research(f_test_, g_test_);
}



double Newton::Vsolve(Array1D<double>& Vx_initial){

	Vx_optimal = Array1D<double>(Vx_initial.dim(), 0.0);
	f_optimal = 0.0;

	Array1D<double> Vx_current = Vx_initial;
	Array1D<double> Vg_current = Vg_test(Vx_current);
//	cout<<"当前梯度" << Vg_current;
	H_current = H_test(Vx_current);//当前向量的海森矩阵

	LU<double> inverse_H(H_current);
	Array2D<double> E(Vx_initial.dim(), Vx_initial.dim(), 0.0);
	for (int i = 0; i < Vx_initial.dim(); i++)
	{
		E[i][i] = 1.0;
	}

	Array2D<double> tempV = inverse_H.solve(E);
//	cout << tempV;
	Array1D<double> Vd_current = tempV*Vg_current;
	Vd_current = -Vd_current;
//	cout<<"当前梯度 \n" << Vd_current;
	double tolerance = 1e-5;


	wolf.Vsolve(Vx_current, Vd_current);
	int k = 0;
	while ((wolf.Vx_next - Vx_current).norm() > tolerance)
	{
		k++;
		//temp = (wolf.Vx_next - Vx_current).norm();

		Vx_current = wolf.Vx_next;
		Vg_current = Vg_test(Vx_current);
		if (Vg_current.norm()<tolerance) break;//梯度为0是退出

		H_current = H_test(Vx_current);

		Vd_current = inverse_H.solve(E)*Vg_current;
		Vd_current = -Vd_current;


		wolf.Vsolve(Vx_current, Vd_current);
		if (Vg_test(wolf.Vx_next).norm()<tolerance) break;//梯度为0是退出
	}

	Vx_optimal = wolf.Vx_next;
	f_optimal = wolf.f_next;


	return f_optimal;



}


























#endif