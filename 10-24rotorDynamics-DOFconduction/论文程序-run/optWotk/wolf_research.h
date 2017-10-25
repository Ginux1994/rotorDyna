#ifndef WOLF_RESEARCH_H
#define WOLF_RESEARCH_H

#include <iostream>
#include <vector>

#include "TNT/tnt_array1d.h"
#include "TNT/tnt_array1d_utils.h"

using namespace TNT;
using namespace std;
class wolf_research
{
public:

	wolf_research();
	wolf_research(double(*f_test)(double&), double(*g_test)(double&));
	wolf_research(double(*f_test)(Array1D<double>&), Array1D<double>(*g_test)(Array1D<double>&));
	//FOR LAGRANG
	wolf_research(double(*f_test)(Array1D<double>&, Array1D<double>&), Array1D<double>(*g_test)(Array1D<double>&, Array1D<double>&));

	void solve(double x_current, double d_current, double rho, double sigma);
	double Vsolve(Array1D<double>& Vx_current, Array1D<double>& Vd_current, double rho, double sigma);
	double LagrangeSolve(Array1D<double>& Vx_current, Array1D<double>& Vd_current, Array1D<double>& Vniu_initial, double rho, double sigma);
protected:
public:

	int dim;
	double x_current;//x在向量空间中的当前点
	Array1D<double> Vx_current;
	double d_current;//函数f_test在x_current下的下降搜索方向，作为已知输入
	Array1D<double> Vd_current;

	/************************************************************************/
	/* 函数参数                                                                     */
	/************************************************************************/
	double(*f_test)(double&);
	double(*g_test)(double&);
	double(*Vf_test)(Array1D<double>&);
	Array1D<double>(*Vg_test)(Array1D<double>&);

	/************************************************************************/
	/*			FLR LAGRANGE                                                                     */
	/************************************************************************/
	double (*Lf_test)(Array1D<double>&, Array1D<double>&);
	Array1D<double> (*Lg_test)(Array1D<double>&, Array1D<double>&);

	double f_current;//当前点的函数值，及其梯度值
	double g_current;
	Array1D<double> Vg_current;//梯度向量

	/************************************************************************/
	/* output param                                                                     */
	/************************************************************************/
	double alpha_acceptable;//return alpha，可接受的输出步长
	double x_next, f_next;
	Array1D<double> Vx_next;
};

wolf_research::wolf_research(){

}

wolf_research::wolf_research(double(*f_test_)(double&), double(*g_test_)(double&)){

	f_test = f_test_;
	g_test = g_test_;

	/************************************************************************/
	/* 函数参数                                                                     */
	/************************************************************************/
	f_current=0.0;//当前点的函数值，及其梯度值
	g_current=0.0;

	alpha_acceptable = 0.0;
	x_next = 0.0;
	f_next = 0.0;

}
wolf_research::wolf_research(double(*f_test_)(Array1D<double>&), Array1D<double>(*g_test_)(Array1D<double>&)){
	
	Vf_test = f_test_;
	Vg_test = g_test_;

	f_current = 0.0;//当前点的函数值，及其梯度值


	alpha_acceptable = 0.0;

	f_next = 0.0;

}
wolf_research::wolf_research(double(*f_test_)(Array1D<double>&, Array1D<double>&), Array1D<double>(*g_test_)(Array1D<double>&, Array1D<double>&)){

	Lf_test = f_test_;
	Lg_test = g_test_;

	f_current = 0.0;//当前点的函数值，及其梯度值


	alpha_acceptable = 0.0;

	f_next = 0.0;

}

void wolf_research::solve(double x_current, double d_current, double rho=0.1, double sigma=0.11){//rho: 可接受系数[0,1]; sigma: 可接受点处切线斜率大雨初始点处切线斜率的倍数, 0<rho<sigma<1

	
	f_current = f_test(x_current);//当前点的函数值，及其梯度值
	g_current = g_test(x_current);

	double x_alpha_k, f_alpha_k;//进步后的x值及其函数值，梯度值，以及梯度与搜索方向之积
	double g_alpha_k, df_alpha_k;

	double f_alpha_lower_k = f_current;//经过两点外插后更新步长后，迭代到x_alpha_k时的函数值，及其梯度，现在作为初始值
	double g_alpha_lower_k = g_current;
	double f_alpha_lower_0 = f_alpha_lower_k;
	//迭代到k次时，目标函数梯度与方向之积，用于计算
	double df_alpha_lower_k = d_current * g_alpha_lower_k;//d_current.transpose();向量与常数的点乘，

	double df_alpha_lower_0 = df_alpha_lower_k;

	double alpha_initial;//初始步长

	/************************************************************************/
	/* 中间变量                                                                     */
	/************************************************************************/
	double delta_alpha_k, alpha_k_temp;

	double tolerance = 1e-15;

	/************************************************************************/
	/* 估计alpha_initial初值，防止过大过小                                                                     */
	/************************************************************************/
	if (abs(df_alpha_lower_k) > tolerance)
	{
		alpha_initial = -2 * f_alpha_lower_k / df_alpha_lower_k;

	}
	else{
		alpha_initial = 1;
	}

	if (alpha_initial < tolerance)
	{
		alpha_initial = 1;
	}

	/************************************************************************/
	/*				步长参数                                                                     */
	/************************************************************************/

	double alpha_lower_k = 0.0;//使用两点外插公式更新后的左步长
	double alpha_upper_k = 1e8;//使用两点内插公式更新后的右步长

	double alpha_k = alpha_initial;//经过k次搜索后，获得的步长




	/************************************************************************/
	/* WOLF条件参数                                                                     */
	/************************************************************************/
	double wolfe_condition1;
	double wolfe_condition2;


	/************************************************************************/
	/* begin iteria                                                                     */
	/************************************************************************/

	int k_max = 1000;
	int k ;
	for (k = 1; k < k_max; k++)
	{
		x_alpha_k = x_current + alpha_k*d_current;
		/************************************************************************/
		/* erro                                                                     */
		/************************************************************************/
		f_alpha_k = f_test(x_alpha_k);//vector of x_alpha_k
		g_alpha_k = g_test(x_alpha_k);

		df_alpha_k = d_current*g_alpha_k;

		//wolf condition
		wolfe_condition1 = f_alpha_k - f_alpha_lower_0 - rho*alpha_k*df_alpha_lower_0;
		wolfe_condition2 = sigma*df_alpha_lower_0 - df_alpha_k;

		if (wolfe_condition1 <= 0)
		{
			if (wolfe_condition2 <= 0)
			{
				alpha_acceptable = alpha_initial;
				x_next = x_alpha_k;
				f_next = f_alpha_k;
				break;
			}
			else{//only sastify condition1
				delta_alpha_k = (alpha_k - alpha_lower_k)*df_alpha_k / (df_alpha_lower_k - df_alpha_k);//对alpha_lower_k, aloha_k插值

				if (delta_alpha_k <= 0)
				{
					alpha_k_temp = 2 * alpha_k;//外插特殊情况
				}
				else{
					alpha_k_temp = alpha_k + delta_alpha_k;
				}

				alpha_lower_k = alpha_k;//更新可接受不长的下边界
				f_alpha_lower_k = f_alpha_k;
				df_alpha_lower_k = df_alpha_k;
				alpha_k = alpha_k_temp;//更新其他值，换位
			}
		}
		else{
			//condition1 not satisfied, 对alpha_lower_k, aloha_k插值

			if (alpha_k < alpha_upper_k)
			{
				alpha_upper_k = alpha_k;//更新可接受不长上边界
			}

			alpha_k_temp = alpha_lower_k + 0.5*(((alpha_k - alpha_lower_k)*(alpha_k - alpha_lower_k)*df_alpha_lower_k) / (f_alpha_lower_k - f_alpha_k + (alpha_k - alpha_lower_k)*df_alpha_lower_k));

			alpha_k = alpha_k_temp;

		}


		if ((alpha_upper_k - alpha_lower_k) < tolerance)//防止在过小区间搜索
		{
			alpha_acceptable = alpha_k;
			x_next = x_alpha_k;
			f_next = f_alpha_k;
		}

	}//k<k_max

	if (wolfe_condition1 > 0 || wolfe_condition2 > 0)
	{
		std::cout << "failed \n";

	}

}

double wolf_research::Vsolve(Array1D<double>& Vx_current, Array1D<double>& Vd_current, double rho=0.1, double sigma=0.11){
	int k_max = 1000;
	int k = 0;
	dim = Vx_current.dim();

	Vx_next = Array1D<double>(dim, 0.0);
	double f_current = Vf_test(Vx_current);
	Array1D<double> Vg_current = Vg_test(Vx_current);

	double f_alpha_lower_k = f_current;
	Array1D<double> Vg_alpha_lower_k = Vg_current;

	double df_alpha_lower_k = Vd_current * Vg_alpha_lower_k;//d_current.transpose();向量与常数的点乘

	double f_alpha_lower_0 = f_alpha_lower_k;
	double df_alpha_lower_0 = df_alpha_lower_k;
	double alpha_initial = 0.0;

	double tolerance = 1e-15;
	/************************************************************************/
	/* 估计alpha_initial初值，防止过大过小                                                                     */
	/************************************************************************/
	if (abs(df_alpha_lower_k) > tolerance)
	{
		alpha_initial = -2 * f_alpha_lower_k / df_alpha_lower_k;

	}
	else{
		alpha_initial = 1;
	}

	if (alpha_initial < tolerance)
	{
		alpha_initial = 1;
	}

	double alpha_lower_k = 0.0;
	double alpha_upper_k = 1e8;
	double alpha_k = alpha_initial;
	/************************************************************************/
	/* output param      alpha_acceptable    Vx_next      f_next                                                      */
	/************************************************************************/
//	double alpha_acceptable = 0.0;//return alpha
//	Array1D<double> Vx_next(dim, 0.0);
//	double f_next = 0.0;
	double delta_alpha_k = 0.0;
	double alpha_k_temp=0.0;
	/************************************************************************/
	/* begin search                                                                     */
	/************************************************************************/
	double wolfe_condition1 = 0.0;
	double wolfe_condition2 = 0.0;
	Array1D<double> Vx_alpha_k(dim, 0.0);
	double f_alpha_k = 0.0;
	Array1D<double> Vg_alpha_k(dim, 0.0);
	double df_alpha_k = 0.0;

	for (k = 1; k < k_max; k++)
	{
		/************************************************************************/
		/* 常量与向量点乘                                                                     */
		/************************************************************************/
		for (int i = 0; i < dim;i++)
		{
			Vx_alpha_k[i] = Vx_current[i] + Vd_current[i] * alpha_k;
		}
//		Vx_alpha_k = Vx_current + Vd_current * alpha_k ;
		f_alpha_k = Vf_test(Vx_alpha_k);//vector of x_alpha_k
		Vg_alpha_k = Vg_test(Vx_alpha_k);
		/************************************************************************/
		/* 向量点乘                                                                     */
		/************************************************************************/
		df_alpha_k = Vd_current * Vg_alpha_k;

		//wolf condition
		wolfe_condition1 = f_alpha_k - f_alpha_lower_0 - rho*alpha_k*df_alpha_lower_0;
		wolfe_condition2 = sigma*df_alpha_lower_0 - df_alpha_k;

		if (wolfe_condition1 <= 0)
		{
			if (wolfe_condition2 <= 0)
			{
				alpha_acceptable = alpha_initial;
				Vx_next = Vx_alpha_k;
				f_next = f_alpha_k;
				break;
			}
			else{//only sastify condition1
				delta_alpha_k = (alpha_k - alpha_lower_k)*df_alpha_k / (df_alpha_lower_k - df_alpha_k);//对alpha_lower_k, aloha_k插值

				if (delta_alpha_k <= 0)
				{
					alpha_k_temp = 2 * alpha_k;//外插特殊情况
				}
				else{
					alpha_k_temp = alpha_k + delta_alpha_k;
				}

				alpha_lower_k = alpha_k;//更新可接受不长的下边界
				f_alpha_lower_k = f_alpha_k;
				df_alpha_lower_k = df_alpha_k;
				alpha_k = alpha_k_temp;//更新其他值，换位
			}
		}
		else{
			//condition1 not satisfied, 对alpha_lower_k, aloha_k插值

			if (alpha_k < alpha_upper_k)
			{
				alpha_upper_k = alpha_k;//更新可接受不长上边界
			}

			alpha_k_temp = alpha_lower_k + 0.5*(((alpha_k - alpha_lower_k)*(alpha_k - alpha_lower_k)*df_alpha_lower_k) / (f_alpha_lower_k - f_alpha_k + (alpha_k - alpha_lower_k)*df_alpha_lower_k));

			alpha_k = alpha_k_temp;

		}


		if ((alpha_upper_k - alpha_lower_k) < tolerance)//防止在过小区间搜索
		{
			alpha_acceptable = alpha_k;
			Vx_next = Vx_alpha_k;
			f_next = f_alpha_k;
		}

	}//k<k_max

	if (wolfe_condition1 > 0 || wolfe_condition2 > 0)
	{
		cout << "failed \n";

	}

	//return alpha_acceptable;
	return alpha_acceptable;




}

double wolf_research::LagrangeSolve(Array1D<double>& Vx_current, Array1D<double>& Vd_current, Array1D<double>& Vniu_initial, double rho = 0.1, double sigma = 0.11){
	int k_max = 1000;
	int k = 0;
	dim = Vx_current.dim();

	Vx_next = Array1D<double>(dim, 0.0);
	double f_current = Lf_test(Vx_current, Vniu_initial);
	Array1D<double> Vg_current = Lg_test(Vx_current, Vniu_initial);

	double f_alpha_lower_k = f_current;
	Array1D<double> Vg_alpha_lower_k = Vg_current;

	double df_alpha_lower_k = Vd_current * Vg_alpha_lower_k;//d_current.transpose();向量与常数的点乘

	double f_alpha_lower_0 = f_alpha_lower_k;
	double df_alpha_lower_0 = df_alpha_lower_k;
	double alpha_initial = 0.0;

	double tolerance = 1e-8;
	/************************************************************************/
	/* 估计alpha_initial初值，防止过大过小                                                                     */
	/************************************************************************/
	if (abs(df_alpha_lower_k) > tolerance)
	{
		alpha_initial = -2 * f_alpha_lower_k / df_alpha_lower_k;

	}
	else{
		alpha_initial = 1;
	}

	if (alpha_initial < tolerance)
	{
		alpha_initial = 1;
	}

	double alpha_lower_k = 0.0;
	double alpha_upper_k = 1e8;
	double alpha_k = alpha_initial;
	/************************************************************************/
	/* output param      alpha_acceptable    Vx_next      f_next                                                      */
	/************************************************************************/
	//	double alpha_acceptable = 0.0;//return alpha
	//	Array1D<double> Vx_next(dim, 0.0);
	//	double f_next = 0.0;
	double delta_alpha_k = 0.0;
	double alpha_k_temp = 0.0;
	/************************************************************************/
	/* begin search                                                                     */
	/************************************************************************/
	double wolfe_condition1 = 0.0;
	double wolfe_condition2 = 0.0;
	Array1D<double> Vx_alpha_k(dim, 0.0);
	double f_alpha_k = 0.0;
	Array1D<double> Vg_alpha_k(dim, 0.0);
	double df_alpha_k = 0.0;

	for (k = 1; k < k_max; k++)
	{
		/************************************************************************/
		/* 常量与向量点乘                                                                     */
		/************************************************************************/
		for (int i = 0; i < dim; i++)
		{
			Vx_alpha_k[i] = Vx_current[i] + Vd_current[i] * alpha_k;
		}
//		Vx_alpha_k = Vx_current + Vd_current * alpha_k;
		f_alpha_k = Lf_test(Vx_alpha_k, Vniu_initial);//vector of x_alpha_k
		Vg_alpha_k = Lg_test(Vx_alpha_k, Vniu_initial);
		/************************************************************************/
		/* 向量点乘                                                                     */
		/************************************************************************/
		df_alpha_k = Vd_current * Vg_alpha_k;

		//wolf condition
		wolfe_condition1 = f_alpha_k - f_alpha_lower_0 - rho*alpha_k*df_alpha_lower_0;
		wolfe_condition2 = sigma*df_alpha_lower_0 - df_alpha_k;

		if (wolfe_condition1 <= 0)
		{
			if (wolfe_condition2 <= 0)
			{
				alpha_acceptable = alpha_initial;
				Vx_next = Vx_alpha_k;
				f_next = f_alpha_k;
				break;
			}
			else{//only sastify condition1
				delta_alpha_k = (alpha_k - alpha_lower_k)*df_alpha_k / (df_alpha_lower_k - df_alpha_k);//对alpha_lower_k, aloha_k插值

				if (delta_alpha_k <= 0)
				{
					alpha_k_temp = 2 * alpha_k;//外插特殊情况
				}
				else{
					alpha_k_temp = alpha_k + delta_alpha_k;
				}

				alpha_lower_k = alpha_k;//更新可接受不长的下边界
				f_alpha_lower_k = f_alpha_k;
				df_alpha_lower_k = df_alpha_k;
				alpha_k = alpha_k_temp;//更新其他值，换位
			}
		}
		else{
			//condition1 not satisfied, 对alpha_lower_k, aloha_k插值

			if (alpha_k < alpha_upper_k)
			{
				alpha_upper_k = alpha_k;//更新可接受不长上边界
			}

			alpha_k_temp = alpha_lower_k + 0.5*(((alpha_k - alpha_lower_k)*(alpha_k - alpha_lower_k)*df_alpha_lower_k) / (f_alpha_lower_k - f_alpha_k + (alpha_k - alpha_lower_k)*df_alpha_lower_k));

			alpha_k = alpha_k_temp;

		}


		if ((alpha_upper_k - alpha_lower_k) < tolerance)//防止在过小区间搜索
		{
			alpha_acceptable = alpha_k;
			Vx_next = Vx_alpha_k;
			f_next = f_alpha_k;
			break;
		}

	}//k<k_max

	if (wolfe_condition1 > 0 || wolfe_condition2 > 0)
	{
//		cout << "failed \n";

	}

	//return alpha_acceptable;
	return alpha_acceptable;




}

























#endif