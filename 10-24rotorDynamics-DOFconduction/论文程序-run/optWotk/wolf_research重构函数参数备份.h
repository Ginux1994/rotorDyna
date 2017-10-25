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
	wolf_research(double x_current, double d_current, double rho, double sigma);
	wolf_research(Array1D<double>& x_current, Array1D<double>& d_currentV, double rho, double sigma);
	void solve(double(*f_test)(double&), double(*g_test)(double&));


	double Vsolve(double(*f_test)(Array1D<double>&), Array1D<double>(*g_test)(Array1D<double>&));
protected:
public:

	int dim;
	double x_current;//x�������ռ��еĵ�ǰ��
	Array1D<double> Vx_current;
	double d_current;//����f_test��x_current�µ��½�����������Ϊ��֪����
	Array1D<double> Vd_current;
	double rho;
	double sigma;//rho: �ɽ���ϵ��[0,1]; sigma: �ɽ��ܵ㴦����б�ʴ����ʼ�㴦����б�ʵı���, 0<rho<sigma<1

	/************************************************************************/
	/* ��������                                                                     */
	/************************************************************************/

	double f_current;//��ǰ��ĺ���ֵ�������ݶ�ֵ
	double g_current;
	Array1D<double> Vg_current;
//	Array1D<double> Vf_current;
	/************************************************************************/
	/* output param                                                                     */
	/************************************************************************/
	double alpha_acceptable;//return alpha���ɽ��ܵ��������
	double x_next, f_next;
	Array1D<double> Vx_next;
};
wolf_research::wolf_research(double x_current_, double d_current_, double rho_, double sigma_){

	x_current = x_current_;
	d_current = d_current_;
	sigma = sigma_;
	rho = rho_;

	/************************************************************************/
	/* ��������                                                                     */
	/************************************************************************/
	f_current=0.0;//��ǰ��ĺ���ֵ�������ݶ�ֵ
	g_current=0.0;

	alpha_acceptable = 0.0;
	x_next = 0.0;
	f_next = 0.0;

}
wolf_research::wolf_research(Array1D<double>& Vx_current_, Array1D<double>& Vd_current_, double rho_, double sigma_){
	dim = Vx_current_.dim();

	Vx_current = Vx_current_;
	Vd_current = Vd_current_;
	rho= rho_;
	sigma=sigma_;

	f_current = 0.0;//��ǰ��ĺ���ֵ�������ݶ�ֵ
	Vg_current = Array1D<double>(dim,0.0);

	alpha_acceptable = 0.0;
	Vx_next = Array1D<double>(dim, 0.0);
	f_next = 0.0;

}

void wolf_research::solve(double(*f_test)(double&), double(*g_test)(double&)){


	f_current = f_test(this->x_current);//��ǰ��ĺ���ֵ�������ݶ�ֵ
	g_current = g_test(this->x_current);

	double x_alpha_k, f_alpha_k;//�������xֵ���亯��ֵ���ݶ�ֵ���Լ��ݶ�����������֮��
	double g_alpha_k, df_alpha_k;

	double f_alpha_lower_k = f_current;//��������������²����󣬵�����x_alpha_kʱ�ĺ���ֵ�������ݶȣ�������Ϊ��ʼֵ
	double g_alpha_lower_k = g_current;
	double f_alpha_lower_0 = f_alpha_lower_k;
	//������k��ʱ��Ŀ�꺯���ݶ��뷽��֮�������ڼ���
	double df_alpha_lower_k = d_current * g_alpha_lower_k;//d_current.transpose();�����볣���ĵ�ˣ�

	double df_alpha_lower_0 = df_alpha_lower_k;

	double alpha_initial;//��ʼ����

	/************************************************************************/
	/* �м����                                                                     */
	/************************************************************************/
	double delta_alpha_k, alpha_k_temp;

	double tolerance = 1e-15;

	/************************************************************************/
	/* ����alpha_initial��ֵ����ֹ�����С                                                                     */
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
	/*				��������                                                                     */
	/************************************************************************/

	double alpha_lower_k = 0.0;//ʹ��������幫ʽ���º���󲽳�
	double alpha_upper_k = 1e8;//ʹ�������ڲ幫ʽ���º���Ҳ���

	double alpha_k = alpha_initial;//����k�������󣬻�õĲ���




	/************************************************************************/
	/* WOLF��������                                                                     */
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
				delta_alpha_k = (alpha_k - alpha_lower_k)*df_alpha_k / (df_alpha_lower_k - df_alpha_k);//��alpha_lower_k, aloha_k��ֵ

				if (delta_alpha_k <= 0)
				{
					alpha_k_temp = 2 * alpha_k;//����������
				}
				else{
					alpha_k_temp = alpha_k + delta_alpha_k;
				}

				alpha_lower_k = alpha_k;//���¿ɽ��ܲ������±߽�
				f_alpha_lower_k = f_alpha_k;
				df_alpha_lower_k = df_alpha_k;
				alpha_k = alpha_k_temp;//��������ֵ����λ
			}
		}
		else{
			//condition1 not satisfied, ��alpha_lower_k, aloha_k��ֵ

			if (alpha_k < alpha_upper_k)
			{
				alpha_upper_k = alpha_k;//���¿ɽ��ܲ����ϱ߽�
			}

			alpha_k_temp = alpha_lower_k + 0.5*(((alpha_k - alpha_lower_k)*(alpha_k - alpha_lower_k)*df_alpha_lower_k) / (f_alpha_lower_k - f_alpha_k + (alpha_k - alpha_lower_k)*df_alpha_lower_k));

			alpha_k = alpha_k_temp;

		}


		if ((alpha_upper_k - alpha_lower_k) < tolerance)//��ֹ�ڹ�С��������
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

double wolf_research::Vsolve(double(*f_test)(Array1D<double>&), Array1D<double>(*g_test)(Array1D<double>&)){
	int k_max = 1000;
	int k = 0;

	double f_current = f_test(Vx_current);
	Array1D<double> Vg_current = g_test(Vx_current);

	double f_alpha_lower_k = f_current;
	Array1D<double> Vg_alpha_lower_k = Vg_current;

	double df_alpha_lower_k = Vd_current * Vg_alpha_lower_k;//d_current.transpose();�����볣���ĵ��

	double f_alpha_lower_0 = f_alpha_lower_k;
	double df_alpha_lower_0 = df_alpha_lower_k;
	double alpha_initial = 0.0;

	double tolerance = 1e-15;
	/************************************************************************/
	/* ����alpha_initial��ֵ����ֹ�����С                                                                     */
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
		/* �������������                                                                     */
		/************************************************************************/
		Vx_alpha_k = Vx_current + Vd_current * alpha_k ;
		f_alpha_k = f_test(Vx_alpha_k);//vector of x_alpha_k
		Vg_alpha_k = g_test(Vx_alpha_k);
		/************************************************************************/
		/* �������                                                                     */
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
				delta_alpha_k = (alpha_k - alpha_lower_k)*df_alpha_k / (df_alpha_lower_k - df_alpha_k);//��alpha_lower_k, aloha_k��ֵ

				if (delta_alpha_k <= 0)
				{
					alpha_k_temp = 2 * alpha_k;//����������
				}
				else{
					alpha_k_temp = alpha_k + delta_alpha_k;
				}

				alpha_lower_k = alpha_k;//���¿ɽ��ܲ������±߽�
				f_alpha_lower_k = f_alpha_k;
				df_alpha_lower_k = df_alpha_k;
				alpha_k = alpha_k_temp;//��������ֵ����λ
			}
		}
		else{
			//condition1 not satisfied, ��alpha_lower_k, aloha_k��ֵ

			if (alpha_k < alpha_upper_k)
			{
				alpha_upper_k = alpha_k;//���¿ɽ��ܲ����ϱ߽�
			}

			alpha_k_temp = alpha_lower_k + 0.5*(((alpha_k - alpha_lower_k)*(alpha_k - alpha_lower_k)*df_alpha_lower_k) / (f_alpha_lower_k - f_alpha_k + (alpha_k - alpha_lower_k)*df_alpha_lower_k));

			alpha_k = alpha_k_temp;

		}


		if ((alpha_upper_k - alpha_lower_k) < tolerance)//��ֹ�ڹ�С��������
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


























#endif