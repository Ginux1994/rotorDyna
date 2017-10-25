#ifndef GUGMENTEDLAGRANGE_H
#define GUGMENTEDLAGRANGE_H

#include "Gradient.h"

class Gradient;

class AugmentedLagrange
{
public:
	AugmentedLagrange(double(*f_test_)(Array1D<double>&, Array1D<double>&), Array1D<double>(*g_test_)(Array1D<double>&, Array1D<double>&), Array1D<double>(*niuUpdate_)(Array1D<double>&, Array1D<double>&));

	double Vsolve(Array1D<double>& Vx_current, Array1D<double>& Vniu_initial);
	double objFunction(Array1D<double>& X){


		double a = X[0];
		double b = X[1];
		double ans = (a - 2)*(a - 2)*(a - 2)*(a - 2) + (a - 2 * b)*(a - 2 * b);//Ŀ�꺯��

		return ans;

	}

protected:
public:
	/************************************************************************/
	/* WOLFһά����                                                                     */
	/************************************************************************/
	Gradient gradient;
	/************************************************************************/
	/* ��������                                                                     */
	/************************************************************************/

	bool stopOrNot(Array1D<double>& X){
		double x1 = X[0];
		double x2 = X[1];
		double x3 = X[2];
		double d = X[3];//��Ʊ���
		double tolerance = 1e-5;

		double restraintF1 = d - 0.3;//��һ������ʽԼ���������Դ���
		//Լ�������� X ���ݶ�

		double restraintF2 =1.8 - d;//��һ������ʽԼ���������Դ���

		double restraintF3 = x1 - d - 0.3;//��һ������ʽԼ���������Դ���

		double restraintF4 = x2 - d - 0.3;//��һ������ʽԼ���������Դ���


		double restraintF5 = x3 - 0.3;//��һ������ʽԼ���������Դ���

		double restraintF6 = asin(x3 / sqrt(x1*x1 + x2*x2 + x3*x3)) - M_PI / 9.0;//��һ������ʽԼ���������Դ���

		double restraintF7 = M_PI / 4.0 - asin(x3 / sqrt(x1*x1 + x2*x2 + x3*x3));//��һ������ʽԼ���������Դ���

		if (sqrt(restraintF1*restraintF1 + restraintF2*restraintF2 + restraintF3*restraintF3 + restraintF4*restraintF4 + restraintF5*restraintF5 + restraintF6*restraintF6 + restraintF7*restraintF7) < tolerance)
		{
			return false;
		}
		else return true;
	}

	/************************************************************************/
	/* WOLFһά����                                                                     */
	/************************************************************************/
	wolf_research wolf;
	/************************************************************************/
	/* ��������                                                                     */
	/************************************************************************/
	double(*Vf_test)(Array1D<double>&, Array1D<double>&);
	Array1D<double>(*Vg_test)(Array1D<double>&, Array1D<double>&);
	Array1D<double>(*Vniu_update)(Array1D<double>&, Array1D<double>&);

	//�������
	Array1D<double> Vx_optimal;
	double x_optimal;
	double f_optimal;


};

AugmentedLagrange::AugmentedLagrange(double(*f_test_)(Array1D<double>&, Array1D<double>&), Array1D<double>(*g_test_)(Array1D<double>&, Array1D<double>&), Array1D<double>(*niuUpdate_)(Array1D<double>&, Array1D<double>&)){
	Vf_test = f_test_;
	Vg_test = g_test_;
	Vniu_update = niuUpdate_;

	gradient = Gradient(Vf_test, Vg_test);
	
}





double AugmentedLagrange::Vsolve(Array1D<double>& Vx_initial, Array1D<double>& Vniu_initial){

	
	f_optimal = gradient.LangrangeSolve(Vx_initial, Vniu_initial);
	Array1D<double> Vx_current = gradient.Vx_optimal;
	Array1D<double> Vniu_current = Vniu_initial;
	Array1D<double> Vniu_previous = Vniu_current;
	int k = 0;
	double f_optimal_last;
	while (stopOrNot(gradient.Vx_optimal))//��ֹ���������ݵ�ǰx�Ƿ���
	{
		k++;
		Vniu_previous = Vniu_current;
		//���������½������������������պ��������ŵ�
		f_optimal_last = f_optimal;
		f_optimal = gradient.LangrangeSolve(Vx_current, Vniu_current);
		Vx_current = gradient.Vx_optimal;

		//������һ���������x ����niu
		Vniu_current = Vniu_update(gradient.Vx_optimal, Vniu_initial);
		


		double temp = abs(f_optimal - f_optimal_last);
		if (temp<1e-1)
		{
			break;
		}
		cout << "=======�� " << k << " ������" << endl;
		cout << "vx_current = " << Vx_current;
		cout << "��ʱ����ܶ�Ϊ " << f_optimal << endl;
		cout << temp;
	}
	
	f_optimal = objFunction(Vx_current);//����Ŀ�꺯��������ֵ�϶���ͨ��Ŀ�꺯��������ã������ݶ��㷨���ں��˱߽�Լ�������ĸ��Ϻ���f_next���������Ҫ����������X���ǶԵ�
	cout << "********����XΪ L1��L2��L3��D \n";
	cout << Vx_current << endl;
	return f_optimal;


}
























#endif
