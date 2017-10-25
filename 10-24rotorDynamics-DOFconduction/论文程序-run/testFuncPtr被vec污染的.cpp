
#include <math.h>
#include <iostream>
#include "wolf_research.h"

using namespace std;
using namespace TNT;
int testFunc(int& a, int& b){//�������������ĺ���
	return a + b;
}

void printfunc(int(*objFunc)(int&, int&), int x, int y){//ֻ����Ϊһ��������ת����ָ�����
	cout << objFunc(x, y);
}

//**************************************************************************x_current d_current shell be a vector
double Wolf_Search(double(*f_test)(Array1D<double>&), Array1D<double>(*g_test)(Array1D<double>&), Array1D<double> Vx_current, Array1D<double> Vd_current, double rho, double sigma){
	int k_max = 1000;
	int k = 0;
	int dim = Vx_current.dim();

	double f_current = f_test(Vx_current);
	Array1D<double> Vg_current = g_test(Vx_current);

	double f_alpha_lower_k = f_current;
	Array1D<double> Vg_alpha_lower_k = Vg_current;

	double df_alpha_lower_k = Vd_current * Vg_alpha_lower_k;//d_current.transpose();�����볣���ĵ��

	double f_alpha_lower_0 = f_alpha_lower_k;
	double df_alpha_lower_0 = df_alpha_lower_k;
	double alpha_initial;

	double tolerance = 1e-15;
	/************************************************************************/
	/* ����alpha_initial��ֵ����ֹ�����С                                                                     */
	/************************************************************************/
	if (abs(df_alpha_lower_k)>tolerance)
	{
		alpha_initial = -2 * f_alpha_lower_k / df_alpha_lower_k;

	}
	else{
		alpha_initial = 1;
	}

	if (alpha_initial<tolerance)
	{
		alpha_initial = 1;
	}
	
	double alpha_lower_k = 0.0;
	double alpha_upper_k = 1e8;
	double alpha_k = alpha_initial;
	double alpha_acceptable;//return alpha
	Array1D<double> Vx_next(dim,0.0);
	double f_next;
	double delta_alpha_k, alpha_k_temp;
	/************************************************************************/
	/* begin search                                                                     */
	/************************************************************************/
	double wolfe_condition1;
	double wolfe_condition2;
	Array1D<double> Vx_alpha_k(dim, 0.0);
	double f_alpha_k;
	Array1D<double> Vg_alpha_k(dim, 0.0);
	double df_alpha_k;
	
	for (k = 1; k < k_max; k++)
	{
		/************************************************************************/
		/* �������������                                                                     */
		/************************************************************************/
		Vx_alpha_k = Vx_current +  Vd_current*alpha_k;
		f_alpha_k = f_test(Vx_alpha_k);//vector of x_alpha_k
		Vg_alpha_k = g_test(Vx_alpha_k);
		/************************************************************************/
		/* �������                                                                     */
		/************************************************************************/
		df_alpha_k = Vd_current * Vg_alpha_k;

		//wolf condition
		wolfe_condition1 = f_alpha_k - f_alpha_lower_0 - rho*alpha_k*df_alpha_lower_0;
		wolfe_condition2 = sigma*df_alpha_lower_0 - df_alpha_k;

		if (wolfe_condition1<=0)
		{
			if (wolfe_condition2<=0)
			{
				alpha_acceptable = alpha_initial;
				Vx_next = Vx_alpha_k;
				f_next = f_alpha_k;
				break;
			}
			else{//only sastify condition1
				delta_alpha_k = (alpha_k - alpha_lower_k)*df_alpha_k / (df_alpha_lower_k - df_alpha_k);//��alpha_lower_k, aloha_k��ֵ

				if (delta_alpha_k<=0)
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

			if (alpha_k<alpha_upper_k)
			{
				alpha_upper_k = alpha_k;//���¿ɽ��ܲ����ϱ߽�
			}

			alpha_k_temp = alpha_lower_k + 0.5*(((alpha_k - alpha_lower_k)*(alpha_k - alpha_lower_k)*df_alpha_lower_k) / (f_alpha_lower_k - f_alpha_k + (alpha_k - alpha_lower_k)*df_alpha_lower_k));

			alpha_k = alpha_k_temp;

		}


		if ((alpha_upper_k-alpha_lower_k)<tolerance)//��ֹ�ڹ�С��������
		{
			alpha_acceptable = alpha_k;
			Vx_next = Vx_alpha_k;
			f_next = f_alpha_k;
		}

	}//k<k_max

	if (wolfe_condition1>0||wolfe_condition2>0)
	{
		cout << "failed \n";
		
	}

	//return alpha_acceptable;
	return alpha_acceptable;


}

double f_test(double& x){
	return -3 * x*sin(0.75*x) + exp(-2 * x);
}
double g_test(double& x){
	return -2 / exp(2 * x) - 3 * sin((3 * x) / 4) - (9 * x*cos((3 * x) / 4)) / 4;
}


int main(){

	int(*p)(int &, int &);
	p = testFunc;

	printfunc(testFunc, 1, 1);

	printfunc(&testFunc, 1, 1);//����ֱ��ʹ�ô����õĺ�����Ҳ������ǰ��������ת

	double x_current = -2, d_current = 1, rho = 0.1, sigma = 0.11;

	double alpha;

//	alpha = Wolf_Search(f_test, g_test, x_current, d_current, rho, sigma);

//	cout << "***************alpha = " << alpha<<endl;

	wolf_research wolf(x_current, d_current, rho, sigma);
	wolf.solve(f_test, g_test);

	cout << wolf.alpha_acceptable;


	system("pause");
	return 0;
}