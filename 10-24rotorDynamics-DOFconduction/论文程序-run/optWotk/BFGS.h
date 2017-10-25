#ifndef BFGS_H
#define BFGS_H

class BFGS
{
public:
	BFGS(double(*f_test_)(Array1D<double>&), Array1D<double>(*g_test_)(Array1D<double>&));

	Array2D<double> VtoMatrix(const Array1D<double> &A, const Array1D<double> &B);

	double Vsolve(Array1D<double>& Vx_initial);
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
};

BFGS::BFGS(double(*f_test_)(Array1D<double>&), Array1D<double>(*g_test_)(Array1D<double>&)){
	Vf_test = f_test_;
	Vg_test = g_test_;


	wolf = wolf_research(f_test_, g_test_);
}


double BFGS::Vsolve(Array1D<double>& Vx_initial){

	int dim = Vx_initial.dim();

	Vx_optimal = Array1D<double>(dim, 0.0);
	f_optimal = 0.0;

	Array1D<double> Vx_current = Vx_initial;
	Array1D<double> Vg_current = Vg_test(Vx_current);

	Array2D<double> B_current(Vx_initial.dim(), Vx_initial.dim(), 0.0);
	Array2D<double> B_previous(Vx_initial.dim(), Vx_initial.dim(), 0.0);
	for (int i = 0; i < Vx_initial.dim(); i++)
	{
		B_current[i][i] = 1.0;
	}

	Array1D<double> Vd_current = -(B_current*Vg_current);

	double tolerance = 1e-4;


	wolf.Vsolve(Vx_current, Vd_current);
	int k = 0;
	//临时中转变量
	Array1D<double> Vx_previous(dim, 0.0);
	Array1D<double> Vg_previous(dim, 0.0);
	Array1D<double> Vd_previous(dim, 0.0);
	Array1D<double> deltaX_current(dim, 0.0);//g_k+1 - g_k
	Array1D<double> deltaG_current(dim, 0.0);//x_k+1 - x_k

	double beta_current;
	while ((wolf.Vx_next - Vx_current).norm() > tolerance)
	{
		k++;
		Vx_previous = Vx_current;
		Vg_previous = Vg_current;
		B_previous = B_current;

		Vx_current = wolf.Vx_next;
		Vg_current = Vg_test(Vx_current);

		deltaX_current = Vx_current - Vx_previous;
		deltaG_current = Vg_current - Vg_previous;

		if (deltaX_current*deltaG_current <= 0)
		{
			for (int i = 0; i < Vx_initial.dim(); i++)
			{
				for (int j = 0; j < Vx_initial.dim(); j++)
				{
					B_current[i][j] = (i == j) ? 1.0 : 0.0;
				}

			}
		}
		else{
			double tempC1 = deltaG_current*deltaX_current;
			double tempC2 = 1 + (deltaG_current*B_previous*deltaG_current) / tempC1;
			Array2D<double> tempM3 = VtoMatrix(deltaX_current, deltaX_current)/tempC1;
			Array2D<double> tempM4 = tempM3*tempC2;
			Array2D<double> tempM5 = (VtoMatrix(deltaX_current, deltaG_current)*B_previous + B_previous* VtoMatrix(deltaG_current, deltaX_current)) / tempC1;
			B_current = B_previous + tempM4 - tempM5;
		}

		Vd_current = -(B_current*Vg_current);

		wolf.Vsolve(Vx_current, Vd_current);
		if (Vg_test(wolf.Vx_next).norm() < tolerance) break;//梯度为0是退出
	}

	Vx_optimal = wolf.Vx_next;
	f_optimal = wolf.f_next;


	return f_optimal;



}



Array2D<double> BFGS::VtoMatrix(const Array1D<double> &A, const Array1D<double> &B)
{
	int m = A.dim1();
	int n = B.dim1();

	Array2D<double> C(m, n, 0.0);

	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			C[i][j] = A[i] * B[j];
		}

	}
	return C;

}





























#endif