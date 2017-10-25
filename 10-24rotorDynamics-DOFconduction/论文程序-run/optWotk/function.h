double relative_density(Array1D<double>& X){
	double x1 = X[0];
	double x2 = X[1];
	double x3 = X[2];
	double d = X[3];


	return M_PI*d*d*sqrt(x1*x1 + x2*x2 + x3*x3) / x1*x2*x3;

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

double modulus(Array1D<double>& X){
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

	return E;
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

double plastic_strength(Array1D<double>& X){
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

	return 4.0*x1*x2*temp_length*A0 / (M_PI*d*d*d);

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