#include "stdafx.h"
//#include "spline.h"
#include "numeric.h"
#include "bem.h"
#include <omp.h>

double timer;
void tic() { printf("tic\n"); timer = omp_get_wtime(); };
void toc() { printf("toc ... %.ef sec\n" ,omp_get_wtime()-timer); };

int main() {
	int n = 11;
	Eigen::MatrixX2d xy(n,2);
	xy.setZero();
	for (int i = 0; i < n; i++) {
		double t = (double)i / (n - 1.);
		double theta = M_PI *sqrt(t);
		xy(i, 0) = sin(theta) *(1 + 0.1 *cos(6*theta));
		xy(i, 1) = cos(theta) *(1 + 0.1 *cos(6*theta));
	}

	Bem bem;
	bem.settings.indexShift(0);
	//bem.settings.yBC.end.set(Spline::BC::Mix, 0., );
	
	tic();
	bem.initialize(xy);
	toc();
	
	
	bem.settings.print();
	printf("uint--:  %16.16f\n  ", Numeric::N[2][2](0.0));

	//Spline sp0;
	//sp0.node(xy);	
	//sp0.x(Spline::BC::Odd, Spline::BC::Odd);
	//sp0.y(Spline::BC::Even, Spline::BC::Even);
	Spline sp(bem.sp());
	
	
	

	
	Eigen::IOFormat fmt(Eigen::FullPrecision, 0, "\t", "\n", "", "", "", "");
	std::ofstream file("./Output/test.txt");
	file << sp.x().format(fmt) << '\n' << sp.y().format(fmt) << '\n';	
	file << sp.h().format(fmt) << '\n';	
	
	
	double out = 0;
	int o = 4;
	for (int i = 0; i < o * 5; i++) {		
		out -= Numeric::lqd[o][2 * i + 1 ] * cos(Numeric::lqd[o][2 * i]);
	}
	printf("%16.16f = \n", out);
	double K, E;
	Numeric::KEPQ(0.99, K, E);
	printf("K = %16.16f\nE = %16.16f\n", K,E);




	
	return 0;
}

