#include "stdafx.h"
#include "spline.h"
#include "numeric.h"
#define M_PI   3.14159265358979323846264338327950288

int main()
{
	int n = 21;
	Eigen::MatrixX2d xy(n,2);
	xy.setZero();
	for (int i = 0; i < n; i++) {
		double t = (double)i / (n - 1.);
		double theta = M_PI *sqrt(t);
		xy(i, 0) = sin(theta) *(1 + 0.1 *cos(6*theta));
		xy(i, 1) = cos(theta) *(1 + 0.1 *cos(6*theta));
	}

	Spline sp0;
	sp0.node(xy);	
	sp0.x(Spline::BC::Odd, Spline::BC::Odd);
	sp0.y(Spline::BC::Even, Spline::BC::Even);
	Spline sp(sp0);
	
	//sp = sp;
	//sp.localArc(5,1.0);
	for (int i = 0; i < n - 1; i++) {
		double tmp = sp.arc2t(i, sp.localArc(i)*0.5);
		std::cout << sp.localArc(i)*0.5 - sp.localArc(i, tmp) << "\n";		
	}
	//sp.x(Spline::BC::Mix, Spline::BC::Mix,0,0,0,0);
	//sp.y(Spline::BC::Mix, Spline::BC::Mix, 0, 0, 0, 0);
	//std::cout << sp.x() <<std::endl;
	//std::cout << sp.y() << std::endl;
	Eigen::IOFormat fmt(Eigen::FullPrecision, 0, "\t", "\n", "", "", "", "");
	std::ofstream file("./Output/test.txt");
	file << sp.x().format(fmt) << '\n' << sp.y().format(fmt) << '\n';	
	file << sp.h().format(fmt) << '\n';
	
	//sp.y(Spline::BC::Even, Spline::BC::Even);
	
	double out = 0;
	int o = 4;
	for (int i = 0; i < o * 5; i++) {		
		out -= Numeric::lqd[o][2 * i + 1 ] * cos(Numeric::lqd[o][2 * i]);
	}
	printf("\n\n\n\n%16.16f = \n", out);
	double K, E;
	Numeric::KEPQ(0.99, K, E);
	printf("\n\n\n\n%16.16f\t%16.16f = \n", K,E);

	return 0;
}

