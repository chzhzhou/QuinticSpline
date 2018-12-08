#include "stdafx.h"
//#include "spline.h"
#include "numeric.h"
#include "bem.h"
#include <omp.h>
#include <Eigen/Dense>

Eigen::IOFormat fmt(Eigen::FullPrecision, 0, "\t", "\n", "", "", "", "");

int main() {
	int n = 11;
	int ntest = 2;

	for (int global = 5; global < 24; global++) {
		n = (int)pow(sqrt(2.0), global) + 1;
		//n = (int)pow(2.0,4)+1;
		Eigen::MatrixX2d xy(n, 2);
		xy.setZero();
		for (int i = 0; i < n; i++) {
			double t = ((double)i) / (n - 1.);
			double theta = M_PI * t;
			xy(i, 0) = 1. * sin(theta) * (1+ 0.0 * cos(4.* theta)) ;
			xy(i, 1) = 1. * cos(theta) * (1 + 0.0 *cos(4.*theta));
		}
		Bem bem;
		bem.settings.indexShift(0);
		bem.settings.order(1);
		bem.settings.qdOrder(10);
		bem.initialize(xy);
		std::ofstream file("./Output/test.txt");
		file << bem.sp().x().format(fmt) << '\n'
			<< bem.sp().y().format(fmt) << '\n';
		file << bem.sp().h().format(fmt) << '\n';
		file.close();
		std::ofstream file2("./Output/tt.txt");
		file2 << bem.node().r.format(fmt) << '\n'
			<< bem.node().z.format(fmt) << '\n';
		file2.close();
		//bem.node().r;
		//bem.settings.yBC.end.set(Spline::BC::Mix, 0., );
		//bem.settings.print();		
		Eigen::MatrixXd S, D, B;	
		bem.assembly(S, D);
		//std::cout << (D).row(D.rows() - 2).sum()+0.5 << std::endl;

		
		
		B.resize(D.rows(), D.cols());
		B.setIdentity();
		B *= 0.5;
		


	

		
		const Eigen::VectorXd &z = bem.node().z.col(0);
		const Eigen::VectorXd &dr = bem.node().r.col(1);
		const Eigen::VectorXd &dz = bem.node().z.col(1);
		Eigen::VectorXd nr = -dz.array() *(dr.array().square() + dz.array().square()).rsqrt();
		Eigen::VectorXd nz = dr.array() *(dr.array().square() + dz.array().square()).rsqrt();


		/*const Eigen::VectorXd lhs = (B+D) * z;
		printf("{%04d, %16.16f},\n", n - 1, (S.fullPivLu().solve(lhs) +		
			 - nz
			).cwiseAbs().maxCoeff());*/
		
		const Eigen::VectorXd lhs = (S) * (z);
		//std::cout << ((B - D).fullPivLu().solve(lhs)) ;
		printf("{%04d, %16.16f},\n", n - 1, ((B - D).fullPivLu().solve(lhs) -0.5 * z).cwiseAbs().maxCoeff());

		
		
	}



	//std::ofstream file("./Output/test.txt");
	//file << bem.sp().x().format(fmt) << '\n' 
	//	 << bem.sp().y().format(fmt) << '\n';
	//file << bem.sp().h().format(fmt) << '\n';
	
	
	




	
	return 0;
}

