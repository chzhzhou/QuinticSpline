#include "stdafx.h"
//#include "spline.h"
#include "numeric.h"
#include "bem.h"
#include <omp.h>
#include <Eigen/Dense>

Eigen::IOFormat fmt(Eigen::FullPrecision, 0, "\t", "\n", "", "", "", "");

int main() {
	int n = 11;
	int ntest = 1;

	for (int global = 3; global < 21	; global++) {
		n = 8 * (int)pow(1.414, global) + 1;
		//n = (int)pow(2.0,4)+1;
		Eigen::MatrixX2d xy(n, 2);
		xy.setZero();
		for (int i = 0; i < n; i++) {
			double t = ((double)i) / (n - 1.);
			double theta = M_PI * t;
			xy(i, 0) = 50. * sin(theta) *(1 + 0.3 *cos(8 * theta));
			xy(i, 1) = 50. * cos(theta) *(1 + 0.3 *cos(8 * theta));
		}
		Bem bem;
		bem.settings.indexShift(0);
		bem.settings.order(1);
		bem.settings.qdOrder(6);
		//bem.node().r;
		//bem.settings.yBC.end.set(Spline::BC::Mix, 0., );
		bem.initialize(xy);

		

		
		

		//bem.settings.print();

		
		Eigen::MatrixXd S, D, B;
		Bem::Properties::tic();
		bem.assembly(S, D);
		//std::cout << Eigen::VectorXd(D.rowwise().sum().array()+0.5).format(fmt) <<std::endl;
		Bem::Properties::toc();
		
		D.diagonal() -= D.rowwise().sum();

	


		const Eigen::VectorXd &z = bem.node().z.col(0);
		const Eigen::VectorXd &dr = bem.node().r.col(1);
		const Eigen::VectorXd &dz = bem.node().z.col(1);
		Eigen::VectorXd nr = -dz.array() *(dr.array().square() + dz.array().square()).rsqrt();
		Eigen::VectorXd nz = dr.array() *(dr.array().square() + dz.array().square()).rsqrt();


		const Eigen::VectorXd lhs = D * z;
		//std::cout << lhs;
		//std::cout <<( nz);
		Bem::Properties::tic();
		printf("{%04d, %16.16f},\n", n - 1, (S.fullPivLu().solve(lhs) - nz).norm() / sqrt(n - 1.));
		Bem::Properties::toc();


		std::ofstream file("./Output/test.txt");
		file << bem.sp().x().format(fmt) << '\n'
			<< bem.sp().y().format(fmt) << '\n';
		file << bem.sp().h().format(fmt) << '\n';
		file.close();
		std::ofstream file2("./Output/tt.txt");
		file2 << bem.node().r.format(fmt) << '\n'
			<< bem.node().z.format(fmt) << '\n';
		file2.close();
		
		
	}



	//std::ofstream file("./Output/test.txt");
	//file << bem.sp().x().format(fmt) << '\n' 
	//	 << bem.sp().y().format(fmt) << '\n';
	//file << bem.sp().h().format(fmt) << '\n';
	
	
	




	
	return 0;
}

