#include "stdafx.h"
//#include "spline.h"
#include "numeric.h"
#include "bem.h"
#include <omp.h>
#include <Eigen/Dense>
//#define COMPUTE
//#define TEST1

Eigen::IOFormat fmt(Eigen::FullPrecision, 0, "\t", "\n", "", "", "", "");
Eigen::MatrixX2d circle(double angle, int n) {
	Eigen::MatrixX2d xy(n, 2);
	xy.setZero();
	for (int i = 0; i < n; i++) {
		double t = ((double)i) / (n - 1.);
		double theta = angle * t;// pow(t, 0.8);
		xy(i, 0) = 2. * sin(theta) * (1 + 0.25 * cos(8. * theta - M_PI));
		xy(i, 1) = 2. * cos(theta) * (1 + 0.25 * cos(8. * theta - M_PI));
	}
	return xy;
}
Eigen::MatrixX2d line(double x0, double y0, double x1, double y1, int n) {
	Eigen::MatrixX2d xy(n, 2);
	xy.setZero();
	for (int i = 0; i < n; i++) {
		double t = ((double)i) / (n - 1.);
		t = pow(t, 1.0);
		double xt = x0 + (x1 - x0) * t;
		double yt = y0 + (y1 - y0) * t;
		xy(i, 0) = xt;
		xy(i, 1) = yt;
	}
	return xy;

}

double phi(double r, double z) {
	double xx = z * z /(r*r + z * z);
	return  0.5 * ( 3 * z * z  - (r*r + z * z) *1.0);
}
double phin(double nr, double nz, double r, double z) {
	//double x = z / sqrt(r*r + z * z);
	return nr * (-r) + nz * (2. * z);
}

//nr.array() * (-r.array())  + nz.array() * (2. * z.array()))

void printSpline(const Spline &sp, const std::string &name) {
	std::ofstream file(name);	
	file << sp.x().format(fmt) << '\n'	<< sp.y().format(fmt) << '\n' << sp.h().format(fmt) << '\n';
	file.close();
}

int main() {
	int n = 11;
	int ntest = 2;

	for (int global = 7; global < 20; global++) {
		n = (int)pow(sqrt(2.0), global) + 1;
		//n = (int)pow(2.0, 9) + 1;
		Eigen::MatrixX2d xy0 = circle(M_PI/2, n);
		Bem bem0;
		bem0.settings.indexShift(0);
		bem0.settings.order(2);
		bem0.settings.qdOrder(20);
		bem0.settings.xBC.end.set(Spline::BC::Even,0,0);		
		bem0.settings.yBC.end.set(Spline::BC::Odd,0,0);
		bem0.initialize(xy0);
		printSpline(bem0.sp(),"./Output/sp0.txt");
		int n0 = bem0.node().r.rows();
		//bem0.settings.print();

		Eigen::MatrixX2d xy1 = line(xy0(xy0.rows() - 1, 0 ), xy0(xy0.rows() - 1, 1), 0.0, 0.0,  n );
		Bem bem1;
		bem1.settings.indexShift(bem0.node().r.rows() );
		bem1.settings.order(2);
		bem1.settings.qdOrder(20);
		bem1.settings.xBC.end.set(Spline::BC::Odd, 0, 0);		
		bem1.settings.yBC.end.set(Spline::BC::Even, 0, 0);
		bem1.initialize(xy1);
		printSpline(bem1.sp(), "./Output/sp1.txt");
		int n1 = bem1.node().r.rows();

		Eigen::MatrixXd S, D, B, L,R;
		int nTotal = n0  + n1;
//		std::cout << nTotal;
		S.setZero(nTotal, nTotal);
		B.setZero(nTotal, nTotal);
		D.setZero(nTotal, nTotal);
		R.setZero(nTotal, nTotal);
		L.setZero(nTotal, nTotal);
		
		B.setIdentity();
		B *= 0.5;
		B(n0-1 , n0-1) = 0.25;
		B(n1, n0) = 0.25;


		std::cout.precision(15);
		Bem::assembly(bem0, bem0, S, D);						
		//std::cout << "\n--------------\n";
		Bem::assembly(bem0, bem1, S, D);				
		//std::cout << "\n--------------\n";
		Bem::assembly(bem1, bem0, S, D);
		//std::cout << "\n--------------\n";
	    Bem::assembly(bem1, bem1, S, D);
		//std::cout << "\n--------------\n";
		
		D = D + B;		
		//std::cout << (D).rowwise().sum() << std::endl;
		std::cout.precision(3);
		//std::cout << D.row(11)- D.row(12) <<std::endl;
		//std::cout <<  <<std::endl;
		Eigen::VectorXd rhs, lhs;
		rhs.setZero(nTotal);
		lhs.setZero(nTotal);
		for (int k = 0; k < nTotal; k++) {
			if (k < n0) {
				double r = bem0.node().r(k, 0), dr = bem0.node().r(k, 1);
				double z = bem0.node().z(k, 0), dz = bem0.node().z(k, 1);
				double nr = -dz / sqrt(dr * dr + dz * dz);
				double nz = dr / sqrt(dr * dr + dz * dz);				
				rhs(k) = phin(nr,nz,r,z);
				lhs(k) = phi(r, z);
			}
			else {
				int kk = k - n0;
				double r = bem1.node().r(kk, 0), dr = bem1.node().r(kk, 1);
				double z = bem1.node().z(kk, 0), dz = bem1.node().z(kk, 1);
				double nr = -dz / sqrt(dr * dr + dz * dz);
				double nz = dr / sqrt(dr * dr + dz * dz);
				rhs(k) = phi(r, z);
				lhs(k) = phin(nr, nz, r, z);

			}
		}
		//std::cout << lhs;
		D.row(n0 - 1) *= 0.;		
		D(n0 -1 , n0 - 1) = 1.0;
		D(n0 -1,  n0 ) = -1.0;
		S.row(n0 - 1) *= 0.;
		

		L.topLeftCorner(n0, n0) = D.topLeftCorner(n0, n0);
		L.topRightCorner(n0, n1) = -S.topRightCorner(n0, n1);
		L.bottomLeftCorner(n1, n0) = D.bottomLeftCorner(n1, n0);
		L.bottomRightCorner(n1, n1) = -S.bottomRightCorner(n1, n1);

		R.topLeftCorner(n0, n0) = S.topLeftCorner(n0, n0);
		R.topRightCorner(n0, n1) = -D.topRightCorner(n0, n1);
		R.bottomLeftCorner(n1, n0) = S.bottomLeftCorner(n1, n0);
		R.bottomRightCorner(n1, n1) = -D.bottomRightCorner(n1, n1);

		
		//std::cout << L.row(n0 - 1);

		//std::cout << "\n--------------\n";
		//std::cout << R.row(n0 - 1);

		//std::cout << "\n--------------\n";

		Eigen::VectorXd answer =  R.fullPivLu().solve(L*lhs);
		//std::cout << answer(n0 - 1) << "\t ";
		//std::cout << answer(n0  ) << "\t ";
		//std::cout << answer(n0 + 1) << "\t \n";
		//std::cout << rhs(n0 - 1) << "\t ";
		//std::cout << rhs(n0) << "\t ";
		//std::cout << rhs(n0 + 1) << "\t \n";
		Eigen::VectorXd errorDirchlet(n0);
		for (int kk = 0; kk < n0; kk++) { 
			errorDirchlet(kk) = answer(kk) - rhs(kk);
		}
		Eigen::VectorXd errorNeumann(n1);
		for (int kk = 0; kk < n1; kk++) {
			errorNeumann(kk) = answer(n0 + kk) - rhs(n0 + kk);
		}

		printf("{%04d, %16.16f},\n", 2*n, (errorDirchlet	).cwiseAbs().maxCoeff());

#ifdef TEST1

		Eigen::MatrixX2d xy0 = circle(M_PI, n);
		Bem bem0;
		bem0.settings.indexShift(0);
		bem0.settings.order(2);
		bem0.settings.qdOrder(20);		
		bem0.initialize(xy0);
		printSpline(bem0.sp(), "./Output/sp0.txt");

		Eigen::MatrixXd S, D, B;
		int nTotal = bem0.node().r.rows();
		std::cout << nTotal;
		S.setZero(nTotal, nTotal);
		B.setZero(nTotal, nTotal);
		D.setZero(nTotal, nTotal);
		Bem::assembly(bem0, bem0, S, D);

		for (int i = 0; i < D.rows(); i++) {
			D(i, i) = -(D.row(i).sum() - D(i, i));
		}
			   
		B.resize(D.rows(), D.cols());
		B.setIdentity();
		B *= 0.5;

		const Eigen::VectorXd &z = bem0.node().z.col(0);
		const Eigen::VectorXd &r = bem0.node().r.col(0);
		const Eigen::VectorXd &dr = bem0.node().r.col(1);
		const Eigen::VectorXd &dz = bem0.node().z.col(1);
		const Eigen::VectorXd &cosTh = z.array() *(r.array().square() + z.array().square()).rsqrt();
		const Eigen::VectorXd &R2 = (r.array().square() + z.array().square());

		Eigen::VectorXd nr = -dz.array() *(dr.array().square() + dz.array().square()).rsqrt();
		Eigen::VectorXd nz = dr.array() *(dr.array().square() + dz.array().square()).rsqrt();
		Eigen::VectorXd phi = R2.array() * (0.5 * (3. * cosTh.array().square() - 1.0));

		const Eigen::VectorXd lhs = D * phi;
		
		printf("{%04d, %16.16f},\n", 2 * (n - 1) + 1, (S.fullPivLu().solve(lhs) -
			Eigen::VectorXd(nr.array() * (-r.array()) + nz.array() * (2. * z.array()))
			).cwiseAbs().maxCoeff());

#endif // TEST1





	

#ifdef COMPUTE



		Eigen::MatrixXd S, D, B;	
		Bem::assembly(bem0, bem0,S, D);
		//std::cout << (D).row(D.rows() - 2).sum()+0.5 << std::endl;
		//std::cout << (D).rowwise().sum().array() + 0.5  << std::endl;

		for (int i = 0; i < D.rows(); i++){

			D(i, i) = -(D.row(i).sum() - D(i, i));
		}

		
		
		B.resize(D.rows(), D.cols());
		B.setIdentity();
		B *= 0.5;
			
		const Eigen::VectorXd &z = bem.node().z.col(0);
		const Eigen::VectorXd &r = bem.node().r.col(0);
		const Eigen::VectorXd &dr = bem.node().r.col(1);
		const Eigen::VectorXd &dz = bem.node().z.col(1);
		const Eigen::VectorXd &cosTh = z.array() *(r.array().square() + z.array().square()).rsqrt();
		const Eigen::VectorXd &R2 = (r.array().square() + z.array().square());

		
		Eigen::VectorXd nr = -dz.array() *(dr.array().square() + dz.array().square()).rsqrt();
		Eigen::VectorXd nz = dr.array() *(dr.array().square() + dz.array().square()).rsqrt();
		Eigen::VectorXd phi = R2.array() * (0.5 * (3. * cosTh.array().square() - 1.0));

		const Eigen::VectorXd lhs = D * phi;

		//std::cout << lhs;
		printf("{%04d, %16.16f},\n", 2 * (n -1) + 1, (S.fullPivLu().solve(lhs) -		
			  Eigen::VectorXd (nr.array() * (-r.array())  + nz.array() * (2. * z.array()))
			).cwiseAbs().maxCoeff());
		
		//const Eigen::VectorXd lhs = (S) * (z);
		//std::cout << ((B - D).fullPivLu().solve(lhs)) ;
		//printf("{%04d, %16.16f},\n", n - 1, ((B - D).fullPivLu().solve(lhs) -0.5 * z).cwiseAbs().maxCoeff());

#endif // COMPUTE
		
		
	}



	//std::ofstream file("./Output/test.txt");
	//file << bem.sp().x().format(fmt) << '\n' 
	//	 << bem.sp().y().format(fmt) << '\n';
	//file << bem.sp().h().format(fmt) << '\n';
	
	
	




	
	return 0;
}

