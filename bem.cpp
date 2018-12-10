#include "stdafx.h"
#include "bem.h"
#include "numeric.h"
#include <omp.h>
#ifndef M_PI
#define M_PI   3.14159265358979323846264338327950288
#endif // !M_PI
void Bem::initialize(const Eigen::MatrixX2d &xy) {
	_sp.node(xy);
	const Properties::bc &xbc = settings.xBC;
	_sp.x(xbc.begin.type, xbc.end.type, xbc.begin.a, xbc.begin.b, xbc.end.a, xbc.end.b);
	const Properties::bc &ybc = settings.yBC;
	_sp.y(ybc.begin.type, ybc.end.type, ybc.begin.a, ybc.begin.b, ybc.end.a, ybc.end.b);
	settings.nElm(_sp.node().rows() - 1);
	_e.resize(settings.nElm());
	for (std::size_t i = 0; i < _e.size(); i++) {
		_e[i].init(_sp, i, settings.order(), settings.qdOrder());
	}
	setrz();

	

	
}

void Bem::setrz() {
	int nNode = settings.nElm() * settings.order() + 1;
	_node.r.setZero(nNode, 3);
	_node.z.setZero(nNode, 3);
	for (int idElement = 0; idElement < settings.nElm(); idElement++) {
		const std::vector<double> &t = e()[idElement].t();
		for (int k = 0; k < settings.order(); k++) {
			int idNode = idElement * settings.order() + k;
			_node.r.row(idNode) = sp().d(sp().x(), idElement, t[k]);
			_node.z.row(idNode) = sp().d(sp().y(), idElement, t[k]);
		}
	}
	int idNode = _node.r.rows() - 1;
	int idElement = settings.nElm() - 1;
	const std::vector<double> &t = e()[idElement].t();
	_node.r.row(idNode) = sp().d(sp().x(), idElement + 1, 0.);
	_node.z.row(idNode) = sp().d(sp().y(), idElement + 1, 0.);
	//printf("ff%15.15f\n", _node.r.row(idNode)(0));
};

int Bem::checkBounbdaryPosition(const Bem &bem0, const Bem &bem1) {
	int s0 = bem0.settings.indexShift();
	int s1 = bem1.settings.indexShift();
	int ds = s1-s0;

	if (ds == 0) { return 0; }
	if (ds > 0 && ds == bem0.node().r.rows()) {
		return 1;
	}
	if (ds < 0 && abs(ds) == bem1.node().r.rows()) {
		return -1;
	}
	return 2;

};

void Bem::assembly(const Bem &bem0, const Bem &bem1, Eigen::MatrixXd &S, Eigen::MatrixXd &D) {

	const int nNode = bem0.node().r.rows();
	const int shift0 = bem0.settings.indexShift();
	const Eigen::VectorXd &r = bem0.node().r.col(0);
	const Eigen::VectorXd &z = bem0.node().z.col(0);

	const int shift1 = bem1.settings.indexShift();	
	const int o = bem1.settings.order();
	const int nElem = bem1.settings.nElm();
	const int check = checkBounbdaryPosition(bem0, bem1);
	
	//std::cout << "check " << checkBounbdaryPosition(bem0, bem1) << "\n";

#pragma omp parallel for 
	for (int i = 0; i < nNode; i++) {		//std::cout << i << "\n";
		if (abs(r(i)) < 1e-13) {			//printf("Axis %03d %16.16f: \n", i, abs(r(i)));
			for (int j = 0; j < nElem; j++) {
				const std::vector<double > axisIntegral = bem1.axis(z(i), j);
				for (int k = 0; k <= o; k++) {
					S(i + shift0, shift1 + o * j + k) += axisIntegral[k];
					D(i + shift0, shift1 + o * j + k) += axisIntegral[o + 1 + k];
					
				}
			}
		}
		else {
			switch (check) {
			case 0: {
				for (int j = 0; j < i / o - 1 + (i % o); j++) {
					const std::vector<double > regularIntegral = bem1.regular(r(i), z(i), j);
					for (int k = 0; k <= o; k++) {
						S(i + shift0, o * j + k + shift1) += regularIntegral[k];
						D(i + shift0, o * j + k + shift1) += regularIntegral[o + 1 + k];
					}
				}
				if (i % o == 0) {
					int j = i / o - 1;

					if (j >= 0 && j < nElem) {
						std::vector<double > singularIntegral = bem1.singular(1, j);
						for (int k = 0; k <= o; k++) {
							S(i + shift0, o * j + k + shift1) += singularIntegral[k];
							D(i + shift0, o * j + k + shift1) += singularIntegral[o + 1 + k];
						}
					}

					j = i / o;
					if (j >= 0 && j < nElem) {
						std::vector<double > singularIntegral = bem1.singular(0, j);
						for (int k = 0; k <= o; k++) {
							S(i + shift0, o * j + k + shift1) += singularIntegral[k];
							D(i + shift0, o * j + k + shift1) += singularIntegral[o + 1 + k];
						}
					}
				}
				else {
					int j = i / o;
					std::vector<double > singularIntegral = bem1.singular(bem1.e()[j].t()[1], j);
					for (int k = 0; k <= o; k++) {
						S(i + shift0, o * j + k + shift1) += singularIntegral[k];
						D(i + shift0, o * j + k + shift1) += singularIntegral[o + 1 + k];
					}
				}
				for (int j = i / o + 1; j < nElem; j++) {
					const std::vector<double > regularIntegral = bem1.regular(r(i), z(i), j);
					for (int k = 0; k <= o; k++) {
						S(i + shift0, o * j + k + shift1) += regularIntegral[k];
						D(i + shift0, o * j + k + shift1) += regularIntegral[o + 1 + k];
					}
				}
			
				break;
			}

			case -1: {
				for (int j = 0; j < nElem; j++) {
					if (i == 0 && j == nElem - 1) {
						std::vector<double > singularIntegral = bem1.singular(1.0, j);
						for (int k = 0; k <= o; k++) {
							S(i + shift0, o * j + k + shift1) += singularIntegral[k];
							D(i + shift0, o * j + k + shift1) += singularIntegral[o + 1 + k];
						}

					}
					else {

						const std::vector<double > regularIntegral = bem1.regular(r(i), z(i), j);
						for (int k = 0; k <= o; k++) {
							S(i + shift0, o * j + k + shift1) += regularIntegral[k];
							D(i + shift0, o * j + k + shift1) += regularIntegral[o + 1 + k];
						}
					}
				}
				break;
			}

			case 1: {
				for (int j = 0; j < nElem; j++) {
					if (i == nNode - 1 && j==0) {
						std::vector<double > singularIntegral = bem1.singular(0, j);
						for (int k = 0; k <= o; k++) {
							S(i + shift0, o * j + k + shift1) += singularIntegral[k];
							D(i + shift0, o * j + k + shift1) += singularIntegral[o + 1 + k];
						}

					}
					else {

						const std::vector<double > regularIntegral = bem1.regular(r(i), z(i), j);
						for (int k = 0; k <= o; k++) {
							S(i + shift0, o * j + k + shift1) += regularIntegral[k];
							D(i + shift0, o * j + k + shift1) += regularIntegral[o + 1 + k];
						}
					}
				}
				break;
			}

			default: {
				for (int j = 0; j < nElem; j++) {
					const std::vector<double > regularIntegral = bem1.regular(r(i), z(i), j);
					for (int k = 0; k <= o; k++) {
						S(i + shift0, o * j + k + shift1) += regularIntegral[k];
						D(i + shift0, o * j + k + shift1) += regularIntegral[o + 1 + k];
					}
				}
				break;
			}
			
			}




















			/*
			for (int j = 0; j < i / o - 1 + (i % o); j++) {
				const std::vector<double > regularIntegral = bem1.regular(r(i), z(i), j);
				for (int k = 0; k <= o; k++) {
					S(i + shift0, o * j + k + shift1) += regularIntegral[k];
					D(i + shift0, o * j + k + shift1) += regularIntegral[o + 1 + k];
				}
			}
			if (i % o == 0) {
				int j = i / o - 1;

				if (j >= 0 && j < nElem) {
					std::vector<double > singularIntegral = bem1.singular(1, j);
					for (int k = 0; k <= o; k++) {
						S(i + shift0, o * j + k + shift1) += singularIntegral[k];
						D(i + shift0, o * j + k + shift1) += singularIntegral[o + 1 + k];
					}
				}
				
				j = i / o;
				if (j >= 0 && j < nElem) {
					std::vector<double > singularIntegral = bem1.singular(0, j);
					for (int k = 0; k <= o; k++) {
						S(i + shift0, o * j + k + shift1) += singularIntegral[k];
						D(i + shift0, o * j + k + shift1) += singularIntegral[o + 1 + k];
					}
				}
			}
			else {
				int j = i / o;
				std::vector<double > singularIntegral = bem1.singular(bem1.e()[j].t()[1], j);
				for (int k = 0; k <= o; k++) {
					S(i + shift0, o * j + k + shift1) += singularIntegral[k];
					D(i + shift0, o * j + k + shift1) += singularIntegral[o + 1 + k];
				}
			}
			for (int j = i / o + 1; j < nElem; j++) {
				const std::vector<double > regularIntegral = bem1.regular(r(i), z(i), j);
				for (int k = 0; k <= o; k++) {
					S(i + shift0, o * j + k + shift1) += regularIntegral[k];
					D(i + shift0, o * j + k + shift1) += regularIntegral[o + 1 + k];
				}
			}
		*/
		}

		

	
		
	}
}


//------------local integration -------------

const std::vector<double > Bem::regular(double rp, double zp, int idElement) const  {
	const Element &ee = _e[idElement];
	const int nqd = ee.r().size();
	const int &o = ee.order();
	const std::vector<std::vector<double>> &basis = ee.basis();
	std::vector<double > output(2 * (o + 1));	

	for (int k = 0; k < nqd; k++) {
		const double &r = ee.r()[k], &z = ee.z()[k], &dr = ee.dr()[k], &dz = ee.dz()[k], &J = ee.J()[k];
		const double &ab = Numeric::qd_GL_x[nqd][k], &wt = Numeric::qd_GL_w[nqd][k];

		double a, b, m;
		double K, E;
		double f_single_K, f_double_K, f_double_E;
		
		abm(rp, zp, r, z, a, b, m);		
		Numeric::KEPQ(m, K, E);		
		fKE(rp, zp, r, z, dr, dz, J, a, b, f_single_K, f_double_K, f_double_E);		
		for (int j = 0; j <= o; j++) {
			const double &N = basis[j][k];
			output[j] += f_single_K * K * N * wt;
			output[o + 1 + j] += (f_double_K * K + f_double_E * E) * N * wt;
		}
	}
	return output;
};

const std::vector<double > Bem::axis(double zp, int idElement) const  {
	const Element &e = _e[idElement];
	const int nqd = e.r().size();	
	const int &o = e.order();
	const std::vector<std::vector<double>> &basis = e.basis();
	std::vector<double > output(2 * (o + 1));

	for (int k = 0; k < nqd; k++) {
		const double &r = e.r()[k];
		const double &z = e.z()[k];
		const double &dr = e.dr()[k];
		const double &dz = e.dz()[k];
		const double &J = e.J()[k];

		const double ab = Numeric::qd_GL_x[nqd][k];
		const double wt = Numeric::qd_GL_w[nqd][k];

		double f_single_axis = J * r / 2. / sqrt(r * r + (z - zp) * (z - zp));
		double f_double_axis = r / 2. *(dz * r + dr *(zp -z))/ pow(r * r + (z - zp) * (z - zp), 1.5);		

		for (int j = 0; j <= o; j++) {
			const double &N = basis[j][k];
			output[j] += f_single_axis * N *  wt;
			output[o + 1 + j] += f_double_axis * N * wt;
		}
	}
	return output;
};

const std::vector<double > Bem::singular(double tau, int idElement) const {


	const Element &ee = _e[idElement];	
	const int &o = ee.order();
	std::vector<double > output(2 * (o + 1));
	double rp = sp().d(sp().x(), idElement, tau)(0);
	double zp = sp().d(sp().y(), idElement, tau)(0);
	if (abs(tau) < 1e-12) {
		rp = node().r(idElement * o, 0);
		zp = node().z(idElement * o, 0);
	}
	if (abs(1. - tau) < 1e-12) {
		rp = node().r((idElement + 1) * o , 0);
		zp = node().z((idElement + 1) * o, 0);
	}
	
	
	const int nqd_regular = 20;//settings.qdOrder() * 2;
	const int nqd_singular = 20;//settings.qdOrder() * 2;
	
	switch (o)	{
	case 1: {	
		double * qdx = NULL;
		qdx = new double[nqd_regular];
		for (int l = 0; l < nqd_regular; l++) {
			const double &ab = Numeric::qd_GL_x[nqd_regular][l];
			qdx[l] = ab;
		}
		Element e(ee);
		e.init(_sp, idElement, nqd_regular, qdx);
		for (int k = 0; k < nqd_regular; k++) {			
			const double &r = e.r()[k], &z = e.z()[k], &dr = e.dr()[k], &dz = e.dz()[k], &J = e.J()[k];
			const double &ab = Numeric::qd_GL_x[nqd_regular][k],  &wt = Numeric::qd_GL_w[nqd_regular][k];
			double a, b, m, K, E, PK, QK, PE, QE, RK,RE, f_single_K, f_double_K, f_double_E;			
			abm(rp, zp, r, z, a, b, m);			
			Numeric::KEPQ(m, K, E, PK, QK, PE, QE);
			RK = RKE(PK, QK, m, ab, tau);
			RE = RKE(PE, QE, m, ab, tau);			
			fKE(rp, zp, r, z, dr, dz, J, a, b, f_single_K, f_double_K, f_double_E);
			for (int j = 0; j <= o; j++) {
				const double &N = e.basis()[j][k];
				output[j] += f_single_K * RK * N * wt;
				output[o + 1 + j] += (f_double_K * RK + f_double_E * RE) * N * wt;
			}
		}		
		delete[] qdx;
		qdx = NULL;
		
		qdx = new double[nqd_singular];		
		for (int l = 0; l < nqd_singular; l++) {
			const double &ab = Numeric::qd_LOG_x[nqd_singular][l];
			qdx[l] = (1. - tau) * ab + tau * (1. -ab);
		}
		//Element e(ee);
		e.init(_sp, idElement, nqd_singular, qdx);
		for (int k = 0; k < nqd_singular; k++) {			
			const double &r = e.r()[k], &z = e.z()[k], &dr = e.dr()[k], &dz = e.dz()[k], &J = e.J()[k];
			const double &wt = Numeric::qd_LOG_w[nqd_singular][k];
			double a, b, m, K, E, PK, QK, PE, QE, f_single_K, f_double_K, f_double_E;			
			abm(rp, zp, r, z, a, b, m);			
			Numeric::KEPQ(m, K, E, PK, QK, PE, QE);			
			fKE(rp, zp, r, z, dr, dz, J, a, b, f_single_K, f_double_K, f_double_E);
			for (int j = 0; j <= o; j++) {
				const double &N = e.basis()[j][k];
				output[j] += 2. * f_single_K * QK * N * wt;
				output[o + 1 + j] += 2. * (f_double_K * QK + f_double_E * QE) * N * wt;
			}
		}
		delete[] qdx;  
		qdx = NULL;

		break;
	}
	case 2: {		
		if ((abs(tau) < 5.e-15) || (abs(1. - tau) < 5.e-15)) {
			double * qdx = NULL;
			qdx = new double[nqd_regular];
			for (int l = 0; l < nqd_regular; l++) {
				const double &ab = Numeric::qd_GL_x[nqd_regular][l];
				qdx[l] = ab;
			}
			Element e(ee);
			e.init(_sp, idElement, nqd_regular, qdx);
			for (int k = 0; k < nqd_regular; k++) {
				const double &r = e.r()[k], &z = e.z()[k], &dr = e.dr()[k], &dz = e.dz()[k], &J = e.J()[k];
				const double &ab = Numeric::qd_GL_x[nqd_regular][k], &wt = Numeric::qd_GL_w[nqd_regular][k];
				double a, b, m, K, E, PK, QK, PE, QE, RK, RE, f_single_K, f_double_K, f_double_E;
				abm(rp, zp, r, z, a, b, m);
				Numeric::KEPQ(m, K, E, PK, QK, PE, QE);
				RK = RKE(PK, QK, m, ab, tau);
				RE = RKE(PE, QE, m, ab, tau);
				fKE(rp, zp, r, z, dr, dz, J, a, b, f_single_K, f_double_K, f_double_E);
				for (int j = 0; j <= o; j++) {
					const double &N = e.basis()[j][k];
					output[j] += f_single_K * RK * N * wt;
					output[o + 1 + j] += (f_double_K * RK + f_double_E * RE) * N * wt;
				}
			}
			/*============= switch a new pointer  =============*/
			delete[] qdx;
			qdx = NULL;
			/*============= switch a new pointer  =============*/
			qdx = new double[nqd_singular];
			for (int l = 0; l < nqd_singular; l++) {
				const double &ab = Numeric::qd_LOG_x[nqd_singular][l];
				qdx[l] = (1. - tau) * ab + tau * (1. - ab);
			}
			//Element e(ee);
			e.init(_sp, idElement, nqd_singular, qdx);
			for (int k = 0; k < nqd_singular; k++) {
				const double &r = e.r()[k], &z = e.z()[k], &dr = e.dr()[k], &dz = e.dz()[k], &J = e.J()[k];
				const double &wt = Numeric::qd_LOG_w[nqd_singular][k];
				double a, b, m, K, E, PK, QK, PE, QE, f_single_K, f_double_K, f_double_E;
				abm(rp, zp, r, z, a, b, m);
				Numeric::KEPQ(m, K, E, PK, QK, PE, QE);
				fKE(rp, zp, r, z, dr, dz, J, a, b, f_single_K, f_double_K, f_double_E);
				for (int j = 0; j <= o; j++) {
					const double &N = e.basis()[j][k];
					output[j] += 2. * f_single_K * QK * N * wt;
					output[o + 1 + j] += 2. * (f_double_K * QK + f_double_E * QE) * N * wt;
				}
			}
			delete[] qdx;
			qdx = NULL;

		}
		else {
			//printf("%d %f\n" , idElement , tau);
			double * qdx = NULL;
			Element e(ee);
			qdx = new double[nqd_regular];
			/*============= regular 0 to tau =============*/				
			for (int l = 0; l < nqd_regular; l++) {
				const double &ab = Numeric::qd_GL_x[nqd_regular][l];
				qdx[l] = tau * ab;
			}
			e.init(_sp, idElement, nqd_regular, qdx);
			for (int k = 0; k < nqd_regular; k++) {
				const double &r = e.r()[k], &z = e.z()[k], &dr = e.dr()[k], &dz = e.dz()[k], &J = e.J()[k];
				const double &ab = qdx[k], &wt = Numeric::qd_GL_w[nqd_regular][k];
				double a, b, m, K, E, PK, QK, PE, QE, RK, RE, f_single_K, f_double_K, f_double_E;
				abm(rp, zp, r, z, a, b, m);
				Numeric::KEPQ(m, K, E, PK, QK, PE, QE);
				RK = RKE(PK, QK, m, ab, tau);
				RE = RKE(PE, QE, m, ab, tau);
				fKE(rp, zp, r, z, dr, dz, J, a, b, f_single_K, f_double_K, f_double_E);
				for (int j = 0; j <= o; j++) {
					const double &N = e.basis()[j][k];
					output[j] += f_single_K * RK * N * wt * tau;
					output[j] += -2. * f_single_K * QK * N * wt * tau * log(tau);
					output[o + 1 + j] += (f_double_K * RK + f_double_E * RE) * N * wt * tau;					
					output[o + 1 + j] += -2. * (f_double_K * QK + f_double_E * QE) * N * wt * tau * log(tau);

				}
			}
			/*============= regular tau to 1 =============*/
			for (int l = 0; l < nqd_regular; l++) {
				const double &ab = Numeric::qd_GL_x[nqd_regular][l];
				qdx[l] = tau + (1. - tau) * ab;
			}
			e.init(_sp, idElement, nqd_regular, qdx);
			for (int k = 0; k < nqd_regular; k++) {
				const double &r = e.r()[k], &z = e.z()[k], &dr = e.dr()[k], &dz = e.dz()[k], &J = e.J()[k];
				const double &ab = qdx[k], &wt = Numeric::qd_GL_w[nqd_regular][k];
				double a, b, m, K, E, PK, QK, PE, QE, RK, RE, f_single_K, f_double_K, f_double_E;
				abm(rp, zp, r, z, a, b, m);
				Numeric::KEPQ(m, K, E, PK, QK, PE, QE);
				RK = RKE(PK, QK, m, ab, tau);
				RE = RKE(PE, QE, m, ab, tau);
				fKE(rp, zp, r, z, dr, dz, J, a, b, f_single_K, f_double_K, f_double_E);
				for (int j = 0; j <= o; j++) {
					const double &N = e.basis()[j][k];
					output[j] += f_single_K * RK * N * wt * (1. - tau);
					output[j] += -2. * f_single_K * QK * N * wt * (1. - tau) * log(1. - tau);
					output[o + 1 + j] += (f_double_K * RK + f_double_E * RE) * N * wt * (1. - tau) ;
					output[o + 1 + j] += -2. * (f_double_K * QK + f_double_E * QE) * N * wt * (1. - tau) * log(1. - tau);

				}
			}
			/*============= switch a new pointer  =============*/
			delete[] qdx;
			qdx = NULL;
			qdx = new double[nqd_singular];
			/*============= singular 0 to tau  =============*/			
			for (int l = 0; l < nqd_singular; l++) {
				const double &ab = Numeric::qd_LOG_x[nqd_singular][l];
				qdx[l] = tau * (1. - ab);
			}
			e.init(_sp, idElement, nqd_singular, qdx);
			for (int k = 0; k < nqd_singular; k++) {
				const double &r = e.r()[k], &z = e.z()[k], &dr = e.dr()[k], &dz = e.dz()[k], &J = e.J()[k];
				const double &wt = Numeric::qd_LOG_w[nqd_singular][k];
				double a, b, m, K, E, PK, QK, PE, QE, f_single_K, f_double_K, f_double_E;
				abm(rp, zp, r, z, a, b, m);
				Numeric::KEPQ(m, K, E, PK, QK, PE, QE);
				fKE(rp, zp, r, z, dr, dz, J, a, b, f_single_K, f_double_K, f_double_E);
				for (int j = 0; j <= o; j++) {
					const double &N = e.basis()[j][k];
					output[j] += 2. * f_single_K * QK * N * wt * tau;
					output[o + 1 + j] += 2. * (f_double_K * QK + f_double_E * QE) * N * wt * tau;
				}
			}
			/*============= singular tau to 1  =============*/
			for (int l = 0; l < nqd_singular; l++) {
				const double &ab = Numeric::qd_LOG_x[nqd_singular][l];
				qdx[l] = tau + (1. - tau) * ab;
			}
			e.init(_sp, idElement, nqd_singular, qdx);
			for (int k = 0; k < nqd_singular; k++) {
				const double &r = e.r()[k], &z = e.z()[k], &dr = e.dr()[k], &dz = e.dz()[k], &J = e.J()[k];
				const double &wt = Numeric::qd_LOG_w[nqd_singular][k];
				double a, b, m, K, E, PK, QK, PE, QE, f_single_K, f_double_K, f_double_E;
				abm(rp, zp, r, z, a, b, m);
				Numeric::KEPQ(m, K, E, PK, QK, PE, QE);
				fKE(rp, zp, r, z, dr, dz, J, a, b, f_single_K, f_double_K, f_double_E);
				for (int j = 0; j <= o; j++) {
					const double &N = e.basis()[j][k];
					output[j] += 2. * f_single_K * QK * N * wt * ( 1. - tau);
					output[o + 1 + j] += 2. * (f_double_K * QK + f_double_E * QE) * N * wt * (1. - tau);
				}
			}





			
			//*============ free memory 
			delete[] qdx;
			qdx = NULL;

			

		}
		
		
		
		
		
		
		
		
		
		
		
	



		break;


	}





	/*case 2: {//printf("tau %1.5f\t %04d\n", tau, idElement);
		if (abs(tau) < 1e-14 || abs(1.0 - tau) < 1e-14) {
			qd = Element::getGLQuad(nqd_regular);
			e.init(_sp, idElement, qd);
			sKdKdE(rp, zp, e, sK, dK, dE, m);
			RKRE(m, Eigen::VectorXd(qd.row(0)), tau, RK, RE);
			basis = e.basis().transpose();
			for (int k = 0; k <= settings.order(); k++) {
				output(k) += qd.row(1).dot(Eigen::VectorXd(
					RK.array() * sK.array()* basis.col(k).array()));
				output(settings.order() + 1 + k) += qd.row(1).dot(Eigen::VectorXd(
					(RK.array() * dK.array() + RE.array() * dE.array())	* basis.col(k).array()));
			}
			//std::cout << tau << " m\t" << m.transpose() << "\n";
			qd = Element::getLogQuad(nqd_singular);
			qd.row(0).array() = tau + qd.row(0).array() * (1. - 2.*tau);
			e.init(_sp, idElement, qd);
			sKdKdE(rp, zp, e, sK, dK, dE, m);
			KEPKQKPEQE(m, K, E, PK, QK, PE, QE);
			basis = e.basis().transpose();
			for (int k = 0; k <= settings.order(); k++) {
				output(k) += 2 * qd.row(1).dot(
					Eigen::VectorXd(sK.array()* QK.array()* basis.col(k).array()));
				output(settings.order() + 1 + k) += 2 * qd.row(1).dot(
					Eigen::VectorXd((dE.array()* QE.array() + dK.array()* QK.array())* basis.col(k).array()
					));
			}

		}
		else {
			//printf("shit");
			qd = Element::getGLQuad(nqd_regular);
			qd.row(0).array() *= tau;
			qd.row(1).array() *= tau;
			e.init(_sp, idElement, qd);
			sKdKdE(rp, zp, e, sK, dK, dE, m);
			RKRE(m, Eigen::VectorXd(qd.row(0)), tau, RK, RE);
			basis = e.basis().transpose();
			for (int k = 0; k <= settings.order(); k++) {
				output(k) += qd.row(1).dot(
					Eigen::VectorXd(RK.array() * sK.array()* basis.col(k).array()));
				output(settings.order() + 1 + k) += qd.row(1).dot(
					Eigen::VectorXd((RK.array() * dK.array() + RE.array() * dE.array())	* basis.col(k).array()));
			}
			qd = Element::getGLQuad(nqd_regular);
			qd.row(0).array() = tau + (1. - tau) * qd.row(0).array();
			qd.row(1).array() *= 1. - tau;
			e.init(_sp, idElement, qd);
			sKdKdE(rp, zp, e, sK, dK, dE, m);
			RKRE(m, Eigen::VectorXd(qd.row(0)), tau, RK, RE);
			basis = e.basis().transpose();
			for (int k = 0; k <= settings.order(); k++) {
				output(k) += qd.row(1).dot(
					Eigen::VectorXd(RK.array() * sK.array()* basis.col(k).array()));
				output(settings.order() + 1 + k) += qd.row(1).dot(
					Eigen::VectorXd((RK.array() * dK.array() + RE.array() * dE.array())	* basis.col(k).array()));
			}
			//
			//printf("sss");		
			qd = Element::getGLQuad(nqd_regular);
			qd.row(0).array() *= tau;
			qd.row(1).array() *= tau * log(tau);
			e.init(_sp, idElement, qd);
			sKdKdE(rp, zp, e, sK, dK, dE, m);
			KEPKQKPEQE(m, K, E, PK, QK, PE, QE);
			basis = e.basis().transpose();
			for (int k = 0; k <= settings.order(); k++) {
				output(k) += -2 * qd.row(1).dot(
					Eigen::VectorXd(sK.array()* QK.array()* basis.col(k).array()));
				output(settings.order() + 1 + k) += -2 * qd.row(1).dot(
					Eigen::VectorXd((dE.array()* QE.array() + dK.array()* QK.array())* basis.col(k).array()
					));
			}
			qd = Element::getGLQuad(nqd_regular);
			qd.row(0).array() = (1. - tau) * qd.row(0).array() + tau;
			qd.row(1).array() *= (1. - tau) * log(1. - tau);
			e.init(_sp, idElement, qd);
			sKdKdE(rp, zp, e, sK, dK, dE, m);
			KEPKQKPEQE(m, K, E, PK, QK, PE, QE);
			basis = e.basis().transpose();
			for (int k = 0; k <= settings.order(); k++) {
				output(k) += -2 * qd.row(1).dot(
					Eigen::VectorXd(sK.array()* QK.array()* basis.col(k).array()));
				output(settings.order() + 1 + k) += -2 * qd.row(1).dot(
					Eigen::VectorXd((dE.array()* QE.array() + dK.array()* QK.array())* basis.col(k).array()
					));
			}

			qd = Element::getLogQuad(nqd_singular);
			qd.row(0).array() = tau * (1. -  qd.row(0).array());
			qd.row(1).array() *= tau;
			e.init(_sp, idElement, qd);
			sKdKdE(rp, zp, e, sK, dK, dE, m);
			KEPKQKPEQE(m, K, E, PK, QK, PE, QE);
			basis = e.basis().transpose();
			for (int k = 0; k <= settings.order(); k++) {
				output(k) += 2 * qd.row(1).dot(
					Eigen::VectorXd(sK.array()* QK.array()* basis.col(k).array()));
				output(settings.order() + 1 + k) += 2 * qd.row(1).dot(
					Eigen::VectorXd((dE.array()* QE.array() + dK.array()* QK.array())* basis.col(k).array()
					));
			}
			qd = Element::getLogQuad(nqd_singular);
			qd.row(0).array() = tau + (1. - tau) * qd.row(0).array();
			qd.row(1) = qd.row(1) * ( 1. - tau );
			e.init(_sp, idElement, qd);
			sKdKdE(rp, zp, e, sK, dK, dE, m);
			KEPKQKPEQE(m, K, E, PK, QK, PE, QE);
			basis = e.basis().transpose();
			for (int k = 0; k <= settings.order(); k++) {
				output(k) += 2 * qd.row(1).dot(
					Eigen::VectorXd(sK.array()* QK.array()* basis.col(k).array()));
				output(settings.order() + 1 + k) += 2 * qd.row(1).dot(
					Eigen::VectorXd((dE.array()* QE.array() + dK.array()* QK.array())* basis.col(k).array()
					));
			}

		}







		break;
	}
	*/
	
	}
		
	
	/*qd = Element::getLogQuad(nqd_singular);
	qd.row(0).array() = tau + (1. - 2.*tau) *qd.row(0).array();
	e.init(_sp, idElement, qd);
	sKdKdE(rp, zp, e, sK, dK, dE, m);
	KEPKQKPEQE(m, K, E, PK, QK, PE, QE);
	basis = e.basis().transpose();
	for (int k = 0; k <= settings.order(); k++) {
		output(k) += 2 * qd.row(1).dot(
			Eigen::VectorXd(sK.array()* QK.array()* basis.col(k).array()));
		output(settings.order() + 1 + k) += 2 * qd.row(1).dot(
			Eigen::VectorXd((dE.array()* QE.array() + dK.array()* QK.array())* basis.col(k).array()
			));
	}*/

	//std::cout << "singular\n" << output <<"\n----------------";
	return output;
}

//------------helper functions -------------

void Bem::abm(double rp, double zp, double r, double z, double &a, double &b, double &m) {
	a = r * r + rp * rp + (zp - z) * (zp - z);
	b = 2. * rp * r;
	m = 2. * b / (a + b);
};

double Bem::RKE(double P, double Q, double m, double t, double tau) {
	return P - Q * log((1. - m) / (t * t  - 2 * tau * t + tau * tau));
};

void Bem::fKE(double rp, double zp, double r, double z, double dr, double dz, double J, double a, double b,
	double &f_single_K, double &f_double_K, double &f_double_E) {
	f_single_K = r * J / M_PI / sqrt(a + b);
	f_double_K = dz / 2. / M_PI / sqrt(a + b);
	f_double_E = ((dr * (zp - z) - dz * (rp - r)) * r / (a - b) - dz / 2.) / M_PI / sqrt(a + b);
};

//------------properties -------------

void Bem::Properties::bc::_bc::print() const {
	switch (type) {
	case Spline::BC::Even:
		printf("Even ");
		break;
	case Spline::BC::Odd:
		printf("Odd ");
		break;
	case Spline::BC::Mix:
		printf("Mix a(%.3e) b(%.3e) ", a, b);
		break;
	default:
		break;
	}
}

void Bem::Properties::print() const {
	printf("Index shift: %d\n", _indexShift);
	printf("Elements#: %d\n", _nElm);
	printf("Element order: %d\n", _order);
	printf("Preprocess G-L quadrature order: %d\n", _qdOrder);
	printf("xBC: ");
	xBC.begin.print();
	xBC.end.print();
	printf("\n");
	printf("yBC: ");
	yBC.begin.print();
	yBC.end.print();
	printf("\n");
}

void Bem::Properties::tic() { printf("tic\n"); _timer = omp_get_wtime(); };
void Bem::Properties::toc() { printf("toc ... %.2e sec\n", omp_get_wtime() - _timer); };

double Bem::Properties::_timer = 0;
const double Bem::eps = 1E-14;