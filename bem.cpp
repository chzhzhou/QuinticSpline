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
	for (unsigned int i = 0; i < _e.size(); i++) {
		_e[i].init(_sp, i, settings.order(), settings.qdOrder());
	}
	
	//int idSingular = 0;
	//Element ee(_e[idSingular]);
	//const Eigen::Matrix2Xd qd = Element::getGLQuad(14);
	//ee.init(_sp, 1, qd);
	//double rp = sp().d(sp().x(), idSingular, 0.0)(0);
	//double zp = sp().d(sp().y(), idSingular,0.0)(0);
	//std::cout.precision(15);
	//std::cout <<regular(rp, zp,10) <<"\n";
	//printf("tmid = %15.15f\n",_e[idSingular].t()(1));
	//settings.tic();
	//singular(_e[idSingular].t()(1), idSingular);
	//settings.toc();
	//singular(0., idSingular);
	//singular(1., idSingular);
	setrz();
	//std::cout << node().r;
	//Eigen::VectorXd out = singular(e()[0].t()(1), 0);
	//out = regular(node().r(1,0), node().z(1, 0),15);
	//std::cout << node().r(1, 0) << "abc" << node().z(1, 0) << "abd";
	

	//std::cout << out;
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

void Bem::assembly(Eigen::MatrixXd &S, Eigen::MatrixXd &D) {
	const int nNode = settings.nElm() * settings.order() + 1;
	const int nElem = settings.nElm();
	const int nOrder = settings.order();
	const Eigen::VectorXd &r = node().r.col(0);
	const Eigen::VectorXd &z = node().z.col(0);
	S.setZero(nNode, nNode);
	D.setZero(nNode, nNode);

	
	//std::cout << "asdfasdfasdf" <<nNode;
#pragma omp parallel for 
	for (int idSource = 0; idSource < nNode; idSource++) {
		//printf("wtf %02d %16.16f: \n", idSource, abs(r(idSource)));
		//printf("source %02d: %02d:,%02d:\n ", idSource, idElementHasSource_Low, idElementHasSource_Up);
		if (abs(r(idSource)) < 1e-13) {
			//printf("Axis %03d %16.16f: \n", idSource, abs(r(idSource)));
			for (int idElement = 0; idElement < nElem; idElement++) {
				const std::vector<double > axisIntegral = axis(z(idSource), idElement);
				for (int local = 0; local <= nOrder; local++) {
					int idReceiver = nOrder * idElement + local;
					S(idSource, idReceiver) += axisIntegral[local];
					D(idSource, idReceiver) += axisIntegral[nOrder + 1 + local];
					
					
					//printf("{%02d,%02d} ", idSource, idReceiver);

				}					
				//printf("\n");

			}
		}
		else {
			int idElementHasSource_Low, idElementHasSource_Up;
			if (idSource % nOrder == 0) {
				idElementHasSource_Up = Numeric::clamp(idSource / nOrder, 0, nElem - 1);
				idElementHasSource_Low = Numeric::clamp(idSource / nOrder - 1, 0, nElem - 1);
			}
			else {
				idElementHasSource_Up = Numeric::clamp(idSource / nOrder, 0, nElem - 1);
				idElementHasSource_Low = idElementHasSource_Up;
			}
		
			for (int idElement = 0; idElement < nElem; idElement++) {
				
				if (idElement < idElementHasSource_Low || idElement > idElementHasSource_Up) {	//regular
					const std::vector<double > regularIntegral = regular(r(idSource), z(idSource), idElement);
					for (int local = 0; local <= nOrder; local++) {
						int idReceiver = nOrder * idElement + local;
						S(idSource, idReceiver) += regularIntegral[local];
						D(idSource, idReceiver) += regularIntegral[nOrder + 1 + local];
						//printf("(%02d,%02d) ", idSource, idReceiver);
					}
					//printf("\n");

				}
				else {		//singular
					int idSourcelocal = idSource - nOrder * idElement;

					//printf("source %03d element%03d local%d tau %1.2f",	idSource, idElement, idSourcelocal, e()[idElement].t()[idSourcelocal]);
					const std::vector<double > singularIntegral = singular(e()[idElement].t()[idSourcelocal], idElement);
					


					for (int local = 0; local <= nOrder; local++) {
						int idReceiver = nOrder * idElement + local;
						S(idSource, idReceiver) += singularIntegral[local];
						D(idSource, idReceiver) += singularIntegral[nOrder + 1 + local];
						//printf("[%d: %02d,%02d] ", idSourcelocal, idSource, idReceiver);

					}
					//printf("\n");

				}
			}
			//printf("\n");
		}
		 //printf("\n"); 
	}
}

//------------local integration -------------

const std::vector<double > Bem::regular(double rp, double zp, int idElement) {
	const Element &e = _e[idElement];

	const int nqd = e.r().size();
//	Element::getGLQuad(nqd, qdx, qdw);	
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


		double a = r * r + rp * rp + (zp - z) * (zp - z);
		double b = 2. * rp * r;
		double m = 2. * b / (a + b);
		double K, E;
		Numeric::KEPQ(m, K, E);

		double f_single_K = r * J / M_PI / sqrt(a + b);
		double f_double_K = dz / 2. / M_PI / sqrt(a + b);
		double f_double_E = ((dr * (zp - z) - dz * (rp - r)) * r / (a - b) - dz / 2.) / M_PI / sqrt(a + b);

		for (int j = 0; j <= o; j++) {
			const double &N = basis[j][k];
			output[j] += f_single_K * K * N * wt;
			output[o + 1 + j] += (f_double_K * K + f_double_E * E) * N * wt;
		}
	}
	return output;
};


const std::vector<double > Bem::axis(double zp, int idElement) {
	const Element &e = _e[idElement];
	const int nqd = e.r().size();
	//	Element::getGLQuad(nqd, qdx, qdw);	
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
		double f_double_axis = r / 2. *(dz * r + dr *(zp -z))/ pow(r * r + (z - zp) * (z - zp),1.5);		

		for (int j = 0; j <= o; j++) {
			const double &N = basis[j][k];
			output[j] += f_single_axis * N *  wt;
			output[o + 1 + j] += f_double_axis * N * wt;
		}
	}
	return output;
};


const std::vector<double > Bem::singular(double tau, int idElement) {
	double rp = sp().d(sp().x(), idElement, tau)(0);
	double zp = sp().d(sp().y(), idElement, tau)(0);
	const Element &ee = _e[idElement];
	Element e(ee);	
	const int &o = e.order();
	std::vector<double > output(2 * (o + 1));
	std::vector<double> qqd0;
	std::vector<double> qqd1;

	const int nqd_regular = 6;
	const int nqd_singular = 6;	
	
	switch (o)	{
	case 1: {
		
		qqd0.resize(nqd_regular);
		for (int l = 0; l < qqd0.size(); l++) {
			qqd0[l] = Numeric::qd_GL_x[nqd_regular][l];
		}		
		Element e0(ee);
		e0.init(_sp, idElement,nqd_regular, qqd0);
		
		for (int k = 0; k < nqd_regular; k++) {
			
			const double &r = e0.r()[k];
			const double &z = e0.z()[k];
			const double &dr = e0.dr()[k];
			const double &dz = e0.dz()[k];
			const double &J = e0.J()[k];
			const double &ab = Numeric::qd_GL_x[nqd_regular][k];
			const double &wt = Numeric::qd_GL_w[nqd_regular][k];

			double a = r * r + rp * rp + (zp - z) * (zp - z);
			double b = 2. * rp * r;
			double m = 2. * b / (a + b);
			double K, E, PK, QK,PE,QE;
			Numeric::KEPQ(m, K, E,PK,QK,PE,QE);
			double RK = PK - QK * 2.0 * log(sqrt(1. - m) / abs(ab - tau));
			double RE = PE - QE * 2.0 * log(sqrt(1. - m) / abs(ab - tau));;

			double f_single_K = r * J / M_PI / sqrt(a + b);
			double f_double_K = dz / 2. / M_PI / sqrt(a + b);
			double f_double_E = ((dr * (zp - z) - dz * (rp - r)) * r / (a - b) - dz / 2.) / M_PI / sqrt(a + b);

			for (int j = 0; j <= o; j++) {
				const double &N = e0.basis()[j][k];
				output[j] += f_single_K * RK * N * wt;
				output[o + 1 + j] += (f_double_K * RK + f_double_E * RE) * N * wt;
			}
		}

		
		qqd1.resize(nqd_singular);
		for (int l = 0; l < qqd1.size(); l++) {
			double ab = Numeric::qd_LOG_x[nqd_singular][l];
			qqd1[l] = (1. - tau) * ab + tau * (1. -ab);
		}
		Element e1(ee);
		e1.init(_sp, idElement, nqd_singular, qqd1);
		//std::cout << e1.r().size() << " abc";
		for (int k = 0; k < nqd_singular; k++) {
			
			const double &r = e1.r()[k];
			const double &z = e1.z()[k];
			const double &dr = e1.dr()[k];
			const double &dz = e1.dz()[k];
			const double &J = e1.J()[k];
			const double &wt = Numeric::qd_LOG_w[nqd_singular][k];
		
			double a = r * r + rp * rp + (zp - z) * (zp - z);
			double b = 2. * rp * r;
			double m = 2. * b / (a + b);
			double K, E, PK, QK, PE, QE;
			Numeric::KEPQ(m, K, E, PK, QK, PE, QE);
			
			
			double f_single_K = r * J / M_PI / sqrt(a + b);
			double f_double_K = dz / 2. / M_PI / sqrt(a + b);
			double f_double_E = ((dr * (zp - z) - dz * (rp - r)) * r / (a - b) - dz / 2.) / M_PI / sqrt(a + b);
			for (int j = 0; j <= o; j++) {
				const double &N = e1.basis()[j][k];
				output[j] += 2. * f_single_K * QK * N * wt;
				output[o + 1 + j] += 2. * (f_double_K * QK + f_double_E * QE) * N * wt;
				

			}
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

double Bem::Properties::_timer = 0;
const double Bem::eps = 1E-12;