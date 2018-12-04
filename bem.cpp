#include "stdafx.h"
#include "bem.h"
#include "numeric.h"
#include <omp.h>

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

	int idSingular = 0;
	//Element ee(_e[idSingular]);
	//const Eigen::Matrix2Xd qd = Element::getGLQuad(14);
	//ee.init(_sp, 1, qd);
	double rp = sp().d(sp().x(), idSingular, 0.0)(0);
	double zp = sp().d(sp().y(), idSingular,0.0)(0);
	std::cout.precision(15);
	std::cout <<regular(rp, zp,10) <<"\n";
	printf("tmid = %15.15f\n",_e[idSingular].t()(1));
	//settings.tic();
	//singular(_e[idSingular].t()(1), idSingular);
	//settings.toc();
	//singular(0., idSingular);
	//singular(1., idSingular);
}


void Bem::assembly(Eigen::MatrixXd &S, Eigen::MatrixXd &D) {

	int nNode = settings.nElm() * settings.order() + 1;
	Eigen::VectorXd r; r.setZero(nNode);
	Eigen::VectorXd z; z.setZero(nNode);
	for (int idElement = 0; idElement < settings.nElm(); idElement++) {
		const Eigen::VectorXd &t = e()[idElement].t();
		for (int k = 0; k < t.size(); k++) {
			int idNode = idElement * settings.order() + k;			
			r(idNode) = sp().d(sp().x(), idElement, t(k))(0);
			z(idNode) = sp().d(sp().y(), idElement, t(k))(0);
		}
	}
	S.setZero(nNode, nNode);
	D.setZero(nNode, nNode);

	for (int idSource = 0; idSource < nNode; idSource++) {
		for (int k = 0; k <= settings.order(); k++) {

		}

	}
	


}

//------------local integration -------------

Eigen::VectorXd Bem::regular(double rp, double zp, int idElement) {
	const Element &e = _e[idElement];
	const Eigen::Matrix2Xd qd = Element::getGLQuad(e.r().size());	
	Eigen::VectorXd sK, dK, dE, K, E;
	sKdKdE(rp,zp,e, sK, dK, dE, K, E);
	const Eigen::MatrixXd &basis = e.basis().transpose();		

	Eigen::VectorXd output(2 * (settings.order() + 1));
	for (int k = 0; k <= settings.order(); k++) {
		output(k) =	qd.row(1).dot(
			Eigen::VectorXd(K.array() * sK.array()* basis.col(k).array())
		);

		output(settings.order() + 1 + k) = qd.row(1).dot(Eigen::VectorXd(
				(K.array() * dK.array() + E.array() * dE.array())* basis.col(k).array()
			));
	}
	return output;
};

Eigen::VectorXd Bem::singular(double tau, int idElement) {
	double rp = sp().d(sp().x(), idElement, tau)(0);
	double zp = sp().d(sp().y(), idElement, tau)(0);
	Eigen::VectorXd output(2 * (settings.order() + 1)); output.setZero();
	const Element &ee = _e[idElement];
	Element e(ee);

	Eigen::Matrix2Xd qd;
	Eigen::MatrixXd basis;
	Eigen::VectorXd sK, dK, dE, m, RK, RE, QK, QE, PK, PE,K,E;

	
	int nqd_regular = 10;
	
	if (abs(tau) < eps || abs(1.0 - tau) < eps) {
		qd = Element::getGLQuad(nqd_regular);
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
	}
	else {
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
		qd.row(0).array() = tau + (1.-tau) * qd.row(0).array();
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

	}

	int nqd_singular = 2;

	// singular 
	if (abs(tau) < eps || abs(1.0 - tau) < eps) {
		qd = Element::getLogQuad(nqd_singular);
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
		}
	}	
	else {
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
		qd.row(0).array() =  tau * (1. - qd.row(0).array());
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
		qd.row(1).array() *= 1. - tau;
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

void Bem::sKdKdE(double rp, double zp, const Element &e,
	Eigen::VectorXd &sK,	Eigen::VectorXd &dK,	Eigen::VectorXd &dE,
	Eigen::VectorXd &K,	Eigen::VectorXd &E) {
	const Eigen::VectorXd &r = e.r();
	const Eigen::VectorXd &dr = e.dr();
	const Eigen::VectorXd &z = e.z();
	const Eigen::VectorXd &dz = e.dz();
	const Eigen::VectorXd &J = e.J();
	const Eigen::VectorXd &xi = e.xi();
	

	if (abs(rp) > eps ) {	
		Eigen::VectorXd a = r.array().square() + rp * rp + (z.array() - zp).square();
		Eigen::VectorXd b = r.array() * 2 * rp;
		Eigen::VectorXd m = 2. * b.array() / (a + b).array();
		Eigen::VectorXd pi_sqrt_a_plus_b = M_PI * (a + b).array().sqrt();
		sK = r.array() * J.array() / pi_sqrt_a_plus_b.array();
		dK = 0.5 * dz.array() / pi_sqrt_a_plus_b.array();
		dE = (
			(dr.array() *(zp - z.array()) - dz.array() *(rp - r.array()))
			* r.array() / (a - b).array() - 0.5 * dz.array()) / pi_sqrt_a_plus_b.array();
		KEPKQKPEQE(m, K, E);
	}
	else {
		E.setZero(r.size());
		dE.setZero(r.size());
		Eigen::VectorXd tmp = (1.0 + ((z.array() - zp) / r.array()).square()).rsqrt();
		K = 0.5 * tmp;
		sK = J;
		dK = (dz.array()/r.array() + dr.array() * (zp - z.array())/(r.array().square()))*
			tmp.array().square()	;
	}
};

void Bem::sKdKdE(double rp, double zp, const Element &e,
	Eigen::VectorXd &sK, Eigen::VectorXd &dK, Eigen::VectorXd &dE,
	Eigen::VectorXd &m) {
	const Eigen::VectorXd &r = e.r();
	const Eigen::VectorXd &dr = e.dr();
	const Eigen::VectorXd &z = e.z();
	const Eigen::VectorXd &dz = e.dz();
	const Eigen::VectorXd &J = e.J();
	const Eigen::VectorXd &xi = e.xi();

	Eigen::VectorXd a = r.array().square() + rp * rp + (z.array() - zp).square();
	Eigen::VectorXd b = r.array() * 2 * rp;
	m = 2. * b.array() / (a + b).array();
	Eigen::VectorXd pi_sqrt_a_plus_b = M_PI * (a + b).array().sqrt();
	sK = r.array() * J.array() / pi_sqrt_a_plus_b.array();
	dK = 0.5 * dz.array() / pi_sqrt_a_plus_b.array();
	dE = (
		(dr.array() *(zp - z.array()) - dz.array() *(rp - r.array()))
		* r.array() / (a - b).array() - 0.5 * dz.array()) / pi_sqrt_a_plus_b.array();

	
};

void Bem::RKRE(const Eigen::VectorXd &m, const Eigen::VectorXd &t, double tau, Eigen::VectorXd &RK, Eigen::VectorXd &RE) {
	Eigen::VectorXd K, E, PK, QK, PE, QE;	
	KEPKQKPEQE(m, K, E, PK, QK, PE, QE);
	RK.setZero(m.size());
	RE.setZero(m.size());
	Eigen::VectorXd LOG(((1. - m.array()) / (t.array() - tau).square()).log());
	RK = PK.array() - QK.array() *LOG.array();
	RE = PE.array() - QE.array() *LOG.array();
};

void Bem::KEPKQKPEQE( const Eigen::VectorXd &m,
	Eigen::VectorXd &K,  Eigen::VectorXd &E,
	Eigen::VectorXd &PK, Eigen::VectorXd &QK,
	Eigen::VectorXd &PE, Eigen::VectorXd &QE) {
	K.resize(m.size());
	E.resize(m.size());
	PK.resize(m.size());
	QK.resize(m.size());
	PE.resize(m.size());
	QE.resize(m.size());
	for (int k = 0; k < m.size(); k++) {
		Numeric::KEPQ(m(k), K(k),E(k), PK(k), QK(k), PE(k), QE(k));
	}
};

void Bem::KEPKQKPEQE(const Eigen::VectorXd &m, Eigen::VectorXd &K, Eigen::VectorXd &E) {
	K.resize(m.size());
	E.resize(m.size());
	for (int k = 0; k < m.size(); k++) {
		Numeric::KEPQ(m(k), K(k), E(k));
	}
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

double Bem::Properties::_timer = 0;
const double Bem::eps = 1e-14;