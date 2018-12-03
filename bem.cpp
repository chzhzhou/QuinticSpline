#include "stdafx.h"
#include "bem.h"
#include "numeric.h"
void Bem::Element::init(const Spline &sp, int i, int o, int nqd) {
	_arc = sp.localArc(i);
	_t.setZero(o + 1);
	_t(o) = 1.0;
	_t(1) = sp.arc2t(i, 0.5 * _arc);
	_r.resize(nqd);
	_dr.resize(nqd);
	_ddr.resize(nqd);
	_z.resize(nqd);
	_dz.resize(nqd);
	_ddz.resize(nqd);
	_J.resize(nqd);
	_xi.resize(nqd);
	_basis.resize(o + 1, nqd);
	//Eigen::Map<const Eigen::Matrix2Xd> qd(Numeric::qd[nqd], 2, nqd);
	const Eigen::Matrix2Xd qd = getGLQuad(nqd);	
	for (int k = 0; k < nqd; k++) {
		double t = qd(0, k);
		const Eigen::Vector3d d_r = sp.d(sp.x(), i, t);
		const Eigen::Vector3d d_z = sp.d(sp.y(), i, t);
		_r(k) = d_r(0);
		_dr(k) = d_r(1);
		_ddr(k) = d_r(2);
		_z(k) = d_z(0);
		_dz(k) = d_z(1);
		_ddz(k) = d_z(2);
		_J(k) = sqrt(d_r(1) * d_r(1) + d_z(1) * d_z(1));
		_xi(k) = sp.localArc(i, t) / _arc;
	}
	

}

void Bem::initialize(const Eigen::MatrixX2d &xy) {
	_sp.node(xy);
	const Properties::bc &xbc = settings.xBC;
	_sp.x(xbc.begin.type, xbc.end.type, xbc.begin.a, xbc.begin.b, xbc.end.a, xbc.end.b);
	const Properties::bc &ybc = settings.yBC;
	_sp.y(ybc.begin.type, ybc.end.type, ybc.begin.a, ybc.begin.b, ybc.end.a, ybc.end.b);
	settings.nElm(_sp.node().rows() - 1);

	e.resize(settings.nElm());
	for (unsigned int i = 0; i < e.size(); i++) {
		e[i].init(_sp, i, settings.order(),settings.qdOrder());
		//std::cout << e[i].xi().transpose() <<std::endl;
	}
}



Eigen::VectorXd Bem::regular(double rp, double zp, const Element &e) {
	const Eigen::Matrix2Xd qd = Element::getGLQuad(settings.qdOrder());
	
	const Eigen::VectorXd &r = e.r();
	const Eigen::VectorXd &dr = e.dr();
	const Eigen::VectorXd &z = e.z();
	const Eigen::VectorXd &dz = e.dz();
	const Eigen::VectorXd &J = e.J();
	const Eigen::VectorXd &xi = e.xi();
	

	Eigen::VectorXd a = r.array().square() + rp * rp + (z.array() - zp).square();
	Eigen::VectorXd b = r.array() *2 * rp;
	Eigen::VectorXd m = 2. * b.array() / (a+b).array();
	Eigen::VectorXd pi_sqrt_a_plus_b = M_PI * (a + b).array().sqrt();
	Eigen::VectorXd single_fK = r.array() * J.array() / pi_sqrt_a_plus_b.array();
	Eigen::VectorXd double_fK = 0.5 * dz.array() / pi_sqrt_a_plus_b.array();
	Eigen::VectorXd double_fE = (
		(dr.array() *(zp - z.array()) - dz.array() *(rp - r.array()))
		* r.array() / (a - b).array() - 0.5 * dz.array()) / pi_sqrt_a_plus_b.array();
	Eigen::VectorXd K, E;
	Eigen::MatrixXd N(settings.order() + 1, K.size());
	KEPKQKPEQE(m, K, E);
	
		

	



	

	return Eigen::VectorXd::Ones(1);
};




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

const Eigen::Matrix2Xd Bem::Element::getGLQuad(int nqd) {	
	return Eigen::Map<const Eigen::Matrix2Xd>(Numeric::qd[nqd], 2, nqd);
};

const Eigen::Matrix2Xd Bem::Element::getLogQuad(int nqd) {	
	return Eigen::Map<const Eigen::Matrix2Xd>(Numeric::lqd[nqd], 2, nqd * 5);
};

void Bem::KEPKQKPEQE(
	const Eigen::VectorXd &m,
	Eigen::VectorXd K,
	Eigen::VectorXd E,
	Eigen::VectorXd PK,
	Eigen::VectorXd QK,
	Eigen::VectorXd PE,
	Eigen::VectorXd QE) {
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
void Bem::KEPKQKPEQE(
	const Eigen::VectorXd &m,
	Eigen::VectorXd K,
	Eigen::VectorXd E) {
	K.resize(m.size());
	E.resize(m.size());
	for (int k = 0; k < m.size(); k++) {
		Numeric::KEPQ(m(k), K(k), E(k));
	}
};
