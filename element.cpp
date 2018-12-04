#include "stdafx.h"
#include "element.h"
#include "numeric.h"

void Element::init(const Spline &sp, int i, int o, int nqd) {	
	_arc = sp.localArc(i);
	_order = o;	
	_t.setZero(o + 1);	
	_t(1) = sp.arc2t(i, 0.5 * _arc);
	_t(o) = 1.0;
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
		for (int j = 0; j <= o; j++) {
			_basis(j, k) = Numeric::N[o][j](_xi(k));
		}
	}

}

Element::Element(const Element &e) {	
	_arc = e.arc();
	_t = e.t();	
	_order = e.order();
};

void Element::init(const Spline &sp, int i, const Eigen::Matrix2Xd &qd) {
	int nqd = qd.cols(); 	
	int o = _order; 
	_r.resize(nqd);
	_dr.resize(nqd);
	_ddr.resize(nqd);
	_z.resize(nqd);
	_dz.resize(nqd);
	_ddz.resize(nqd); 
	_J.resize(nqd);
	_xi.resize(nqd);
	_basis.resize( o + 1, nqd);
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
		for (int j = 0; j <= o; j++) {
			_basis(j, k) = Numeric::N[o][j](_xi(k));			
		}
	
	}

}


const Eigen::Matrix2Xd Element::getGLQuad(int nqd) {
	return Eigen::Map<const Eigen::Matrix2Xd>(Numeric::qd[nqd], 2, nqd);
};

const Eigen::Matrix2Xd Element::getLogQuad(int nqd) {
	return Eigen::Map<const Eigen::Matrix2Xd>(Numeric::lqd[nqd], 2, nqd * 5);
};