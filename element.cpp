#include "stdafx.h"
#include "element.h"
#include "numeric.h"

void Element::init(const Spline &sp, int i, int o, int nqd) {	
	_arc = sp.localArc(i,1.0);	
	_order = o;	
	_t.resize(o + 1);	
	_t[1] = sp.arc2t(i, 0.5 * _arc);
	_t[o] = 1.0;
	_r.resize(nqd);
	_dr.resize(nqd);	
	_z.resize(nqd);
	_dz.resize(nqd);	
	_J.resize(nqd);
	_xi.resize(nqd);
	_basis.resize(o + 1);
	for (auto & member : _basis) { 
		member.resize(nqd); 
	};
	
	for (int k = 0; k < nqd; k++) {
		double tt = Numeric::qd_GL_x[nqd][k];
		const Eigen::Vector3d d_r = sp.d(sp.x(), i, tt);
		const Eigen::Vector3d d_z = sp.d(sp.y(), i, tt);
		_r[k] = d_r(0);
		_z[k] = d_z(0);	
		_dr[k] = d_r(1);
		_dz[k] = d_z(1);		
		_J[k] = sqrt(d_r(1) * d_r(1) + d_z(1) * d_z(1));
		_xi[k] = sp.localArc(i, tt) / _arc;
		for (int j = 0; j <= o; j++) {
			_basis[j][k] = Numeric::N[o][j](_xi[k]);
		}
	}

}

Element::Element(const Element &e) {	
	_arc = e.arc();
	_t = e.t();	
	_order = e.order();
};

void Element::init(const Spline &sp, int i, int nqd, const std::vector<double> & qdx) {
	int o = order(); 
	_r.resize(nqd);
	_dr.resize(nqd);
	_z.resize(nqd);
	_dz.resize(nqd);
	_J.resize(nqd);
	_xi.resize(nqd);
	_basis.resize(o + 1);
	for (auto & member : _basis) {
		member.resize(nqd);
	};

	for (int k = 0; k < nqd; k++) {
		double tt = qdx[k];
		const Eigen::Vector3d d_r = sp.d(sp.x(), i, tt);
		const Eigen::Vector3d d_z = sp.d(sp.y(), i, tt);
		_r[k] = d_r(0);
		_z[k] = d_z(0);
		_dr[k] = d_r(1);
		_dz[k] = d_z(1);
		_J[k] = sqrt(d_r(1) * d_r(1) + d_z(1) * d_z(1));
		_xi[k] = sp.localArc(i, tt) / _arc;
		for (int j = 0; j <= o; j++) {
			_basis[j][k] = Numeric::N[o][j](_xi[k]);
		}
	}

}


void Element::getGLQuad(int nqd, double const * qdx, double const * qdw) {
	qdx = Numeric::qd_GL_x[nqd];
	qdw = Numeric::qd_GL_w[nqd];
	//return Eigen::Map<const Eigen::Matrix2Xd>(Numeric::qd[nqd], 2, nqd);
};

void Element::getGLQuad(int nqd, double const * qdx) {
	qdx = Numeric::qd_GL_x[nqd];
	//qdw = Numeric::qd_GL_w[nqd];
	//return Eigen::Map<const Eigen::Matrix2Xd>(Numeric::qd[nqd], 2, nqd);
};

void Element::getLOGQuad(int nqd, double const * qdx, double const * qdw) {
	qdx = Numeric::qd_LOG_x[nqd];
	qdw = Numeric::qd_LOG_w[nqd];
	//return Eigen::Map<const Eigen::Matrix2Xd>(Numeric::qd[nqd], 2, nqd);
};

void Element::getLOGQuad(int nqd, double const * qdx) {
	qdx = Numeric::qd_LOG_x[nqd];
	//qdw = Numeric::qd_LOG_w[nqd];
	//return Eigen::Map<const Eigen::Matrix2Xd>(Numeric::qd[nqd], 2, nqd);
};

