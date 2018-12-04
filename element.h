#pragma once
#ifndef ELEMENT_H
#define ELEMENT_H
#include "spline.h"
#include <Eigen\Core>

class Element {
		double _arc = 0;
		int _order = 2;
		Eigen::VectorXd _t;
		Eigen::VectorXd _r, _dr, _ddr, _z, _dz, _ddz, _J, _xi;
		Eigen::MatrixXd _basis;
	public:
		Element() {};
		Element(const Element &e);
		
		void init(const Spline &sp, int i, int o, int nqd);
		void init(const Spline &sp, int i, const Eigen::Matrix2Xd &qd);


		const double &arc() const { return _arc; };
		const int &order() const { return _order; };
		const Eigen::VectorXd &t() const { return _t; };
		const Eigen::VectorXd &r() const { return _r; };
		const Eigen::VectorXd &dr() const { return _dr; };
		const Eigen::VectorXd &ddr() const { return _ddr; };
		const Eigen::VectorXd &z() const { return _z; };
		const Eigen::VectorXd &dz() const { return _dz; };
		const Eigen::VectorXd &ddz() const { return _ddz; };
		const Eigen::VectorXd &J() const { return _J; };
		const Eigen::VectorXd &xi() const { return _xi; };
		const Eigen::MatrixXd &basis() const { return _basis; };
		static const Eigen::Matrix2Xd getGLQuad(int nqd);
		static const Eigen::Matrix2Xd getLogQuad(int nqd);
	};



#endif