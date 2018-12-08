#pragma once
#ifndef ELEMENT_H
#define ELEMENT_H
#include "spline.h"
//#include <Eigen\Core>
#include <vector>

class Element {
		double _arc = 0;
		int _order = 2;
		//int _nqd = 6;
		std::vector<double> _t;
		std::vector<double> _r, _dr,  _z, _dz, _J, _xi;
		std::vector<std::vector<double>> _basis;		
	public:
		Element() {};
		Element(const Element &e);		
		void init(const Spline &sp, int i, int o, int nqd);
		void init(const Spline &sp, int i, int nqd, double const * qdx);

		const double &arc() const { return _arc; };
		const int &order() const { return _order; };
		//const int &nqd() const { return _nqd; };
		const std::vector<double> &t() const { return _t; };
		const std::vector<double> &r() const { return _r; };
		const std::vector<double> &dr() const { return _dr; };		
		const std::vector<double> &z() const { return _z; };
		const std::vector<double> &dz() const { return _dz; };
		
		const std::vector<double> &J() const { return _J; };
		const std::vector<double> &xi() const { return _xi; };
		const std::vector<std::vector<double>> &basis() const { return _basis; };
		static void getGLQuad(int nqd, double const * qdx, double const * qdw);
		static void getLOGQuad(int nqd, double  const * qdx,  double const * qdw);
		static void getGLQuad(int nqd, double const * qdx);
		static void getLOGQuad(int nqd, double const * qdx);

	};



#endif