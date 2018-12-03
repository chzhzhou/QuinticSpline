#pragma once
#ifndef BEM_H
#define BEM_H
#ifndef M_PI
#define M_PI   3.14159265358979323846264338327950288
#endif // !M_PI
#include <Eigen/Core>
#include "spline.h"
#include <vector>

class Bem {	
	class Element {
		double _arc = 0;
		Eigen::VectorXd _t;
		Eigen::VectorXd _r, _dr, _ddr, _z, _dz, _ddz, _J, _xi;
		Eigen::MatrixXd _basis;
	public:
		void init(const Spline &sp, int i, int o, int nqd);
		const double &arc() { return _arc; };		
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
		static const Eigen::Matrix2Xd getGLQuad (int nqd);
		static const Eigen::Matrix2Xd getLogQuad(int nqd);

	};


public:
	Bem() {};
	std::vector<Element> e;
	void initialize(const Eigen::MatrixX2d &xy);
	const Spline &sp() { return _sp; };
	Eigen::VectorXd regular(double rp, double zp, const Element &e);


private:
	Spline _sp;
	//static const Eigen::VectorXd
	void KEPKQKPEQE(
		const Eigen::VectorXd &m,
		Eigen::VectorXd K,
		Eigen::VectorXd E,
		Eigen::VectorXd PK,
		Eigen::VectorXd QK,
		Eigen::VectorXd PE,
		Eigen::VectorXd QE);
	void KEPKQKPEQE(
		const Eigen::VectorXd &m,
		Eigen::VectorXd K,
		Eigen::VectorXd E);

public:	
	struct Properties {		
		struct bc {
			struct _bc {
				Spline::BC type = Spline::BC::Odd;
				double a = 0.0, b = 0.0;
				_bc() {};				
				void set(Spline::BC tmp, double aa, double bb) { type = tmp; a = aa; b = bb; };				
				void print() const;
			};
			_bc begin, end;			
		};		
		Properties() { 
			xBC.begin.set(Spline::BC::Odd,0,0); 
			xBC.end.set(Spline::BC::Odd, 0, 0);
			yBC.begin.set(Spline::BC::Even, 0, 0);
			yBC.end.set(Spline::BC::Even, 0, 0);
		};
		bc xBC , yBC;
		const int &nElm() { return _nElm; };
		const int &order() { return _order; };
		const int &qdOrder() { return _qdOrder; };
		const int &indexShift() { return _indexShift; };
		void indexShift(int i) { _indexShift = i; };			
		void nElm(int i) { _nElm = i; };
		void print() const;
	private:
		int _nElm = 0;
		//default
		int _order = 2;
		int _indexShift = 0;
		int _qdOrder = 6;
	};
	Properties settings;
};

#endif

