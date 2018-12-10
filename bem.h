#pragma once
#ifndef BEM_H
#define BEM_H
#ifndef M_PI
#define M_PI   3.14159265358979323846264338327950288
#endif // !M_PI
#include <Eigen/Core>
#include "spline.h"
#include "element.h"
#include <vector>

class Bem {
private:
	typedef struct Node {
		Eigen::MatrixX3d r, z;
	} Node;
public:
	Bem() {};	
	void initialize(const Eigen::MatrixX2d &xy);
	const std::vector<Element> &e() const { return _e; };
	const Spline &sp() const { return _sp; };
	//const Eigen::MatrixX3d &r() const { return _r; };
	//const Eigen::MatrixX3d &z() const { return _z; };
	const Node &node() const { return _node; };
	const std::vector<double > regular (double rp, double zp, int idElement) const;
	const std::vector<double > axis (double zp, int idElement) const;
	const std::vector<double > singular (double tau, int idElement) const;
	static void assembly(const Bem &source, const Bem &reciever,  Eigen::MatrixXd &S, Eigen::MatrixXd &D);	



private:
	// the real data
	std::vector<Element> _e;	
	Node  _node;
	Spline _sp;
	// under the hood
	static const double eps;	
	void setrz();
	static void abm(double rp, double zp, double r, double z,	double &a, double &b, double &m	);
	static double RKE(double P, double Q, double m, double t, double tau);
	static void fKE(
		double rp, double zp, double r, double z, double dr, double dz, double J, double a, double b,
		double &f_single_K, double &f_double_K, double &f_double_E
	);
	static int checkBounbdaryPosition(const Bem &bem0, const Bem &bem1);


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
		const int &nElm() const{ return _nElm; };
		const int &order() const { return _order; };
		const int &qdOrder() const { return _qdOrder; };
		const int &indexShift() const { return _indexShift; };
		void indexShift(int i) { _indexShift = i; };			
		void qdOrder(int i) { _qdOrder = i; };
		void order(int i) { _order = i; };
		void nElm(int i) { _nElm = i; };
		void print() const;
		static void tic();
		static void toc(); 
	private:
		int _nElm = 0;
		//default
		int _order = 2;
		int _indexShift = 0;
		int _qdOrder = 6;
		static double _timer;
	};
	Properties settings;
};

#endif

