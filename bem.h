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
public:
	Bem() {};
	
	void initialize(const Eigen::MatrixX2d &xy);
	const std::vector<Element> &e() const { return _e; };
	const Spline &sp() const { return _sp; };
	Eigen::VectorXd regular(double rp, double zp, int idElement);	
	Eigen::VectorXd singular(double tau, int idElement);
	void assembly(Eigen::MatrixXd &S, Eigen::MatrixXd &D);
	
	
	
private:
	std::vector<Element> _e;
	Spline _sp;	
	static const double eps;
	void RKRE(const Eigen::VectorXd &m, const Eigen::VectorXd &t, double tau, Eigen::VectorXd &RK, Eigen::VectorXd &RE);
	void sKdKdE(double rp, double zp, const Element &e,	Eigen::VectorXd &sK, Eigen::VectorXd &dK, Eigen::VectorXd &dE, Eigen::VectorXd &K,	Eigen::VectorXd &E);
	void sKdKdE(double rp, double zp, const Element &e,	Eigen::VectorXd &sK,Eigen::VectorXd &dK,	Eigen::VectorXd &dE,Eigen::VectorXd &m);
	void KEPKQKPEQE(const Eigen::VectorXd &m, Eigen::VectorXd &K , Eigen::VectorXd &E,	Eigen::VectorXd &PK, Eigen::VectorXd &QK,	Eigen::VectorXd &PE, Eigen::VectorXd &QE);
	void KEPKQKPEQE(const Eigen::VectorXd &m, Eigen::VectorXd &K, Eigen::VectorXd &E);
	

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
		void qdOrder(int i) { _qdOrder = i; };
		void nElm(int i) { _nElm = i; };
		void print() const;
		static void tic() { printf("tic\n"); _timer = omp_get_wtime(); };
		static void toc() { printf("toc ... %.ef sec\n", omp_get_wtime() - _timer); };
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

