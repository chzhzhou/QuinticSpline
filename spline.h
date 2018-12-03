#pragma once
#ifndef SPLINE_H
#define SPLINE_H
#include <Eigen/Core>

class Spline {
	typedef Eigen::Matrix<double, Eigen::Dynamic, 6> Coef;	
public:
	enum class BC { Even, Odd, Mix };
	Spline& operator = (const Spline &sp);
	Spline() {};
	Spline(const Spline &sp) { *this = sp; }
	//getter
	const Coef &x() const { return _x; };
	const Coef &y() const { return _y; };	
	const Eigen::VectorXd &h() const { return _h; }; //get rid of later
	const Eigen::MatrixX2d &node() const { return _node; };
	//setter
	void node(const Eigen::MatrixX2d &xy);
	void x(BC bc0, BC bc1, double a0 = 0, double b0 = 0, double a1 = 0, double b1 = 0);
	void y(BC bc0, BC bc1, double a0 = 0, double b0 = 0, double a1 = 0, double b1 = 0);
	//API
	double localArc(int i, double t = 1.0, int nqd = 10) const;
	double arc2t(int i, double arc, int nqd = 10) const;		
private:
	Coef _x, _y;	
	Eigen::VectorXd _h;
	Eigen::MatrixX2d _node;

	void h(const Eigen::MatrixX2d &node);
	void setComponent(int i, Coef &x);
	void computeCoef(Coef &x,BC bc0, BC bc1, double a0 = 0, double b0 = 0, double a1 = 0, double b1 = 0);	
	Eigen::Vector3d d(const Coef &x, int i, double t) const;
};





#endif