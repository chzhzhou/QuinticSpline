#pragma once
#ifndef TAYLORCONE_H
#define TAYLORCONE_H
#include "bem.h"
#include <Eigen/Dense>

class TaylorCone {

public:
	int _nSp0, _nSp1, _nSp2;
	Eigen::MatrixX2d _xy0, _xy1, _xy2;
	Bem bem0, bem1, bem2;
	
	double a[5];
	double b[5];
	double c[5];

	TaylorCone(double c1, double b0) { computeCoefabc(c1, b0); };
	void init(double c1, double b0);
	void computeCoefabc(double c1, double b0);
	void prepareBem(int type, const Eigen::MatrixX2d &xy, int shift, Bem &bem);
	void setFluidBC(const Bem &bemCone, const Bem &bemPatch, Eigen::VectorXd &fluidBC) const;
	void computeResidue(const Bem &bemCone, const Eigen::VectorXd &phi,   Eigen::VectorXd &residue) const;


public:
	static Eigen::MatrixX2d generateCircle(double angle0, double angle1, double radius, int n);
	static Eigen::MatrixX2d generateCircle(double r0, double z0, int end, int n);
	static Eigen::MatrixX2d generateCone(double rc, double rstar, double c[5], int n);


	static double c3Cone(double r, double rc, double c[5]);
	static double fD2(int order, double h0, double h1, double y0, double y1, double y2, int location);
	static double fD3(int order, double h0, double h1, double h2, double y0, double y1, double y2, double y3, int location);
	static double fD4(int order, double h0, double h1, double h2, double h3, double y0, double y1, double y2, double y3, double y4);
	static void circleDerivativeBegin(const Eigen::MatrixX2d &xy,	double &dx, double &ddx, double &dy, double &ddy);
	static void coneDerivativeEnd(const Eigen::MatrixX2d &xy, double c[5], double &dx, double &ddx, double &dy, double &ddy);
	static double harmonicGrow(double r, double z, int l, int divide = 0);
	static double harmonicDecay(double r, double z, int l, int divide = 0);
	static double velocityPotentialFarField(double r, double z, const double (&a)[5]);

	static void SD2LR(const Eigen::MatrixXd &S, const Eigen::MatrixXd &D, int nSwap, Eigen::MatrixXd &L, Eigen::MatrixXd &R);

	


};

#endif