#include "stdafx.h"
#include "taylorCone.h"

void TaylorCone::init(double c1, double b0){
	computeCoefabc(c1, b0);




}



double TaylorCone::fD2(int order, double h0, double h1, double y0, double y1, double y2, int location) {
	switch (order)	{
	case 1: {
		switch (location) {
		case 0: {
			return (-1. / h0 - 1. / (h0 + h1))*y0 + (1. / h0 + 1 / h1)*y1 + (-1. / h1 + 1 / (h0 + h1))*y2;
			break;
		}
		case 1: {
			return (-1. / h0 + 1. / (h0 + h1))*y0 + (1. / h0 - 1. / h1)*y1 + (1. / h1 - 1. / (h0 + h1))*y2;
			break;
		}
		case 2: {
			return (1. / h0 - 1. / (h0 + h1))*y0 + (-1. / h0 - 1. / h1)*y1 + (1. / h1 + 1. / (h0 + h1))*y2;
			break;
		}
		}
		break;
	}
	case 2: {
		return (2.*y0) / (h0*(h0 + h1)) - (2.*y1) / (h0*h1) + (2. / (h0*h1) - 2. / (h0*(h0 + h1)))*y2;
		break;
	}
	default:
		return 0;
		break;
	}
};

double TaylorCone::fD3(int order, double h0, double h1, double h2, double y0, double y1, double y2, double y3, int location) {
	switch (order) {
	case 1: {
		switch (location) {
		case 0: {	
			return ((-1.*(2.*h0 + h1)) / (h0*(h0 + h1)) - 1. / (h0 + h1 + h2))*y0
				+ ((h0 + h1) / (h0*h1) + (h0 + h1) / (h1*(h1 + h2)))*y1
				+ ((-1.*h0) / (h1*(h0 + h1)) - (1.*h0) / (h1*h2))*y2
				+ (h0 / (h1*h2) + (-1.*h0 - 1.*h1) / (h1*(h1 + h2)) + 1 / (h0 + h1 + h2))*y3;
			break;
		}
		case 3: {	
			return (1 / (h0 + h1) - (1.*h2) / (h0*(h0 + h1)) - 1. / (h0 + h1 + h2))*y0
				+ (1 / h1 + h2 / (h0*h1) - 1. / (h1 + h2))*y1
				+ ((-1.*(h0 + 2.*h1)) / (h1*(h0 + h1)) - 1. / h2 - (1.*h2) / (h1*(h0 + h1)))*y2
				+ (1 / h2 + 1 / (h1 + h2) + 1 / (h0 + h1 + h2))*y3;
			break;
		}		
		}
		break;
	}
	case 2: {	
		switch (location) {
		case 0: {
			return (2. / (h0*(h0 + h1)) + (2.*(2.*h0 + h1)) / (h0*(h0 + h1)*(h0 + h1 + h2)))*y0
				+ (-2. / (h0*h1) - (2.*(2.*h0 + h1)) / (h0*h1*(h1 + h2)))*y1
				+ (2. / (h1*(h0 + h1)) + (2.*(2.*h0 + h1)) / (h1*(h0 + h1)*h2))*y2
				+ ((-2.*(2.*h0 + h1)) / (h1*(h0 + h1)*h2) + (2.*(2.*h0 + h1)) / (h0*h1*(h1 + h2)) - (2.*(2.*h0 + h1)) / (h0*(h0 + h1)*(h0 + h1 + h2)))*y3;
			
			break;
		}
		case 3: {
			return (-4. / (h0*(h0 + h1)) + (2.*(2.*h0 + h1)) / (h0*(h0 + h1)*(h0 + h1 + h2)))*y0
				+ (4. / (h0*h1) + (2.*(h0 - 1.*h1)) / (h0*h1*(h1 + h2)))*y1
				+ (-4. / (h1*(h0 + h1)) - (2.*(h0 + 2.*h1)) / (h1*(h0 + h1)*h2))*y2
				+ ((2.*(h0 + 2.*h1)) / (h1*(h0 + h1)*h2) - (2.*(h0 - 1.*h1)) / (h0*h1*(h1 + h2)) - (2.*(2.*h0 + h1)) / (h0*(h0 + h1)*(h0 + h1 + h2)))*y3;
			break;
		}
		}		
		return 0;
		break;
	}
	default:
		
		return 0;

		break;
	}	
	return 0;
};
double TaylorCone::fD4(int order, double h0, double h1, double h2, double h3, double y0, double y1, double y2, double y3, double y4) {
	switch (order)
	{
	case 1:
		return (h3*(h2 + h3)*(h1 + h2 + h3)*y0
			) / (h0*(h0 + h1)*(h0 + h1 + h2)*(h0 + h1 + h2 + h3)) - (1.*h3*(h2 + h3)*(h0 + h1 + h2 + h3)*y1
				) / (h0*h1*(h1 + h2)*(h1 + h2 + h3)) + (h3*(h1 + h2 + h3)*(h0 + h1 + h2 + h3)*y2
					) / (h1*(h0 + h1)*h2*(h2 + h3)) - (1.*(h2 + h3)*(h1 + h2 + h3)*(h0 + h1 + h2 + h3)*y3
						) / (h2*(h1 + h2)*(h0 + h1 + h2)*h3) + ((pow(h1, 2)*(h2 + 2.*h3) + pow(h2 + h3, 2)*(h2 + 4.*h3) + 2.*h1*(pow(h2, 2) + 4.*h2*h3 + 3.*pow(h3, 2)) + h0 * (pow(h2, 2) + 4.*h2*h3 + 3.*pow(h3, 2) + h1 * (h2 + 2.*h3)))*y4) / (h3*(h2 + h3)*(h1 + h2 + h3)*(h0 + h1 + h2 + h3));
	case 2:
		return (2.*(pow(h2, 2) + 4.*h2*h3 + 3.*pow(h3, 2) + h1 * (h2 + 2.*h3))*y0
			) / (h0*(h0 + h1)*(h0 + h1 + h2)*(h0 + h1 + h2 + h3)) - (2.*(pow(h2, 2) + 4.*h2*h3 + 3.*pow(h3, 2) + h0 * (h2 + 2.*h3) + h1 * (h2 + 2.*h3))*y1
				) / (h0*h1*(h1 + h2)*(h1 + h2 + h3)) + (2.*(pow(h1, 2) + pow(h2, 2) + 4.*h2*h3 + 3.*pow(h3, 2) + 2.*h1*(h2 + 2.*h3) + h0 * (h1 + h2 + 2.*h3))*y2
					) / (h1*(h0 + h1)*h2*(h2 + h3)) - (2.*(pow(h1, 2) + 4.*h1*(h2 + h3) + 3.*pow(h2 + h3, 2) + h0 * (h1 + 2.*(h2 + h3)))*y3
						) / (h2*(h1 + h2)*(h0 + h1 + h2)*h3) + (2.*(pow(h1, 2) + 4.*h1*h2 + 3.*pow(h2, 2) + 6.*h1*h3 + 9.*h2*h3 + 6.*pow(h3, 2) + h0 * (h1 + 2.*h2 + 3.*h3))*y4) / (h3*(h2 + h3)*(h1 + h2 + h3)*(h0 + h1 + h2 + h3));
	default:
		break;
	}
	
};
void TaylorCone::circleDerivativeBegin(const Eigen::MatrixX2d &xy, double &dx, double &ddx, double &dy, double &ddy) {	
	double x0 = xy(0, 0), y0 = xy(0, 1);
	double x1 = xy(1, 0), y1 = xy(1, 1);
	double x2 = xy(2, 0), y2 = xy(2, 1);
	double x3 = xy(3, 0), y3 = xy(3, 1);
	double h0 = sqrt((x1 - x0)*(x1 - x0) + (y1 - y0)*(y1 - y0));
	double h1 = sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1));
	double h2 = sqrt((x3 - x2)*(x3 - x2) + (y3 - y2)*(y3 - y2));

	double curvature = -2.0 / sqrt(x0 * x0 + y0 * y0);
	double slope = abs(x0) / abs(y0);
	if (xy(xy.rows() - 1, 1) > 0) {	curvature = curvature * -1.;}
	
	dx = fD3(1, h0, h1, h2, x0, x1, x2, x3, 0);
	ddx = fD3(2, h0, h1, h2, x0, x1, x2, x3, 0);	
	dy = dx * slope;
	ddy = (ddx * dy + curvature * pow(dx * dx + dy * dy,1.5)) / dx - dy * (dx * dx + dy * dy) / dx / x0;
};

void TaylorCone::coneDerivativeEnd(const Eigen::MatrixX2d &xy, double c[5], double &dx, double &ddx, double &dy, double &ddy) {
	const int n = xy.rows();
	double xn1 = xy(n - 5, 0), yn1 = xy(n - 5, 1);
	double x0 = xy(n - 4, 0), y0 = xy(n - 4, 1);
	double x1 = xy(n - 3, 0), y1 = xy(n - 3, 1);
	double x2 = xy(n - 2, 0), y2 = xy(n - 2, 1);
	double x3 = xy(n - 1, 0), y3 = xy(n - 1, 1);


	double hn1 = sqrt((x0 - xn1)*(x0 - xn1) + (y0 - yn1)*(y0 - yn1));
	double h0 = sqrt((x1 - x0)*(x1 - x0) + (y1 - y0)*(y1 - y0));
	double h1 = sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1));
	double h2 = sqrt((x3 - x2)*(x3 - x2) + (y3 - y2)*(y3 - y2));

	double r = x3;

	double curvature = (0.75*pow(r, 4.5)*c[1] + 15.75*pow(r, 1.5)*c[3] + pow(r, 6)*(1. + pow(c[0] + (-0.5*pow(r, 4.5)*c[1] - 3.5*pow(r, 1.5)*c[3] - 5.*c[4]) / pow(r, 6), 2))*(c[0] + (-0.5*pow(r, 4.5)*c[1] - 3.5*pow(r, 1.5)*c[3] - 5.*c[4]) / pow(r, 6)) + 30.*c[4]) / (pow(r, 7)*pow(1. + pow(c[0] + (-0.5*pow(r, 4.5)*c[1] - 3.5*pow(r, 1.5)*c[3] - 5.*c[4]) / pow(r, 6), 2), 1.5));
	double slope = c[0] - c[1] / (2.*pow(r, 1.5)) - (7. * c[3]) / (2.*pow(r, 4.5)) - (5. * c[4]) / pow(r, 6);
	
	//dx  = fD3(1, h0, h1, h2, x0, x1, x2, x3, 3);
	dx  = fD4(1, hn1, h0, h1, h2, xn1, x0, x1, x2, x3);
	//ddx = fD3(2, h0, h1, h2, x0, x1, x2, x3, 3);
	ddx = fD4(2, hn1, h0, h1, h2, xn1, x0, x1, x2, x3);
	dy = dx * slope;
	ddy = (ddx * dy + curvature * pow(dx * dx + dy * dy, 1.5)) / dx - dy * (dx * dx + dy * dy) / dx / x3;
}



void TaylorCone::computeCoefabc(double c1, double b0) {
	a[0] = 2.7142175397111330 * c1;
	b[0] = b0;
	c[0] = -0.8604366861256783;

	a[1] = 0.1442586135181731 * a[0] * a[0] - 0.4749808176761397 * b0 * b0 - c[0];
	b[1] = 0;
	c[1] = c1;

	c[2] = 0.0;
}

Eigen::MatrixX2d TaylorCone::generateCircle(double angle0, double angle1, double radius, int n) {
	Eigen::MatrixX2d xy(n, 2);
	xy.setZero();
	for (int i = 0; i < n; i++) {
		double t = ((double)i) / (n - 1.);
		double theta = angle0 + (angle1 - angle0) * t;// pow(t, 0.8);
		xy(i, 0) = radius * sin(theta) ;
		xy(i, 1) = radius * cos(theta) ;
	}
	if (abs(angle0) < 1e-12) { xy(0, 0) = 0.0; }
	if (abs(angle1 - M_PI) < 1e-12) { xy(n - 1, 0) = 0.0; }

	return xy;
}

Eigen::MatrixX2d TaylorCone::generateCircle(double r0, double z0, int end, int n) {
	Eigen::MatrixX2d xy(n, 2);
	xy.setZero();
	double radius = sqrt(r0* r0 + z0 * z0);
	double angle0 = acos(z0 / radius);
	double angle1 = 0.;
	
	xy(0, 0) = r0;
	xy(0, 1) = z0;
	if (end != 0) { angle1 = M_PI; };	
	
	for (int i = 1; i < n; i++) {
		double t = ((double)i) / (n - 1.);
		double theta = angle0 + (angle1 - angle0) * t;// pow(t, 0.8);
		xy(i, 0) = radius * sin(theta);
		xy(i, 1) = radius * cos(theta);
	}	
	xy(n - 1, 0) = 0;
	if (end == 0) { xy(n - 1, 1) = radius; }
	else { xy(n - 1, 1) = -radius; }
	return xy;
}

Eigen::MatrixX2d TaylorCone::generateCone(double rc, double rstar, double c[5], int n) {
	Eigen::MatrixX2d xy(n, 2);
	for (int i = 0; i < n; i++) {
		double t = ((double)i) / (n - 1.);
		double r = rstar * t;
		xy(i, 0) = r;
		xy(i, 1) = c3Cone(r,rc,c);
		
	}
	
	return xy;
};

double TaylorCone::c3Cone(double r, double rc, double c[5]) {
	double f0, f2, f4, f6;
	f0 = (231 * c[4]) / (16.*pow(rc, 5)) + (1045 * c[3]) / (128.*pow(rc, 3.5)) + (195 * c[1]) / (128.*sqrt(rc)) + (5 * c[0] * rc) / 16.;
	f2 = (-495 * c[4]) / (16.*pow(rc, 7)) - (1995 * c[3]) / (128.*pow(rc, 5.5)) - (117 * c[1]) / (128.*pow(rc, 2.5)) + (15 * c[0]) / (16.*rc);
	f4 = (385 * c[4]) / (16.*pow(rc, 9)) + (1463 * c[3]) / (128.*pow(rc, 7.5)) + (65 * c[1]) / (128.*pow(rc, 4.5)) - (5 * c[0]) / (16.*pow(rc, 3));
	f6 = (-105 * c[4]) / (16.*pow(rc, 11)) - (385 * c[3]) / (128.*pow(rc, 9.5)) - (15 * c[1]) / (128.*pow(rc, 6.5)) + c[0] / (16.*pow(rc, 5));
	if (r < rc) {
		return f0 + f2 * r * r + f4 * r * r * r * r + f6 * r * r * r * r * r * r;
	}
	else {
		return c[0] * r + c[1] / sqrt(r) + c[3] / pow(r, 3.5) + c[4] / pow(r, 5.0);
	}
};