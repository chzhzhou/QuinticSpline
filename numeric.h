#pragma once
#ifndef NUMERIC_H
#define NUMERIC_H


class Numeric {
public:
	static const double* qd[21];
	static const double* lqd[21];
	static void KEPQ(double m, double &K, double &E, double &KP, double &KQ, double &EP, double &EQ);
	static void KEPQ(double m, double &K, double &E);
	static double(**N[3])(double);
	template<class T>
	static constexpr const T& clamp(const T& v, const T& lo, const T& hi) {
		return v<= lo ? lo : v>= hi ? hi : v;
	}
private:
	static const double qd1[2], qd2[4], qd3[6], qd4[8], qd5[10],
		qd6[12], qd7[14], qd8[16], qd9[18], qd10[20],
		qd11[22], qd12[24], qd13[26], qd14[28], qd15[30],
		qd16[32], qd17[34], qd18[36], qd19[38], qd20[40];
	//static const double lqd1[10], lqd2[20], lqd3[30], lqd4[40];
	static const double lqd1[2], lqd2[4], lqd3[6], lqd4[8], lqd5[10],
		lqd6[12], lqd7[14], lqd8[16], lqd9[18], lqd10[20],
		lqd11[22], lqd12[24], lqd13[26], lqd14[28], lqd15[30],
		lqd16[32], lqd17[34], lqd18[36], lqd19[38], lqd20[40];
	static const double ellipKP[11], ellipKQ[11], ellipEP[11], ellipEQ[11];
	static double N00(double x);
	static double N10(double x);
	static double N11(double x);
	static double N20(double x);
	static double N21(double x);
	static double N22(double x);
	static double(*N0[1])(double);
	static double(*N1[2])(double);
	static double(*N2[3])(double);
};

#endif