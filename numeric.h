#pragma once
#ifndef NUMERIC_H
#define NUMERIC_H
#ifndef M_PI
#define M_PI   3.14159265358979323846264338327950288
#endif // !M_PI

class Numeric {
public:
	static const double* qd_GL_x[21];
	static const double* qd_GL_w[21];
	static const double* qd_LOG_x[21];
	static const double* qd_LOG_w[21];
	static void KEPQ(double m, double &K, double &E, double &KP, double &KQ, double &EP, double &EQ);
	static void KEPQ(double m, double &K, double &E);
	static double(**N[3])(double);
	static double legendreP(int l, double x);
	static double legendreP(int l, int d, double x);

	template<class T>
	static constexpr const T& clamp(const T& v, const T& lo, const T& hi) {
		return v<= lo ? lo : v>= hi ? hi : v;
	}
private:
	static const double qd_GL_x01[1], qd_GL_x02[2], qd_GL_x03[3], qd_GL_x04[4], qd_GL_x05[5], qd_GL_x06[6], qd_GL_x07[7], qd_GL_x08[8], qd_GL_x09[9], qd_GL_x10[10], qd_GL_x11[11], qd_GL_x12[12], qd_GL_x13[13], qd_GL_x14[14], qd_GL_x15[15], qd_GL_x16[16], qd_GL_x17[17], qd_GL_x18[18], qd_GL_x19[19], qd_GL_x20[20];
	static const double qd_GL_w01[1], qd_GL_w02[2], qd_GL_w03[3], qd_GL_w04[4], qd_GL_w05[5], qd_GL_w06[6], qd_GL_w07[7], qd_GL_w08[8], qd_GL_w09[9], qd_GL_w10[10], qd_GL_w11[11], qd_GL_w12[12], qd_GL_w13[13], qd_GL_w14[14], qd_GL_w15[15], qd_GL_w16[16], qd_GL_w17[17], qd_GL_w18[18], qd_GL_w19[19], qd_GL_w20[20];
	static const double qd_LOG_x01[1], qd_LOG_x02[2], qd_LOG_x03[3], qd_LOG_x04[4], qd_LOG_x05[5], qd_LOG_x06[6], qd_LOG_x07[7], qd_LOG_x08[8], qd_LOG_x09[9], qd_LOG_x10[10], qd_LOG_x11[11], qd_LOG_x12[12], qd_LOG_x13[13], qd_LOG_x14[14], qd_LOG_x15[15], qd_LOG_x16[16], qd_LOG_x17[17], qd_LOG_x18[18], qd_LOG_x19[19], qd_LOG_x20[20];
	static const double qd_LOG_w01[1], qd_LOG_w02[2], qd_LOG_w03[3], qd_LOG_w04[4], qd_LOG_w05[5], qd_LOG_w06[6], qd_LOG_w07[7], qd_LOG_w08[8], qd_LOG_w09[9], qd_LOG_w10[10], qd_LOG_w11[11], qd_LOG_w12[12], qd_LOG_w13[13], qd_LOG_w14[14], qd_LOG_w15[15], qd_LOG_w16[16], qd_LOG_w17[17], qd_LOG_w18[18], qd_LOG_w19[19], qd_LOG_w20[20];

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