#pragma once
#include <assert.h>
#include <cmath>
#include <cstdint>
#include <climits>
#include <memory>

// =========== others ============

inline double Lerp(double t, double x, double y) { return (1 - t) * x + t * y; }
inline double Clamp(double x, double a, double b) { if (x < a) return a; if (x > b) return b; return x; }
inline double PowerHeuristic(int nf, double fPdf, int ng, double gPdf) {
	double f = nf * fPdf, g = ng * gPdf;
	return (f * f) / (f * f + g * g);
}

#define MachineEpsilon (std::numeric_limits<double>::epsilon() * 0.5)
inline double gamma(int n) {
	return (n * MachineEpsilon) / (1 - n * MachineEpsilon);
}

inline uint64_t DoubleToBits(double f) {
	uint64_t ui;
	memcpy(&ui, &f, sizeof(double));
	return ui;
}

inline double BitsToDouble(uint64_t ui) {
	double f;
	memcpy(&f, &ui, sizeof(uint64_t));
	return f;
}

inline double NextDoubleUp(double v) {
	if (std::isinf(v) && v > 0.0) return v;
	if (v == -0.f) v = 0.f;
	uint64_t ui = DoubleToBits(v);
	if (v >= 0) ++ui;
	else --ui;
	return BitsToDouble(ui);
}

inline double NextDoubleDown(double v) {
	if (std::isinf(v) && v < 0.0) return v;
	if (v == 0.f) v = -0.f;
	uint64_t ui = DoubleToBits(v);
	if (v > 0) --ui;
	else ++ui;
	return BitsToDouble(ui);
}

inline double AddRoundUp(double a, double b) { return NextDoubleUp(a + b); }
inline double AddRoundDown(double a, double b) { return NextDoubleDown(a + b); }
inline double SubRoundUp(double a, double b) { return AddRoundUp(a, -b); }
inline double SubRoundDown(double a, double b) { return AddRoundDown(a, -b); }
inline double MulRoundUp(double a, double b) { return NextDoubleUp(a * b); }
inline double MulRoundDown(double a, double b) { return NextDoubleDown(a * b); }
inline double DivRoundUp(double a, double b) { return NextDoubleUp(a / b); }
inline double DivRoundDown(double a, double b) { return NextDoubleDown(a / b); }

// =========== complex ============

class complex {
public:
	double re, im;
	complex() : re(0), im(0) {}
	complex(double _re) : re(_re), im(0) {}
	complex(double _re, double _im) : re(_re), im(_im) {}

	double norm() const { return re * re + im * im; }
	complex operator -() const { return { -re, -im }; }
	complex operator +(complex x) const { return { re + x.re, im + x.im }; }
	complex operator -(complex x) const { return { re - x.re, im - x.im }; }
	complex operator *(complex x) const { return { re * x.re - im * x.im, re * x.im + im * x.re }; }
	complex operator /(complex x) const {
		double scale = 1.0 / x.norm();
		return { scale * (re * x.re + im * x.im), scale * (im * x.re - re * x.im) };
	}

	complex sqrt() const {
		double n = std::sqrt(norm());
		if (n == 0) return { 0 };
		double t1 = std::sqrt(0.5 * (n + std::abs(re)));
		double t2 = 0.5 * im / t1;
		if (re >= 0) return { t1, t2 };
		return { std::abs(t2), im > 0 ? t1 : -t1 };
	}
};

// =========== vector ============

const double pi = 3.141592653589793238463;

class vec2
{
public:
	double a[2];

	vec2(double x = 0) { a[0] = a[1] = x; }
	vec2(double x, double y) { a[0] = x; a[1] = y; }

	double x() { return a[0]; }
	const double x() const { return a[0]; }
	double y() { return a[1]; }
	const double y() const { return a[1]; }

	vec2 operator -() const { return vec2(-a[0], -a[1]); }
	double operator [](int i) const { return a[i]; }
	double& operator [](int i) { return a[i]; }
	vec2& operator +=(const vec2& v) { a[0] += v.a[0]; a[1] += v.a[1]; return *this; }
	vec2& operator *=(const double& x) { a[0] *= x; a[1] *= x; return *this; }
	vec2& operator /=(const double& x) { return *this *= (1 / x); }

	double length() const { return std::sqrt(lengthSquared()); }
	double lengthSquared() const { return a[0] * a[0] + a[1] * a[1]; }

	vec2 xy() const { return vec2(a[0], a[1]); }
};

class vec3
{
public : 
	double a[3];

	vec3(double x = 0) { a[0] = a[1] = a[2] = x; }
	vec3(double x, double y, double z) { a[0] = x; a[1] = y; a[2] = z; }

	double x() { return a[0]; }
	double y() { return a[1]; }
	double z() { return a[2]; }
	double x() const { return a[0]; }
	double y() const { return a[1]; }
	double z() const { return a[2]; }

	vec3 operator -() const { return vec3(-a[0], -a[1], -a[2]); }
	double operator [](int i) const { return a[i]; }
	double& operator [](int i) { return a[i]; }
	vec3& operator +=(const vec3& v) { a[0] += v.a[0]; a[1] += v.a[1]; a[2] += v.a[2]; return *this; }
	vec3& operator *=(const double& x) { a[0] *= x; a[1] *= x; a[2] *= x; return *this; }
	vec3& operator *=(const vec3& x) { a[0] *= x.a[0]; a[1] *= x.a[1]; a[2] *= x.a[2]; return *this; }
	vec3& operator /=(const double& x) { return *this *= (1 / x); }

	double length() const { return std::sqrt(lengthSquared()); }
	double lengthSquared() const { return a[0] * a[0] + a[1] * a[1] + a[2] * a[2]; }

	vec2 xy() const { return vec2(a[0], a[1]); }
	bool near_zero() const {
		auto eps = 1e-8;
		return (fabs(a[0]) < eps) && (fabs(a[1]) < eps) && (fabs(a[2]) < eps);
	}
};

inline vec2 sin(const vec2& x) { return vec2(sin(x.a[0]), sin(x.a[1])); }
inline vec3 sin(const vec3& x) { return vec3(sin(x.a[0]), sin(x.a[1]), sin(x.a[2])); }

inline vec2 operator +(const vec2& x, const vec2& y) { return vec2(x.a[0] + y.a[0], x.a[1] + y.a[1]); }
inline vec3 operator +(const vec3& x, const vec3& y) { return vec3(x.a[0] + y.a[0], x.a[1] + y.a[1], x.a[2] + y.a[2]); }

inline vec2 operator -(const vec2& x, const vec2& y) { return vec2(x.a[0] - y.a[0], x.a[1] - y.a[1]); }
inline vec3 operator -(const vec3& x, const vec3& y) { return vec3(x.a[0] - y.a[0], x.a[1] - y.a[1], x.a[2] - y.a[2]); }

inline vec2 operator *(const double y, const vec2& x) { return vec2(x.a[0] * y, x.a[1] * y); }
inline vec2 operator *(const vec2& x, const double y) { return vec2(x.a[0] * y, x.a[1] * y); }
inline vec2 operator *(const vec2& x, const vec2& y) { return vec2(x.a[0] * y.a[0], x.a[1] * y.a[1]); }
inline vec3 operator *(const double y, const vec3& x) { return vec3(x.a[0] * y, x.a[1] * y, x.a[2] * y); }
inline vec3 operator *(const vec3& x, const double y) { return vec3(x.a[0] * y, x.a[1] * y, x.a[2] * y); }
inline vec3 operator *(const vec3& x, const vec3& y) { return vec3(x.a[0] * y.a[0], x.a[1] * y.a[1], x.a[2] * y.a[2]); }

inline vec2 operator /(const vec2& x, const double y) { return vec2(x.a[0] / y, x.a[1] / y); }
inline vec3 operator /(const vec3& x, const double y) { return vec3(x.a[0] / y, x.a[1] / y, x.a[2] / y); }

inline vec3 Abs(vec3 x) { return vec3(std::abs(x.a[0]), std::abs(x.a[1]), std::abs(x.a[2])); }

inline double Sqr(double x) { return x * x; }

inline double CosTheta(vec3 x) { return x.z(); }
inline double Cos2Theta(vec3 x) { return x.z() * x.z(); }
inline double Sin2Theta(vec3 x) { return std::fmax(0.0, 1 - Cos2Theta(x)); }
inline double SinTheta(vec3 x) { return std::sqrt(Sin2Theta(x)); }
inline double TanTheta(vec3 x) { return SinTheta(x) / CosTheta(x); }
inline double Tan2Theta(vec3 x) { return Sin2Theta(x) / Cos2Theta(x); }

inline double CosPhi(vec3 x) { double sinTheta = SinTheta(x); return (sinTheta == 0) ? 1 : Clamp(x.x() / sinTheta, -1, 1); }
inline double SinPhi(vec3 x) { double sinTheta = SinTheta(x); return (sinTheta == 0) ? 0 : Clamp(x.y() / sinTheta, -1, 1); }

inline double AbsCosTheta(vec3 x) { return std::abs(x.z()); }
inline vec3 normalize(const vec3& x) { return x / x.length(); }
inline double dot(const vec2& x, const vec2& y) { return x.a[0] * y.a[0] + x.a[1] * y.a[1]; }
inline double dot(const vec3& x, const vec3& y) { return x.a[0] * y.a[0] + x.a[1] * y.a[1] + x.a[2] * y.a[2]; }
inline vec3 cross(const vec3& x, const vec3& y) { return vec3(x.a[1] * y.a[2] - x.a[2] * y.a[1], x.a[2] * y.a[0] - x.a[0] * y.a[2], x.a[0] * y.a[1] - x.a[1] * y.a[0]); };
inline vec3 Reflect(vec3 x, vec3 y) { return -x + 2 * dot(x, y) * y; }
inline bool Refract(vec3 wi, vec3 n, double eta, double* etap, vec3* wt) {
	double cosThetai = dot(n, wi);
	if (cosThetai < 0) { eta = 1 / eta; cosThetai = -cosThetai; n = -n; }

	double sin2Thetai = std::fmax(0.0, 1 - cosThetai * cosThetai);
	double sin2Thetat = sin2Thetai / (eta * eta);
	if (sin2Thetat >= 1) return false;
	double cosThetat = std::sqrt(1.0 - sin2Thetat);

	*wt = -wi / eta + (cosThetai / eta - cosThetat) * n;
	if (etap) *etap = eta;
	return true;
}
inline double FrDielectric(double cosThetai, double eta) {
	cosThetai = Clamp(cosThetai, -1, 1);
	if (cosThetai < 0) { eta = 1 / eta; cosThetai = -cosThetai; }

	double sin2Thetai = 1 - cosThetai * cosThetai;
	double sin2Thetat = sin2Thetai / (eta * eta);
	if (sin2Thetat >= 1) return 1.0;
	double cosThetat = std::sqrt(1.0 - sin2Thetat);

	double r_parl = (eta * cosThetai - cosThetat) / (eta * cosThetai + cosThetat);
	double r_perp = (cosThetai - eta * cosThetat) / (cosThetai + eta * cosThetat);
	return (r_parl * r_parl + r_perp * r_perp) / 2.0;
}
inline double FrComplex(double cosThetai, complex eta) {
	cosThetai = Clamp(cosThetai, 0, 1);

	double sin2Thetai = 1 - cosThetai * cosThetai;
	complex sin2Thetat = complex(sin2Thetai) / (eta * eta);
	complex cosThetat = (complex(1) - sin2Thetat).sqrt();

	complex r_parl = (eta * complex(cosThetai) - cosThetat) / (eta * complex(cosThetai) + cosThetat);
	complex r_perp = (complex(cosThetai) - eta * cosThetat) / (complex(cosThetai) + eta * cosThetat);
	return (r_parl.norm() + r_perp.norm()) / 2;
}

// ======== interval ==========
class Interval
{
public:
	double low, high;
	Interval(double _low, double _high) : low(std::fmin(_low, _high)), high(std::fmax(_low, _high)) {}
	Interval(double x) : low(x), high(x) {}
	Interval() : low(0), high(0) {}

	Interval operator +(const Interval& other) const { return Interval(AddRoundDown(low, other.low), AddRoundUp(high, other.high)); }
	Interval operator -(const Interval& other) const { return Interval(SubRoundDown(low, other.low), SubRoundUp(high, other.high)); }
	Interval operator *(const Interval& other) const {
		double a = MulRoundDown(low, other.low), b = MulRoundDown(low, other.high), c = MulRoundDown(high, other.low), d = MulRoundDown(high, other.high);
		double A = MulRoundUp(low, other.low), B = MulRoundUp(low, other.high), C = MulRoundUp(high, other.low), D = MulRoundUp(high, other.high);
		return Interval(std::min(std::min(a, b), std::min(c, d)), std::max(std::max(A, B), std::max(C, D)));
	}
	Interval operator /(const Interval& other) const {
		if (other.low <= 0 && other.high >= 0) return Interval(-INFINITY, INFINITY);
		double a = DivRoundDown(low, other.low), b = DivRoundDown(low, other.high), c = DivRoundDown(high, other.low), d = DivRoundDown(high, other.high);
		double A = DivRoundUp(low, other.low), B = DivRoundUp(low, other.high), C = DivRoundUp(high, other.low), D = DivRoundUp(high, other.high);
		return Interval(std::min(std::min(a, b), std::min(c, d)), std::max(std::max(A, B), std::max(C, D)));
	}
	Interval operator -() const { return Interval(-high, -low); }

	Interval& operator +=(Interval other) { return *this = *this + other; }
	Interval& operator -=(Interval other) { return *this = *this - other; }
	Interval& operator *=(Interval other) { return *this = *this * other; }
	Interval& operator /=(Interval other) { return *this = *this / other; }

	Interval& operator +=(double f) { return *this += Interval(f); }
	Interval& operator -=(double f) { return *this -= Interval(f); }
	Interval& operator *=(double f) {
		if (f > 0) *this = Interval(MulRoundDown(low, f), MulRoundUp(high, f));
		else *this = Interval(MulRoundDown(high, f), MulRoundUp(low, f));
		return *this;
	}
	Interval& operator /=(double f) {
		if (f > 0) *this = Interval(DivRoundDown(low, f), DivRoundUp(high, f));
		else *this = Interval(DivRoundDown(high, f), DivRoundUp(low, f));
		return *this;
	}

	bool operator ==(const Interval& other) const { return low == other.low && high == other.high; }
	bool operator !=(const Interval& other) const { return low != other.low || high != other.high; }

	double Midpoint() const { return (low + high) / 2; }
	double Width() const { return high - low; }
};

inline Interval operator +(Interval i, double f) { return i + Interval(f); }
inline Interval operator -(Interval i, double f) { return i - Interval(f); }
inline Interval operator *(Interval i, double f) { 
	if (f > 0) return Interval(MulRoundDown(i.low, f), MulRoundUp(i.high, f));
	return Interval(MulRoundDown(i.high, f), MulRoundUp(i.low, f));
}
inline Interval operator /(Interval i, double f) { 
	if (f > 0) return Interval(DivRoundDown(i.low, f), DivRoundUp(i.high, f));
	return Interval(DivRoundDown(i.high, f), DivRoundUp(i.low, f));
}

static Interval FromValueAndError(double v, double err) {
	if (err == 0) return Interval(v);
	return Interval(SubRoundDown(v, err), AddRoundUp(v, err));
}

// =========== vector3fi ============
class Vector3fi
{
public:
	Interval x, y, z;
	Vector3fi(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}
	Vector3fi(Interval _x, Interval _y, Interval _z) : x(_x), y(_y), z(_z) {}
	Vector3fi(vec3 p) : x(p.x()), y(p.y()), z(p.z()) {}
	Vector3fi(vec3 v, vec3 e) : x(FromValueAndError(v.x(), e.x())), y(FromValueAndError(v.y(), e.y())), z(FromValueAndError(v.z(), e.z())) {}
	Vector3fi() : x(0), y(0), z(0) {}

	vec3 Error() const { return vec3(x.Width() / 2, y.Width() / 2, z.Width() / 2); }
	bool IsExact() const { return x.Width() == 0 && y.Width() == 0 && z.Width() == 0; }

	Vector3fi operator +(const Vector3fi& other) const { return Vector3fi(x + other.x, y + other.y, z + other.z); }
	Vector3fi operator -(const Vector3fi& other) const { return Vector3fi(x - other.x, y - other.y, z - other.z); }
	Vector3fi operator *(const Vector3fi& other) const { return Vector3fi(x * other.x, y * other.y, z * other.z); }
	Vector3fi operator /(const Vector3fi& other) const { return Vector3fi(x / other.x, y / other.y, z / other.z); }
	Vector3fi operator +(const double& f) const { return Vector3fi(x + f, y + f, z + f); }
	Vector3fi operator -(const double& f) const { return Vector3fi(x - f, y - f, z - f); }
	Vector3fi operator *(const double& f) const { return Vector3fi(x * f, y * f, z * f); }
	Vector3fi operator /(const double& f) const { return Vector3fi(x / f, y / f, z / f); }
};

// =========== matrix ============

template <int N> class SquareMatrix
{
public:
	SquareMatrix()
	{
		for (int i = 0; i < N; ++i)
			for (int j = 0; j < N; ++j)
				a[i][j] = (i == j) ? 1 : 0;
	}
	SquareMatrix(const double _a[N][N]) {
		for (int i = 0; i < N; ++i)
			for (int j = 0; j < N; ++j)
				a[i][j] = _a[i][j];
	}
	template <typename... Args> SquareMatrix(double x, Args... args) {
		static_assert(1 + sizeof...(Args) == N * N,
			"Incorrect number of values provided to SquareMatrix constructor");
		init(0, 0, x, args...);
	}

	double Determinant() const
	{
		if (N == 2) return a[0][0] * a[1][1] - a[0][1] * a[1][0];
		if (N == 3) return a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) -
			a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]) +
			a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
		if (N == 4)
		{
			double det = 0.0;
			for (int i = 0; i < 4; ++i) {
				double submat[3][3];
				for (int j = 0; j < 3; j++) {
					int k = 0;
					for (int l = 0; l < 4; l++) {
						if (l == i) continue;
						submat[j][k] = a[j + 1][l];
						k++;
					}
				}
				double subDet = submat[0][0] * (submat[1][1] * submat[2][2] - submat[1][2] * submat[2][1]) -
					submat[0][1] * (submat[1][0] * submat[2][2] - submat[1][2] * submat[2][0]) +
					submat[0][2] * (submat[1][0] * submat[2][1] - submat[1][1] * submat[2][0]);
				if (i % 2 == 0) det += subDet;
				else det -= subDet;
			}
			return det;
		}
		return 0;
	}
	SquareMatrix<N> Invert() const
	{
		SquareMatrix<N> result;
		double det = Determinant();
		assert(det != 0);
		double invDet = 1.0 / det;
		if (N == 2) {
			result[0][0] = a[1][1] * invDet;
			result[0][1] = -a[0][1] * invDet;
			result[1][0] = -a[1][0] * invDet;
			result[1][1] = a[0][0] * invDet;
		}
		if (N == 3)
		{
			result[0][0] = (a[1][1] * a[2][2] - a[1][2] * a[2][1]) * invDet;
			result[0][1] = (a[0][2] * a[2][1] - a[0][1] * a[2][2]) * invDet;
			result[0][2] = (a[0][1] * a[1][2] - a[0][2] * a[1][1]) * invDet;
			result[1][0] = (a[1][2] * a[2][0] - a[1][0] * a[2][2]) * invDet;
			result[1][1] = (a[0][0] * a[2][2] - a[0][2] * a[2][0]) * invDet;
			result[1][2] = (a[0][2] * a[1][0] - a[0][0] * a[1][2]) * invDet;
			result[2][0] = (a[1][0] * a[2][1] - a[1][1] * a[2][0]) * invDet;
			result[2][1] = (a[0][1] * a[2][0] - a[0][0] * a[2][1]) * invDet;
			result[2][2] = (a[0][0] * a[1][1] - a[0][1] * a[1][0]) * invDet;
		}
		if (N == 4)
		{
			double A00 = a[0][0], A01 = a[0][1], A02 = a[0][2], A03 = a[0][3];
			double A10 = a[1][0], A11 = a[1][1], A12 = a[1][2], A13 = a[1][3];
			double A20 = a[2][0], A21 = a[2][1], A22 = a[2][2], A23 = a[2][3];
			double A30 = a[3][0], A31 = a[3][1], A32 = a[3][2], A33 = a[3][3];

			double cofactor00 = A11 * (A22 * A33 - A23 * A32) - A12 * (A21 * A33 - A23 * A31) + A13 * (A21 * A32 - A22 * A31);
			double cofactor01 = -(A10 * (A22 * A33 - A23 * A32) - A12 * (A20 * A33 - A23 * A30) + A13 * (A20 * A32 - A22 * A30));
			double cofactor02 = A10 * (A21 * A33 - A23 * A31) - A11 * (A20 * A33 - A23 * A30) + A13 * (A20 * A31 - A21 * A30);
			double cofactor03 = -(A10 * (A21 * A32 - A22 * A31) - A11 * (A20 * A32 - A22 * A30) + A12 * (A20 * A31 - A21 * A30));

			double cofactor10 = -(A01 * (A22 * A33 - A23 * A32) - A02 * (A21 * A33 - A23 * A31) + A03 * (A21 * A32 - A22 * A31));
			double cofactor11 = A00 * (A22 * A33 - A23 * A32) - A02 * (A20 * A33 - A23 * A30) + A03 * (A20 * A32 - A22 * A30);
			double cofactor12 = -(A00 * (A21 * A33 - A23 * A31) - A01 * (A20 * A33 - A23 * A30) + A03 * (A20 * A31 - A21 * A30));
			double cofactor13 = A00 * (A21 * A32 - A22 * A31) - A01 * (A20 * A32 - A22 * A30) + A02 * (A20 * A31 - A21 * A30);

			double cofactor20 = A01 * (A12 * A33 - A13 * A32) - A02 * (A11 * A33 - A13 * A31) + A03 * (A11 * A32 - A12 * A31);
			double cofactor21 = -(A00 * (A12 * A33 - A13 * A32) - A02 * (A10 * A33 - A13 * A30) + A03 * (A10 * A32 - A12 * A30));
			double cofactor22 = A00 * (A11 * A33 - A13 * A31) - A01 * (A10 * A33 - A13 * A30) + A03 * (A10 * A31 - A11 * A30);
			double cofactor23 = -(A00 * (A11 * A32 - A12 * A31) - A01 * (A10 * A32 - A12 * A30) + A02 * (A10 * A31 - A11 * A30));

			double cofactor30 = -(A01 * (A12 * A23 - A13 * A22) - A02 * (A11 * A23 - A13 * A21) + A03 * (A11 * A22 - A12 * A21));
			double cofactor31 = A00 * (A12 * A23 - A13 * A22) - A02 * (A10 * A23 - A13 * A20) + A03 * (A10 * A22 - A12 * A20);
			double cofactor32 = -(A00 * (A11 * A23 - A13 * A21) - A01 * (A10 * A23 - A13 * A20) + A03 * (A10 * A21 - A11 * A20));
			double cofactor33 = A00 * (A11 * A22 - A12 * A21) - A01 * (A10 * A22 - A12 * A20) + A02 * (A10 * A21 - A11 * A20);

			result[0][0] = cofactor00 * invDet;
			result[0][1] = cofactor01 * invDet;
			result[0][2] = cofactor02 * invDet;
			result[0][3] = cofactor03 * invDet;
			result[1][0] = cofactor10 * invDet;
			result[1][1] = cofactor11 * invDet;
			result[1][2] = cofactor12 * invDet;
			result[1][3] = cofactor13 * invDet;
			result[2][0] = cofactor20 * invDet;
			result[2][1] = cofactor21 * invDet;
			result[2][2] = cofactor22 * invDet;
			result[2][3] = cofactor23 * invDet;
			result[3][0] = cofactor30 * invDet;
			result[3][1] = cofactor31 * invDet;
			result[3][2] = cofactor32 * invDet;
			result[3][3] = cofactor33 * invDet;
		}
		return result;
	}
	const double* operator [](int i) const { return a[i]; }
	double* operator [](int i) { return a[i]; }

	SquareMatrix operator +(const SquareMatrix<N>& other) const {
		SquareMatrix r = *this;
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
				r.a[i][j] += other.a[i][j];
		return r;
	}

private:
	double a[N][N];

	void init(int, int) {}
	template <typename... Args> void init(int i, int j, double x, Args... args) {
		a[i][j] = x;
		if (j + 1 == N) init(i + 1, 0, args...);
		else init(i, j + 1, args...);
	}
};

vec3 operator *(const SquareMatrix<3>& m, const vec3& v);

template<int N> SquareMatrix<N> operator *(const SquareMatrix<N>& m1, const SquareMatrix<N>& m2) {
	SquareMatrix<N> ans;
	for (int i = 0; i < N; ++i)
		for (int j = 0; j < N; ++j) {
			ans[i][j] = 0;
			for (int k = 0; k < N; ++k)
				ans[i][j] += m1[i][k] * m2[k][j];
		}
	return ans;
}

// =========== transform ============
class Transform
{
public : 
	Transform(const SquareMatrix<4> _mat) : mat(_mat) { inv = _mat.Invert(); }
	Transform(const SquareMatrix<4> _mat, const SquareMatrix<4> _inv) : mat(_mat), inv(_inv) {}
	Transform(const double mat[4][4]) : Transform(SquareMatrix<4>(mat)) {}
	Transform() : Transform(SquareMatrix<4>()) {}

	Transform Inverse() {
		return Transform(inv, mat);
	}

	vec3 operator()(const vec3& p) const;
	Vector3fi operator()(const Vector3fi& p) const;
	static Transform Translate(const vec3& delta);
	static Transform Scale(double x, double y, double z);
	static Transform RotateX(double theta);
	static Transform RotateY(double theta);
	static Transform RotateZ(double theta);
	static Transform Rotate(double theta, const vec3& axis);

private : 
	SquareMatrix<4> mat;
	SquareMatrix<4> inv;
};

// =========== frame ================
class Frame {
public : 
	vec3 x, y, z;
	Frame() : x(1, 0, 0), y(0, 1, 0), z(0, 0, 1) {}
	Frame(vec3 _x, vec3 _y, vec3 _z) : x(_x), y(_y), z(_z) {}
	vec3 ToLocal(vec3 v) const { return vec3(dot(v, x), dot(v, y), dot(v, z)); }
	vec3 FromLocal(vec3 v) const { return v.x() * x + v.y() * y + v.z() * z; }
};
static Frame FromXZ(vec3 x, vec3 z) { return Frame(x, cross(x, z), z); }

// =========== others ============
inline vec3 FaceForward(vec3 n, vec3 v) { return (dot(n, v) < 0) ? -n : n; }
