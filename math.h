#pragma once
#include <cmath>
#include <assert.h>

// =========== vector ============

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

using point3 = vec3;

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

inline vec3 normalize(const vec3& x) { return x / x.length(); }
inline double dot(const vec2& x, const vec2& y) { return x.a[0] * y.a[0] + x.a[1] * y.a[1]; }
inline double dot(const vec3& x, const vec3& y) { return x.a[0] * y.a[0] + x.a[1] * y.a[1] + x.a[2] * y.a[2]; }
inline vec3 cross(const vec3& x, const vec3& y) { return vec3(x.a[1] * y.a[2] - x.a[2] * y.a[1], x.a[2] * y.a[0] - x.a[0] * y.a[2], x.a[0] * y.a[1] - x.a[1] * y.a[0]); };
inline vec3 reflect(const vec3& x, const vec3& y)
{
	return x - 2 * dot(x, y) * y;
}
inline vec3 refract(const vec3& uv, const vec3& n, double etai_over_etat) {
	auto cos_theta = fmin(dot(-uv, n), 1.0);
	vec3 r_out_perp = etai_over_etat * (uv + cos_theta * n);
	vec3 r_out_parallel = -sqrt(fabs(1.0 - r_out_perp.lengthSquared())) * n;
	return r_out_perp + r_out_parallel;
}

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

	const double* operator [](int i) const { return a[i]; }
	double* operator [](int i) { return a[i]; }

	SquareMatrix operator +(const SquareMatrix<N>& other) const {
		SquareMatrix r = *this;
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
				r.a[i][j] += other.a[i][j];
		return r;
	}

	double Determinant() const {
		if (N == 2) return a[0][0] * a[1][1] - a[0][1] * a[1][0];
		if (N == 3) return a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) -
			a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]) +
			a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
		return 0;
	}

	SquareMatrix Invert() const {
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
		return result;
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

// =========== others ============

inline double Lerp(double t, double x, double y)
{
	return (1 - t) * x + t * y;
}
