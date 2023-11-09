#pragma once
#include <cmath>
#include <assert.h>

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

// =========== others ============

inline double Lerp(double t, double x, double y)
{
	return (1 - t) * x + t * y;
}

class Frame {
public : 
	vec3 x, y, z;
	Frame() : x(1, 0, 0), y(0, 1, 0), z(0, 0, 1) {}
	Frame(vec3 _x, vec3 _y, vec3 _z) : x(_x), y(_y), z(_z) {}
	vec3 ToLocal(vec3 v) const { return vec3(dot(v, x), dot(v, y), dot(v, z)); }
	vec3 FromLocal(vec3 v) const { return v.x() * x + v.y() * y + v.z() * z; }
};
static Frame FromXZ(vec3 x, vec3 z) { return Frame(x, cross(x, z), z); }
