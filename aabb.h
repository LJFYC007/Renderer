#pragma once
#include "math.h"
#include "interval.h"
#include "ray.h"

class AABB 
{
public:
	interval x, y, z;

	AABB() {}
	AABB(const interval& _x, const interval& _y, const interval& _z) : x(_x), y(_y), z(_z) {}
	AABB(const vec3& a, const vec3& b) {
		x = interval(fmin(a[0], b[0]), fmax(a[0], b[0]));
		y = interval(fmin(a[1], b[1]), fmax(a[1], b[1]));
		z = interval(fmin(a[2], b[2]), fmax(a[2], b[2]));
	}
	AABB(const AABB& a, const AABB& b) {
		x = interval(a.x, b.x);
		y = interval(a.y, b.y);
		z = interval(a.z, b.z);
	}

	AABB Union(vec3 a)
	{
		x = interval(fmin(x.Min, a[0]), fmax(x.Max, a[0]));
		y = interval(fmin(y.Min, a[1]), fmax(y.Max, a[1]));
		z = interval(fmin(z.Min, a[2]), fmax(z.Max, a[2]));
		return *this;
	}

	const interval& axis(int n) const {
		if (n == 0) return x;
		if (n == 1) return y;
		if (n == 2) return z;
		return x;
	}

	bool Intersect(const ray& r, interval t, double& hitt0, double& hitt1) const {
		for (int a = 0; a < 3; ++a) {
			auto invD = 1 / r.rd[a];
			auto orig = r.ro[a];

			auto t0 = (axis(a).Min - orig) * invD;
			auto t1 = (axis(a).Max - orig) * invD;

			if (invD < 0)
				std::swap(t0, t1);

			if (t0 > t.Min) t.Min = t0;
			if (t1 < t.Max) t.Max = t1;

			if (t.Max <= t.Min)
				return false;
		}
		hitt0 = t.Min; hitt1 = t.Max;
		return true;
	}
};


