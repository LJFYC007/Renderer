#pragma once
#include "math.h"
#include "interval.h"
#include "ray.h"

class bvh
{
public : 
	interval x, y, z;

	bvh() {}
	bvh(const interval& _x, const interval& _y, const interval& _z) : x(_x), y(_y), z(_z) {}
	bvh(const vec3& a, const vec3& b) {
		x = interval(fmin(a[0], b[0]), fmax(a[0], b[0]));
		y = interval(fmin(a[1], b[1]), fmax(a[1], b[1]));
		z = interval(fmin(a[2], b[2]), fmax(a[2], b[2]));
	}
	bvh(const bvh& a, const bvh& b) {
		x = interval(a.x, b.x);
		y = interval(a.y, b.y);
		z = interval(a.z, b.z);
	}

	const interval& axis(int n) const {
		if (n == 0) return x;
		if (n == 1) return y;
		if (n == 2) return z;
		return x;
	}

	bool hit(const ray& r, interval t) const {
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
		return true;
	}
};
