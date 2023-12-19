#pragma once
#include "math.h"
#include "interval.h"
#include "ray.h"

class AABB 
{
public:
	vec3 pMin, pMax;

	AABB() {}
	AABB(const interval& _x, const interval& _y, const interval& _z) {
		pMin = vec3(_x.Min, _y.Min, _z.Min);
		pMax = vec3(_x.Max, _y.Max, _z.Max);
	}
	AABB(const vec3& a, const vec3& b) {
		pMin = vec3(fmin(a[0], b[0]), fmin(a[1], b[1]), fmin(a[2], b[2]));
		pMax = vec3(fmax(a[0], b[0]), fmax(a[1], b[1]), fmax(a[2], b[2]));
	}
	AABB(const vec3& a, const vec3& b, const vec3& c) {
		pMin = vec3(fmin(a[0], fmin(b[0], c[0])), fmin(a[1], fmin(b[1], c[1])), fmin(a[2], fmin(b[2], c[2])));
		pMax = vec3(fmax(a[0], fmax(b[0], c[0])), fmax(a[1], fmax(b[1], c[1])), fmax(a[2], fmax(b[2], c[2])));
	}
	AABB(const AABB& a, const AABB& b) {
		pMin = vec3(fmin(a.pMin[0], b.pMin[0]), fmin(a.pMin[1], b.pMin[1]), fmin(a.pMin[2], b.pMin[2]));
		pMax = vec3(fmax(a.pMax[0], b.pMax[0]), fmax(a.pMax[1], b.pMax[1]), fmax(a.pMax[2], b.pMax[2]));
	}
	AABB(const vec3& a) { pMin = pMax = a; }

	vec3 Diagonal() const {
		return pMax - pMin;
	}

	vec3 Offset(vec3 p) const {
		vec3 o = p - pMin;
		if (pMax[0] > pMin[0]) o[0] /= pMax[0] - pMin[0];
		if (pMax[1] > pMin[1]) o[1] /= pMax[1] - pMin[1];
		if (pMax[2] > pMin[2]) o[2] /= pMax[2] - pMin[2];
		return o;
	}

	double SurfaceArea() const {
		vec3 d = Diagonal();
		return 2 * (d[0] * d[1] + d[0] * d[2] + d[1] * d[2]);
	}

	int MaxDimension() const {
		vec3 d = Diagonal();
		if (d[0] > d[1] && d[0] > d[2]) return 0;
		if (d[1] > d[2]) return 1;
		return 2;
	}

	vec3 operator[](int i) const {
		return i == 0 ? pMin : pMax;
	}

	bool Intersect(const vec3& ro, const vec3& rd, double raytMax, const vec3& invDir, const int dirIsNeg[3]) const {
		const AABB& bounds = *this;
		double tMin = (bounds[dirIsNeg[0]].x() - ro.x()) * invDir.x();
		double tMax = (bounds[1 - dirIsNeg[0]].x() - ro.x()) * invDir.x();
		double tyMin = (bounds[dirIsNeg[1]].y() - ro.y()) * invDir.y();
		double tyMax = (bounds[1 - dirIsNeg[1]].y() - ro.y()) * invDir.y();
		tMax *= 1 + 2 * gamma(3);
		tyMax *= 1 + 2 * gamma(3);
		
		if (tMin > tyMax || tyMin > tMax) return false;
		if (tyMin > tMin) tMin = tyMin;
		if (tyMax < tMax) tMax = tyMax;

		double tzMin = (bounds[dirIsNeg[2]].z() - ro.z()) * invDir.z();
		double tzMax = (bounds[1 - dirIsNeg[2]].z() - ro.z()) * invDir.z();
		tzMax *= 1 + 2 * gamma(3);
		if (tMin > tzMax || tzMin > tMax) return false;
		if (tzMin > tMin) tMin = tzMin;
		if (tzMax < tMax) tMax = tzMax;
		
		return (tMin < raytMax) && (tMax > 0);
	}
};


