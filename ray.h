#pragma once
#include "math.h"

class Ray
{
public: 
	vec3 ro, rd;
	double time;
	Ray(vec3 _ro, vec3 _rd, double _time = 0) : ro(_ro), rd(_rd), time(_time) {}
	vec3 operator()(double t) const { return ro + t * rd; }
};

inline vec3 OffsetRayOrigin(Vector3fi pi, vec3 n, vec3 w)
{
	double d = dot(Abs(n), pi.Error());
	vec3 offset = d * n;
	if (dot(w, n) < 0) offset = -offset;
	vec3 po = vec3(pi.x.Midpoint(), pi.y.Midpoint(), pi.z.Midpoint()) + offset;
	for (int i = 0; i < 3; ++i)
		if (offset[i] > 0) po[i] = NextDoubleUp(po[i]);
		else if (offset[i] < 0) po[i] = NextDoubleDown(po[i]);
	return po;
}

inline Ray SpawnRay(Vector3fi pi, vec3 n, double time, vec3 d) {
	return Ray(OffsetRayOrigin(pi, n, d), d, time);
}

inline Ray SpawnRayTo(Vector3fi pFrom, vec3 n, double time, vec3 pTo) {
	vec3 d = pTo - vec3(pFrom);
	return SpawnRay(pFrom, n, time, d);
}

inline Ray SpawnRayTo(Vector3fi pFrom, vec3 nFrom, Vector3fi pTo, vec3 nTo) {
	vec3 pf = OffsetRayOrigin(pFrom, nFrom, vec3(pTo) - vec3(pFrom));
	vec3 pt = OffsetRayOrigin(pTo, nTo, pf - vec3(pTo));
	return Ray(pf, pt - pf);
}

class RayDifferential : public Ray {
public:
	RayDifferential(vec3 ro, vec3 rd, double time = 0) : Ray(ro, rd, time) {}
	bool hasDifferentials = false;
	vec3 rxOrigin, ryOrigin, rxDirection, ryDirection;
};
