#pragma once
#include "math.h"

class ray
{
public : 
	vec3 ro;
	vec3 rd;
	double tMax, time;

	ray() {}
	ray(const vec3& origin, const vec3& direction) { ro = origin; rd = direction; }
	vec3 at(double t) const { return ro + t * rd; }
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

inline ray SpawnRay(Vector3fi pi, vec3 n, vec3 d) {
	return ray(OffsetRayOrigin(pi, n, d), d);
}

inline ray SpawnRayTo(Vector3fi pFrom, vec3 n, vec3 pTo) {
	vec3 d = pTo - vec3(pFrom);
	return SpawnRay(pFrom, n, d);
}

inline ray SpawnRayTo(Vector3fi pFrom, vec3 nFrom, Vector3fi pTo, vec3 nTo) {
	vec3 pf = OffsetRayOrigin(pFrom, nFrom, vec3(pTo) - vec3(pFrom));
	vec3 pt = OffsetRayOrigin(pTo, nTo, pf - vec3(pTo));
	return ray(pf, pt - pf);
}
