#pragma once
#include "math.h"

class ray
{
public : 
	point3 ro;
	vec3 rd;
	double tMax, time;

	ray() {}
	ray(const point3& origin, const vec3& direction) { ro = origin; rd = direction; }
	vec3 at(double t) const { return ro + t * rd; }
};
