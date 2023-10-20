#pragma once
#include "math.h"

class ray
{
public : 
	vec3 ro, rd;

	ray() {}
	ray(const vec3& origin, const vec3& direction) { ro = origin; rd = direction; }
	vec3 at(double t) const { return ro + t * rd; }
};
