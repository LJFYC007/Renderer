#pragma once
#include "math.h"

vec3 operator *(const SquareMatrix<3>& m, const vec3& v)
{
	return vec3(m[0][0] * v.a[0] + m[0][1] * v.a[1] + m[0][2] * v.a[2],
		m[1][0] * v.a[0] + m[1][1] * v.a[1] + m[1][2] * v.a[2],
		m[2][0] * v.a[0] + m[2][1] * v.a[1] + m[2][2] * v.a[2]
	);
}
