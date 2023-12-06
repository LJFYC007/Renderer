#pragma once

class Interaction
{
public: 
	Interaction(vec3 _p) : p(_p) {}
	Interaction(vec3 _p, vec3 _n, vec2 uv) : p(_p), n(_n), uv(uv) {}

	vec3 p, wo, n;
	double t;
	vec2 uv;
};
