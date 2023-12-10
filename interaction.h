#pragma once
#include <memory>
using std::shared_ptr;

class Material;
class Light;

class Interaction
{
public: 
	Interaction(vec3 _p) : p(_p) {}
	Interaction(vec3 _p, vec3 _n, vec2 _uv) : p(_p), n(_n), uv(_uv) {}
	Interaction(vec3 _p, vec3 _n, vec2 _uv, vec3 _wo, double _t) : p(_p), n(_n), uv(_uv), wo(_wo), t(_t) {}

	vec3 p, wo, n;
	double t;
	vec2 uv;
};

class SurfaceInteraction : public Interaction
{
public:
	SurfaceInteraction(vec3 p, vec2 uv, vec3 wo, vec3 dpdu, vec3 dpdv, vec3 dndu, vec3 dndv, double t, bool flipNormal) :
		Interaction(p, normalize(cross(dpdu, dpdv)), uv, wo, t), dpdu(dpdu), dpdv(dpdv), dndu(dndu), dndv(dndv) {
		if (flipNormal) n = -n;
		shading.n = n;
		shading.dpdu = dpdu;
		shading.dpdv = dpdv;
		shading.dndu = dndu;
		shading.dndv = dndv;
	}

	vec3 dpdu, dpdv, dndu, dndv;
	struct {
		vec3 n, dpdu, dpdv, dndu, dndv;
	} shading;
	shared_ptr<Material> material;
	shared_ptr<Light> areaLight;
	vec3 dpdx, dpdy;
	double dudx, dvdx, dudy, dvdy;
};
