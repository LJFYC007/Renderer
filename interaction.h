#pragma once
#include "math.h"
#include "spectrum.h"
#include "ray.h"

#include <memory>
using std::shared_ptr;

class Material;
class Light;

class Interaction
{
public: 
	Interaction(Vector3fi _pi) : pi(_pi) {}
	Interaction(Vector3fi _pi, vec3 _n, vec2 _uv) : pi(_pi), n(_n), uv(_uv) {}
	Interaction(Vector3fi _pi, vec3 _n, vec2 _uv, vec3 _wo, double _t) : pi(_pi), n(_n), uv(_uv), wo(_wo), t(_t) {}
	vec3 p() const { return vec3(pi); }

	Vector3fi pi;
	vec3 wo, n;
	double t;
	vec2 uv;
};

class SurfaceInteraction : public Interaction
{
public:
	SurfaceInteraction(Vector3fi pi, vec2 uv, vec3 wo, vec3 dpdu, vec3 dpdv, vec3 dndu, vec3 dndv, double t, bool flipNormal) :
		Interaction(pi, normalize(cross(dpdu, dpdv)), uv, wo, t), dpdu(dpdu), dpdv(dpdv), dndu(dndu), dndv(dndv) {
		if (flipNormal) n = -n;
		shading.n = n;
		shading.dpdu = dpdu;
		shading.dpdv = dpdv;
		shading.dndu = dndu;
		shading.dndv = dndv;
	}

	void SetIntersectionProperties(shared_ptr<Material> _material, shared_ptr<Light> _areaLight) { material = _material; areaLight = _areaLight; }

	void SetShadingGeometry(vec3 ns, vec3 dpdus, vec3 dpdvs, vec3 dndus, vec3 dndvs) {
		shading.n = FaceForward(ns, n);
		shading.dpdu = dpdus;
		shading.dpdv = dpdvs;
		shading.dndu = dndus;
		shading.dndv = dndvs;
	}

	SampledSpectrum Le(vec3 w, const SampledWaveLengths& lambda) const;

	vec3 dpdu, dpdv, dndu, dndv;
	struct {
		vec3 n, dpdu, dpdv, dndu, dndv;
	} shading;
	shared_ptr<Material> material;
	shared_ptr<Light> areaLight;
	vec3 dpdx, dpdy;
	double dudx, dvdx, dudy, dvdy;
	int faceIndex = 0;
};
