#pragma once
#include "math.h"
#include "spectrum.h"
#include "interaction.h"
#include "color.h"
#include "colorspace.h"

#include <memory>
using std::shared_ptr;
using std::make_shared;

struct TextureEvalContext 
{
	TextureEvalContext(const Interaction& intr) : p(intr.pi), uv(intr.uv) {}
	TextureEvalContext(const SurfaceInteraction& isect) 
		: p(isect.p()), dpdx(isect.dpdx), dpdy(isect.dpdy), n(isect.n), uv(isect.uv), dudx(isect.dudx), dudy(isect.dudy), dvdx(isect.dvdx), dvdy(isect.dvdy), faceIndex(isect.faceIndex) {}
	TextureEvalContext(vec3 p, vec3 dpdx, vec3 dpdy, vec3 n, vec2 uv, double dudx, double dudy, double dvdx, double dvdy, int faceIndex = 0)
		: p(p), dpdx(dpdx), dpdy(dpdy), n(n), uv(uv), dudx(dudx), dudy(dudy), dvdx(dvdx), dvdy(dvdy), faceIndex(faceIndex) {}
	vec3 p, dpdx, dpdy, n;
	vec2 uv;
	double dudx = 0, dudy = 0, dvdx = 0, dvdy = 0;
	int faceIndex = 0;
};

class DoubleTexture
{
public:
	virtual double Evaluate(TextureEvalContext ctx) const = 0;
};

class SpectrumTexture
{
public:
	virtual SampledSpectrum Evaluate(TextureEvalContext ctx, SampledWaveLengths lambda) const = 0;
	virtual vec3 Evaluate(TextureEvalContext ctx) const {
		return vec3(0.0);
	}
};
