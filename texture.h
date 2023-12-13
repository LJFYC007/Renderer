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
};
