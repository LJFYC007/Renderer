#pragma once
#include "bvh.h"
#include "sample.h"
#include "math.h"
#include "interval.h"
#include "spectrum.h"
#include "color.h"
#include "colorspace.h"
#include "bsdf.h"
#include "lights.h"
#include "lightsampler.h"

class Camera;

static bool Unoccluded(const BVHAggregate& bvh, const Interaction& p0, const Interaction& p1)
{
	Ray r = SpawnRayTo(p0.pi, p0.n, p1.pi, p1.n);
	std::optional<ShapeIntersection> isect = bvh.Intersect(r, interval(0, 0.9999999));
	return !isect;
}

SampledSpectrum SampleLd(const SurfaceInteraction& intr, const BSDF& bsdf, Camera* camera, const SampledWaveLengths& lambda, const BVHAggregate& bvh);

SampledSpectrum Li(RayDifferential r, Camera* camera, SampledWaveLengths& lambda, const BVHAggregate& bvh);

