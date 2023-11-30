#pragma once
#include "spectrum.h"
#include "shape.h"
#include "spectrum.h"
#include "ray.h"

#include <optional>

enum class LightType {
	DeltaPosition,
	DeltaDirection,
	Area,
	Infinite
};

class LightSampleContext {
public :
	vec3 p, n, ns;
	LightSampleContext(const hitRecord& rec) : p(rec.p), n(rec.normal), ns(rec.normal) {};
};

class LightLiSample {
public : 
	SampledSpectrum L;
	vec3 wi;
	double pdf;
	vec3 p;
	LightLiSample(SampledSpectrum _L, vec3 _wi, double _pdf, vec3 _p) : L(_L), wi(_wi), pdf(_pdf), p(_p) {}
};

class Light
{
public : 
	Light(LightType _type, const Transform& _renderFromLight) : type(_type), renderFromLight(_renderFromLight) {}

	LightType Type() const { return type; }

	virtual SampledSpectrum L(vec3 p, vec3 n, vec2 uv, vec3 w, const SampledWaveLengths& lambda) const {
		return SampledSpectrum(0.0);
	}

	virtual SampledSpectrum Le(ray r, const SampledWaveLengths& lambda) const {
		return SampledSpectrum(0.0);
	}

	virtual	SampledSpectrum Phi(SampledWaveLengths lambda) const = 0;
	virtual std::optional<LightLiSample> SampleLi(LightSampleContext sample, vec2 u, SampledWaveLengths lambda, bool allowIncompletePDF = false) const = 0;
	virtual double PDF_Li(LightSampleContext sample, vec3 wi, bool allowIncompletePDF = false) const = 0;
	// virtual void PDF_Le(ray r, double* pdfPos, double* pdfDir) const = 0;

public : 
	LightType type;
	Transform renderFromLight;
};
