#pragma once
#include "spectrum.h"
#include "shape.h"
#include "ray.h"
#include "interaction.h"

#include <optional>

enum class LightType {
	DeltaPosition,
	DeltaDirection,
	Area,
	Infinite
};

inline bool IsDeltaLight(LightType type) {
	return (type == LightType::DeltaPosition) || (type == LightType::DeltaDirection);
}

class LightSampleContext {
public :
	Vector3fi pi;
	vec3 n, ns;
	LightSampleContext() = default;
	LightSampleContext(const SurfaceInteraction& intr) : pi(intr.pi), n(intr.n), ns(intr.shading.n) {}
	vec3 p() const { return vec3(pi); }
};

class LightLiSample {
public : 
	SampledSpectrum L;
	vec3 wi;
	double pdf;
	Interaction pLight;
	LightLiSample(SampledSpectrum _L, vec3 _wi, double _pdf, const Interaction &_pLight) : L(_L), wi(_wi), pdf(_pdf), pLight(_pLight) {}
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
