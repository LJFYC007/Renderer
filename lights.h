#pragma once
#include "light.h"
#include "math.h"

class PointLight : public Light
{
public : 
	PointLight(const Transform &_renderFromLight, const Spectrum &_intensity, double _scale) 
		: Light(LightType::DeltaPosition, _renderFromLight), intensity(_intensity), scale(_scale) {}

	std::optional<LightLiSample> SampleLi(LightSampleContext sample, vec2 u, SampledWaveLengths lambda, bool allowIncompletePDF = false) const override {
		vec3 p = renderFromLight(vec3(0.0));
		vec3 wi = normalize(p - sample.p);
		SampledSpectrum Li = intensity.Sample(lambda) * scale; // TODO: fix this, need add / (p - sample.p).lengthSquared();
		return LightLiSample(Li, wi, 1.0, p);
	}

	double PDF_Li(LightSampleContext sample, vec3 wi, bool allowIncompletePDF = false) const override {
		return 0.0;
	}

	SampledSpectrum Phi(SampledWaveLengths lambda) const override {
		return intensity.Sample(lambda) * scale * 4 * pi;
	}
private : 
	DenselySampledSpectrum intensity;
	double scale;
};

class DistantLight : public Light
{
public:  
	DistantLight(const Transform &_renderFromLight, const Spectrum &_intensity, double _scale) 
		: Light(LightType::DeltaDirection, _renderFromLight), intensity(_intensity), scale(_scale) {}

	std::optional<LightLiSample> SampleLi(LightSampleContext sample, vec2 u, SampledWaveLengths lambda, bool allowIncompletePDF = false) const override {
		vec3 wi = normalize(renderFromLight(vec3(0, 0, 1)));
		vec3 pOutside = sample.p + wi * 1000000.0; // TODO: fix this, need to multiply R
		SampledSpectrum Li = intensity.Sample(lambda) * scale;
		return LightLiSample(Li, wi, 1.0, pOutside);
	}

	double PDF_Li(LightSampleContext sample, vec3 wi, bool allowIncompletePDF = false) const override {
		return 0.0;
	}

	SampledSpectrum Phi(SampledWaveLengths lambda) const override {
		return intensity.Sample(lambda) * scale * pi; // TODO: fix this, need to multiply R^2 
	}
private:
	DenselySampledSpectrum intensity;
	double scale;
};
