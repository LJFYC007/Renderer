#pragma once
#include "light.h"
#include "math.h"
#include "shape.h"

class PointLight : public Light
{
public : 
	PointLight(const Transform &_renderFromLight, const Spectrum &_intensity, double _scale) 
		: Light(LightType::DeltaPosition, _renderFromLight), intensity(_intensity), scale(_scale) {}

	std::optional<LightLiSample> SampleLi(LightSampleContext sample, vec2 u, SampledWaveLengths lambda, bool allowIncompletePDF = false) const override {
		vec3 p = renderFromLight(vec3(0.0));
		vec3 wi = normalize(p - sample.p);
		SampledSpectrum Li = intensity.Sample(lambda) * scale; // TODO: fix this, need add / (p - sample.p).lengthSquared();
		return LightLiSample(Li, wi, 1.0, Interaction(p));
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
		return LightLiSample(Li, wi, 1.0, Interaction(pOutside));
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

class DiffuseAreaLight : public Light
{
public: 
	DiffuseAreaLight(const Transform& _renderFromLight, const Spectrum& _intensity, double _scale, shared_ptr<Shape> _shape) :
		Light(LightType::Area, _renderFromLight), intensity(_intensity), scale(_scale), shape(_shape), area(_shape->Area()) {}

	SampledSpectrum L(vec3 p, vec3 n, vec2 uv, vec3 w, const SampledWaveLengths& lambda) const override {
		return intensity.Sample(lambda) * scale;
	}

	std::optional<LightLiSample> SampleLi(LightSampleContext sample, vec2 u, SampledWaveLengths lambda, bool allowIncompletePDF = false) const override {
		ShapeSampleContext shapeSample(sample.p, sample.n, sample.ns, 0);
		std::optional<ShapeSample> ss = shape->Sample(shapeSample, u);
		if(!ss || ss->pdf == 0.0 || (ss->intr.p - sample.p).lengthSquared() == 0.0) return {};

		vec3 wi = normalize(ss->intr.p - sample.p);
		SampledSpectrum Li = L(ss->intr.p, ss->intr.n, ss->intr.uv, -wi, lambda);
		if (!Li) return {};
		return LightLiSample(Li, wi, ss->pdf, ss->intr);
	}

	double PDF_Li(LightSampleContext sample, vec3 wi, bool allowIncompletePDF = false) const override {
		ShapeSampleContext shapeSample(sample.p, sample.n, sample.ns, 0);
		return shape->PDF(shapeSample, wi);
	}

	SampledSpectrum Phi(SampledWaveLengths lambda) const override {
		return intensity.Sample(lambda) * scale * pi * area;
	}
private: 
	DenselySampledSpectrum intensity;
	double scale, area;
	shared_ptr<Shape> shape;
};
