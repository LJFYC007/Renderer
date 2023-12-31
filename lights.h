#pragma once
#include "light.h"
#include "math.h"
#include "shape.h"
#include "image.h"

class PointLight : public Light
{
public : 
	PointLight(const Transform &_renderFromLight, const Spectrum &_intensity, double _scale) 
		: Light(LightType::DeltaPosition, _renderFromLight), intensity(_intensity), scale(_scale) {}

	std::optional<LightLiSample> SampleLi(LightSampleContext sample, vec2 u, SampledWaveLengths lambda, bool allowIncompletePDF) const override {
		vec3 p = renderFromLight(vec3(0.0));
		vec3 wi = normalize(p - sample.p());
		SampledSpectrum Li = intensity.Sample(lambda) * scale; // TODO: fix this, need add / (p - sample.p).lengthSquared();
		return LightLiSample(Li, wi, 1.0, Interaction(p));
	}

	double PDF_Li(LightSampleContext sample, vec3 wi, bool allowIncompletePDF) const override {
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

	std::optional<LightLiSample> SampleLi(LightSampleContext sample, vec2 u, SampledWaveLengths lambda, bool allowIncompletePDF) const override {
		vec3 wi = normalize(renderFromLight(vec3(0, 0, 1)));
		vec3 pOutside = sample.p() + wi * 100000.0; // TODO: fix this, need to multiply R
		SampledSpectrum Li = intensity.Sample(lambda) * scale;
		return LightLiSample(Li, wi, 1.0, Interaction(pOutside));
	}

	double PDF_Li(LightSampleContext sample, vec3 wi, bool allowIncompletePDF) const override {
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

	std::optional<LightLiSample> SampleLi(LightSampleContext sample, vec2 u, SampledWaveLengths lambda, bool allowIncompletePDF) const override {
		ShapeSampleContext shapeSample(sample.p(), sample.n, sample.ns, 0);
		std::optional<ShapeSample> ss = shape->Sample(shapeSample, u);
		if(!ss || ss->pdf == 0.0 || (ss->intr.p() - sample.p()).lengthSquared() == 0.0) return {};

		vec3 wi = normalize(ss->intr.p() - sample.p());
		SampledSpectrum Li = L(ss->intr.p(), ss->intr.n, ss->intr.uv, -wi, lambda);
		if (!Li) return {};
		return LightLiSample(Li, wi, ss->pdf, ss->intr);
	}

	double PDF_Li(LightSampleContext sample, vec3 wi, bool allowIncompletePDF) const override {
		ShapeSampleContext shapeSample(sample.p(), sample.n, sample.ns, 0);
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

class UniformInfiniteLight : public Light
{
public:
	UniformInfiniteLight(const Transform& _renderFromLight, const Spectrum& _intensity, double _scale) :
		Light(LightType::Infinite, _renderFromLight), intensity(_intensity), scale(_scale) {}

	SampledSpectrum Le(const Ray& ray, const SampledWaveLengths& lambda) const override {
		return intensity.Sample(lambda) * scale;
	}

	std::optional<LightLiSample> SampleLi(LightSampleContext sample, vec2 u, SampledWaveLengths lambda, bool allowIncompletePDF) const override {
		if (allowIncompletePDF) return {};
		vec3 wi = SampleUniformSphere(u);
		double pdf = UniformSpherePDF();
		vec3 pOutside = sample.p() + wi * 100000.0; // TODO: fix this, need to multiply R
		SampledSpectrum Li = intensity.Sample(lambda) * scale;
		return LightLiSample(Li, wi, pdf, Interaction(pOutside));
	}

	double PDF_Li(LightSampleContext sample, vec3 wi, bool allowIncompletePDF) const override {
		if (allowIncompletePDF) return 0.0;
		return UniformSpherePDF();
	}

	SampledSpectrum Phi(SampledWaveLengths lambda) const override {
		return intensity.Sample(lambda) * scale * 4 * pi * pi; // TODO: fix this, need to multiply R * R
	}

private:
	DenselySampledSpectrum intensity;
	double scale;
	vec3 sceneCenter, sceneRadius;
};

class ImageInfiniteLight : public Light
{
public:
	ImageInfiniteLight(const Transform& _renderFromLight, double _scale, std::string filename) :
		Light(LightType::Infinite, _renderFromLight), scale(_scale), image(make_shared<HDRImage>(filename)){}

	SampledSpectrum ImageLe(vec2 uv, const SampledWaveLengths& lambda) const {
		return image->LookUp(uv, lambda, false) * scale;
	}

	SampledSpectrum Le(const Ray& ray, const SampledWaveLengths& lambda) const override {
		vec3 wLight = normalize(renderFromLight.ApplyInverse(ray.rd));	
		vec2 uv(SphericalPhi(wLight) / 2 / pi, SphericalTheta(wLight) / pi);
		return ImageLe(uv, lambda);
	}

	std::optional<LightLiSample> SampleLi(LightSampleContext sample, vec2 u, SampledWaveLengths lambda, bool allowIncompletePDF) const override {
		double mapPDF = 1.0;// / (image->width * image->height);
		vec2 uv = u;
		double theta = uv[1] * pi, phi = uv[0] * 2 * pi;
		double cosTheta = std::cos(theta), sinTheta = std::sin(theta);
		double sinPhi = std::sin(phi), cosPhi = std::cos(phi);
		vec3 wLight(sinTheta * cosPhi, sinTheta * sinPhi, cosTheta);
		
		vec3 wi = renderFromLight(wLight);
		double pdf = mapPDF / (2 * pi * pi * sinTheta);
		vec3 pOutside = sample.p() + wi * 100000.0; // TODO: fix this, need to multiply R
		return LightLiSample(ImageLe(uv, lambda), wi, pdf, Interaction(pOutside));
	}

	double PDF_Li(LightSampleContext sample, vec3 wi, bool allowIncompletePDF) const override {
		vec3 wLight = normalize(renderFromLight.ApplyInverse(wi));
		double theta = SphericalTheta(wLight), phi = SphericalPhi(wLight);
		double sinTheta = std::sin(theta);
		if (sinTheta == 0) return 0;
		double pdf = 1.0;// / (image->width * image->height);
		return pdf / (2 * pi * pi * sinTheta);
	}

	SampledSpectrum Phi(SampledWaveLengths lambda) const override {
		return scale * 4 * pi * pi; // TODO: fix this, need to multiply R * R
	}

private:
	shared_ptr<HDRImage> image;
	double scale;
	vec3 sceneCenter, sceneRadius;
};
