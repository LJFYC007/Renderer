#pragma once
#include "spectrum.h"
#include "math.h"

#include <optional>

enum class BxDFReflTransFlags {
	Unset = 0,
	Reflection = 1 << 0,
	Transmission = 1 << 1,
	All = Reflection | Transmission
};

enum BxDFFlags {
	Unset = 0,
	Reflection = 1 << 0,
	Transmission = 1 << 1,
	Diffuse = 1 << 2,
	Glossy = 1 << 3,
	Specular = 1 << 4,
	DiffuseReflection = Diffuse | Reflection,
	DiffuseTransmission = Diffuse | Transmission,
	GlossyReflection = Glossy | Reflection,
	GlossyTransmission = Glossy | Transmission,
	SpecularReflection = Specular | Reflection,
	SpecularTransmission = Specular | Transmission,
	All = Diffuse | Glossy | Specular | Reflection | Transmission
};

struct BSDFSample {
	SampledSpectrum f;
	vec3 wi;
	double pdf = 0;
	BxDFFlags flags;
	BSDFSample(SampledSpectrum _f, vec3 _wi, double _pdf, BxDFFlags _flags) : f(_f), wi(_wi), pdf(_pdf), flags(_flags) {}
};

class BxDF
{
public : 
	SampledSpectrum f(vec3 wo, vec3 wi) const;
	std::optional<BSDFSample> Sample_f(vec3 wo, double uc, vec2 u, BxDFReflTransFlags sampleFlags = BxDFReflTransFlags::All) const;
	double PDF(vec3 wo, vec3 wi, BxDFReflTransFlags sampleFlags = BxDFReflTransFlags::All) const;
	SampledSpectrum rho(vec3 wo, int nSamples, double* uc, vec2* u2) const;
	SampledSpectrum rho(int nSamples, vec2* u1, double* uc, vec2* u2) const;
};
