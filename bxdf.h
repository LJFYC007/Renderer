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
inline BxDFReflTransFlags operator |(BxDFReflTransFlags a, BxDFReflTransFlags b) { return BxDFReflTransFlags((int)a | (int)b); }
inline BxDFReflTransFlags operator &(BxDFReflTransFlags a, BxDFReflTransFlags b) { return BxDFReflTransFlags((int)a & (int)b); }

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
inline BxDFFlags operator |(BxDFFlags a, BxDFFlags b) { return BxDFFlags((int)a | (int)b); }
inline BxDFFlags operator &(BxDFFlags a, BxDFFlags b) { return BxDFFlags((int)a & (int)b); }
inline bool IsNonSpecular(BxDFFlags f) { return f & (BxDFFlags::Diffuse | BxDFFlags::Glossy); }

struct BSDFSample {
	SampledSpectrum f;
	vec3 wi;
	double pdf = 0;
	BxDFFlags flags;
	double eta = 1;
	BSDFSample(SampledSpectrum _f, vec3 _wi, double _pdf, BxDFFlags _flags, double _eta = 1) : f(_f), wi(_wi), pdf(_pdf), flags(_flags), eta(_eta) {}

	bool IsTransmission() const { return flags & BxDFFlags::Transmission; }
};

class BxDF
{
public : 
	virtual SampledSpectrum f(vec3 wo, vec3 wi) const = 0;
	virtual std::optional<BSDFSample> Sample_f(vec3 wo, double uc, vec2 u, BxDFReflTransFlags sampleFlags = BxDFReflTransFlags::All) const = 0;
	virtual double PDF(vec3 wo, vec3 wi, BxDFReflTransFlags sampleFlags = BxDFReflTransFlags::All) const = 0;
	virtual SampledSpectrum rho(vec3 wo, int nSamples, double* uc, vec2* u2) const;
	virtual SampledSpectrum rho(int nSamples, vec2* u1, double* uc, vec2* u2) const;
	virtual BxDFFlags Flags() const = 0;
};
