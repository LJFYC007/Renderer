#pragma once
#include "bxdf.h"
#include "math.h"
#include "sample.h"

class DiffuseBxDF : public BxDF
{
public:
	DiffuseBxDF(SampledSpectrum _R) : R(_R) {}

	SampledSpectrum f(vec3 wo, vec3 wi) const override
	{
		if (!SameHemisphere(wo, wi))
			return SampledSpectrum(0.0);
		return R / pi;
	}

	std::optional<BSDFSample> Sample_f(vec3 wo, double uc, vec2 u, BxDFReflTransFlags sampleFlags = BxDFReflTransFlags::All) const override {
		if ((sampleFlags & BxDFReflTransFlags::Reflection) == BxDFReflTransFlags::Unset) return {};
		vec3 wi = SampleCosineHemisphere(u);
		if (wo.z() < 0.0) wi.a[2] *= -1;
		double pdf = CosineHemispherePDF(abs(wi.z()));
		return BSDFSample(R / pi, wi, pdf, BxDFFlags::DiffuseReflection);
	}

	double PDF(vec3 wo, vec3 wi, BxDFReflTransFlags sampleFlags = BxDFReflTransFlags::All) const override {
		return CosineHemispherePDF(abs(wi.z()));
	}

private:
	SampledSpectrum R;
};
