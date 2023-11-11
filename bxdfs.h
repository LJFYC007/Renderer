#pragma once
#include "bxdf.h"
#include "math.h"
#include "sample.h"

class DiffuseBxDF : public BxDF
{
public:
	DiffuseBxDF() = default;
	DiffuseBxDF(SampledSpectrum _R) : R(_R) {}

	SampledSpectrum f(vec3 wo, vec3 wi) const override {
		if (!SameHemisphere(wo, wi)) return SampledSpectrum(0.0);
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
		if ((sampleFlags & BxDFReflTransFlags::Reflection) == BxDFReflTransFlags::Unset) return 0.0;
		if (!SameHemisphere(wo, wi)) return 0.0;
		return CosineHemispherePDF(AbsCosTheta(wi));
	}

private:
	SampledSpectrum R;
};

SampledSpectrum FrComplex(double cosThetai, SampledSpectrum eta, SampledSpectrum k) {
	SampledSpectrum result;
	for (int i = 0; i < NSpectrumSamples; ++i)
		result[i] = FrComplex(cosThetai, complex(eta[i], k[i]));
	return result;
}

class TrowbridgeReitzDistribution
{
public: 
	TrowbridgeReitzDistribution() = default;
	TrowbridgeReitzDistribution(double _alphax, double _alphay) : alphax(_alphax), alphay(_alphay) {}

	bool EffectivelySmooth() const { return std::fmax(alphax, alphay) < 1e-3f; }

	double D(vec3 wm) const {
	}

private:
	double alphax, alphay;
};

class ConductorBxDF : public BxDF
{
public:
	ConductorBxDF() = default;
	ConductorBxDF(const TrowbridgeReitzDistribution& _mfDistrib, SampledSpectrum _eta, SampledSpectrum _k) : mfDistrib(_mfDistrib), eta(_eta), k(_k) {}

	BxDFFlags Flags() const { return mfDistrib.EffectivelySmooth() ? BxDFFlags::SpecularReflection : BxDFFlags::GlossyReflection; }

	SampledSpectrum f(vec3 wo, vec3 wi) const override {
		if (!SameHemisphere(wo, wi)) return SampledSpectrum(0.0);
		if (mfDistrib.EffectivelySmooth()) return {};
		return {};
	}

	std::optional<BSDFSample> Sample_f(vec3 wo, double uc, vec2 u, BxDFReflTransFlags sampleFlags = BxDFReflTransFlags::All) const override {
		if ((sampleFlags & BxDFReflTransFlags::Reflection) == BxDFReflTransFlags::Unset) return {};
		if (mfDistrib.EffectivelySmooth()) {
			vec3 wi(-wo.x(), -wo.y(), wo.z());
			SampledSpectrum f = FrComplex(AbsCosTheta(wi), eta, k) / AbsCosTheta(wi);
			return BSDFSample(f, wi, 1, BxDFFlags::SpecularReflection);
		}
		return {};
	}

	double PDF(vec3 wo, vec3 wi, BxDFReflTransFlags sampleFlags = BxDFReflTransFlags::All) const override {
		if ((sampleFlags & BxDFReflTransFlags::Reflection) == BxDFReflTransFlags::Unset || !SameHemisphere(wo, wi)) return 0.0;
		if (!SameHemisphere(wo, wi)) return 0.0;
		if (mfDistrib.EffectivelySmooth()) return 0.0;
		return 0.0;
	}

private:
	TrowbridgeReitzDistribution mfDistrib;
	SampledSpectrum eta, k;
};
