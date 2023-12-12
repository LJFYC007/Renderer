#pragma once
#include "bxdf.h"
#include "math.h"
#include "sample.h"

class DiffuseBxDF : public BxDF
{
public:
	DiffuseBxDF() = default;
	DiffuseBxDF(SampledSpectrum _R) : R(_R) {}

	BxDFFlags Flags() const override { return R ? BxDFFlags::DiffuseReflection : BxDFFlags::Unset; }

	SampledSpectrum f(vec3 wo, vec3 wi) const override {
		if (!SameHemisphere(wo, wi)) return	{};
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

inline SampledSpectrum FrComplex(double cosThetai, SampledSpectrum eta, SampledSpectrum k) {
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
		double tan2Theta = Tan2Theta(wm);
		if (tan2Theta > 1e16 || tan2Theta < -1e16) return 0;
		double cos4Theta = Sqr(Cos2Theta(wm));
		if (cos4Theta < 1e-16) return 0;
		double e = tan2Theta * (Sqr(CosPhi(wm) / alphax) + Sqr(SinPhi(wm) / alphay));
		return 1.0 / (pi * alphax * alphay * cos4Theta * Sqr(1 + e));
	}

	double Lambda(vec3 w) const {
		double tan2Theta = Tan2Theta(w);
		if (tan2Theta > 1e16 || tan2Theta < -1e16) return 0;
		double alpha2 = Sqr(CosPhi(w) * alphax) + Sqr(SinPhi(w) * alphay);
		return (std::sqrt(1 + alpha2 * tan2Theta) - 1) / 2.0;
	}

	double G1(vec3 w) const { return 1.0 / (1 + Lambda(w)); }
	double G(vec3 wo, vec3 wi) const { return 1.0 / (1 + Lambda(wo) + Lambda(wi)); }
	double D(vec3 w, vec3 wm) const { return G1(w) / AbsCosTheta(w) * D(wm) * std::abs(dot(w, wm)); }
	double PDF(vec3 w, vec3 wm) const { return D(w, wm); }

	vec3 Sample_wm(vec3 w, vec2 u) const {
		vec3 wh = normalize(vec3(alphax * w.x(), alphay * w.y(), w.z()));
		if (wh.z() < 0) wh = -wh;
		vec3 T1 = (wh.z() < 0.99999) ? normalize(cross(vec3(0, 0, 1), wh)) : vec3(1, 0, 0);
		vec3 T2 = cross(wh, T1);
		vec3 T3 = wh;

		vec2 p = SampleUniformDiskPolar(u);
		double h = std::sqrt(1 - p.x() * p.x());
		p.a[1] = Lerp((1 + wh.z() / 2), h, p.y()); // ºÃÇÉÃîµÄ±ä»»

		double pz = std::sqrt(std::fmax(0, 1 - p.lengthSquared()));
		vec3 nh = p.x() * T1 + p.y() * T2 + pz * T3;
		return normalize(vec3(alphax * nh.x(), alphay * nh.y(), std::fmax(1e-6, nh.z())));
	}

private:
	double alphax, alphay;
};

class ConductorBxDF : public BxDF
{
public:
	ConductorBxDF() = default;
	ConductorBxDF(const TrowbridgeReitzDistribution& _mfDistrib, SampledSpectrum _eta, SampledSpectrum _k) : mfDistrib(_mfDistrib), eta(_eta), k(_k) {}

	BxDFFlags Flags() const override { return mfDistrib.EffectivelySmooth() ? BxDFFlags::SpecularReflection : BxDFFlags::GlossyReflection; }

	SampledSpectrum f(vec3 wo, vec3 wi) const override {
		if (!SameHemisphere(wo, wi)) return {};
		if (mfDistrib.EffectivelySmooth()) return {};

		double cosThetao = AbsCosTheta(wo), cosThetai = AbsCosTheta(wi);
		if (cosThetai == 0 || cosThetao == 0) return {};
		vec3 wm = wi + wo;
		if (wm.lengthSquared() == 0) return {};
		wm = normalize(wm);
		SampledSpectrum F = FrComplex(std::abs(dot(wo, wm)), eta, k);
		return F * mfDistrib.D(wm) * mfDistrib.G(wo, wi) / (4 * cosThetai * cosThetao);
	}

	std::optional<BSDFSample> Sample_f(vec3 wo, double uc, vec2 u, BxDFReflTransFlags sampleFlags = BxDFReflTransFlags::All) const override {
		if ((sampleFlags & BxDFReflTransFlags::Reflection) == BxDFReflTransFlags::Unset) return {};
		if (mfDistrib.EffectivelySmooth()) {
			vec3 wi = vec3(-wo.x(), -wo.y(), wo.z());
			SampledSpectrum f = SampledSpectrum(FrComplex(AbsCosTheta(wi), eta, k) / AbsCosTheta(wi));
			return BSDFSample(f, wi, 1, BxDFFlags::SpecularReflection);
		}

		vec3 wm = mfDistrib.Sample_wm(wo, u);
		vec3 wi = Reflect(wo, wm);
		if (!SameHemisphere(wo, wi)) return {};
		double pdf = mfDistrib.PDF(wo, wm) / (4 * std::abs(dot(wo, wm)));
		double cosThetao = AbsCosTheta(wo), cosThetai = AbsCosTheta(wi);
		if (cosThetai == 0 || cosThetao == 0) return {};
		SampledSpectrum F = FrComplex(std::abs(dot(wo, wm)), eta, k);
		SampledSpectrum f = F * mfDistrib.D(wm) * mfDistrib.G(wo, wi) / (4 * cosThetai * cosThetao);
		return BSDFSample(f, wi, pdf, BxDFFlags::GlossyReflection);
	}

	double PDF(vec3 wo, vec3 wi, BxDFReflTransFlags sampleFlags = BxDFReflTransFlags::All) const override {
		if ((sampleFlags & BxDFReflTransFlags::Reflection) == BxDFReflTransFlags::Unset || !SameHemisphere(wo, wi)) return 0.0;
		if (!SameHemisphere(wo, wi)) return 0.0;
		if (mfDistrib.EffectivelySmooth()) return 0.0;
		vec3 wm = wo + wi;
		if (wm.lengthSquared() == 0) return 0.0;
		wm = normalize(wm);
		if (dot(wm, vec3(0, 0, 1)) < 0) wm = -wm;
		return mfDistrib.PDF(wo, wm) / (4 * std::abs(dot(wo, wm)));
	}

private:
	TrowbridgeReitzDistribution mfDistrib;
	SampledSpectrum eta, k;
};

class DielectricBxDF : public BxDF
{
public:
	DielectricBxDF() = default;
	DielectricBxDF(const TrowbridgeReitzDistribution& _mfDistrib, double _eta) : mfDistrib(_mfDistrib), eta(_eta) {}

	BxDFFlags Flags() const override { 
		BxDFFlags flags = (eta == 1) ? BxDFFlags::Transmission : (BxDFFlags::Transmission | BxDFFlags::Reflection);
		return flags | (mfDistrib.EffectivelySmooth() ? BxDFFlags::Specular : BxDFFlags::Glossy); 
	}

	SampledSpectrum f(vec3 wo, vec3 wi) const override {
		if (eta == 1 || mfDistrib.EffectivelySmooth()) return SampledSpectrum(0.0);

		double cosThetao = CosTheta(wo), cosThetai = CosTheta(wi);
		if (cosThetai == 0 || cosThetao == 0) return SampledSpectrum(0.0);
		bool reflect = cosThetai * cosThetao > 0;
		double etap = 1;
		if (!reflect) etap = cosThetao > 0 ? eta : (1 / eta);

		vec3 wm = wi * etap + wo;
		if (wm.lengthSquared() == 0.0 ) return SampledSpectrum(0.0);
		wm = normalize(wm);
		if (dot(wm, vec3(0, 0, 1)) < 0) wm = -wm;
		if (dot(wm, wi) * cosThetai < 0 || dot(wm, wo) * cosThetao < 0) return SampledSpectrum(0.0);

		double F = FrDielectric(dot(wo, wm), eta);
		if (reflect)
			return SampledSpectrum(mfDistrib.D(wm) * mfDistrib.G(wo, wi) * F / std::abs(4 * cosThetai * cosThetao));
		else {
			double denom = Sqr(dot(wi, wm) + dot(wo, wm) / etap) * cosThetai * cosThetao;
			return SampledSpectrum(mfDistrib.D(wm) * mfDistrib.G(wo, wi) * (1 - F) * std::abs(dot(wi, wm) * dot(wo, wm) / denom));
		}
	}

	std::optional<BSDFSample> Sample_f(vec3 wo, double uc, vec2 u, BxDFReflTransFlags sampleFlags = BxDFReflTransFlags::All) const override {
		if (eta == 1 || mfDistrib.EffectivelySmooth()) {
			double R = FrDielectric(CosTheta(wo), eta);

			double pr = R, pt = 1 - R;
			if ((sampleFlags & BxDFReflTransFlags::Reflection) == BxDFReflTransFlags::Unset) pr = 0;
			if ((sampleFlags & BxDFReflTransFlags::Transmission) == BxDFReflTransFlags::Unset) pt = 0;
			if (pr == 0 && pt == 0) return {};

			if (uc < pr / (pr + pt)) {
				vec3 wi = vec3(-wo.x(), -wo.y(), wo.z());
				SampledSpectrum fr = SampledSpectrum(R / AbsCosTheta(wi));
				return BSDFSample(fr, wi, pr / (pr + pt), BxDFFlags::SpecularReflection);
			}
			else {
				vec3 wi; double etap;
				if (!Refract(wo, vec3(0, 0, 1), eta, &etap, &wi)) return {};
				SampledSpectrum ft = SampledSpectrum((1 - R) / AbsCosTheta(wi));
				return BSDFSample(ft, wi, pt / (pr + pt), BxDFFlags::SpecularTransmission, etap);
			}
		}

		vec3 wm = mfDistrib.Sample_wm(wo, u);
		double R = FrDielectric(dot(wo, wm), eta);
		double pr = R, pt = 1 - R;
		if ((sampleFlags & BxDFReflTransFlags::Reflection) == BxDFReflTransFlags::Unset) pr = 0;
		if ((sampleFlags & BxDFReflTransFlags::Transmission) == BxDFReflTransFlags::Unset) pt = 0;
		if (pr == 0 && pt == 0) return {};

		if (uc < pr / (pr + pt)) {
			vec3 wi = Reflect(wo, wm);
			if (!SameHemisphere(wo, wi)) return {};
			double pdf =  mfDistrib.PDF(wo, wm) / (4 * std::abs(dot(wo, wm))) * pr / (pr + pt);
			SampledSpectrum f = SampledSpectrum(mfDistrib.D(wm) * mfDistrib.G(wo, wi) * R / (4 * CosTheta(wi) * CosTheta(wo)));
			return BSDFSample(f, wi, pdf, BxDFFlags::GlossyReflection);
		}
		else {
			double etap; vec3 wi;
			if (!Refract(wo, wm, eta, &etap, &wi)) return {};
			if (SameHemisphere(wo, wi) || wi.z() == 0) return {};
			double denom = Sqr(dot(wi, wm) + dot(wo, wm) / etap);
			double pdf = mfDistrib.PDF(wo, wm) * std::abs(dot(wi, wm)) / denom * pt / (pr + pt);
			SampledSpectrum f = SampledSpectrum(mfDistrib.D(wm) * mfDistrib.G(wo, wi) * (1 - R) * std::abs(dot(wi, wm) * dot(wo, wm) / (CosTheta(wi) * CosTheta(wo) * denom)));
			return BSDFSample(f, wi, pdf, BxDFFlags::GlossyTransmission, etap);
		}
	}

	double PDF(vec3 wo, vec3 wi, BxDFReflTransFlags sampleFlags = BxDFReflTransFlags::All) const override {
		if (eta == 1 || mfDistrib.EffectivelySmooth()) return 0.0;

		double cosThetao = CosTheta(wo), cosThetai = CosTheta(wi);
		if (cosThetai == 0 || cosThetao == 0) return 0.0;
		bool reflect = cosThetai * cosThetao > 0;
		double etap = 1;
		if (!reflect) etap = cosThetao > 0 ? eta : (1 / eta);

		vec3 wm = wi * etap + wo;
		if (wm.lengthSquared() == 0.0) return 0.0;
		wm = normalize(wm);
		if (dot(wm, vec3(0, 0, 1)) < 0) wm = -wm;
		if (dot(wm, wi) * cosThetai < 0 || dot(wm, wo) * cosThetao < 0) return 0.0;

		double R = FrDielectric(dot(wo, wm), eta);
		double pr = R, pt = 1 - R;
		if ((sampleFlags & BxDFReflTransFlags::Reflection) == BxDFReflTransFlags::Unset) pr = 0;
		if ((sampleFlags & BxDFReflTransFlags::Transmission) == BxDFReflTransFlags::Unset) pt = 0;
		if (pr == 0 && pt == 0) return 0.0;

		if (reflect) 
			return mfDistrib.PDF(wo, wm) / (4 * std::abs(dot(wo, wm))) * pr / (pr + pt);
		else {
			double denom = Sqr(dot(wi, wm) + dot(wo, wm) / etap);
			return mfDistrib.PDF(wo, wm) * std::abs(dot(wi, wm)) / denom * pt / (pr + pt);
		}
	}

private:
	TrowbridgeReitzDistribution mfDistrib;
	double eta;
};
