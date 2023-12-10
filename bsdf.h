#pragma once
#include "bxdf.h"
#include "math.h"

#include <memory>
using std::shared_ptr;

class BSDF {
public:
	BSDF() = default;
	BSDF(vec3 ns, vec3 dpdus, shared_ptr<BxDF> _bxdf) : bxdf(_bxdf), shadingFrame(FromXZ(normalize(dpdus), ns)) {}

	vec3 RenderToLocal(vec3 v) const { return shadingFrame.ToLocal(v); }
	vec3 LocalToRender(vec3 v) const { return shadingFrame.FromLocal(v); }

	SampledSpectrum f(vec3 woRender, vec3 wiRender) const {
		vec3 wo = RenderToLocal(woRender), wi = RenderToLocal(wiRender);
		if (wo.z() == 0) return {};
		return bxdf->f(wo, wi);
	}

	std::optional<BSDFSample> Sample_f(vec3 woRender, double u, vec2 u2, BxDFReflTransFlags sampleFlags = BxDFReflTransFlags::All) const {
		vec3 wo = RenderToLocal(woRender);
		if (wo.z() == 0) return{};
		std::optional<BSDFSample> sample = bxdf->Sample_f(wo, u, u2, sampleFlags);
		if (!sample || !sample->f || sample->pdf == 0 || sample->wi.z() == 0) return {};
		sample->wi = LocalToRender(sample->wi);
		return sample;
	}

	double PDF(vec3 woRender, vec3 wiRender, BxDFReflTransFlags sampleFlags = BxDFReflTransFlags::All) const {
		vec3 wo = RenderToLocal(woRender), wi = RenderToLocal(wiRender);
		if (wo.z() == 0) return 0.0;
		return bxdf->PDF(wo, wi, sampleFlags);
	}

	SampledSpectrum rho(vec3 woRender, int nSamples, double* uc, vec2* u2) const {
		vec3 wo = RenderToLocal(woRender);
		return bxdf->rho(wo, nSamples, uc, u2);
	}

	SampledSpectrum rho(int nSamples, vec2* u1, double* uc, vec2* u2) const {
		return bxdf->rho(nSamples, u1, uc, u2);
	}

	BxDFFlags Flags() const { return bxdf->Flags(); }

private:
	shared_ptr<BxDF> bxdf;
	Frame shadingFrame;
};
