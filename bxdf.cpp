#include "bxdf.h"
#include "sample.h"

SampledSpectrum BxDF::rho(vec3 wo, int nSamples, double* uc, vec2* u2) const
{
    SampledSpectrum r(0.0);
    for (size_t i = 0; i < nSamples; ++i) {
        std::optional<BSDFSample> sample = Sample_f(wo, uc[i], u2[i]);
        if (sample) r = r + sample->f * AbsCosTheta(sample->wi) / sample->pdf;
    }
    return r / nSamples;
}

SampledSpectrum BxDF::rho(int nSamples, vec2* u1, double* uc, vec2* u2) const
{
    SampledSpectrum r(0.0);
    for (size_t i = 0; i < nSamples; ++i) {
        vec3 wo = SampleUniformHemisphere(u1[i]);
        if (wo.z() == 0) continue;
        double pdfo = UniformHemispherePDF();
        std::optional<BSDFSample> sample = Sample_f(wo, uc[i], u2[i]);
        if (sample) r = r + sample->f * AbsCosTheta(sample->wi) * AbsCosTheta(wo) / (pdfo * sample->pdf);
    }
    return r / (pi * nSamples);
}
