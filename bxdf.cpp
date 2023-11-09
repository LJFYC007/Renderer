#include "bxdf.h"
#include "sample.h"

SampledSpectrum BxDF::f(vec3 wo, vec3 wi) const
{
    return SampledSpectrum(0.0);
}

std::optional<BSDFSample> BxDF::Sample_f(vec3 wo, double uc, vec2 u, BxDFReflTransFlags sampleFlags) const
{
    return std::optional<BSDFSample>();
}

double BxDF::PDF(vec3 wo, vec3 wi, BxDFReflTransFlags sampleFlags) const
{
    return 0.0;
}

SampledSpectrum BxDF::rho(vec3 wo, int nSamples, double* uc, vec2* u2) const
{
    SampledSpectrum r(0.0);
    for (size_t i = 0; i < nSamples; ++i) {
        std::optional<BSDFSample> sample = Sample_f(wo, uc[i], u2[i]);
        if (sample) r = r + sample->f * std::abs(sample->wi.z()) / sample->pdf;
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
        if (sample) r = r + sample->f * std::abs(sample->wi.z()) * std::abs(wo.z()) / (pdfo * sample->pdf);
    }
    return r / (pi * nSamples);
}
