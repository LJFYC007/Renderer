#pragma once
#include "math.h"
#include "textures.h"
#include "interaction.h"
#include "bsdf.h"
#include "bxdfs.h"

#include <memory>
using std::shared_ptr;

class Material
{
public:
    virtual ~Material() = default;
    virtual BSDF GetBSDF(SurfaceInteraction intr, const SampledWaveLengths& lambda) const = 0;
};

class DiffuseMaterial : public Material {
public:
    DiffuseMaterial(shared_ptr<SpectrumTexture> _albedo) : albedo(_albedo) {}
    BSDF GetBSDF(SurfaceInteraction intr, const SampledWaveLengths& lambda) const override {
        SampledSpectrum r = albedo->Evaluate(intr, lambda);
        return BSDF(intr.n, intr.shading.dpdu, make_shared<DiffuseBxDF>(r));
    }

private:
    shared_ptr<SpectrumTexture> albedo;
};


class ConductorMaterial : public Material {
public:
    ConductorMaterial(double _alphax, double _alphay, double _eta, double _k) : alphax(_alphax), alphay(_alphay), eta(_eta), k(_k) {}

    BSDF GetBSDF(SurfaceInteraction intr, const SampledWaveLengths& lambda) const override {
        return BSDF(intr.n, intr.shading.dpdu, make_shared<ConductorBxDF>(TrowbridgeReitzDistribution(alphax, alphay), SampledSpectrum(eta), SampledSpectrum(k)));
    }

private:
    double alphax, alphay, eta, k;
};

class DielectricMaterial : public Material {
public:
    DielectricMaterial(double _alphax, double _alphay, double _eta) : alphax(_alphax), alphay(_alphay), eta(_eta) {}

    BSDF GetBSDF(SurfaceInteraction intr, const SampledWaveLengths& lambda) const override {
        return BSDF(intr.n, intr.shading.dpdu, make_shared<DielectricBxDF>(TrowbridgeReitzDistribution(alphax, alphay), eta));
    }

private:
    double alphax, alphay, eta;
};
