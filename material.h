#pragma once
#include "math.h"
#include "texture.h"
#include "bsdf.h"
#include "bxdfs.h"

#include <memory>
using std::shared_ptr;

class Material
{
public:
    virtual ~Material() = default;
    virtual BSDF GetBSDF(vec3 n, vec3 dpduv, const SampledWaveLengths& sample) const = 0;
};

class DiffuseMaterial : public Material {
public:
    DiffuseMaterial(const vec3& a) : albedo(make_shared<solidColor>(a)) {}
    DiffuseMaterial(shared_ptr<texture> _albedo) : albedo(_albedo) {}
    BSDF GetBSDF(vec3 n, vec3 dpduv, const SampledWaveLengths& sample) const override {
        SampledSpectrum r = albedo->value(0, 0, vec3(0)).Sample(sample);
        return BSDF(n, dpduv, make_shared<DiffuseBxDF>(r));
    }

private:
    shared_ptr<texture> albedo;
};


class ConductorMaterial : public Material {
public:
    ConductorMaterial(double _alphax, double _alphay, double _eta, double _k) : alphax(_alphax), alphay(_alphay), eta(_eta), k(_k) {}

    BSDF GetBSDF(vec3 n, vec3 dpduv, const SampledWaveLengths& sample) const override {
        return BSDF(n, dpduv, make_shared<ConductorBxDF>(TrowbridgeReitzDistribution(alphax, alphay), SampledSpectrum(eta), SampledSpectrum(k)));
    }

private:
    double alphax, alphay, eta, k;
};
