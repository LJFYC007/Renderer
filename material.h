#pragma once
#include "math.h"
#include "textures.h"
#include "interaction.h"
#include "bsdf.h"
#include "bxdfs.h"

#include <memory>
using std::shared_ptr;

struct NormalBumpEvalContext {
    NormalBumpEvalContext(const SurfaceInteraction &intr) 
        : p(intr.p()), uv(intr.uv), n(intr.n), dudx(intr.dudx), dudy(intr.dudy), dvdx(intr.dvdx), dvdy(intr.dvdy), dpdx(intr.dpdx), dpdy(intr.dpdy) {
        shading.n = intr.shading.n;
        shading.dpdu = intr.shading.dpdu;
        shading.dpdv = intr.shading.dpdv;
        shading.dndu = intr.shading.dndu;
        shading.dndv = intr.shading.dndv;
    }
    operator TextureEvalContext() const {
        return TextureEvalContext(p, dpdx, dpdy, n, uv, dudx, dudy, dvdx, dvdy);
    }

    vec3 p, n; vec2 uv;
    struct {
        vec3 n, dpdu, dpdv, dndu, dndv;
    } shading;
    double dudx = 0, dudy = 0, dvdx = 0, dvdy = 0;
    vec3 dpdx, dpdy;
};

inline void NormalMap(shared_ptr<ImageTexture> normalMap, const NormalBumpEvalContext& ctx, vec3* dpdu, vec3* dpdv) {
    vec3 ns = 2 * normalMap->Evaluate(ctx) - vec3(1.0);
    ns = normalize(ns);
    Frame frame = FromXZ(normalize(ctx.shading.dpdu), ctx.shading.n);
    ns = frame.FromLocal(ns);
    double ulen = ctx.shading.dpdu.length(), vlen = ctx.shading.dpdv.length();
    *dpdu = normalize(GramSchmidt(ctx.shading.dpdu, ns)) * ulen;
    *dpdv = normalize(cross(ns, *dpdu)) * vlen;
}

class Material
{
public:
    virtual ~Material() = default;
    virtual BSDF GetBSDF(SurfaceInteraction& intr, const SampledWaveLengths& lambda) const = 0;
};

class DiffuseMaterial : public Material {
public:
    DiffuseMaterial(shared_ptr<SpectrumTexture> _albedo) : albedo(_albedo) {}
    BSDF GetBSDF(SurfaceInteraction& intr, const SampledWaveLengths& lambda) const override {
        SampledSpectrum r = albedo->Evaluate(intr, lambda);
        if (normalMap) {
            vec3 dpdu, dpdv;
            NormalMap(normalMap, NormalBumpEvalContext(intr), &dpdu, &dpdv);
            vec3 ns = normalize(cross(dpdu, dpdv));
            intr.SetShadingGeometry(ns, dpdu, dpdv, intr.shading.dndu, intr.shading.dndv);
        }
        return BSDF(intr.shading.n, intr.shading.dpdu, make_shared<DiffuseBxDF>(r));
    }

    void SetNormalMap(shared_ptr<ImageTexture> _normalMap) {
        normalMap = _normalMap;
    }

private:
    shared_ptr<SpectrumTexture> albedo;
    shared_ptr<ImageTexture> normalMap;
};


class ConductorMaterial : public Material {
public:
    ConductorMaterial(double _alphax, double _alphay, double _eta, double _k) : alphax(_alphax), alphay(_alphay), eta(_eta), k(_k) {}

    BSDF GetBSDF(SurfaceInteraction& intr, const SampledWaveLengths& lambda) const override {
        return BSDF(intr.shading.n, intr.shading.dpdu, make_shared<ConductorBxDF>(TrowbridgeReitzDistribution(alphax, alphay), SampledSpectrum(eta), SampledSpectrum(k)));
    }

private:
    double alphax, alphay, eta, k;
};

class DielectricMaterial : public Material {
public:
    DielectricMaterial(double _alphax, double _alphay, double _eta) : alphax(_alphax), alphay(_alphay), eta(_eta) {}

    BSDF GetBSDF(SurfaceInteraction& intr, const SampledWaveLengths& lambda) const override {
        return BSDF(intr.shading.n, intr.shading.dpdu, make_shared<DielectricBxDF>(TrowbridgeReitzDistribution(alphax, alphay), eta));
    }

private:
    double alphax, alphay, eta;
};
