#pragma once
#include "math.h"
#include "textures.h"
#include "interaction.h"
#include "bsdf.h"
#include "bxdfs.h"

#include <memory>
using std::shared_ptr;
using std::make_shared;

struct MaterialEvalContext : public TextureEvalContext
{
    MaterialEvalContext(const SurfaceInteraction& intr) : TextureEvalContext(intr), wo(intr.wo), ns(intr.shading.n), dpdus(intr.shading.dpdu) {}
    vec3 wo, ns, dpdus;
};

struct NormalBumpEvalContext {
    NormalBumpEvalContext(const SurfaceInteraction& intr)
        : p(intr.p()), uv(intr.uv), UV(intr.UV), n(intr.n), dudx(intr.dudx), dudy(intr.dudy), dvdx(intr.dvdx), dvdy(intr.dvdy), dpdx(intr.dpdx), dpdy(intr.dpdy) {
        shading.n = intr.shading.n;
        shading.dpdu = intr.shading.dpdu;
        shading.dpdv = intr.shading.dpdv;
        shading.dpdU = intr.shading.dpdU;
        shading.dpdV = intr.shading.dpdV;
        shading.dndu = intr.shading.dndu;
        shading.dndv = intr.shading.dndv;
    }
    operator TextureEvalContext() const {
        return TextureEvalContext(p, dpdx, dpdy, n, uv, UV, dudx, dudy, dvdx, dvdy);
    }

    vec3 p, n; vec2 uv, UV;
    struct {
        vec3 n, dpdu, dpdv, dpdU, dpdV, dndu, dndv;
    } shading;
    double dudx = 0, dudy = 0, dvdx = 0, dvdy = 0;
    vec3 dpdx, dpdy;
};

inline void NormalMap(shared_ptr<SpectrumTexture> normalMap, const NormalBumpEvalContext& ctx, vec3* dpdu, vec3* dpdv, vec3* dpdU, vec3* dpdV) {
    vec3 ns = 2 * normalMap->Evaluate(ctx) - vec3(1.0);
    ns = normalize(ns);
    Frame frame = FromXZ(normalize(ctx.shading.dpdu), ctx.shading.n);
    ns = frame.FromLocal(ns);
    double ulen = ctx.shading.dpdu.length(), vlen = ctx.shading.dpdv.length();
    *dpdu = normalize(GramSchmidt(ctx.shading.dpdu, ns)) * ulen;
    *dpdv = normalize(cross(ns, *dpdu)) * vlen;

    frame = FromXZ(normalize(ctx.shading.dpdU), ctx.shading.n);
    ns = frame.FromLocal(ns);
    ulen = ctx.shading.dpdU.length(), vlen = ctx.shading.dpdV.length();
    *dpdU = normalize(GramSchmidt(ctx.shading.dpdU, ns)) * ulen;
    *dpdV = normalize(cross(ns, *dpdU)) * vlen;
}

class Material
{
public:
    Material() { normalMap = bumpMap = nullptr; }
    virtual bool IsMixMaterial() { return false; }
    virtual shared_ptr<Material> GetMaterial(MaterialEvalContext ctx, double u) {
        std::cerr << "Should not call GetMaterial on non-mix material";
        return nullptr;
    }
    void SetNormalMap(shared_ptr<SpectrumTexture> _normalMap) { normalMap = _normalMap; }
    void SetBumpMap(shared_ptr<SpectrumTexture> _bumpMap) { bumpMap = _bumpMap; }
    shared_ptr<SpectrumTexture> GetNormalMap() const { return normalMap; }
    shared_ptr<SpectrumTexture> GetBumpMap() const { return bumpMap; }

    virtual shared_ptr<BxDF> GetBxDF(const MaterialEvalContext& ctx, SampledWaveLengths& lambda) const = 0;

    BSDF GetBSDF(MaterialEvalContext ctx, SampledWaveLengths& lambda) {
        shared_ptr<BxDF> bxdf = GetBxDF(ctx, lambda);
        return BSDF(ctx.ns, ctx.dpdus, bxdf);
    }

private:
    shared_ptr<SpectrumTexture> normalMap = nullptr;
    shared_ptr<SpectrumTexture> bumpMap = nullptr;
};

class DiffuseMaterial : public Material {
public:
    DiffuseMaterial(shared_ptr<SpectrumTexture> _albedo) : albedo(_albedo) {}

    shared_ptr<BxDF> GetBxDF(const MaterialEvalContext& ctx, SampledWaveLengths& lambda) const override {
        SampledSpectrum reflectance = albedo->Evaluate(ctx, lambda);
        return make_shared<DiffuseBxDF>(reflectance);
    }

private:
    shared_ptr<SpectrumTexture> albedo;
};


class ConductorMaterial : public Material {
public:
    ConductorMaterial(shared_ptr<SpectrumTexture> _reflectance, shared_ptr<SpectrumTexture> _metallicRoughness) : reflectance(_reflectance), metallicRoughness(_metallicRoughness) {}
    ConductorMaterial(shared_ptr<Spectrum> _eta, shared_ptr<Spectrum> _k, shared_ptr<SpectrumTexture> _metallicRoughness) : eta(_eta), k(_k), metallicRoughness(_metallicRoughness) {}

    shared_ptr<BxDF> GetBxDF(const MaterialEvalContext& ctx, SampledWaveLengths& lambda) const override {
        vec3 x = metallicRoughness->Evaluate(ctx);
        double alpha = RoughnessToAlpha(x[1]);
        SampledSpectrum r = reflectance->Evaluate(ctx, lambda);
        SampledSpectrum etas, ks;
        if (!eta) {
            etas = SampledSpectrum(1.0);
            ks = Sqrt(r) / Sqrt(ClampZero(SampledSpectrum(1.0) - r)) * 2;
        }
        else {
            etas = eta->Sample(lambda);
            ks = k->Sample(lambda);
        }
        return make_shared<ConductorBxDF>(TrowbridgeReitzDistribution(alpha, alpha), etas, ks);
    }

private:
    shared_ptr<SpectrumTexture> reflectance;
    shared_ptr<SpectrumTexture> metallicRoughness;
    shared_ptr<Spectrum> eta, k;
};

class DielectricMaterial : public Material {
public:
    DielectricMaterial(shared_ptr<SpectrumTexture> _metallicRoughness, shared_ptr<Spectrum> _eta) 
        : metallicRoughness(_metallicRoughness), eta(_eta) {}

    shared_ptr<BxDF> GetBxDF(const MaterialEvalContext& ctx, SampledWaveLengths& lambda) const override {
        vec3 x = metallicRoughness->Evaluate(ctx);
        double alpha = RoughnessToAlpha(x[1]);
        double sampledEta = eta->operator()(lambda[0]);
        if (!eta->IsConstant()) lambda.Terminate();
        if (sampledEta == 0) sampledEta = 1;
        return make_shared<DielectricBxDF>(TrowbridgeReitzDistribution(alpha, alpha), sampledEta);
    }

private:
    shared_ptr<SpectrumTexture> metallicRoughness;
    shared_ptr<Spectrum> eta;
};

class MixMaterial : public Material {
public:
    MixMaterial(shared_ptr<Material> _m1, shared_ptr<Material> _m2, shared_ptr<SpectrumTexture> _metallicRoughness) 
        : materials{ _m1, _m2 }, metallicRoughness(_metallicRoughness) {}
	bool IsMixMaterial() override { return true; }

    shared_ptr<Material> GetMaterial(MaterialEvalContext ctx, double u) override {
        if (metallicRoughness != nullptr) {
			vec3 x = metallicRoughness->Evaluate(ctx);
            double metallic = x[2];
			if (u < metallic ) return materials[0];
			else return materials[1];
        }
		double cosTheta = std::abs(dot(ctx.ns, ctx.wo));
		double f0 = 0.04;
		double fr = f0 + (1 - f0) * pow(1 - cosTheta, 5);
		if (u < 1 - fr) return materials[0];
		else return materials[1];
    }

    shared_ptr<BxDF> GetBxDF(const MaterialEvalContext& ctx, SampledWaveLengths& lambda) const override {
        std::cerr << "Should not call GetBxDF on MixMaterial";
        return nullptr;
    }

private:
    shared_ptr<Material> materials[2];
    shared_ptr<SpectrumTexture> metallicRoughness;
};
