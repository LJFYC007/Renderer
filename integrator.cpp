#include "integrator.h"
#include "camera.h"

SampledSpectrum SampleLd(const SurfaceInteraction& intr, const BSDF& bsdf, Camera* camera, const SampledWaveLengths& lambda, const BVHAggregate& bvh)
{
	LightSampleContext ctx(intr);
	BxDFFlags flags = bsdf.Flags();
	if (IsReflective(flags) && !IsTransmissive(flags)) ctx.pi = intr.OffsetRayOrigin(intr.wo);
	else if (IsTransmissive(flags) && !IsReflective(flags)) ctx.pi = intr.OffsetRayOrigin(-intr.wo);

	std::optional<SampledLight> sampledLight = camera->lightSampler->Sample(randomDouble());
	if (!sampledLight) return {};
	shared_ptr<Light> light = sampledLight->light;
	vec2 u = vec2Random();
	std::optional<LightLiSample> ls = light->SampleLi(ctx, u, lambda, true);
	if (!ls || !ls->L || ls->pdf == 0.0) return { 0.0 };

	vec3 wo = intr.wo, wi = ls->wi;
	SampledSpectrum f = bsdf.f(wo, wi) * std::abs(dot(wi, intr.shading.n));
	if (!f || !Unoccluded(bvh, intr, ls->pLight)) return { 0.0 };

	double p_l = sampledLight->p * ls->pdf;
	if (IsDeltaLight(light->Type()))
		return ls->L * f / p_l;
	double p_b = bsdf.PDF(wo, wi);
	double w_l = PowerHeuristic(1, p_l, 1, p_b);
	return ls->L * w_l * f / p_l;
}

SampledSpectrum Li(Ray r, Camera* camera, SampledWaveLengths& lambda, const BVHAggregate& bvh)
{
	SampledSpectrum L(0.0), beta(1.0);
	int depth = 0;
	double p_b = 1.0, eta_scale = 1.0;
	bool specularBounce = false;
	LightSampleContext prevIntrCtx;

	while (beta) {
		r.rd = normalize(r.rd);
		std::optional<ShapeIntersection> isect = bvh.Intersect(r, interval(0, infinity));
		if (!isect)
		{
			for (const auto& light : camera->infiniteLights)
			{
				SampledSpectrum Le = light->Le(r, lambda);
				if (depth == 0 || specularBounce)
					L = L + beta * Le;
				else {
					double p_l = camera->lightSampler->PMF(prevIntrCtx, light) * light->PDF_Li(prevIntrCtx, r.rd, true);
					double w_b = PowerHeuristic(1, p_b, 1, p_l);
					L = L + beta * Le * w_b;
				}
			}
			break;
		}

		SampledSpectrum Le = isect->intr.Le(-r.rd, lambda);
		if (Le) {
			if (depth == 0 || specularBounce)
				L = L + beta * Le;
			else {
				shared_ptr<Light> areaLight = isect->intr.areaLight;
				double p_l = camera->lightSampler->PMF(prevIntrCtx, areaLight) * areaLight->PDF_Li(prevIntrCtx, r.rd, true);
				double w_l = PowerHeuristic(1, p_b, 1, p_l);
				L = L + beta * Le * w_l;
			}
		}

		SurfaceInteraction intr = isect->intr;
		if (depth++ == camera->maxDepth) break;
		BSDF bsdf = intr.GetBSDF(r, lambda, camera);
		if (IsNonSpecular(bsdf.Flags())) {
			SampledSpectrum Ld = SampleLd(intr, bsdf, camera, lambda, bvh);
			L = L + beta * Ld;
		}

		std::optional<BSDFSample> bs = bsdf.Sample_f(-r.rd, randomDouble(), vec2Random());
		if (!bs) break;
		beta = beta * bs->f * std::abs(dot(bs->wi, intr.shading.n)) / bs->pdf;
		p_b = bs->pdf;
		specularBounce = bs->IsSpecular();

		if (bs->IsTransmission())
			eta_scale *= Sqr(bs->eta);
		prevIntrCtx = intr;
		r = intr.SpawnRay(r, bsdf, bs->wi, bs->flags, bs->eta);

		SampledSpectrum rrBeta = beta * eta_scale;
		if (rrBeta.MaxComponentValue() < 1.0 && depth > 1)
		{
			double q = std::fmax(0.0, 1 - rrBeta.MaxComponentValue());
			if (randomDouble() < q) break;
			beta = beta / (1 - q);
		}

		// -------------------- naive brdf ---------------------
		/*
		double pdf;
		vec3 wi;
		BxDFFlags flags = bsdf.Flags();
		if (IsReflective(flags) && IsTransmissive(flags)) {
		wi = SampleUniformSphere(vec2Random());
		pdf = UniformSpherePDF();
		}
		else {
		wi = SampleUniformHemisphere(vec2Random());
		pdf = UniformHemispherePDF();
		if (IsReflective(flags) && dot(wo, intr.n) * dot(wi, intr.n) < 0)
		wi = -wi;
		else if (IsTransmissive(flags) && dot(wo, intr.n) * dot(wi, intr.n) > 0)
		wi = -wi;
		}
		beta = beta * bsdf.f(wo, wi) * std::abs(dot(wi, intr.shading.n)) / pdf;
		r = Ray(intr.p, wi);
		*/
	}
	return L;
}
