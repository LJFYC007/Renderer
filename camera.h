#pragma once
#include "bvh.h"
#include "sample.h"
#include "math.h"
#include "interval.h"
#include "spectrum.h"
#include "color.h"
#include "colorspace.h"
#include "bsdf.h"
#include "lights.h"
#include "lightsampler.h"

#include <Windows.h>
#include <omp.h>
#include <optional>

static vec3 ans[3010][2210];

class camera
{
public:
	const int ImageWidth = 1200;
	const int ImageHeight = 1200;
	double fov = 40.0;
	vec3 lookfrom = vec3(0.0, 0.0, -800.0);
	vec3 lookat = vec3(0.0, 0.0, 0.0);
	vec3 vup = vec3(0.0, 1.0, 0.0);
	double defocusAngle = 0.0;
	double focusDist = 10.0;
	int samplePixel = 1024;
	int maxDepth = 10;

	void render(const BVHAggregate& bvh, const std::vector<shared_ptr<Light>>& _lights)
	{
		initialize(_lights);
#pragma omp parallel for schedule(dynamic) 
		for (int j = 0; j < ImageHeight; ++j)
		{
			std::clog << "\rScanlines remaining: " << (ImageHeight - j) << ' ' << std::flush;
			for (int i = 0; i < ImageWidth; ++i)
			{
				XYZ xyz(0.0);
				SampledSpectrum spec(0.0);
				vec3 col = vec3(0.0);
				for (int oi = 0; oi < sqrtSPP; ++oi)
					for (int oj = 0; oj < sqrtSPP; ++oj)
					{
						vec3 pixelCenter = pixel00Location + i * pixelDeltaU + j * pixelDeltaV;
						vec3 pixel = pixelCenter + (-0.5 + (1.0 / sqrtSPP) * (oi + randomDouble())) * pixelDeltaU + (-0.5 + (1.0 / sqrtSPP) * (oj + randomDouble())) * pixelDeltaV;
						vec3 ro = (defocusAngle <= 0.0) ? cameraCenter : cameraCenter + defocusDiskSample(defocusDiskU, defocusDiskV);
						vec3 rd = pixel - ro;
						SampledWaveLengths sample(randomDouble());
						RGBColor rgb = sRGB.ToRGB((Li(ray(ro, rd), maxDepth, sample, bvh)).ToXYZ(sample));
						col = col + vec3(rgb.r, rgb.g, rgb.b) / samplePixel;
					}

				if (col.a[0] < 0) col.a[0] = 0;
				if (col.a[1] < 0) col.a[1] = 0;
				if (col.a[2] < 0) col.a[2] = 0;
				col.a[0] = pow(col.a[0], 1.0 / 2.2);
				col.a[1] = pow(col.a[1], 1.0 / 2.2);
				col.a[2] = pow(col.a[2], 1.0 / 2.2);

				static const interval intensity(0.000, 0.999);
				col.a[0] = static_cast<int>(256 * intensity.clamp(col.a[0]));
				col.a[1] = static_cast<int>(256 * intensity.clamp(col.a[1]));
				col.a[2] = static_cast<int>(256 * intensity.clamp(col.a[2]));
				ans[i][j] = col;
			}
		}

		std::cout << "P3\n" << ImageWidth << ' ' << ImageHeight << "\n255\n";
		for (int j = ImageHeight - 1; j >= 0; --j)
			for (int i = 0; i < ImageWidth; ++i)
				std::cout << ans[i][j].x() << ' ' << ans[i][j].y() << ' ' << ans[i][j].z() << '\n';
		std::clog << "\rDone.                 \n";
	}
private:
	int sqrtSPP;
	vec3 defocusDiskU;
	vec3 defocusDiskV;
	vec3 pixel00Location;
	vec3 cameraCenter;
	vec3 pixelDeltaU;
	vec3 pixelDeltaV;

	std::vector<shared_ptr<Light>> lights;
	std::vector<shared_ptr<Light>> infiniteLights;
	shared_ptr<LightSampler> lightSampler;

	void initialize(const std::vector<shared_ptr<Light>>& _lights)
	{
		lights = _lights;
		for (auto& light : lights)
			if (light->Type() == LightType::Infinite)
				infiniteLights.push_back(light);	
		lightSampler = std::make_shared<UniformLightSampler>(lights);

		sqrtSPP = static_cast<int>(sqrt(samplePixel));
		assert(sqrtSPP * sqrtSPP == samplePixel);

		cameraCenter = lookfrom;
		double h = tan(radians(fov) / 2.0);
		double viewportHeight = 2.0 * h * focusDist;
		double viewportWidth = viewportHeight * (ImageWidth * 1.0 / ImageHeight);

		vec3 w = normalize(lookfrom - lookat);
		vec3 u = normalize(cross(vup, w));
		vec3 v = normalize(cross(w, u));

		vec3 viewportU = viewportWidth * u;
		vec3 viewportV = viewportHeight * v;
		pixelDeltaU = viewportU / ImageWidth;
		pixelDeltaV = viewportV / ImageHeight;

		vec3 viewportLowerLeft = cameraCenter - focusDist * w - viewportU / 2.0 - viewportV / 2.0;
		pixel00Location = viewportLowerLeft + (pixelDeltaU + pixelDeltaV) / 2.0;

		double defocusRadius = focusDist * tan(radians(defocusAngle / 2.0));
		defocusDiskU = u * defocusRadius;
		defocusDiskV = v * defocusRadius;
	}

	bool Unoccluded(const BVHAggregate& bvh, const Interaction& p0, const Interaction& p1) const {
		ray r = SpawnRayTo(p0.pi, p0.n, p1.pi, p1.n);
		std::optional<ShapeIntersection> isect = bvh.Intersect(r, interval(0, 0.9999999));
		return !isect;
	}

	SampledSpectrum SampleLd(vec3 wo, const SurfaceInteraction& intr, const BSDF& bsdf, SampledWaveLengths& lambda, const BVHAggregate& bvh) {
		std::optional<SampledLight> sampledLight = lightSampler->Sample(randomDouble());
		if (!sampledLight) return {};
		shared_ptr<Light> light = sampledLight->light;
		vec2 u = vec2Random();
		std::optional<LightLiSample> ls = light->SampleLi(LightSampleContext(intr), u, lambda);
		if (!ls || !ls->L || ls->pdf == 0.0) return SampledSpectrum(0.0);

		vec3 wi = ls->wi;
		SampledSpectrum f = bsdf.f(wo, wi) * std::abs(dot(wi, intr.shading.n));
		if(!f || !Unoccluded(bvh, intr, ls->pLight)) return SampledSpectrum(0.0);

		double p_l = sampledLight->p * ls -> pdf;
		if (IsDeltaLight(light->Type()))
			return ls->L * f / p_l;
		double p_b = bsdf.PDF(wo, wi);
		double w_l = PowerHeuristic(1, p_l , 1, p_b);
		return ls->L * w_l * f / p_l;
	}

	SampledSpectrum Li(ray r, const int maxDepth, SampledWaveLengths& lambda, const BVHAggregate& bvh) {
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
				for (const auto &light : infiniteLights)
				{ 
					SampledSpectrum Le = light->Le(r, lambda);
					if (depth == 0 || specularBounce)
						L = L + beta * Le;
					else {
						double p_l = lightSampler->PMF(prevIntrCtx, light) * light->PDF_Li(prevIntrCtx, r.rd, true);
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
					double p_l = lightSampler->PMF(prevIntrCtx, areaLight) * areaLight->PDF_Li(prevIntrCtx, r.rd, true);
					double w_l = PowerHeuristic(1, p_b, 1, p_l);
					L = L + beta * Le * w_l;	
				}
			}

			SurfaceInteraction intr = isect->intr;
			if (depth++ == maxDepth) break;
			BSDF bsdf = intr.material->GetBSDF(intr, lambda);
			if (IsNonSpecular(bsdf.Flags())) {
				SampledSpectrum Ld = SampleLd(-r.rd, intr, bsdf, lambda, bvh);
				L = L + beta * Ld;	
			}

			vec3 wo = -r.rd;
			std::optional<BSDFSample> bs = bsdf.Sample_f(-r.rd, randomDouble(), vec2Random());
			if (!bs) break;
			beta = beta * bs->f * std::abs(dot(bs->wi, intr.shading.n)) / bs->pdf;
			p_b = bs->pdf;
			specularBounce = bs->IsSpecular();

			if (bs->IsTransmission())
				eta_scale *= Sqr(bs->eta);
			prevIntrCtx = intr;
			r = SpawnRay(intr.pi, intr.shading.n, bs->wi);

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
			r = ray(intr.p, wi);
			*/
		}
		return L;
	}

	vec3 defocusDiskSample(vec3 u, vec3 v)
	{
		vec3 p = randInDisk();
		return p.x() * u + p.y() * v;
	}
};
