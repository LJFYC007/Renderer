#include "interaction.h"
#include "lights.h"
#include "camera.h"

RayDifferential SurfaceInteraction::SpawnRay(const RayDifferential& rayi, const BSDF& bsdf, vec3 wi, int flags, double eta) const
{
	RayDifferential rd(Interaction::SpawnRay(wi));
	if (rayi.hasDifferentials) {
		vec3 n = shading.n;
		vec3 dndx = shading.dndu * dudx + shading.dndv * dvdx;
		vec3 dndy = shading.dndu * dudy + shading.dndv * dvdy;
		vec3 dwodx = -rayi.rxDirection - wo, dwody = -rayi.ryDirection - wo;
		if (flags == BxDFFlags::SpecularReflection) {
			rd.hasDifferentials = true;
			rd.rxOrigin = p() + dpdx;
			rd.ryOrigin = p() + dpdy;
			double dwoDotn_dx = dot(dwodx, n) + dot(wo, dndx);
			double dwoDotn_dy = dot(dwody, n) + dot(wo, dndy);
			rd.rxDirection = wi - dwodx + 2 * vec3(dot(wo, n) * dndx + dwoDotn_dx * n);
			rd.ryDirection = wi - dwody + 2 * vec3(dot(wo, n) * dndy + dwoDotn_dy * n);
		}
		else if (flags == BxDFFlags::SpecularTransmission) {
			rd.hasDifferentials = true;
			rd.rxOrigin = p() + dpdx;
			rd.ryOrigin = p() + dpdy;
			if (dot(wo, n) < 0) { n = -n; dndx = -dndx; dndy = -dndy; }
			double dwoDotn_dx = dot(dwodx, n) + dot(wo, dndx);
			double dwoDotn_dy = dot(dwody, n) + dot(wo, dndy);
			double mu = dot(wo, n) / eta - std::abs(dot(wi, n));
			double dmudx = dwoDotn_dx * (1 / eta + 1 / Sqr(eta) * dot(wo, n) / dot(wi, n));
			double dmudy = dwoDotn_dy * (1 / eta + 1 / Sqr(eta) * dot(wo, n) / dot(wi, n));
			rd.rxDirection = wi - eta * dwodx + vec3(mu * dndx + dmudx * n);
			rd.ryDirection = wi - eta * dwody + vec3(mu * dndy + dmudy * n);
		}
	}
	if (rd.rxDirection.lengthSquared() > 1e16f || rd.ryDirection.lengthSquared() > 1e16f || rd.rxOrigin.lengthSquared() > 1e16f || rd.ryOrigin.lengthSquared() > 1e16f)
		rd.hasDifferentials = false;
	return rd;
}

void SurfaceInteraction::ComputeDifferentials(const RayDifferential& ray, Camera* camera) {
	if (ray.hasDifferentials && dot(n, ray.rxDirection) != 0 && dot(n, ray.ryDirection) != 0) {
		double d = -dot(n, p());
		double tx = (-dot(n, ray.rxOrigin) - d) / dot(n, ray.rxDirection);
		vec3 px = ray.rxOrigin + tx * ray.rxDirection;
		double ty = (-dot(n, ray.ryOrigin) - d) / dot(n, ray.ryDirection);
		vec3 py = ray.ryOrigin + ty * ray.ryDirection;
		dpdx = px - p();
		dpdy = py - p();
	}
	else {
		camera->Approximate_dp_dxy(p(), n, time, &dpdx, &dpdy);
	}

	double ata00 = dot(dpdu, dpdu), ata01 = dot(dpdu, dpdv);
	double ata11 = dot(dpdv, dpdv);
	double invDet = 1.0 / (ata00 * ata11 - ata01 * ata01);

	double atb0x = dot(dpdu, dpdx), atb1x = dot(dpdv, dpdx);
	double atb0y = dot(dpdu, dpdy), atb1y = dot(dpdv, dpdy);

	dudx = (ata11 * atb0x - ata01 * atb1x) * invDet;
	dvdx = (ata00 * atb1x - ata01 * atb0x) * invDet;
	dudy = (ata11 * atb0y - ata01 * atb1y) * invDet;
	dvdy = (ata00 * atb1y - ata01 * atb0y) * invDet;
}

SampledSpectrum SurfaceInteraction::Le(vec3 w, const SampledWaveLengths& lambda) const {
	return areaLight ? areaLight->L(p(), n, uv, w, lambda) : SampledSpectrum(0.0);
}

BSDF SurfaceInteraction::GetBSDF(const RayDifferential &ray, const SampledWaveLengths& lambda, Camera* camera) {
	ComputeDifferentials(ray, camera);
	shared_ptr<SpectrumTexture> normalMap = material->GetNormalMap();
	if (normalMap) {
		vec3 dpdu, dpdv;
		NormalMap(normalMap, NormalBumpEvalContext(*this), &dpdu, &dpdv);
		vec3 ns = normalize(cross(dpdu, dpdv));
		SetShadingGeometry(ns, dpdu, dpdv, shading.dndu, shading.dndv);
	}
	shared_ptr<SpectrumTexture> bumpMap = material->GetBumpMap();
	if (bumpMap) {
		vec3 dpdu, dpdv;
		BumpMap(bumpMap, NormalBumpEvalContext(*this), &dpdu, &dpdv);
		vec3 ns = normalize(cross(dpdu, dpdv));
		SetShadingGeometry(ns, dpdu, dpdv, shading.dndu, shading.dndv);
	}
	BSDF bsdf = material->GetBSDF(*this, lambda);
	return bsdf;
}
