#include "interaction.h"
#include "lights.h"
#include "camera.h"

Ray SurfaceInteraction::SpawnRay(const Ray& rayi, const BSDF& bsdf, vec3 wi, int flags, double eta) const
{
	return Interaction::SpawnRay(wi);
}

SampledSpectrum SurfaceInteraction::Le(vec3 w, const SampledWaveLengths& lambda) const {
	return areaLight ? areaLight->L(p(), n, uv, w, lambda) : SampledSpectrum(0.0);
}

BSDF SurfaceInteraction::GetBSDF(const Ray &ray, SampledWaveLengths& lambda, Camera* camera) {
	shared_ptr<SpectrumTexture> normalMap = material->GetNormalMap();
	if (normalMap) {
		vec3 dpdu, dpdv, dpdU, dpdV;
		NormalMap(normalMap, NormalBumpEvalContext(*this), &dpdu, &dpdv, &dpdU, &dpdV);
		vec3 ns = normalize(cross(dpdu, dpdv));
		SetShadingGeometry(ns, dpdu, dpdv, dpdU, dpdV, shading.dndu, shading.dndv);
	}
	while ( material->IsMixMaterial() == true )
		material = material->GetMaterial(*this, randomDouble());
	BSDF bsdf = material->GetBSDF(*this, lambda);
	return bsdf;
}
