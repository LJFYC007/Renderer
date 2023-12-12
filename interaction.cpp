#include "interaction.h"
#include "lights.h"

SampledSpectrum SurfaceInteraction::Le(vec3 w, const SampledWaveLengths& lambda) const {
	return areaLight ? areaLight->L(p(), n, uv, w, lambda) : SampledSpectrum(0.0);
}
