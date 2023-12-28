#pragma once
#include "math.h"
#include "spectrum.h"
#include "color.h"
#include "colorspace.h"

#include <tiny_gltf.h>
#include <memory>
#include <stb_image.h>
using std::shared_ptr;

class Image
{
public:
	float* data;
	int width, height, component;

	Image(const tinygltf::Image& image) {
		width = image.width;
		height = image.height;
		component = image.component;
		data = new float[width * height * component];
		for (int i = 0; i < width * height * component; ++i)
			data[i] = static_cast<float>(image.image.data()[i]) / 255.0f;
	}

	Image(std::string filename) {
		data = stbi_loadf(filename.c_str(), &width, &height, &component, 0);
	}

	double LookUpAlpha(vec2 uv) {
		int i = static_cast<int>(width * uv[0]);
		int j = static_cast<int>(height * uv[1]);
		i = std::max(0, std::min(i, width - 1));
		j = std::max(0, std::min(j, height - 1));
		int pixelIndex = (j * width + i) * component;
		return static_cast<double>(data[pixelIndex + 3]);
	}

	vec3 LookUp(vec2 uv, bool gammaCorrection) {
		int i = static_cast<int>(width * uv[0]);
		int j = static_cast<int>(height * uv[1]);
		i = std::max(0, std::min(i, width - 1));
		j = std::max(0, std::min(j, height - 1));
		int pixelIndex = (j * width + i) * component;
		double r = static_cast<double>(data[pixelIndex]);
		double g = static_cast<double>(data[pixelIndex + 1]);
		double b = static_cast<double>(data[pixelIndex + 2]);
		if (gammaCorrection) {
			r = pow(r, 2.2);
			g = pow(g, 2.2);
			b = pow(b, 2.2);
		}
		return vec3(r, g, b);
	}

	SampledSpectrum LookUp(vec2 uv, const SampledWaveLengths& lambda, bool gammaCorrection) {
		RGBIlluminantSpectrum spec(sRGB, RGBColor(LookUp(uv, gammaCorrection)));
		return spec.Sample(lambda);
	}
};

