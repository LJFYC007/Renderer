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
	unsigned char* data;
	int width, height, component;

	Image(const tinygltf::Image& image) {
		width = image.width;
		height = image.height;
		component = image.component;
		data = new unsigned char[width * height * component];
		for (int i = 0; i < width * height * component; ++i)
			data[i] = image.image.data()[i];
	}

	double LookUpAlpha(vec2 uv) {
		int i = static_cast<int>(width * uv[0]);
		int j = static_cast<int>(height * uv[1]);
		i = std::max(0, std::min(i, width - 1));
		j = std::max(0, std::min(j, height - 1));
		int pixelIndex = (j * width + i) * component;
		return data[pixelIndex + 3] / 255.0;
	}

	vec3 LookUp(vec2 uv, bool gammaCorrection) {
		int i = static_cast<int>(width * uv[0]);
		int j = static_cast<int>(height * uv[1]);
		i = std::max(0, std::min(i, width - 1));
		j = std::max(0, std::min(j, height - 1));
		int pixelIndex = (j * width + i) * component;
		double r = data[pixelIndex] / 255.0;
		double g = data[pixelIndex + 1] / 255.0;
		double b = data[pixelIndex + 2] / 255.0;
		if (gammaCorrection) {
			r = pow(r, 2.2);
			g = pow(g, 2.2);
			b = pow(b, 2.2);
		}
		return vec3(r, g, b);
	}

	SampledSpectrum LookUp(vec2 uv, const SampledWaveLengths& lambda, bool gammaCorrection) {
		RGBAlbedoSpectrum spec(sRGB, RGBColor(LookUp(uv, gammaCorrection)));
		return spec.Sample(lambda);
	}
};

class HDRImage
{
public:
	float* data;
	int width, height, component;

	HDRImage(std::string filename) {
		data = stbi_loadf(filename.c_str(), &width, &height, &component, 0);
	}
	
	vec3 LookUp(vec2 uv, bool gammaCorrection) {
		uv[0] = 1 - uv[0];
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

