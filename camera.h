#pragma once
#include "primitive.h"
#include "sample.h"
#include "math.h"
#include "interval.h"
#include "spectrum.h"
#include "color.h"
#include "colorspace.h"

#include <Windows.h>
#include <omp.h>

static vec3 ans[3010][2210];


class camera
{
public:
	const int ImageWidth = 200;
	const int ImageHeight = 100;
	double fov = 20.0;
	vec3 lookfrom = vec3(13.0, 2.0, 3.0);
	vec3 lookat = vec3(0.0, 0.0, 0.0);
	vec3 vup = vec3(0.0, 1.0, 0.0);
	double defocusAngle = 0.6;
	double focusDist = 10.0;
	int samplePixel = 1024;
	int maxDepth = 20;
	vec3 background = vec3(0.0);

	void render(const bvhNode& World)
	{
		initialize();
#pragma omp parallel for schedule(dynamic) 
		for (int j = 0; j < ImageHeight; ++j)
		{
			std::clog << "\rScanlines remaining: " << (ImageHeight - j) << ' ' << std::flush;
			for (int i = 0; i < ImageWidth; ++i)
			{
				XYZ xyz(0.0);
				SampledSpectrum spec(0.0);
				for (int oi = 0; oi < sqrtSPP; ++oi)
					for (int oj = 0; oj < sqrtSPP; ++oj)
					{
						vec3 pixelCenter = pixel00Location + i * pixelDeltaU + j * pixelDeltaV;
						vec3 pixel = pixelCenter + (-0.5 + (1.0 / sqrtSPP) * (oi + randomDouble())) * pixelDeltaU + (-0.5 + (1.0 / sqrtSPP) * (oj + randomDouble())) * pixelDeltaV;
						vec3 ro = (defocusAngle <= 0.0) ? cameraCenter : cameraCenter + defocusDiskSample(defocusDiskU, defocusDiskV);
						vec3 rd = pixel - ro;
						SampledWaveLengths sample = (randomDouble());
						xyz = xyz + (renderRay(ray(ro, rd), maxDepth, sample, World)).ToXYZ(sample);
					}
				RGBColor rgb = sRGB.ToRGB(xyz / (double)(samplePixel));
				vec3 col = vec3(rgb.r, rgb.g, rgb.b);

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

	void initialize()
	{
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

	SampledSpectrum renderRay(ray r, const int depth, const SampledWaveLengths& sample, const bvhNode& World)
	{
		if (depth <= 0) return SampledSpectrum(0.0);
		r.rd = normalize(r.rd);

		std::optional<hitRecord> rec = World.Intersect(r, interval(0.001, infinity));
		if (!rec)
		{
			RGBAlbedoSpectrum spec(sRGB, RGBColor(background));
			return spec.Sample(sample);
		}

		ray scattered;
		SampledSpectrum attenuation(0.0);
		SampledSpectrum emissionColor = rec->mat->emitted(rec->u, rec->v, rec->p, sample);
		if (!rec->mat->scatter(r, rec.value(), sample, attenuation, scattered))
			return emissionColor;

		SampledSpectrum scatterColor = attenuation * renderRay(scattered, depth - 1, sample, World);
		return emissionColor + scatterColor;
	}

	vec3 defocusDiskSample(vec3 u, vec3 v)
	{
		vec3 p = randInDisk();
		return p.x() * u + p.y() * v;
	}
};
