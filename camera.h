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
#include "integrator.h"

#include <Windows.h>
#include <omp.h>
#include <optional>

static vec3 ans[3010][2210];

class Camera
{
public:
	int ImageWidth = 500;
	int ImageHeight = 500;
	double fov = 20.0;
	vec3 lookfrom = vec3(-5.5, 1.2, 9.0);
	vec3 lookat = vec3(0.5, 1.2, 0.0);
	vec3 vup = vec3(0.0, 1.0, 0.0);
	double defocusAngle = 0.0;
	double focusDist = 10.0;
	int samplePixel = 256;
	int maxDepth = 10;
	std::vector<shared_ptr<Light>> lights;
	std::vector<shared_ptr<Light>> infiniteLights;
	shared_ptr<LightSampler> lightSampler;

	void Render(const BVHAggregate& bvh, const std::vector<shared_ptr<Light>>& wholeLights)
	{
		lights = wholeLights;
		Initialize();
#pragma omp parallel for schedule(dynamic) 
		for (int j = 0; j < ImageHeight; ++j)
		{
			std::clog << "\rScanlines remaining: " << (ImageHeight - j) << ' ' << std::flush;
			for (int i = 0; i < ImageWidth; ++i)
			{
				XYZ xyz(0.0);
				SampledSpectrum spec(0.0);
				vec3 col(0.0);
				for (int oi = 0; oi < sqrtSpp; ++oi)
					for (int oj = 0; oj < sqrtSpp; ++oj)
					{
						vec3 pixelCenter = pixel00Location + i * pixelDeltaU + j * pixelDeltaV;
						vec3 pixel = pixelCenter + (-0.5 + (1.0 / sqrtSpp) * (oi + randomDouble())) * pixelDeltaU + (-0.5 + (1.0 / sqrtSpp) * (oj + randomDouble())) * pixelDeltaV;
						vec3 ro = (defocusAngle <= 0.0) ? cameraCenter : cameraCenter + DefocusDiskSample(defocusDiskU, defocusDiskV);
						vec3 rd = pixel - ro;
						SampledWaveLengths sample(randomDouble());

						RayDifferential ray(ro, rd, 0.0);
						ray.hasDifferentials = true;
						double scale = focusDist / rd.z(); 
						vec3 dPdx = normalize((pixelCenter + pixelDeltaU) - ro) * scale;
						vec3 dPdy = normalize((pixelCenter + pixelDeltaV) - ro) * scale; 
						ray.rxOrigin = ray.ryOrigin = ro;
						ray.rxDirection = rd + (dPdx - rd);
						ray.ryDirection = rd + (dPdy - rd);

						RGBColor rgb = sRGB.ToRGB((::Li(ray, this, sample, bvh)).ToXYZ(sample));
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

	void Approximate_dp_dxy(const vec3& p, const vec3& n, double time, vec3* dpdx, vec3* dpdy) {
		vec3 pVector = normalize(p - cameraCenter);
		vec3 p_dx = cameraCenter + (pVector + pixelDeltaU * (1.0 / sqrtSpp)) * focusDist;
		vec3 p_dy = cameraCenter + (pVector + pixelDeltaV * (1.0 / sqrtSpp)) * focusDist;
		*dpdx = p_dx - p;
		*dpdy = p_dy - p;
		vec3 normal = normalize(n);
		*dpdx = *dpdx - normal * dot(*dpdx, normal);
		*dpdy = *dpdy - normal * dot(*dpdy, normal);
	}

private:
	int sqrtSpp = 0;
	vec3 defocusDiskU;
	vec3 defocusDiskV;
	vec3 pixel00Location;
	vec3 cameraCenter;
	vec3 pixelDeltaU;
	vec3 pixelDeltaV;

	void Initialize() {
		for (auto& light : lights)
			if (light->Type() == LightType::Infinite)
				infiniteLights.emplace_back(light);	
		lightSampler = std::make_shared<UniformLightSampler>(lights);
		
		sqrtSpp = static_cast<int>(sqrt(samplePixel));
		assert(sqrtSpp * sqrtSpp == samplePixel);

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

	static vec3 DefocusDiskSample(const vec3& u, const vec3& v)
	{
		vec3 p = randInDisk();
		return p.x() * u + p.y() * v;
	}
};

