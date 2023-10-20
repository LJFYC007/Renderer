#pragma once
#include "math.h"
#include "sample.h"
#include "ray.h"
#include "primitive.h"
#include "texture.h"
#include "spectrum.h"

#include <memory>
using std::shared_ptr;
using std::make_shared;

class hitRecord;

class material
{
public : 
	virtual ~material() = default;
	virtual bool scatter(const ray& r, const hitRecord& rec, const SampledWaveLengths& sample, SampledSpectrum& attenuation, ray& scattered) const = 0;
    virtual vec3 emitted(double u, double v, const vec3& p) const { return vec3(0.0); }
    virtual double scatterPDF(const ray& r, const hitRecord& rec, const ray& scattered) const { return 0; }
};

class lambertian : public material {
public:
    lambertian(const vec3& a) : albedo(make_shared<solidColor>(a)) {}
    lambertian(shared_ptr<texture> _albedo) : albedo(_albedo) {}
    bool scatter(const ray& r, const hitRecord& rec, const SampledWaveLengths& sample, SampledSpectrum& attenuation, ray& scattered) const override;
    double scatterPDF(const ray& r, const hitRecord& rec, const ray& scattered) const override;

private:
    shared_ptr<texture> albedo;
};

class dielectric : public material {
public:
    dielectric(double index_of_refraction) : ir(index_of_refraction) {}
    bool scatter(const ray& r, const hitRecord& rec, const SampledWaveLengths& sample, SampledSpectrum& attenuation, ray& scattered) const override;

private:
    double ir;

    static double reflectance(double cosine, double ref_idx) {
        auto r0 = (1 - ref_idx) / (1 + ref_idx);
        r0 = r0 * r0;
        return r0 + (1 - r0) * pow((1 - cosine), 5);
    }
};

class metal : public material {
public:
    metal(const vec3& a, double f) : albedo(RGBAlbedoSpectrum(sRGB, RGBColor(a))), fuzz(f < 1 ? f : 1) {}
    bool scatter(const ray& r, const hitRecord& rec, const SampledWaveLengths& sample, SampledSpectrum& attenuation, ray& scattered) const override;

private:
    RGBAlbedoSpectrum albedo;
    double fuzz;
};
