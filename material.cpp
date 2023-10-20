#pragma once
#include "material.h"

bool lambertian::scatter(const ray& r, const hitRecord& rec, vec3& attenuation, ray& scattered) const 
{
	auto scatter_direction = rec.normal + randUnitVector();
	if (scatter_direction.near_zero())
		scatter_direction = rec.normal;
	scattered = ray(rec.p, scatter_direction);
	attenuation = albedo -> value(rec.u, rec.v, rec.p);
	return true;
}

double lambertian::scatterPDF(const ray& r, const hitRecord& rec, const ray& scattered) const
{
	double cosTheta = dot(rec.normal, normalize(scattered.rd));
	return cosTheta < 0.0 ? 0.0 : cosTheta / pi;
}

bool metal::scatter(const ray& r, const hitRecord& rec, vec3& attenuation, ray& scattered) const
{
	vec3 reflected = reflect(normalize(r.rd), rec.normal);
	scattered = ray(rec.p, reflected + fuzz * randUnitVector());
	attenuation = albedo;
	return (dot(scattered.rd, rec.normal) > 0);
}

bool dielectric::scatter(const ray& r, const hitRecord& rec, vec3& attenuation, ray& scattered) const
{
	attenuation = vec3(1.0, 1.0, 1.0);
	double refraction_ratio = rec.frontFace ? (1.0 / ir) : ir;

	vec3 unit_direction = normalize(r.rd);
	double cos_theta = fmin(dot(-unit_direction, rec.normal), 1.0);
	double sin_theta = sqrt(1.0 - cos_theta * cos_theta);

	bool cannot_refract = refraction_ratio * sin_theta > 1.0;
	vec3 direction;
	if (cannot_refract || reflectance(cos_theta, refraction_ratio) > randomDouble())
		direction = reflect(unit_direction, rec.normal);
	else
		direction = refract(unit_direction, rec.normal, refraction_ratio);

	scattered = ray(rec.p, direction);
	return true;
}

bool diffuseLight::scatter(const ray& r, const hitRecord& rec, vec3& attenuation, ray& scatterd) const
{
	return false;
}
