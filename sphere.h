#pragma once
#include "primitive.h"
#include "math.h"
#include "ray.h"

class sphere : public primitive
{
private : 
    vec3 center;
    double radius;
	shared_ptr<material> mat;
	bvh box;

public : 
	sphere(point3 _center, double _radius, shared_ptr<material> _mat) : center(_center), radius(_radius), mat(_mat) {
		vec3 v(radius);
		box = bvh(center - v, center + v);
	}

	bool hit(const ray& r, interval t, hitRecord& rec) const override
	{
		vec3 oc = r.ro - center;
		double a = dot(r.rd, r.rd);
		double halfb = dot(oc, r.rd);
		double c = dot(oc, oc) - radius * radius;

		double discriminant = halfb * halfb - a * c;
		if (discriminant < 0.0) return false;
		discriminant = sqrt(discriminant);
		double root = (-halfb - discriminant) / a;
		if (!t.surrounds(root))
		{
			root = (-halfb + discriminant) / a;
			if (!t.surrounds(root)) return false;
		}

		rec.t = root; rec.p = r.at(root);
		rec.normal = (rec.p - center) / radius;
		rec.mat = mat;

		rec.frontFace = dot(rec.normal, r.rd) < 0.0;
		rec.normal = rec.frontFace ? rec.normal : -rec.normal;
		return true;
	}

	bvh boundingBox() const override { return box; }
};
