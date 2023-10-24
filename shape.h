#pragma once
#include "aabb.h"
#include "interaction.h"
#include "material.h"

#include <optional>
#include <vector>

using std::vector;

class QuadricIntersection
{
public:
	double tHit;
	vec3 p;
};

class Shape
{
	virtual AABB Bounds() const = 0;
	virtual std::optional<hitRecord> Intersect(const ray& r, interval t) const = 0;
	//virtual bool IntersectP(const ray& r, interval t) const = 0;
	// virtual double Area() const = 0;
};

class Sphere : public Shape
{
public:
	Sphere(point3 _center, double _radius, shared_ptr<material> _mat) : center(_center), radius(_radius), mat(_mat) {}

	AABB Bounds() const override {
		return AABB(center - vec3(radius), center + vec3(radius));
	}

	std::optional<hitRecord> Intersect(const ray& r, interval t) const override {
		std::optional<QuadricIntersection> ints = BasicIntersect(r, t);
		if (!ints) return {};

		hitRecord rec;
		rec.t = ints->tHit; rec.p = ints->p;
		rec.mat = mat;
		rec.normal = (rec.p - center) / radius;

		rec.frontFace = dot(rec.normal, r.rd) < 0.0;
		rec.normal = rec.frontFace ? rec.normal : -rec.normal;
		return rec;
	}

	std::optional<QuadricIntersection> BasicIntersect(const ray& r, interval t) const {
		vec3 oc = r.ro - center;
		double a = dot(r.rd, r.rd);
		double halfb = dot(oc, r.rd);
		double c = dot(oc, oc) - radius * radius;

		double discriminant = halfb * halfb - a * c;
		if (discriminant < 0.0) return {};
		discriminant = sqrt(discriminant);
		double root = (-halfb - discriminant) / a;
		if (!t.surrounds(root))
		{
			root = (-halfb + discriminant) / a;
			if (!t.surrounds(root)) return {};
		}
		return QuadricIntersection{ root, r.at(root) };
	}

private:
	vec3 center;
	double radius;
	shared_ptr<material> mat;
};

class TriangleMesh
{
public:
	int nTriangles, nVertices;
	vector<int> vertexIndices;
	vector<vec3> p, s, n;
	vector<vec2> uv;
	vector<int> faceIndices;

	TriangleMesh(vector<int> _vertexIndices, vector<vec3> _p, vector<vec3> _s, vector<vec3> _n, vector<vec2> _uv, vector<int> _faceIndices) :
		nTriangles(_vertexIndices.size() / 3), nVertices(_p.size()), vertexIndices(_vertexIndices), p(_p), s(_s), n(_n), uv(_uv), faceIndices(_faceIndices) {}
};

class TriangleIntersection 
{
public:
	double t, u, v;
};

class Triangle : public Shape
{
public:
	Triangle(vec3 _p0, vec3 _p1, vec3 _p2, shared_ptr<material> _mat) : p0(_p0), p1(_p1), p2(_p2), mat(_mat) {}

	AABB Bounds() const override {
		/*
		const TriangleMesh* mesh = GetMesh();
		const int* v = &mesh->vertexIndices[3 * triIndex];
		vec3 p0 = mesh->p[v[0]], p1 = mesh->p[v[1]], p2 = mesh->p[v[2]];
		*/
		return AABB(p0, p1).Union(p2);
	}

	std::optional<hitRecord> Intersect(const ray& r, interval t) const override {
		std::optional<TriangleIntersection> ints = BasicIntersect(r, t, p0, p1, p2);
		if (!ints) return {};

		hitRecord rec;
		rec.t = ints->t; rec.p = p0 + (p1 - p0) * ints->u + (p2 - p0) * ints->v;
		rec.mat = mat;
		return rec;
	}

	std::optional<TriangleIntersection> BasicIntersect(const ray& r, interval t, vec3 p0, vec3 p1, vec3 p2) const {
		vec3 e1 = p1 - p0, e2 = p2 - p0, s = r.ro - p0;
		vec3 s1 = cross(r.rd, e2), s2 = cross(s, e1);

		double coeff = 1.0 / dot(s1, e1);
		double T = coeff * dot(s2, e2);
		double U = coeff * dot(s1, s);
		double V = coeff * dot(s2, r.rd);

		if (!t.surrounds(T)) return {};

		if (T >= 0 && U >= 0 && V >= 0 && (1 - U - V) >= 0)
			return TriangleIntersection{ T, U, V };
		return {};
	}

private:
	vec3 p0, p1, p2;
	shared_ptr<material> mat;

	/*
	int meshIndex = -1, triIndex = -1;
	static vector<const TriangleMesh*>* allMeshes;

	const TriangleMesh* GetMesh() const {
		return (*allMeshes)[meshIndex];
	}
	*/
};
