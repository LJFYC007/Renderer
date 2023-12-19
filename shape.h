#pragma once
#include "aabb.h"
#include "interaction.h"
#include "light.h"
#include "material.h"

#include <optional>
#include <vector>
#include <memory>
using std::shared_ptr;
using std::make_shared;
using std::vector;

struct ShapeIntersection
{
	SurfaceInteraction intr;
	double tHit;
};

struct ShapeSample
{
	Interaction intr;
	double pdf;
};

struct ShapeSampleContext
{
	ShapeSampleContext(Vector3fi _pi, vec3 _n, vec3 _ns, double _t) : pi(_pi), n(_n), ns(_ns), t(_t) {}
	vec3 p() const { return vec3(pi); }
	Vector3fi pi;
	vec3 n, ns;
	double t;
};

class Shape
{
public:
	virtual AABB Bounds() const = 0;
	virtual std::optional<ShapeIntersection> Intersect(const ray& r, interval t) const = 0;
	// virtual bool IntersectP(const ray& r, interval t) const = 0;
	virtual double Area() const = 0;
	virtual std::optional<ShapeSample> Sample(vec2 u) const = 0;
	virtual double PDF(ShapeSampleContext sample) const {
		return 1.0 / Area();
	}
	virtual std::optional<ShapeSample> Sample(ShapeSampleContext sample, vec2 u) const {
		return {};
	}
	virtual double PDF(ShapeSampleContext sample, vec3 wi) const {
		return 0.0;
	}
};

class TriangleIntersection
{
public:
	double t, b0, b1, b2;
};

struct Vertex {
	vec3 p, n; vec2 uv;
	Vertex() {}
	Vertex(vec3 _p, vec3 _n) : p(_p), n(_n) {}
	Vertex(vec3 _p, vec3 _n, vec2 _uv) : p(_p), n(_n), uv(_uv) {}
};

class TriangleMesh
{
public:
	int nTriangles, nVertices;
	vector<int> vertexIndices;
	vector<Vertex> vertices;
	bool uvExists, normalExists;

	TriangleMesh(const Transform& ObjectToWorld, vector<int> _vertexIndices, vector<Vertex> _vertices, bool _uvExists, bool _normalExists) :
		nTriangles(static_cast<int>(_vertexIndices.size()) / 3), nVertices(static_cast<int>(_vertices.size())), vertexIndices(_vertexIndices), uvExists(_uvExists), normalExists(_normalExists) {
		vertices.resize(nVertices);
		for (int i = 0; i < nVertices; ++i)
		{
			vertices[i].p = ObjectToWorld(_vertices[i].p);
			vertices[i].n = ObjectToWorld(_vertices[i].n);
			vertices[i].uv = _vertices[i].uv;
		}
	}
};
static vector<TriangleMesh> meshes;

class Triangle : public Shape
{
public:
	Triangle(int _meshIndex, int _triIndex) : meshIndex(_meshIndex), triIndex(_triIndex) {}

	const TriangleMesh* GetMesh() const {
		return &meshes[meshIndex];
	}

	AABB Bounds() const override {
		const TriangleMesh* mesh = GetMesh();
		const int* v = &mesh->vertexIndices[3 * triIndex];
		vec3 p0 = mesh->vertices[v[0]].p, p1 = mesh->vertices[v[1]].p, p2 = mesh->vertices[v[2]].p;
		return AABB(p0, p1, p2);
	}

	std::optional<ShapeIntersection> Intersect(const ray& r, interval t) const override {
		const TriangleMesh* mesh = GetMesh();
		const int* v = &mesh->vertexIndices[3 * triIndex];
		vec3 p0 = mesh->vertices[v[0]].p, p1 = mesh->vertices[v[1]].p, p2 = mesh->vertices[v[2]].p;

		std::optional<TriangleIntersection> ints = BasicIntersect(r, t, p0, p1, p2);
		if (!ints) return {};

		vec2 uv0, uv1, uv2;
		if (!mesh->uvExists) { uv0 = vec2(0, 0); uv1 = vec2(1, 0); uv2 = vec2(1, 1); }
		else { uv0 = mesh->vertices[v[0]].uv; uv1 = mesh->vertices[v[1]].uv; uv2 = mesh->vertices[v[2]].uv; }
		vec3 p = p0 * ints->b0 + p1 * ints->b1 + p2 * ints->b2;
		vec2 uv = uv0 * ints->b0 + uv1 * ints->b1 + uv2 * ints->b2;

		vec2 duv02 = uv0 - uv2, duv12 = uv1 - uv2;
		vec3 dp02 = p0 - p2, dp12 = p1 - p2;
		double determinant = duv02[0] * duv12[1] - duv02[1] * duv12[0];
		double invdet = 1.0 / determinant;
		vec3 dpdu = (duv12[1] * dp02 - duv02[1] * dp12) * invdet;
		vec3 dpdv = (duv02[0] * dp12 - duv12[0] * dp02) * invdet;

		vec3 pAbsSum = Abs(ints->b0 * p0) + Abs(ints->b1 * p1) + Abs(ints->b2 * p2);
		vec3 pError = gamma(7) * pAbsSum;
		SurfaceInteraction intr(Vector3fi(p, pError), uv, -r.rd, dpdu, dpdv, vec3(), vec3(), ints->t, false);
		intr.n = intr.shading.n = normalize(cross(dp02, dp12));

		if (mesh->normalExists) {
			vec3 ns = ints->b0 * mesh->vertices[v[0]].n + ints->b1 * mesh->vertices[v[1]].n + ints->b2 * mesh->vertices[v[2]].n;
			ns = ns.lengthSquared() > 0 ? normalize(ns) : intr.n;
			vec3 ss = intr.dpdu;
			vec3 ts = cross(ns, ss);
			if (ts.lengthSquared() > 0) ss = cross(ts, ns);
			else CoordinateSystem(ns, ss, ts);

			vec3 dn1 = mesh->vertices[v[0]].n - mesh->vertices[v[2]].n, dn2 = mesh->vertices[v[1]].n - mesh->vertices[v[2]].n;
			vec3 dndu, dndv;
			if (determinant < 1e-9) {
				vec3 dn = cross(mesh->vertices[v[2]].n - mesh->vertices[v[0]].n, mesh->vertices[v[1]].n - mesh->vertices[v[0]].n);
				if (dn.lengthSquared() == 0) dndu = dndv = vec3(0);
				else CoordinateSystem(dn, dndu, dndv);
			}
			else {
				dndu = (duv12[1] * dn1 - duv02[1] * dn2) * invdet;
				dndv = (duv02[0] * dn2 - duv12[0] * dn1) * invdet;
			}
			intr.SetShadingGeometry(ns, ss, ts, dndu, dndv);
		}

		return ShapeIntersection{ intr, ints->t };
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
			return TriangleIntersection{ T, 1 - U - V, U, V };
		return {};
	}

	double Area() const override {
		const TriangleMesh* mesh = GetMesh();
		const int* v = &mesh->vertexIndices[3 * triIndex];
		vec3 p0 = mesh->vertices[v[0]].p, p1 = mesh->vertices[v[1]].p, p2 = mesh->vertices[v[2]].p;
		return 0.5 * cross(p1 - p0, p2 - p0).length();
	}

	std::optional<ShapeSample> Sample(vec2 uv) const override {
		const TriangleMesh* mesh = GetMesh();
		const int* v = &mesh->vertexIndices[3 * triIndex];
		vec3 p0 = mesh->vertices[v[0]].p, p1 = mesh->vertices[v[1]].p, p2 = mesh->vertices[v[2]].p;

		vec3 b = SampleUniformTriangle(uv);
		vec3 p = b.x() * p0 + b.y() * p1 + b.z() * p2;

		vec3 n = normalize(cross(p1 - p0, p2 - p0));
		vec3 ns = b.x() * mesh->vertices[v[0]].n + b.y() * mesh->vertices[v[1]].n + b.z() * mesh->vertices[v[2]].n;
		n = FaceForward(n, ns);
		vec2 uvSample = b.x() * mesh->vertices[v[0]].uv + b.y() * mesh->vertices[v[1]].uv + b.z() * mesh->vertices[v[2]].uv;
		vec3 pAbsSum = Abs(b.x() * p0) + Abs(b.y() * p1) + Abs(b.z() * p2);
		vec3 pError = gamma(6) * pAbsSum;
		return ShapeSample{ Interaction(Vector3fi(p, pError), n, uvSample), 1 / Area() };
	}

	std::optional<ShapeSample> Sample(ShapeSampleContext sample, vec2 u) const override {
		std::optional<ShapeSample> ss = Sample(u);
		if (!ss) return {};
		vec3 wi = ss->intr.p() - sample.p();
		if (wi.lengthSquared() == 0.0) return ShapeSample{ ss->intr, 0.0 };
		wi = normalize(wi);
		ss->pdf /= std::abs(dot(ss->intr.n, -wi)) / (sample.p() - ss->intr.p()).lengthSquared();
		if (std::isinf(ss->pdf)) ss->pdf = 0.0;
		return ss;
	}

	double PDF(ShapeSampleContext sample, vec3 wi) const override {
		ray r = SpawnRay(sample.pi, sample.n, wi);
		std::optional<ShapeIntersection> isect = Intersect(r, interval(0, infinity));
		if (!isect) return 0.0;
		SurfaceInteraction intr = isect->intr;
		double pdf = (1.0 / Area()) / (std::abs(dot(intr.n, -wi)) / (sample.p() - intr.p()).lengthSquared());
		if (std::isinf(pdf)) pdf = 0.0;
		return pdf;
	}

private:
	int meshIndex = -1, triIndex = -1;
};
