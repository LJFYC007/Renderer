#pragma once
#include <iostream>
#include <cmath>
#include <assert.h>

#include "math.h"
#include "camera.h"
#include "material.h"
#include "textures.h"
#include "colorspace.h"
#include "shape.h"
#include "model.h"
#include "lights.h"

std::vector<shared_ptr<Primitive>> World;

void addBox(vec3 a, vec3 b, vec3 c, vec3 n, shared_ptr<Material> mat, Transform t = Transform())
{
    vector<Vertex> vertices;
    vector<int> vertexIndices;
    vertices.push_back(Vertex(a, n));
    vertices.push_back(Vertex(a + b, n));
    vertices.push_back(Vertex(a + c, n));
    vertices.push_back(Vertex(a + b + c, n));
    vertexIndices.push_back(0); vertexIndices.push_back(1); vertexIndices.push_back(2);
    vertexIndices.push_back(1); vertexIndices.push_back(2); vertexIndices.push_back(3);
    meshes.push_back(TriangleMesh(t, vertexIndices, vertices, false, true));
    World.push_back(make_shared<SimplePrimitive>(make_shared<Triangle>(static_cast<int>(meshes.size()) - 1, 0), mat));
    World.push_back(make_shared<SimplePrimitive>(make_shared<Triangle>(static_cast<int>(meshes.size()) - 1, 1), mat));
}

void box(vec3 a, vec3 b, shared_ptr<Material> mat, Transform t = Transform())
{
    auto dx = vec3(b[0] - a[0], 0, 0);
    auto dy = vec3(0, b[1] - a[1], 0);
    auto dz = vec3(0, 0, b[2] - a[2]);
    addBox(vec3(a[0], a[1], b[2]), dx, dy, vec3(0, 0, -1), mat, t);
    addBox(vec3(b[0], a[1], b[2]), -dz, dy, vec3(1, 0, 0), mat, t);
    addBox(vec3(b[0], a[1], a[2]), -dx, dy, vec3(0, 0, 1), mat, t);
    addBox(vec3(a[0], a[1], a[2]), dz, dy, vec3(-1, 0, 0), mat, t);
    addBox(vec3(a[0], b[1], b[2]), dx, -dz, vec3(0, 1, 0), mat, t);
    addBox(vec3(a[0], a[1], a[2]), dx, dz, vec3(0, -1, 0), mat, t);
}

int main()
{
    camera cam;

    auto red = make_shared<DiffuseMaterial>(make_shared<SpectrumConstantTexture>(vec3(.65, .05, .05)));
    auto white = make_shared<DiffuseMaterial>(make_shared<SpectrumConstantTexture>(vec3(.73, .73, .73)));
    auto green = make_shared<DiffuseMaterial>(make_shared<SpectrumConstantTexture>(vec3(.12, .45, .15)));

    auto metal = make_shared<ConductorMaterial>(0.1, 0.1, 2.0, 4.0);
    auto dielectric = make_shared<DielectricMaterial>(0.3, 0.2, 1 / 1.5);
    auto mirror = make_shared<ConductorMaterial>(0.0, 0.0, 2.0, 4.0);

    Model(World, "resources/CaveCrystal01.obj");

    addBox(vec3(277, -278, 0), vec3(0, 555, 0), vec3(0, 0, 555), vec3(-1, 0, 0), green);
    addBox(vec3(-278, -278, 0), vec3(0, 555, 0), vec3(0, 0, 555), vec3(1, 0, 0), red);
    addBox(vec3(-278, -278, 0), vec3(555, 0, 0), vec3(0, 0, 555), vec3(0, 1, 0), white);
    addBox(vec3(277, 277, 555), vec3(-555, 0, 0), vec3(0, 0, -555), vec3(0, -1, 0), white);
    addBox(vec3(-278, -278, 555), vec3(555, 0, 0), vec3(0, 555, 0), vec3(0, 0, -1), white);

    // box(vec3(-148, -278, 65), vec3(17, -113, 230), mirror, Transform::RotateY(pi / 8));
    // box(vec3(-13, -278, 295), vec3(152, 52, 460), metal, Transform::RotateY(-pi / 50));

    vec3 a(65, 276, 332), b(-130, 0, 0), c(0, 0, -105);
    vec3 n = (0, -1, 0);
    std::vector<shared_ptr<Light>> lights;

    {
        vector<Vertex> vertices;
        vector<int> vertexIndices;
        vertices.push_back(Vertex(a, n));
        vertices.push_back(Vertex(a + b, n));
        vertices.push_back(Vertex(a + c, n));
        vertexIndices.push_back(0); vertexIndices.push_back(1); vertexIndices.push_back(2);
        meshes.push_back(TriangleMesh(Transform(), vertexIndices, vertices, false, true));

        shared_ptr<Triangle> triangle = make_shared<Triangle>(static_cast<int>(meshes.size()) - 1, 0);
        DiffuseAreaLight areaLight(Transform(), SpectrasRGB, 10.0, triangle);
        lights.push_back(make_shared<DiffuseAreaLight>(areaLight));
        World.push_back(make_shared<GeometricPrimitive>(triangle, white, make_shared<DiffuseAreaLight>(areaLight)));
    }

    {
        vector<Vertex> vertices;
        vector<int> vertexIndices;
        vertices.push_back(Vertex(a + b, n));
        vertices.push_back(Vertex(a + c, n));
        vertices.push_back(Vertex(a + b + c, n));
        vertexIndices.push_back(0); vertexIndices.push_back(1); vertexIndices.push_back(2);
        meshes.push_back(TriangleMesh(Transform(), vertexIndices, vertices, false, true));

        shared_ptr<Triangle> triangle = make_shared<Triangle>(static_cast<int>(meshes.size()) - 1, 0);
        DiffuseAreaLight areaLight(Transform(), SpectrasRGB, 10.0, triangle);
        lights.push_back(make_shared<DiffuseAreaLight>(areaLight));
        World.push_back(make_shared<GeometricPrimitive>(triangle, white, make_shared<DiffuseAreaLight>(areaLight)));
    }

    //lights.push_back(make_shared<PointLight>(Transform::Translate(vec3(278, 550, 278)), SpectrasRGB, 1.0));
    lights.push_back(make_shared<DistantLight>(Transform::RotateX(pi), SpectrasRGB, 0.05));

    BVHAggregate bvh(World, 3);
    cam.render(bvh, lights);
    return 0;
}