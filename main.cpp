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

/*
void addBox(vec3 a, vec3 b, vec3 c, vec3 n, const shared_ptr<Material>& mat, Transform t = Transform())
{
    vector<Vertex> vertices;
    vector<int> vertexIndices;
    vertices.emplace_back(a, n);
    vertices.emplace_back(a + b, n);
    vertices.emplace_back(a + c, n);
    vertices.emplace_back(a + b + c, n);
    vertexIndices.emplace_back(0); vertexIndices.emplace_back(1); vertexIndices.emplace_back(2);
    vertexIndices.emplace_back(1); vertexIndices.emplace_back(2); vertexIndices.emplace_back(3);
    meshes.emplace_back(t, vertexIndices, vertices, false, false, true);
    World.emplace_back(make_shared<SimplePrimitive>(make_shared<Triangle>(static_cast<int>(meshes.size()) - 1, 0), mat));
    World.emplace_back(make_shared<SimplePrimitive>(make_shared<Triangle>(static_cast<int>(meshes.size()) - 1, 1), mat));
}

void box(vec3 a, vec3 b, const shared_ptr<Material>& mat, const Transform& t = Transform())
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
*/

int main()
{
    Camera cam;
	Model model(World, "resources/glb.glb");

    /*
    auto red = make_shared<DiffuseMaterial>(make_shared<SpectrumConstantTexture>(vec3(.65, .05, .05)));
    auto white = make_shared<DiffuseMaterial>(make_shared<SpectrumConstantTexture>(vec3(.73, .73, .73)));
    auto green = make_shared<DiffuseMaterial>(make_shared<SpectrumConstantTexture>(vec3(.12, .45, .15)));

    auto metal = make_shared<ConductorMaterial>(0.1, 0.1, 2.0, 4.0);
    auto dielectric = make_shared<DielectricMaterial>(0.3, 0.2, 1 / 1.5);
    auto mirror = make_shared<ConductorMaterial>(0.0, 0.0, 2.0, 4.0);

    Model model(World, "resources/CaveCrystal01.obj");

    addBox(vec3(277, -278, 0), vec3(0, 555, 0), vec3(0, 0, 555), vec3(-1, 0, 0), green);
    addBox(vec3(-278, -278, 0), vec3(0, 555, 0), vec3(0, 0, 555), vec3(1, 0, 0), red);
    addBox(vec3(-278, -278, 0), vec3(555, 0, 0), vec3(0, 0, 555), vec3(0, 1, 0), white);
    addBox(vec3(277, 277, 555), vec3(-555, 0, 0), vec3(0, 0, -555), vec3(0, -1, 0), white);
    addBox(vec3(-278, -278, 555), vec3(555, 0, 0), vec3(0, 555, 0), vec3(0, 0, -1), white);

    // box(vec3(-148, -278, 65), vec3(17, -113, 230), mirror, Transform::RotateY(pi / 8));
    // box(vec3(-13, -278, 295), vec3(152, 52, 460), metal, Transform::RotateY(-pi / 50));

    vec3 a(65, 276, 332), b(-130, 0, 0), c(0, 0, -105);
    vec3 n(0, -1, 0);
    std::vector<shared_ptr<Light>> lights;

    {
        vector<Vertex> vertices;
        vector<int> vertexIndices;
        vertices.emplace_back(a, n);
        vertices.emplace_back(a + b, n);
        vertices.emplace_back(a + c, n);
        vertexIndices.emplace_back(0); vertexIndices.emplace_back(1); vertexIndices.emplace_back(2);
        meshes.emplace_back(Transform(), vertexIndices, vertices, false, true);

        shared_ptr<Triangle> triangle = make_shared<Triangle>(static_cast<int>(meshes.size()) - 1, 0);
        DiffuseAreaLight areaLight(Transform(), SpectrasRGB, 10.0, triangle);
        lights.emplace_back(make_shared<DiffuseAreaLight>(areaLight));
        World.emplace_back(make_shared<GeometricPrimitive>(triangle, white, make_shared<DiffuseAreaLight>(areaLight)));
    }

    {
        vector<Vertex> vertices;
        vector<int> vertexIndices;
        vertices.emplace_back(a + b, n);
        vertices.emplace_back(a + c, n);
        vertices.emplace_back(a + b + c, n);
        vertexIndices.emplace_back(0); vertexIndices.emplace_back(1); vertexIndices.emplace_back(2);
        meshes.emplace_back(Transform(), vertexIndices, vertices, false, true);

        shared_ptr<Triangle> triangle = make_shared<Triangle>(static_cast<int>(meshes.size()) - 1, 0);
        DiffuseAreaLight areaLight(Transform(), SpectrasRGB, 10.0, triangle);
        lights.emplace_back(make_shared<DiffuseAreaLight>(areaLight));
        World.emplace_back(make_shared<GeometricPrimitive>(triangle, white, make_shared<DiffuseAreaLight>(areaLight)));
    }
    */

    std::vector<shared_ptr<Light>> lights;
    // lights.emplace_back(make_shared<PointLight>(Transform::Translate(vec3(278, 550, 278)), SpectrasRGB, 1.0));
    lights.emplace_back(make_shared<DistantLight>(Transform::Rotate(-0.20693400502204895,
        0.3364157974720001,
        0.08685380965471268,
        0.9145814180374146), SpectrasRGB, 3.3));
    lights.emplace_back(make_shared<DistantLight>(Transform::Rotate(-0.2957703173160553,
        0.6428548693656921,
        0.29644736647605896,
        0.6413864493370056), SpectrasRGB, 2.0));
    lights.emplace_back(make_shared<ImageInfiniteLight>(Transform::RotateX(pi / 2), 0.1, "1.hdr"));

    BVHAggregate bvh(World, 3);
    cam.Render(bvh, lights);
    return 0;
}