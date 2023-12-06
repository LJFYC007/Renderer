#pragma once
#include <iostream>
#include <cmath>
#include <assert.h>

#include "math.h"
#include "camera.h"
#include "material.h"
#include "texture.h"
#include "colorspace.h"
#include "shape.h"
#include "model.h"
#include "lights.h"

primitiveList World;

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
    meshes.push_back(TriangleMesh(t, vertexIndices, vertices, mat));
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

    /*
    auto difflight = make_shared<diffuseLight>(vec3(4, 4, 4));
    World.add(make_shared<Sphere>(point3(0, 3, 0), 1.0, difflight));

    auto checker = make_shared<checkerBoard>(0.82, vec3(.2, .3, .1), vec3(.9, .9, .9));
    World.add(make_shared<Sphere>(point3(0, -1000, 0), 1000, make_shared<lambertian>(checker)));
    auto material1 = make_shared<dielectric>(1.5);
    World.add(make_shared<Sphere>(point3(0, 1, 0), 1.0, material1));
    auto material2 = make_shared<lambertian>(vec3(0.4, 0.2, 0.1));
    World.add(make_shared<Sphere>(point3(-4, 1, 0), 1.0, material2));
    auto material3 = make_shared<metal>(vec3(0.7, 0.6, 0.5), 0.0);
    World.add(make_shared<Sphere>(point3(4, 1, 0), 1.0, material3));

    for (double a = -11.0; a < 11.0; a += 1.0) {
        for (double b = -11.0; b < 11.0; b += 1.0) {
            double choose_mat = randomDouble();
            vec3 center = vec3(a + 0.9 * randomDouble(), 0.2, b + 0.9 * randomDouble());

            if ((center - vec3(4, 0.2, 0)).length() > 0.9) {
                if (choose_mat < 0.8)
                {
                    vec3 albedo = vec3Random() * vec3Random();
                    auto Sphere_material = make_shared<lambertian>(albedo);
                    World.add(make_shared<Sphere>(center, 0.2, Sphere_material));
                }
                else if (choose_mat < 0.95)
                {
                    vec3 albedo = vec3(randomDouble(0.5, 1));
                    double fuzz = randomDouble(0.5, 1);
                    auto Sphere_material = make_shared<metal>(albedo, fuzz);
                    World.add(make_shared<Sphere>(center, 0.2, Sphere_material));
                }
                else
                {
                    auto Sphere_material = make_shared<dielectric>(1.5);
                    World.add(make_shared<Sphere>(center, 0.2, Sphere_material));
                }
            }
        }
    }

    //Model("resources/cyborg.obj", material2);
    */

    auto red = make_shared<DiffuseMaterial>(vec3(.65, .05, .05));
    auto white = make_shared<DiffuseMaterial>(vec3(.73, .73, .73));
    auto green = make_shared<DiffuseMaterial>(vec3(.12, .45, .15));

    auto metal = make_shared<ConductorMaterial>(0.1, 0.2, 2.0, 4.0);
    auto dielectric = make_shared<DielectrivMaterial>(0.3, 0.2, 1 / 1.5);
    // auto light = make_shared<diffuseLight>(vec3(15, 15, 15));

    addBox(point3(555, 0, 0), vec3(0, 555, 0), vec3(0, 0, 555), vec3(-1, 0, 0), green);
    addBox(point3(0, 0, 0), vec3(0, 555, 0), vec3(0, 0, 555), vec3(1, 0, 0), red);
    // addBox(point3(343, 554, 332), vec3(-130, 0, 0), vec3(0, 0, -105), vec3(0, -1, 0), light);
    addBox(point3(0, 0, 0), vec3(555, 0, 0), vec3(0, 0, 555), vec3(0, 1, 0), white);
    addBox(point3(555, 555, 555), vec3(-555, 0, 0), vec3(0, 0, -555), vec3(0, -1, 0), white);
    addBox(point3(0, 0, 555), vec3(555, 0, 0), vec3(0, 555, 0), vec3(0, 0, -1), white);

    box(vec3(130, 0, 65), vec3(295, 165, 230), white, Transform::RotateY(pi / 8));
    box(vec3(265, 0, 295), vec3(430, 330, 460), metal, Transform::RotateY(-pi / 50));

    for (int i = 0; i < meshes.size(); ++i)
        for (int j = 0; j < meshes[i].nTriangles; ++j)
            World.add(make_shared<Triangle>(i, j));


    std::vector<shared_ptr<Light>> lights;
    lights.push_back(make_shared<PointLight>(Transform::Translate(vec3(278, 550, 278)), SpectrasRGB, 1.0));
    // lights.push_back(make_shared<DistantLight>(Transform::RotateX(4 * pi / 3), SpectrasRGB, 1.0));

    cam.render(bvhNode(make_shared<primitiveList>(World)), lights);
    return 0;
}