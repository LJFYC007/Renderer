#pragma once
#include <iostream>
#include <cmath>
#include <assert.h>

#include "math.h"
#include "camera.h"
#include "sphere.h"
#include "bvh.h"
#include "material.h"
#include "texture.h"
#include "colorspace.h"

primitiveList World;

int main()
{
    camera cam;

    auto difflight = make_shared<diffuseLight>(vec3(4, 4, 4));
    World.add(make_shared<sphere>(point3(0, 3, 0), 1.0, difflight));

    auto checker = make_shared<checkerBoard>(0.82, vec3(.2, .3, .1), vec3(.9, .9, .9));
	World.add(make_shared<sphere>(point3(0, -1000, 0), 1000, make_shared<lambertian>(checker)));
    auto material1 = make_shared<dielectric>(1.5);
    World.add(make_shared<sphere>(point3(0, 1, 0), 1.0, material1));
    auto material2 = make_shared<lambertian>(vec3(0.4, 0.2, 0.1));
    World.add(make_shared<sphere>(point3(-4, 1, 0), 1.0, material2));
    auto material3 = make_shared<metal>(vec3(0.7, 0.6, 0.5), 0.0);
    World.add(make_shared<sphere>(point3(4, 1, 0), 1.0, material3));

    for (double a = -11.0; a < 11.0; a += 1.0) {
        for (double b = -11.0; b < 11.0; b += 1.0) {
            double choose_mat = randomDouble();
            vec3 center = vec3(a + 0.9 * randomDouble(), 0.2, b + 0.9 * randomDouble());

            if ((center - vec3(4, 0.2, 0)).length() > 0.9) {
                if (choose_mat < 0.8)
                {
                    vec3 albedo = vec3Random() * vec3Random();
                    auto sphere_material = make_shared<lambertian>(albedo);
                    World.add(make_shared<sphere>(center, 0.2, sphere_material));
                }
                else if (choose_mat < 0.95)
                {
                    vec3 albedo = vec3(randomDouble(0.5, 1));
                    double fuzz = randomDouble(0.5, 1);
                    auto sphere_material = make_shared<metal>(albedo, fuzz);
                    World.add(make_shared<sphere>(center, 0.2, sphere_material));
                }
                else
                {
                    auto sphere_material = make_shared<dielectric>(1.5);
                    World.add(make_shared<sphere>(center, 0.2, sphere_material));
                }
            }
        }
    }

    cam.render(bvhNode(make_shared<primitiveList>(World)));
	return 0;
}