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

    auto checker = make_shared<checkerBoard>(0.82, vec3(.2, .3, .1), vec3(.9, .9, .9));
	World.add(make_shared<sphere>(point3(0, -1000, 0), 1000, make_shared<lambertian>(checker)));
	auto material1 = make_shared<dielectric>(1.5);
	World.add(make_shared<sphere>(point3(0, 1, 0), 1.0, material1));
	auto material2 = make_shared<lambertian>(vec3(0.4, 0.2, 0.1));
	World.add(make_shared<sphere>(point3(-4, 1, 0), 1.0, material2));
	auto material3 = make_shared<metal>(vec3(0.7, 0.6, 0.5), 0.0);
	World.add(make_shared<sphere>(point3(4, 1, 0), 1.0, material3));

    cam.render(bvhNode(make_shared<primitiveList>(World)));
	return 0;
}