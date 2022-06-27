#pragma once
#ifndef _OBJECT_H_
#define _OBJECT_H_
#include "Vector.h"
#include "Color.h"
#include "Ray.h"
#include <iostream>
class Object {
public:
	enum class ReflType {DIFFUSE, SPECULAR, REFRACTIVE};
	// act as refractive index.
	double n;
	ReflType reflType;
	Color color;
	Vector3 emission;

	Object(Vector3 init_color, ReflType init_type, Vector3 init_emission, double init_n = 1) :color(init_color), reflType(init_type), n(init_n), emission(init_emission) {};
	virtual ~Object() = default;
	/// <summary>
	/// Every Object should have a function to judge Intersection.
	/// </summary>
	/// <param name="ray">ray that want to judge intersection with object.</param>
	/// <returns> [success, go_inside, distance]
	/// Specifically, whether a ray and an object intersect; If they intersect, with the distance from ray beginning to the intersection position.
	/// inside is if it's refraction, then whether in or out.
	/// </returns>
	virtual std::tuple<bool, bool, double> Intersect(const Ray& ray) const = 0;
	virtual Vector3 GetNormal(const Vector3& position) const = 0;
};

class Sphere : public Object
{
public:
	Sphere(Vector3 init_position, Vector3 init_emission, double init_radius, Vector3 init_color, ReflType init_type, double init_n = 1) :position(init_position), \
		radius(init_radius), Object(init_color, init_type, init_emission, init_n){};

	Vector3 position;
	double radius;	
	virtual std::tuple<bool, bool, double> Intersect(const Ray& ray) const override {
		Vector3 disVec = position - ray.from;
		double sqrLength = disVec.Magnitude();
		// Here we demand ray.direction is unit.
		// cosLen may be negative when the direction is reversed. This guarantees rightness of following code.
		double cosLen = Vector3::Dot(ray.direction, disVec);
		double sinLenSqr = sqrLength - cosLen * cosLen;
		double radiusSqr = radius * radius;
		if (sinLenSqr > radiusSqr) // also sinLen > radius
			return { false, false, cosLen };
		else
		{
			double temp = std::sqrt(radiusSqr - sinLenSqr);
			double sol1 = cosLen - temp, sol2 = cosLen + temp;
			if (sol1 < 1e-4) // from inside or half-line not intersect or nead border.
			{
				if (sol2 < 1e-4) // near or out of the border and go outside, we regard it as not intersecting.
					return { false, false, 0.f };
				return { true, false, sol2 }; // else inside the sphere and it wants to go outside in the next bounce.
			}
			return { true, true, sol1 }; // sol1 > 1e-4, go inside.
		}
	};

	virtual Vector3 GetNormal(const Vector3& intersectedPos) const override { return Vector3::Normalized(intersectedPos - position); };
};

#endif // !_MATERIAL_H_