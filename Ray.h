#pragma once
#ifndef _RAY_H_
#define _RAY_H_

#include "Vector.h"
struct Ray
{
	Ray(Vector3 init_from, Vector3 init_direction) : from(init_from), direction(Vector3::Normalized(init_direction)) {};
	const Vector3 from;
	const Vector3 direction;
};

#endif // !_RAY_H_