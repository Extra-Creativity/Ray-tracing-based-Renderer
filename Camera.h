#pragma once

#ifndef _CAMERA_H_
#define _CAMERA_H_

#include "Vector.h"
class Camera
{
public:
	Camera(Vector3 init_position, Vector3 init_gaze, Vector3 init_up) : position(init_position), gaze(init_gaze), up(init_up) { gaze.Normalize(); up.Normalize(); };

	Vector3 position;
	Vector3 gaze;
	Vector3 up;
private:

};

#endif // !_CAMERA_H_