#pragma once
#ifndef _SCENE_H_
#define _SCENE_H_

#include "Object.h"
#include "Ray.h"
#include <vector>
#include<iostream>
class Scene
{
public:
	Scene() =default;

	void AddObject(Object* newObjPtr)
	{
		objects.emplace_back(newObjPtr);
	}
	/// <summary>
	/// Search Ray's intersection with object.
	/// </summary>
	/// <returns>If not finding intersection, return negative number; else return object ID and distance.</returns>
	std::tuple<int, double, bool> SearchIntersect(const Ray& ray) const
	{
		double finalDistance = 1e20;
		int finalID = -1;
		bool finalInside = false;
		int size = static_cast<int>(objects.size());
		for(int i = 0; i < size; i++)
		{
			const auto[success, inside, distance] = objects[i]->Intersect(ray);
			if (success && (distance < finalDistance))
			{
				finalDistance = distance;
				finalID = i;
				finalInside = inside;
			}
		}
		return { finalID, finalDistance, finalInside};
	}

	Object& GetObjectByID(int id) { return *objects[id]; }
private:
	std::vector<std::unique_ptr<Object>> objects;
};

#endif // !_SCENE_H_