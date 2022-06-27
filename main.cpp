#include <cmath>
#include <random>
#include <fstream>
#include <iostream>
#include <format>

#include "Vector.h"
#include "Ray.h"
#include "MathExtension.h"
#include "Color.h"
#include "Scene.h"
#include "Object.h"
#include "Camera.h"

Vector3 Shade(Scene& scene, const Ray& ray, int depth) {
    static thread_local std::random_device rd;
    static thread_local std::default_random_engine generator(rd());
    static thread_local std::uniform_real_distribution<double> distribution(0, 1);

    const double stopPossiblity = 0.9;
    const int minDepth = 5;
    const int minRefDepth = 2; // Guarantee as least 2 bounce in refraction, so that we can see obvious effect.

    const auto [ID, distance, into] = scene.SearchIntersect(ray);
    if (ID < 0)
        return Vector3::zero(); // return ambient radiance, here set as 0.

    const Object& intersectedObj = scene.GetObjectByID(ID); 
    Vector3 newFrom = ray.from + ray.direction * distance;
    Vector3 norm = intersectedObj.GetNormal(newFrom), nl =into ? norm : norm * -1;

    Color color = intersectedObj.color;

    if (++depth > minDepth) { 
        if (distribution(generator) < stopPossiblity) // R.R. continue, need to be scaled.
            color = color * (1 / stopPossiblity); 
        else 
            return intersectedObj.emission; // not continue.
    }

    Vector3 radiance = intersectedObj.emission;
    switch (intersectedObj.reflType)
    {
    case Object::ReflType::DIFFUSE:
    {
        // Diffuse, we generate a semisphere coordinate and uniformly sample.
        Vector3 newX, newZ; // norm is newY.
        Vector3 newY = into ? norm : -norm;
        if (iszero(norm.y()) && iszero(norm.z())) [[unlikely]] // then norm is (x, 0, 0); this is rare theoritically.
        {
            newX = Vector3(0, 1, 0);
            newZ = Vector3(0, 0, 1);
        }
        else // we use (1, 0, 0) to generate an orthogonal coordinate; this won't be parallel with norm so we can orthogonalize it.
        {
            newX = Vector3(1, 0, 0);
            newZ = Vector3::Cross(newY, newX);
            newX = Vector3::Cross(newZ, newY);
            newZ.Normalize(); newX.Normalize();
        }

        double theta = distribution(generator) * 2 * std::numbers::pi, phi = distribution(generator) * std::numbers::pi / 2;
        const auto [sinTheta, cosTheta] = sincos(theta);
        const auto [sinPhi, cosPhi] = sincos(phi);

        double x = sinPhi * cosTheta, y = cosPhi, z = sinPhi * sinTheta;
        Vector3 newDir = x * newX + y * newY + z * newZ;
        radiance += color.AbsorbedRadiance(Shade(scene, { newFrom, newDir }, depth));
        break;
    }
    case Object::ReflType::SPECULAR:
    {
        Vector3 newDir = ray.direction - norm * 2 * Vector3::Dot(norm, ray.direction);
        radiance += color.AbsorbedRadiance((Shade(scene ,{ newFrom, newDir }, depth)));
        break;
    }
    case Object::ReflType::REFRACTIVE:
    {
        Vector3 reflectionDir = ray.direction - norm * 2 * Vector3::Dot(norm, ray.direction);

        double inN, outN;
		if (into) {
			inN = 1.0f; outN = intersectedObj.n;
		}
		else { 
			inN = intersectedObj.n; outN = 1.0f; 
		}

		// Test total reflection.
		double relaN = inN / outN;
        double cosInAngle = Vector3::Dot(norm, ray.direction);
		double sinInAngle = std::sqrt(1 - cosInAngle * cosInAngle);
        if (sinInAngle >= 1.0f / relaN) // Total reflection, same to specular.
		{
			radiance += intersectedObj.color.AbsorbedRadiance(Shade(scene, { newFrom, reflectionDir }, depth + 1));
			break;
		}

        double R0 = (inN - outN) / (inN + outN);
		R0 *= R0;
        double R;
        {
            double temp = 1 - std::abs(cosInAngle);
            double tempsqr = temp * temp;
            R = R0 + (1 - R0) * tempsqr * tempsqr * temp;
        }

        // newX and norm is orthogonalized and in the same plane with refracted light.
        Vector3 newX = Vector3::Normalized(into ? ray.direction + cosInAngle * norm: ray.direction - cosInAngle * norm);

        double sinOutAngle = sinInAngle * relaN;
	    double cosOutAngle = sqrt(1 - sinOutAngle * sinOutAngle);
	    if (!into) // norm is inversed.
		    cosOutAngle = -cosOutAngle;
	    Vector3 refractionDir = newX * sinOutAngle - norm * cosOutAngle;

        if (depth > minRefDepth)
        {
            // Choose only one dir in path tracing.
            // We set the stop possibility same to Fresnel, and they will be cancelled in calculation.
            // Note that reflection and refration is quite similar, it's easy to mistake them.
            if (distribution(generator) < R)
            {
                radiance += color.AbsorbedRadiance(Shade(scene, {newFrom, reflectionDir}, depth));
            }
            else
            {
                radiance += color.AbsorbedRadiance(Shade(scene, { newFrom, refractionDir }, depth));
            }
        }
        else
        {
            radiance += color.AbsorbedRadiance(Shade(scene ,{newFrom, reflectionDir}, depth) * R + Shade(scene, { newFrom, refractionDir }, depth) * (1 - R));
        }
        break;
    }
    default:
        break;
    }

    return radiance;
}

Vector3 GetRandomPosition(const Vector3& centerPosition, const Vector3 scale = Vector3(1.0, 1.0, 1.0))
{
    static std::random_device rd;
    static std::default_random_engine generator(rd());
    static std::uniform_real_distribution<double> distribution(-5, 5);

    Vector3 offset{ distribution(generator), distribution(generator), distribution(generator) };
    return centerPosition + Vector3::ElemwiseMult(offset, scale);
}

int main() {
    Scene scene;
    {
        // prevent the lite from being too bright, so y is scaled at 0.01.
        scene.AddObject(new Sphere{ GetRandomPosition({50, 681.6 - 0.27, 81.6}, {1, 0.01, 1}), {12, 12, 12}, 600, {0, 0, 0}, Object::ReflType::DIFFUSE }); // White Lite
        scene.AddObject(new Sphere{ GetRandomPosition({1.01 - 600, 40.8, 81.6}, {0, 1, 1}), {30, 10, 0}, 600, {0, 0, 0}, Object::ReflType::DIFFUSE });
        scene.AddObject(new Sphere{ GetRandomPosition({98.99 + 600, 40.8, 81.6}, {0, 1, 1}), {30, 10, 0}, 600, {0, 0, 0}, Object::ReflType::DIFFUSE });
        scene.AddObject(new Sphere{ GetRandomPosition({27, 16.5, 47}), {0, 0, 0},  16.5, {254, 254, 254}, Object::ReflType::SPECULAR }); // White specular object
        scene.AddObject(new Sphere{ GetRandomPosition({73, 16.5, 78}), {0, 0, 0}, 16.5, {254, 254, 254}, Object::ReflType::REFRACTIVE, 4.0f / 3.0f }); // With water refraction index white object.
        {
            // make a combination geometry as a Calabash.
            Vector3 Offset = GetRandomPosition({ 0, 0, 0 });
            scene.AddObject(new Sphere{ Vector3{60, 30, 15} + Offset, {0, 0, 0}, 25, {64, 192, 64}, Object::ReflType::DIFFUSE }); // Green Diffuse ball.
            scene.AddObject(new Sphere{ Vector3{60, 50, 15} + Offset, {0, 0, 0}, 10, {120, 0, 240}, Object::ReflType::DIFFUSE }); // Purple Diffuse ball.
        }
        // We only fix the position of the wall.
        scene.AddObject(new Sphere{ {99 - 1e5, 40.8, 81.6}, {0, 0, 0}, 1e5, {192, 64, 64}, Object::ReflType::DIFFUSE }); // Right border; huge raidus to make it flat as rectangle.
        scene.AddObject(new Sphere{ {1 + 1e5, 40.8, 81.6}, {0, 0, 0}, 1e5, {64, 64, 192}, Object::ReflType::DIFFUSE }); // Left border
        scene.AddObject(new Sphere{ {50, 81.6 - 1e5, 81.6}, {0, 0, 0}, 1e5, {192, 192, 192}, Object::ReflType::DIFFUSE }); // Bottom
        scene.AddObject(new Sphere{ {50, 1e5, 81.6}, {0, 0, 0}, 1e5, {192, 192, 192}, Object::ReflType::DIFFUSE }); // Top
        scene.AddObject(new Sphere{ {50, 40.8, 1e5}, {0, 0, 0}, 1e5, {192, 192, 192}, Object::ReflType::DIFFUSE }); // Front
        scene.AddObject(new Sphere{ {50, 40.8, 170 - 1e5}, {0, 0, 0}, 1e5, {0, 0, 0}, Object::ReflType::DIFFUSE }); // Behind
    };

    double screenDepth = 140.0;

    int width = 102, height = 76, rayAmount = 1000;
    int edgeSample = 4;
    std::cout << std::format("Configuration : {0} * {0} MSAA, {1} ray samples per subpixel.", edgeSample, rayAmount) << std::endl;
    std::cout << std::format("Totally {0} lines to render.", height) << std::endl;
    
    Camera camera({ 50, 44, 295.6 }, { 0, 0, -1 }, { 0, 1, 0});
    Vector3 yUnit = camera.up * 0.5135, xUnit = Vector3::Normalized(Vector3::Cross(camera.gaze, camera.up)) * (width * 0.5135 / height);
    Vector3 radiance;
    Color* colorMap = new Color[width * height];

#pragma omp parallel for schedule(dynamic, 1) private(radiance) 
    for (int row = 0; row < height; row++) {   // for each row
        std::cerr << std::format("\rFinish Line {:2d}", row + 1);
        std::random_device rd;
        std::default_random_engine generator;
        std::uniform_real_distribution<double> distr(0.0, 1.0);
        for (int col = 0; col < width; col++)   // for each col
        {
            int i = (height - row - 1) * width + col;
            radiance = Vector3::zero();
            for (int suby = 0; suby < edgeSample; suby++) // antialias by MSAA.
            {
                for (int subx = 0; subx < edgeSample; subx++) // for each subpixel
                { 
                    for (int _ = 0; _ < rayAmount; _++)
                    {
                        // uniformly generate sample offset.
                        double xOffset = 2 * distr(generator) - 1, yOffset = 2 * distr(generator) - 1;
                        // move to the center -> give offset -> arrive at rive pixel -> uniform.
                        double currX = ((subx + 0.5 + xOffset) / edgeSample + col) / width - 0.5, currY = ((suby + 0.5 + yOffset) / edgeSample + row) / height - 0.5;
                        Vector3 direction = xUnit * currX + yUnit * currY + camera.gaze;
                        // to avoid overflow, we divide rayAmount here.
                        radiance += Shade(scene, { camera.position + direction * screenDepth, direction }, 0) / rayAmount;
                    }
                }
            }
            radiance /= (edgeSample * edgeSample);
            colorMap[i] = Color::GammaCorrection(radiance);
        }
    }
	std::ofstream fout("image.ppm");
	fout << std::format("P3\n{0} {1}\n{2}\n", width, height, 255);
    for (int i = 0; i < width * height; i++)
    {
        Color& currColor = colorMap[i];
        fout << std::format("{0} {1} {2} ", int(currColor.r()), int(currColor.g()), int(currColor.b()));
    }
	fout.close();
	delete[] colorMap;
    return 0;
}