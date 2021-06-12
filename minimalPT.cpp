// Author: Domenic Zingsheim
// Date: 12.06.2021
// Build:
//		g++ -o minimalPT minimalPT.cpp

#include <iostream>
#include <cmath>
#include <fstream>

//
//	MATH FUNCTIONS
//

struct Vec3
{
	float x = 0,y = 0,z = 0;
	
	Vec3() = default;
	
	Vec3(const float&x, const float& y, const float& z)
		:x(x),y(y),z(z)
		{}
	
	Vec3 inline
	operator+(const Vec3& other)
	{
		return Vec3(x + other.x, y + other.y, z + other.z);
	}
	
	Vec3 inline 
	operator-(const Vec3& other)
	{
		return Vec3(x - other.x, y - other.y, z - other.z);
	}
	
	Vec3 inline 
	operator-()
	{
		return Vec3(-x, -y, -z);
	}
	
	Vec3 inline 
	operator*(const Vec3& other)
	{
		return Vec3(x * other.x, y * other.y, z * other.z);
	}
	
	Vec3 inline
	operator*(const float& l)
	{
		return Vec3(x*l, y*l, z*l);
	}
	
	float inline
	norm()
	{
		return sqrtf(x*x + y*y + z*z);
	}
	
	float inline
	dot(const Vec3& other)
	{
		return x * other.x + y * other.y + z * other.z;
	}
	
};

//
//	SCENE
//

struct Ray
{
	Vec3 origin, direction;
	
	Ray() = default;
	Ray(const Vec3& origin, const Vec3& direction)
		:origin(origin), direction(direction)
	{}
};

enum MaterialType
{
	DIFFUSE,
	MIRROR,
	REFLECT
};

struct Plane
{
	Vec3 position, normal, albedo;
	MaterialType type;
	
	Plane(const Vec3& position, const Vec3& normal, const Vec3& albedo, const MaterialType& type)
		:position(position), normal(normal), albedo(albedo), type(type)
	{}
	
	bool inline
	intersect(Ray& ray, Vec3& output_direction, float& output_depth)
	{
		float PdotN = position.dot(normal);
		float OdotN = ray.origin.dot(normal);
		float DdotN = ray.direction.dot(normal);
		
		if(fabs(DdotN) < 1e-5f) return false;
		
		float t = (PdotN - OdotN) / DdotN;
		if(t < 0)return false;
		
		output_direction = ray.origin + ray.direction*t;
		output_depth = t;
		
		return true;
	}
};

//
//	RENDERING
//
Vec3
estimateRadiance(const uint32_t& x, const uint32_t& y, const uint32_t& width, const uint32_t& height)
{
	float r = (x+y)%256;
	return Vec3(r,r,r);
}

int main(int argc, char* argv[])
{
	const uint32_t width = 512, height = 512;
	
	char* output = new char[width*height*3];
	
	for(uint32_t x = 0; x < width; ++x)
	{
		for(uint32_t y = 0; y < height; ++y)
		{
			Vec3 radiance = estimateRadiance(x,y,width,height);
			output[3*(y * width + x) + 0] = radiance.x%256;
			output[3*(y * width + x) + 1] = radiance.y%256;
			output[3*(y * width + x) + 2] = radiance.z%256;
		}
	}
	
	std::ofstream out_image;
	out_image.open("output.ppm");
	
	out_image << "P6" << "\n"
	<< width << " "
	<< height << "\n"
	<< 255 << "\n";
	
	out_image.write(output, width * height * 3);
	
	out_image.close();
	
	delete[] output;
	return 0;
}