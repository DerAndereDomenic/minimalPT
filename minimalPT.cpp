// Author: Domenic Zingsheim
// Date: 12.06.2021
// Build:
//		g++ -o minimalPT minimalPT.cpp

#include <iostream>
#define PI 3.1415926535f
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
	
	Vec3 inline
	operator/(const float& l)
	{
		return Vec3(x/l, y/l, z/l);
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
//	GEOMETRY
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
	intersect(Ray& ray, Vec3& output_intersection, float& output_depth)
	{
		float PdotN = position.dot(normal);
		float OdotN = ray.origin.dot(normal);
		float DdotN = ray.direction.dot(normal);
		
		if(fabs(DdotN) < 1e-5f) return false;
		
		float t = (PdotN - OdotN) / DdotN;
		if(t < 0)return false;
		
		output_intersection = ray.origin + ray.direction*t;
		output_depth = t;
		
		return true;
	}
};

//
//	SCENE
//
#define NUM_PLANES 5
Plane planes[NUM_PLANES] = 
{
	Plane(Vec3(0,-1,0), Vec3(0,1,0), Vec3(1,1,1), DIFFUSE),
	Plane(Vec3(0,1,0), Vec3(0,-1,0), Vec3(1,1,1), DIFFUSE),
	Plane(Vec3(-1,0,0), Vec3(1,0,0), Vec3(0,1,0), DIFFUSE),
	Plane(Vec3(1,0,0), Vec3(-1,0,0), Vec3(1,0,0), DIFFUSE),
	Plane(Vec3(0,0,3), Vec3(0,0,-1), Vec3(1,1,1), DIFFUSE)
};

//Light
Vec3 light_position(0,0.9f,2.0f);
Vec3 light_intensity(1,1,1);

//
//	RENDERING
//
struct LocalGeometry
{
	Vec3 P;
	Vec3 N;
	Vec3 albedo;
	float depth = INFINITY;
	MaterialType type;
};

LocalGeometry inline 
computeSceneIntersection(Ray& ray)
{
	LocalGeometry geom;
	for(uint32_t i = 0; i < NUM_PLANES; ++i)
	{
		float depth;
		Vec3 intersection;
		if(planes[i].intersect(ray, intersection, depth))
		{
			if(depth < geom.depth)
			{
				geom.depth = depth;
				geom.P = intersection;
				geom.N = planes[i].normal;
				geom.type = planes[i].type;
				geom.albedo = planes[i].albedo;
			}
		}
	}
	
	return geom;
}

Vec3
estimateRadiance(const uint32_t& x, const uint32_t& y, const uint32_t& width, const uint32_t& height)
{
	Ray ray;
	
	ray.origin = Vec3(0,0,0);
	ray.direction = Vec3(2.0f*static_cast<float>(x)/static_cast<float>(width) - 1.0f,
						 2.0f*static_cast<float>(y)/static_cast<float>(height) - 1.0f,
						 1);
	ray.direction = ray.direction * (1.0f/ray.direction.norm());
						 
	//Compute intersection
	LocalGeometry geom = computeSceneIntersection(ray);
	
	Vec3 radiance(0,0,0);
	if(geom.depth != INFINITY)
	{
		Vec3 light_direction = light_position - geom.P;
		float light_dist = light_direction.norm();
		Vec3 light_radiance = light_intensity / (light_dist * light_dist);
		light_direction = light_direction / light_dist;
		
		radiance = geom.albedo * light_radiance / PI * fmaxf(0.0f, geom.N.dot(light_direction));
	}
	
	return radiance;
}

int main(int argc, char* argv[])
{
	const uint32_t width = 512, height = 512;
	
	uint8_t* output = new uint8_t[width*height*3];
	
	for(uint32_t x = 0; x < width; ++x)
	{
		for(uint32_t y = 0; y < height; ++y)
		{
			Vec3 radiance = estimateRadiance(x,width - y - 1,width,height);
			output[3*(y * width + x) + 0] = static_cast<uint8_t>(fminf(radiance.x, 1.0f)*255.0f); //R
			output[3*(y * width + x) + 1] = static_cast<uint8_t>(fminf(radiance.y, 1.0f)*255.0f); //G
			output[3*(y * width + x) + 2] = static_cast<uint8_t>(fminf(radiance.z, 1.0f)*255.0f); //B
		}
	}
	
	std::ofstream out_image;
	out_image.open("output.ppm", std::ios::binary);
	
	out_image << "P6" << "\n"
	<< width << " "
	<< height << "\n"
	<< 255 << "\n";
	
	out_image.write((char*)output, width * height * 3);
	
	out_image.close();
	
	delete[] output;
	return 0;
}