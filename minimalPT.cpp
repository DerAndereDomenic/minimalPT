// Author: Domenic Zingsheim
// Date: 12.06.2021
// Build:
//		g++ -o minimalPT minimalPT.cpp

#include <iostream>
#define PI 3.1415926535f
#include <cmath>
#include <fstream>
#include <ctime>

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

struct Sphere
{
	Vec3 position, albedo;
	float radius;
	MaterialType type;
	
	Sphere(const Vec3& position, const float& radius, const Vec3& albedo, const MaterialType& type)
		:position(position), radius(radius), albedo(albedo), type(type)
	{}
	
	bool inline
	intersect(Ray& ray, Vec3& output_intersection, float& output_depth)
	{
		Vec3 OS = ray.origin - position;
		float b = 2.0f * ray.direction.dot(OS);
		float c = OS.dot(OS) - radius * radius;
		float disc = b * b - 4 * c;
		
		if(disc > 0)
		{
			float distSqrt = sqrtf(disc);
			float q = b < 0 ? (-b - distSqrt) / 2.0f : (-b + distSqrt) / 2.0f;
			float t0 = q;
			float t1 = c/q;
			
			t0 = fminf(t0,t1);
			t1 = fmaxf(t0,t1);
			
			if(t1 >= 0)
			{
				float t = t0 < 0 ? t1 : t0;
				output_intersection = ray.origin + ray.direction*t;
				output_depth = t;
				return true;
			}
		}
		return false;
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

#define NUM_SPHERES 2
Sphere spheres[NUM_SPHERES] = 
{
	Sphere(Vec3(-0.6f,-0.7f,2.6f), 0.3f, Vec3(1,1,1), DIFFUSE),
	Sphere(Vec3(0.6f,-0.7f,1.7f), 0.3f, Vec3(1,1,1), DIFFUSE)
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

float inline 
rnd()
{
	return ((float)rand())/((float)(RAND_MAX));
}

Vec3 sampleBSDF(Vec3& N)
{
	float xi_1 = rnd();
	float xi_2 = rnd();
	
	float r = sqrtf(xi_1);
	float phi = 2.0f * PI * xi_2;
	
	float x = r*cos(phi);
	float y = r*sin(phi);
	float z = sqrtf(fmaxf(0.0f, 1.0f - x * x - y * y));
	
	float sz = (N.z >= 0) ? 1 : -1;
	float a = 1.0f / (sz + N.z);
	float ya = N.y * a;
	float b = N.x * ya;
	float c = N.x * sz;
	
	Vec3 localX(c * N.x * a - 1.0f, sz * b, c);
	Vec3 localY(b, N.y * ya - sz, N.y);
	
	return localX * x + localY * y + N * z;
}

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

	for(uint32_t i = 0; i < NUM_SPHERES; ++i)
	{
		float depth;
		Vec3 intersection;
		if(spheres[i].intersect(ray, intersection, depth))
		{
			if(depth < geom.depth)
			{
				geom.depth = depth;
				geom.P = intersection;
				geom.N = intersection - spheres[i].position;
				geom.N = geom.N / geom.N.norm();
				geom.type = spheres[i].type;
				geom.albedo = spheres[i].albedo;
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
						 1.5f);
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
		
		Ray shadow_ray;
		shadow_ray.origin = geom.P + light_direction*0.01f;
		shadow_ray.direction = light_direction;
		
		LocalGeometry shadow_geom = computeSceneIntersection(shadow_ray);
		float visibility = 1.0f;
		if(shadow_geom.depth < light_dist - 0.01f)
			visibility = 0.0f;
		
		radiance = geom.albedo * light_radiance / PI * fmaxf(0.0f, geom.N.dot(light_direction)) * visibility;
	}
	
	return radiance;
}

int main(int argc, char* argv[])
{
	srand((unsigned int)time(NULL));
	const uint32_t width = 512, height = 512;
	const uint32_t n_samples = 1;
	
	uint8_t* output = new uint8_t[width*height*3];
	
	for(uint32_t x = 0; x < width; ++x)
	{
		for(uint32_t y = 0; y < height; ++y)
		{
			Vec3 radiance(0,0,0);
			for(uint32_t i = 0; i < n_samples; ++i)
			{
				radiance = radiance + estimateRadiance(x,width - y - 1,width,height) / n_samples;
			}
			
			float mapped_x = powf(1.0f - expf(-radiance.x), 1.0f/2.2f);
			float mapped_y = powf(1.0f - expf(-radiance.y), 1.0f/2.2f);
			float mapped_z = powf(1.0f - expf(-radiance.z), 1.0f/2.2f);
			output[3*(y * width + x) + 0] = static_cast<uint8_t>(fmaxf(0,fminf(mapped_x, 1.0f))*255.0f); //R
			output[3*(y * width + x) + 1] = static_cast<uint8_t>(fmaxf(0,fminf(mapped_y, 1.0f))*255.0f); //G
			output[3*(y * width + x) + 2] = static_cast<uint8_t>(fmaxf(0,fminf(mapped_z, 1.0f))*255.0f); //B
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