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

int main(int argc, char* argv[])
{
	const uint32_t width = 512, height = 512;
	
	char* output = new char[width*height*3];
	
	std::ofstream out_image;
	out_image.open("output.ppm");
	
	out_image << "P6" << "\n"
	<< width << " "
	<< height << "\n"
	<< 3 << "\n";
	
	out_image.write(output, width * height * 3);
	
	out_image.close();
	
	delete[] output;
	return 0;
}