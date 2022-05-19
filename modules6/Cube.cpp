using namespace std;

#include "Cube.hpp"
#include <cmath>
#include <iomanip>
#include <fstream>

Cube::Cube(double R, const Vector3D& c) {
	
	Vector3D v[8];
	v[0] = Vector3D( R + c.x(), R + c.y(), R + c.z());
	v[1] = Vector3D(-R + c.x(), R + c.y(), R + c.z());
	v[2] = Vector3D(-R + c.x(),-R + c.y(), R + c.z());
	v[3] = Vector3D( R + c.x(),-R + c.y(), R + c.z());

	v[4] = Vector3D( R + c.x(), R + c.y(),-R + c.z());
	v[5] = Vector3D(-R + c.x(), R + c.y(),-R + c.z());
	v[6] = Vector3D(-R + c.x(),-R + c.y(),-R + c.z());
	v[7] = Vector3D( R + c.x(),-R + c.y(),-R + c.z());

	f.push(v[0], v[1], v[2], v[3]); 
	f.push(v[1], v[2], v[6], v[5]); 
	f.push(v[6], v[5], v[4], v[7]); 
	f.push(v[4], v[7], v[3], v[0]); 
	f.push(v[2], v[3], v[7], v[6]); 
	f.push(v[0], v[1], v[5], v[4]); 

}

Cube::Cube() {}

FacetBox Cube::operator [] (int k) const {
	return f;
}
