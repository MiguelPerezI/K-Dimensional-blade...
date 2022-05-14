using namespace std;

#include "QuaternionBoxBox.hpp"
#include <cmath>
#include <iomanip>
#include <fstream>


QuaternionBoxBox::QuaternionBoxBox() {};
QuaternionBoxBox::QuaternionBoxBox(const Quaternion& v) {
	
	V = (QuaternionBox *) malloc (1 * sizeof(QuaternionBox));
	m = 1;
	V[0] = QuaternionBox(v);

}

//QuaternionBox::QuaternionBox(const Vector3D& v) {
//
//        V = (Quaternion *) malloc (1 * sizeof(Quaternion));
//        n = 1;
//        V[0] = Quaternion(0, v);
//
//        center = Quaternion(0, v);
//
//}

QuaternionBox QuaternionBoxBox::operator [] (int k) const {
   if (k > m - 1)
      	return V[0];
   else {
	return V[k];
   }
}

int QuaternionBoxBox::getM() const {return m;}

void QuaternionBoxBox::push(const Quaternion& v) {

	if (m == 0) {
	
		V[0] = QuaternionBox(v);
		m = 1;
	} else {
	
		V = (QuaternionBox *) realloc(V, sizeof(QuaternionBox) * (m+1));
		m = m + 1;
		V[m-1] = QuaternionBox(v);
	}
}

void QuaternionBoxBox::push(const Vector3D& v) {

        if (m == 0) {

                V[0] = QuaternionBox(v);
                m = 1;
        } else {

                V = (QuaternionBox *) realloc(V, sizeof(QuaternionBox) * (m+1));
                m = m + 1;
                V[m-1] = QuaternionBox(v);
        }
}

//void QuaternionBox::push(const Vector3D& v) {
//
//        if (n == 0) {
//
//                V[0] = Quaternion(0, v);
//                n = 1;
//        } else {
//
//                V = (Quaternion *) realloc(V, sizeof(Quaternion) * (n+1));
//                n = n + 1;
//                V[n-1] = Quaternion(0, v);
//                updateCenter();
//        }
//}
//
//void QuaternionBox::empty() {
//	
//        n = 0;
//	free(V);
//	V = NULL;
//}
//
//void QuaternionBox::push(double a, double b, double c) {
//
//        if (n == 0) {
//
//                V[0] = Quaternion(0, Vector3D(a, b, c));
//                n = 1;
//        } else {
//
//                V = (Quaternion *) realloc(V, sizeof(Quaternion) * (n+1));
//                n = n + 1;
//                V[n-1] = Quaternion(0, Vector3D(a, b, c));
//                updateCenter();
//        }
//}
//
//void QuaternionBox::replace(int i, const Quaternion& A) {
//
//        V[i] = Quaternion(A);
//        updateCenter();
//}
//
//void QuaternionBox::replace(int i, const Vector3D& A) {
//	
//	V[i] = Quaternion(0, A);
//	updateCenter();
//}
//
//void QuaternionBox::replace(int i, double a, double b, double c) {
//
//        V[i] = Quaternion(0, Vector3D(a, b, c));
//        updateCenter();
//}

ostream& operator << (ostream& os, const QuaternionBoxBox& a) {
	
	os << "\n\n{ \n\n\n" ;
	for (int i = 0; i < a.getM(); i++) {
		os << "List-->(" << i << ")" << a[i] << "\n\n";
	}

	os << "}\n\n";

	return os;
}


