//////////////////////////////////////////////////////////////////////////
// Vector and Matrix Class
//
// Created by: Matt Camuto
//
// Modification:
//		Time            Programmer
//		07-09-99        R. Wenger
//		08-20-99        K. Boner
//		07-15-00        J. Gao
// 		05-25-05		H-W Shen, Liya Li
//////////////////////////////////////////////////////////////////////////

#ifndef _VECTOR_MATRIX_H_
#define _VECTOR_MATRIX_H_

#include <math.h>
#include "common.h"
#include "Tuple.h"

namespace edda{

//////////////////////////////////////////////////////////////////////////
// 2d vector used for Point
//////////////////////////////////////////////////////////////////////////
template<typename Real>
class VECTOR2 : public Tuple2<Real> {
public :
    VECTOR2()  { Zero(); }
    inline void Zero() { this->vec[0] = this->vec[1] = 0; }
};

//////////////////////////////////////////////////////////////////////////
//	vector with 3 components
//////////////////////////////////////////////////////////////////////////
template <typename Real>
class VECTOR3 : public Tuple3<Real> {
    using Tuple3<Real>::vec;
public :
	// constructor
	VECTOR3() {Zero();}

    inline void Zero() {vec[0] = vec[1] = vec[2] = 0.0;}					// make zero vector
    inline Real x() const {return vec[0];}
    inline Real y() const {return vec[1];}
    inline Real z() const {return vec[2];}

    bool operator ==(const VECTOR3& v) const
	{
        return ((fabs(vec[0]-v(0))< limits<Real>::epsilon()) &&
                (fabs(vec[1]-v(1))< limits<Real>::epsilon()) &&
                (fabs(vec[2]-v(2))< limits<Real>::epsilon()));
    }
    // get magnitude
    double GetMag() const {
        return sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
    }

    // get the maximum value
    float GetMax() const {
        Real maxval = vec[0];
        if (vec[1] > maxval) maxval = vec[1];
        if (vec[2] > maxval) maxval = vec[2];
        return maxval;
    }

    // make sure all component<=1.0
    void Clamp()
    {
        for (int i = 0; i < N; i++)
            if (vec[i]>1.0) vec[i] = 1.0;
    }
	void Normalize();					                            // normalize vector
	bool IsSame(VECTOR3& a);
	void scale(const float s);
	void minus(VECTOR3& v1, VECTOR3& v2);
	void add(float a, float b, float c) { vec[0]+= a; vec[1] += b; vec[2] += c; }
};

//////////////////////////////////////////////////////////////////////////
// 4d vector class
//////////////////////////////////////////////////////////////////////////
template <typename Real>
class VECTOR4
{
    using Tuple3<Real>::vec;
public :
	// constructor
	VECTOR4() {Zero();}
	// MOD-BY-LEETEN 02/07/2011-FROM:
		// VECTOR4(VECTOR3 v) {vec[0] = v[0]; vec[1] = v[1]; vec[2] = v[2]; vec[3] = 1.0;}
	// TO:
	VECTOR4(const VECTOR3& v) {vec[0] = v[0]; vec[1] = v[1]; vec[2] = v[2]; vec[3] = 1.0;}
	VECTOR4(const VECTOR4& v) {vec[0] = v[0]; vec[1] = v[1]; vec[2] = v[2]; vec[3] = v[3];}
	// MOD-BY-LEETEN 02/07/2011-END
	VECTOR4(const float x0, const float x1, const float x2, const float x3)
	{vec[0] = x0; vec[1] = x1; vec[2] = x2; vec[3] = x3;}

	void Set(const float x0, const float x1, const float x2, const float x3)
	{vec[0] = x0; vec[1] = x1; vec[2] = x2; vec[3] = x3;}

	int Dimension() const {	return 4;}
	float & operator [](const int i) { return(vec[i]);}			// index i'th element
	// ADD-BY-LEETEN 02/07/2011-BEGIN
	const float & operator [](const int i) const { return(vec[i]);}			// index i'th element
	// ADD-BY-LEETEN 02/07/2011-END
	float operator ()(const int i) const {return(vec[i]);}		// return i'th element
	const VECTOR4 & operator =(const VECTOR4 & v0)				// copy vector v0
	{vec[0] = v0(0); vec[1] = v0(1); vec[2] = v0(2); vec[3] = v0(3);return(*this);}
	const VECTOR4 & operator =(const VECTOR3 & v0)				// copy vector v0
	{vec[0] = v0(0); vec[1] = v0(1); vec[2] = v0(2); vec[3] = 1.0; return(*this);}
	VECTOR3 get_vector3() {VECTOR3 temp(vec[0]/vec[3], vec[1]/vec[3], vec[2]/vec[3]); return temp;}
	void Zero() {vec[0] = vec[1] = vec[2] = vec[3] = 0.0;}		// make zero vector

	void Normalize();											// normalize vector
};

//////////////////////////////////////////////////////////////////////////
// 3d matrix
//////////////////////////////////////////////////////////////////////////
class MATRIX3
{
private :
	VECTOR3 mat[3];       // a vector represents each matrix row

public :

	MATRIX3()                                    // constructor
	{	Identity(); };
	MATRIX3(const VECTOR3 & v0, const VECTOR3 & v1, const VECTOR3 & v2)
	{	mat[0] = v0; mat[1] = v1; mat[2] = v2; };  // constructor
	int Dimension() const
	{	return 3; };
	VECTOR3 & operator [](const int i)           // index row i
	{	return(mat[i]); };
	// Note: reference row i, column j of MATRIX3 m0 as m0[i][j] (not m0[i,j])
	VECTOR3 operator()(const int i) const        // return row i
	{	return(mat[i]); };
	float operator ()(const int i, const int j) const
	{	return(mat[i](j)); };                    // return element (i,j)
	MATRIX3 & operator =(const MATRIX3 & m0)     // copy matrix m0
	{	mat[0] = m0(0); mat[1] = m0(1); mat[2] = m0(2);
	return(*this); };
	void Identity();                             // set to identity

	float det();	// determinant of the matrix

	int inverse( MATRIX3& m) ;//added by lijie to handle curvilinear grid
	MATRIX3 transpose();//added by lijie to handle curvilinear grid
};

//////////////////////////////////////////////////////////////////////////
// 4d matrix
//////////////////////////////////////////////////////////////////////////
class MATRIX4
{

private :
	VECTOR4 mat[4];       // a vector represents each matrix row

public :

	MATRIX4()                                    // constructor
	{	Identity(); };
	MATRIX4(const VECTOR4 & v0, const VECTOR4 & v1,
		const VECTOR4 & v2, const VECTOR4 & v3) // constructor
	{	mat[0] = v0; mat[1] = v1; mat[2] = v2; mat[3] = v3; };
	int Dimension() const
	{	return 4; };
	VECTOR4 & operator [](int i)                 // index row i
	{	return(mat[i]); };
	// Note: reference row i, column j of MATRIX4 m0 as m0[i][j] (not m0[i,j])
	VECTOR4 operator()(const int i) const        // return row i
	{	return(mat[i]); };
	float operator ()(const int i, const int j) const
	{	return(mat[i](j)); };                    // return element (i,j)
	MATRIX4 & operator =(const MATRIX4 & m0)     // copy matrix m0
	{	mat[0] = m0(0); mat[1] = m0(1); mat[2] = m0(2); mat[3] = m0(3);
	return(*this); };
	MATRIX4 & operator =(const MATRIX3 & m0)     // copy matrix m0
	{	mat[0] = m0(0); mat[1] = m0(1); mat[2] = m0(2);
	VECTOR4 temp(0.0,0.0,0.0,1.0);
	mat[3] = temp;
	return(*this); };
	void Identity();                             // set to identity
};

//************************
// VECTOR2 operations
//************************

inline VECTOR2 operator +(const VECTOR2 & v0, const VECTOR2 & v1)
// return v0 + v1
{	return(VECTOR2(v0(0) + v1(0), v0(1) + v1(1))); };

inline VECTOR2 operator -(const VECTOR2 & v0, const VECTOR2 & v1)
// return v0 - v1
{	return(VECTOR2(v0(0) - v1(0), v0(1) - v1(1))); };

inline VECTOR2 operator *(float x0, const VECTOR2 & v0)
// return x0*v0
{	return(VECTOR2(x0*v0(0), x0*v0(1))); };

inline VECTOR2 operator *(const VECTOR2 & v0, float x0)
// return v0*x0 (= x0*v0)
{	return(x0*v0); };


//************************
// VECTOR3 operations
//************************

inline float dot(const VECTOR3 & v0, const VECTOR3 & v1)
// return dot product of v0 and v1
{	return(v0(0)*v1(0) + v0(1)*v1(1) + v0(2)*v1(2)); };

inline VECTOR3 cross(const VECTOR3 & v0, const VECTOR3 & v1)
// return cross product of v0 and v1
{	return(VECTOR3(v0(1)*v1(2) - v0(2)*v1(1),
		   v0(2)*v1(0) - v0(0)*v1(2),
		   v0(0)*v1(1) - v0(1)*v1(0))); };

inline VECTOR3 operator +(const VECTOR3 & v0, const VECTOR3 & v1)
// return v0 + v1
{	return(VECTOR3(v0(0) + v1(0), v0(1) + v1(1), v0(2) + v1(2))); };

inline VECTOR3 operator -(const VECTOR3 & v0, const VECTOR3 & v1)
// return v0 - v1
{	return(VECTOR3(v0(0) - v1(0), v0(1) - v1(1), v0(2) - v1(2))); };

inline float operator *(VECTOR3& v0, VECTOR3& v1)
// return v0*v1T
{	return (v0(0) * v1(0) + v0(1) * v1(1) + v0(2) * v0(2)); };

inline VECTOR3 operator *(float x0, const VECTOR3 & v0)
// return x0*v0
{	return(VECTOR3(x0*v0(0), x0*v0(1), x0*v0(2))); };

inline VECTOR3 operator *(const VECTOR3 & v0, float x0)
// return v0*x0 (= x0*v0)
{	return(x0*v0); };

inline VECTOR4 get_vector4(VECTOR3 vec)
{	VECTOR4 temp(vec[0], vec[1], vec[2], 1.0); return temp;};

//************************
// VECTOR4 operations
//************************

inline float dot(const VECTOR4 & v0, const VECTOR4 & v1)
// return dot product of v0 and v1
{	return(v0(0)*v1(0) + v0(1)*v1(1) + v0(2)*v1(2) + v0(3)*v1(3)); };

inline VECTOR4 operator +(const VECTOR4 & v0, const VECTOR4 & v1)
// return v0 + v1
{	return(VECTOR4(v0(0)+v1(0), v0(1)+v1(1), v0(2)+v1(2), v0(3)+v1(3))); };

inline VECTOR4 operator -(const VECTOR4 & v0, const VECTOR4 & v1)
// return v0 - v1
{	return(VECTOR4(v0(0)-v1(0), v0(1)-v1(1), v0(2)-v1(2), v0(3)-v1(3))); };

inline VECTOR4 operator *(float x0, const VECTOR4 & v0)
// return x0*v0
{	return(VECTOR4(x0*v0(0), x0*v0(1), x0*v0(2), x0*v0(3))); };

inline VECTOR4 operator *(const VECTOR4 & v0, float x0)
// return v0*x0 (= x0*v0)
{	return(x0*v0); };

//************************
// MATRIX3 operations
//************************

MATRIX3 operator +(const MATRIX3 & m0, const MATRIX3 & m1); // return m0 + m1
MATRIX3 operator -(const MATRIX3 & m0, const MATRIX3 & m1); // return m0 - m1
MATRIX3 operator *(const MATRIX3 & m0, const MATRIX3 & m1); // return m0 * m1
MATRIX3 operator *(const float x0, const MATRIX3 & m0);    // return x0 * m0
MATRIX3 operator *(const MATRIX3 & m0, const float x0);    // return m0 * x0
VECTOR3 operator *(const MATRIX3 & m0, const VECTOR3 & v0); // return m0 * v0
VECTOR3 operator *(const VECTOR3 & v0, const MATRIX3 & m0); // return v0 * m0

//************************
// MATRIX4 operations
//************************

MATRIX4 operator +(const MATRIX4 & m0, const MATRIX4 & m1); // return m0 + m1
MATRIX4 operator -(const MATRIX4 & m0, const MATRIX4 & m1); // return m0 - m1
MATRIX4 operator *(const MATRIX4 & m0, const MATRIX4 & m1); // return m0 * m1
MATRIX4 operator *(const float x0, const MATRIX4 & m0);    // return x0 * m0
MATRIX4 operator *(const MATRIX4 & m0, const float x0);    // return m0 * x0
VECTOR4 operator *(const MATRIX4 & m0, const VECTOR4 & v0); // return m0 * v0
VECTOR4 operator *(const VECTOR4 & v0, const MATRIX4 & m0); // return v0 * m0
VECTOR3 operator *(const MATRIX4 & m0, const VECTOR3 & v0); // return m0 * v0
VECTOR3 operator *(const VECTOR3 & v0, const MATRIX4 & m0); // return v0 * m0


MATRIX4 inverse(const MATRIX4 & m);  // return inverse of m; return 0 matrix if
// m is singular
MATRIX4 rotate_matrix(int type, float angle); // type: 1:x, 2:y, 3:z
MATRIX4 translate_matrix(float dx, float dy, float dz);
MATRIX4 scale_matrix(float sx, float sy, float sz);

} // namespace edda

#endif
