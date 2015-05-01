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

#include <cmath>
#include "common.h"
#include "core/tuple.h"

namespace edda{

//////////////////////////////////////////////////////////////////////////
/// for perfomrnace we implement each VECTOR with different lengths
//	vector with 3 components
//////////////////////////////////////////////////////////////////////////
template <typename Real>
class Vector3 : public Tuple<Real, 3>
{
    using Tuple<Real, 3>::vec;
public :
    Vector3() : Tuple<Real, 3>((Real)0) {}
    Vector3(Real x, Real y, Real z) { set(x,y,z); }
	// constructor
    inline Real x() const {return vec[0];}
    inline Real y() const {return vec[1];}
    inline Real z() const {return vec[2];}
    inline void set(Real x, Real y, Real z) { vec[0]=x; vec[1]=y; vec[2]=z; }

    bool operator ==(const Vector3& v) const {
        return (fabs(vec[0]-v(0)) < limits<Real>::epsilon() &&
                fabs(vec[1]-v(1)) < limits<Real>::epsilon() &&
                fabs(vec[2]-v(2)) < limits<Real>::epsilon());
    }
    // get magnitude
    double getMag() const { return sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);  }
    // get the maximum value
    Real getMax() const { return std::max(vec[0], std::max(vec[1], vec[2])); }
    // normalize vector
    void normalize()    { double norm = getMag(); if (norm!=0) (*this) *= (1/norm); }
    Vector3<Real> operator+=(double x) { vec[0]+=x; vec[1]+=x; vec[2]+=x; return *this; }
    Vector3<Real> operator*=(double x) { vec[0]*=x; vec[1]*=x; vec[2]*=x; return *this; }

    // make sure all component<=1.0
    void clamp()
    {
        for (int i = 0; i < this->getLen(); i++)
            if (vec[i]>1.0) vec[i] = 1.0;
    }
};


//////////////////////////////////////////////////////////////////////////
/// for perfomrnace we implement each VECTOR with different lengths
//	vector with 4 components
//////////////////////////////////////////////////////////////////////////
template <typename Real>
class Vector4 : public Tuple<Real, 4>
{
    using Tuple<Real, 4>::vec;
public :
    Vector4() : Tuple<Real, 4>((Real)0) {}
    Vector4(Real x, Real y, Real z, Real w) { set(x,y,z,w); }
    // constructor
    inline Real x() const {return vec[0];}
    inline Real y() const {return vec[1];}
    inline Real z() const {return vec[2];}
    inline Real w() const {return vec[3];}
    inline void set(Real x, Real y, Real z, Real w) { vec[0]=x; vec[1]=y; vec[2]=z; vec[3]=w; }

    bool operator ==(const Vector4<Real>& v) const {
        return (fabs(vec[0]-v(0)) < limits<Real>::epsilon() &&
                fabs(vec[1]-v(1)) < limits<Real>::epsilon() &&
                fabs(vec[2]-v(2)) < limits<Real>::epsilon() &&
                fabs(vec[3]-v(3)) < limits<Real>::epsilon());
    }
    // get magnitude
    double getMag() const { return sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2] + vec[3]*vec[3]);  }
    // get the maximum value
    Real getMax() const { return std::max(vec[0], std::max(vec[1], std::max(vec[2], vec[3]))); }
    // normalize vector
    void normalize()    { double norm = getMag(); if (norm!=0) (*this) *= (1/norm); }
    Vector4<Real> operator+=(double x) { vec[0]+=x; vec[1]+=x; vec[2]+=x; vec[3]+=x; return *this; }
    Vector4<Real> operator*=(double x) { vec[0]*=x; vec[1]*=x; vec[2]*=x; vec[3]+=x; return *this; }

    // make sure all component<=1.0
    void clamp()
    {
        for (int i = 0; i < this->getLen(); i++)
            if (vec[i]>1.0) vec[i] = 1.0;
    }
};
//************************
// VECTOR operations
//************************
template<typename Real> inline Real dot(const Vector3<Real> & v0, const Vector3<Real> & v1)
// return dot product of v0 and v1
{    return( v0[0]*v1[0] + v0[1]*v1[1] + v0[2]*v1[2]); }
template<typename Real> inline Vector3<Real> cross(const Vector3<Real> & v0, const Vector3<Real> & v1)
// return cross product of v0 and v1
{	return Vector3<Real>( v0[1]*v1[2] - v0[2]*v1[1],
                          v0[2]*v1[0] - v0[0]*v1[2],
                          v0[0]*v1[1] - v0[1]*v1[0] );
}
template<typename Real> inline Vector3<Real>  operator +(const Vector3<Real>  & v0, const Vector3<Real>  & v1)
// return v0 + v1
{	return(Vector3<Real> (v0[0] + v1[0], v0[1] + v1[1], v0[2] + v1[2])); }
template<typename Real> inline Vector3<Real>  operator -(const Vector3<Real>  & v0, const Vector3<Real> & v1)
// return v0 - v1
{	return(Vector3<Real>(v0[0] - v1[0], v0[1] - v1[1], v0[2] - v1[2])); }
template<typename Real> inline Real operator *(Vector3<Real>& v0, Vector3<Real>& v1)
// return v0*v1T
{	return (v0[0] * v1[0] + v0[1] * v1[1] + v0[2] * v0[2]); }
template<typename Real> inline Vector3<Real> operator *(float x0, const Vector3<Real> & v0)
// return x0*v0
{	return(Vector3<Real>(x0*v0[0], x0*v0[1], x0*v0[2])); }
template<typename Real> inline Vector3<Real> operator *(const Vector3<Real> & v0, float x0)
// return v0*x0 (= x0*v0)
{	return(x0*v0); }

template<typename Real>
inline Vector4<Real> get_Vector4(Vector3<Real> vec)
{	Vector4<Real> temp(vec[0], vec[1], vec[2], 1.0); return temp;}

/// Vector 4 operators
template<typename Real> inline Real dot(const Vector4<Real> & v0, const Vector4<Real> & v1)
// return dot product of v0 and v1
{    return( v0[0]*v1[0] + v0[1]*v1[1] + v0[2]*v1[2] + v0[3]*v1[3]); }
//template<typename Real> inline Vector4<Real> cross(const Vector4<Real> & v0, const Vector4<Real> & v1)
template<typename Real> inline Vector4<Real>  operator +(const Vector4<Real>  & v0, const Vector4<Real>  & v1)
// return v0 + v1
{	return(Vector4<Real> (v0[0] + v1[0], v0[1] + v1[1], v0[2] + v1[2], v0[3] + v1[3])); }
template<typename Real> inline Vector4<Real>  operator -(const Vector4<Real>  & v0, const Vector4<Real> & v1)
// return v0 - v1
{	return(Vector4<Real>(v0[0] - v1[0], v0[1] - v1[1], v0[2] - v1[2], v0[3] - v1[3])); }
template<typename Real> inline Real operator *(Vector4<Real>& v0, Vector4<Real>& v1)
// return v0*v1T
{	return (v0[0] * v1[0] + v0[1] * v1[1] + v0[2] * v0[2] + v0[3] * v1[3]); }
template<typename Real> inline Vector4<Real> operator *(float x0, const Vector4<Real> & v0)
// return x0*v0
{	return(Vector4<Real>(x0*v0[0], x0*v0[1], x0*v0[2], x0*v0[3])); }
template<typename Real> inline Vector4<Real> operator *(const Vector4<Real> & v0, float x0)
// return v0*x0 (= x0*v0)
{	return(x0*v0); }

// backward compatibility
typedef Vector3<float>  VECTOR3;
typedef Vector4<float>  VECTOR4;


#if 0
//////////////////////////////////////////////////////////////////////////
// 3d matrix
//////////////////////////////////////////////////////////////////////////
class MATRIX3
{
private :
    Vector3 mat[3];       // a vector represents each matrix row

public :

	MATRIX3()                                    // constructor
	{	Identity(); };
    MATRIX3(const Vector3 & v0, const Vector3 & v1, const Vector3 & v2)
	{	mat[0] = v0; mat[1] = v1; mat[2] = v2; };  // constructor
	int Dimension() const
	{	return 3; };
    Vector3 & operator [](const int i)           // index row i
	{	return(mat[i]); };
	// Note: reference row i, column j of MATRIX3 m0 as m0[i][j] (not m0[i,j])
    Vector3 operator()(const int i) const        // return row i
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
    Vector4 mat[4];       // a vector represents each matrix row

public :

	MATRIX4()                                    // constructor
	{	Identity(); };
    MATRIX4(const Vector4 & v0, const Vector4 & v1,
        const Vector4 & v2, const Vector4 & v3) // constructor
	{	mat[0] = v0; mat[1] = v1; mat[2] = v2; mat[3] = v3; };
	int Dimension() const
	{	return 4; };
    Vector4 & operator [](int i)                 // index row i
	{	return(mat[i]); };
	// Note: reference row i, column j of MATRIX4 m0 as m0[i][j] (not m0[i,j])
    Vector4 operator()(const int i) const        // return row i
	{	return(mat[i]); };
	float operator ()(const int i, const int j) const
	{	return(mat[i](j)); };                    // return element (i,j)
	MATRIX4 & operator =(const MATRIX4 & m0)     // copy matrix m0
	{	mat[0] = m0(0); mat[1] = m0(1); mat[2] = m0(2); mat[3] = m0(3);
	return(*this); };
	MATRIX4 & operator =(const MATRIX3 & m0)     // copy matrix m0
	{	mat[0] = m0(0); mat[1] = m0(1); mat[2] = m0(2);
    Vector4 temp(0.0,0.0,0.0,1.0);
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
// Vector3 operations
//************************

inline float dot(const Vector3 & v0, const Vector3 & v1)
// return dot product of v0 and v1
{	return(v0(0)*v1(0) + v0(1)*v1(1) + v0(2)*v1(2)); };

inline Vector3 cross(const Vector3 & v0, const Vector3 & v1)
// return cross product of v0 and v1
{	return(Vector3(v0(1)*v1(2) - v0(2)*v1(1),
		   v0(2)*v1(0) - v0(0)*v1(2),
		   v0(0)*v1(1) - v0(1)*v1(0))); };

inline Vector3 operator +(const Vector3 & v0, const Vector3 & v1)
// return v0 + v1
{	return(Vector3(v0(0) + v1(0), v0(1) + v1(1), v0(2) + v1(2))); };

inline Vector3 operator -(const Vector3 & v0, const Vector3 & v1)
// return v0 - v1
{	return(Vector3(v0(0) - v1(0), v0(1) - v1(1), v0(2) - v1(2))); };

inline float operator *(Vector3& v0, Vector3& v1)
// return v0*v1T
{	return (v0(0) * v1(0) + v0(1) * v1(1) + v0(2) * v0(2)); };

inline Vector3 operator *(float x0, const Vector3 & v0)
// return x0*v0
{	return(Vector3(x0*v0(0), x0*v0(1), x0*v0(2))); };

inline Vector3 operator *(const Vector3 & v0, float x0)
// return v0*x0 (= x0*v0)
{	return(x0*v0); };

inline Vector4 get_Vector4(Vector3 vec)
{	Vector4 temp(vec[0], vec[1], vec[2], 1.0); return temp;};

//************************
// Vector4 operations
//************************

inline float dot(const Vector4 & v0, const Vector4 & v1)
// return dot product of v0 and v1
{	return(v0(0)*v1(0) + v0(1)*v1(1) + v0(2)*v1(2) + v0(3)*v1(3)); };

inline Vector4 operator +(const Vector4 & v0, const Vector4 & v1)
// return v0 + v1
{	return(Vector4(v0(0)+v1(0), v0(1)+v1(1), v0(2)+v1(2), v0(3)+v1(3))); };

inline Vector4 operator -(const Vector4 & v0, const Vector4 & v1)
// return v0 - v1
{	return(Vector4(v0(0)-v1(0), v0(1)-v1(1), v0(2)-v1(2), v0(3)-v1(3))); };

inline Vector4 operator *(float x0, const Vector4 & v0)
// return x0*v0
{	return(Vector4(x0*v0(0), x0*v0(1), x0*v0(2), x0*v0(3))); };

inline Vector4 operator *(const Vector4 & v0, float x0)
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
Vector3 operator *(const MATRIX3 & m0, const Vector3 & v0); // return m0 * v0
Vector3 operator *(const Vector3 & v0, const MATRIX3 & m0); // return v0 * m0

//************************
// MATRIX4 operations
//************************

MATRIX4 operator +(const MATRIX4 & m0, const MATRIX4 & m1); // return m0 + m1
MATRIX4 operator -(const MATRIX4 & m0, const MATRIX4 & m1); // return m0 - m1
MATRIX4 operator *(const MATRIX4 & m0, const MATRIX4 & m1); // return m0 * m1
MATRIX4 operator *(const float x0, const MATRIX4 & m0);    // return x0 * m0
MATRIX4 operator *(const MATRIX4 & m0, const float x0);    // return m0 * x0
Vector4 operator *(const MATRIX4 & m0, const Vector4 & v0); // return m0 * v0
Vector4 operator *(const Vector4 & v0, const MATRIX4 & m0); // return v0 * m0
Vector3 operator *(const MATRIX4 & m0, const Vector3 & v0); // return m0 * v0
Vector3 operator *(const Vector3 & v0, const MATRIX4 & m0); // return v0 * m0


MATRIX4 inverse(const MATRIX4 & m);  // return inverse of m; return 0 matrix if
// m is singular
MATRIX4 rotate_matrix(int type, float angle); // type: 1:x, 2:y, 3:z
MATRIX4 translate_matrix(float dx, float dy, float dz);
MATRIX4 scale_matrix(float sx, float sy, float sz);

#endif

} // namespace edda

#endif
