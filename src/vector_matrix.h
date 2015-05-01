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
#include "core/tuple.h"

namespace edda{

//////////////////////////////////////////////////////////////////////////
//	vector with any number of components
//////////////////////////////////////////////////////////////////////////
template <typename Real, int N>
class VECTOR : public Tuple<Real, N>
{
    using Tuple<Real, N>::vec;
public :
    VECTOR() : Tuple<Real, N>(0) {}

    // make sure all component<=1.0
    void clamp()
    {
        for (int i = 0; i < this->getLen(); i++)
            if (vec[i]>1.0) vec[i] = 1.0;
    }
};

//////////////////////////////////////////////////////////////////////////
/// for perfomrnace we implement each VECTOR with different lengths
//	vector with 3 components
//////////////////////////////////////////////////////////////////////////
template <typename Real>
class VECTOR3 : public VECTOR<Real, 3>
{
    using Tuple<Real, 3>::vec;
public :
    VECTOR3(Real x, Real y, Real z) { vec[0]=x; vec[1]=y; vec[2]=z; }
	// constructor
    inline Real x() const {return vec[0];}
    inline Real y() const {return vec[1];}
    inline Real z() const {return vec[2];}

    bool operator ==(const VECTOR3& v) const {
        return (fabs(vec[0]-v(0)) < limits<Real>::epsilon() &&
                fabs(vec[1]-v(1)) < limits<Real>::epsilon() &&
                fabs(vec[2]-v(2)) < limits<Real>::epsilon());
    }
    // get magnitude
    double getMag() const { return sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);  }
    // get the maximum value
    Real getMax() const { return std::max(vec[0], std::max(vec[1], vec[2])); }
    // normalize vector
    void Normalize()    { double norm = getMag(); if (norm!=0) (*this) *= (1/norm); }
    VECTOR3<Real> operator+=(double x) { vec[0]+=x; vec[1]+=x; vec[2]+=x; return *this; }
    VECTOR3<Real> operator*=(double x) { vec[0]*=x; vec[1]*=x; vec[2]*=x; return *this; }
};


//////////////////////////////////////////////////////////////////////////
/// for perfomrnace we implement each VECTOR with different lengths
//	vector with 4 components
//////////////////////////////////////////////////////////////////////////
template <typename Real>
class VECTOR4 : public VECTOR<Real, 4>
{
    using Tuple<Real, 4>::vec;
public :
    VECTOR4(Real x, Real y, Real z, Real w) { vec[0]=x; vec[1]=y; vec[2]=z; vec[3]=w; }
    // constructor
    inline Real x() const {return vec[0];}
    inline Real y() const {return vec[1];}
    inline Real z() const {return vec[2];}
    inline Real w() const {return vec[3];}

    bool operator ==(const VECTOR4<Real>& v) const {
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
    void Normalize()    { double norm = getMag(); if (norm!=0) (*this) *= (1/norm); }
    VECTOR4<Real> operator+=(double x) { vec[0]+=x; vec[1]+=x; vec[2]+=x; vec[3]+=x; return *this; }
    VECTOR4<Real> operator*=(double x) { vec[0]*=x; vec[1]*=x; vec[2]*=x; vec[3]+=x; return *this; }
};
//************************
// VECTOR operations
//************************
template<typename Real> inline Real dot(const VECTOR3<Real> & v0, const VECTOR3<Real> & v1)
// return dot product of v0 and v1
{    return( v0[0]*v1[0] + v0[1]*v1[1] + v0[2]*v1[2]); }
template<typename Real> inline VECTOR3<Real> cross(const VECTOR3<Real> & v0, const VECTOR3<Real> & v1)
// return cross product of v0 and v1
{	return VECTOR3<Real>( v0[1]*v1[2] - v0[2]*v1[1],
                          v0[2]*v1[0] - v0[0]*v1[2],
                          v0[0]*v1[1] - v0[1]*v1[0] );
}
template<typename Real> inline VECTOR3<Real>  operator +(const VECTOR3<Real>  & v0, const VECTOR3<Real>  & v1)
// return v0 + v1
{	return(VECTOR3<Real> (v0[0] + v1[0], v0[1] + v1[1], v0[2] + v1[2])); }
template<typename Real> inline VECTOR3<Real>  operator -(const VECTOR3<Real>  & v0, const VECTOR3<Real> & v1)
// return v0 - v1
{	return(VECTOR3<Real>(v0[0] - v1[0], v0[1] - v1[1], v0[2] - v1[2])); }
template<typename Real> inline Real operator *(VECTOR3<Real>& v0, VECTOR3<Real>& v1)
// return v0*v1T
{	return (v0[0] * v1[0] + v0[1] * v1[1] + v0[2] * v0[2]); }
template<typename Real> inline VECTOR3<Real> operator *(float x0, const VECTOR3<Real> & v0)
// return x0*v0
{	return(VECTOR3<Real>(x0*v0[0], x0*v0[1], x0*v0[2])); }
template<typename Real> inline VECTOR3<Real> operator *(const VECTOR3<Real> & v0, float x0)
// return v0*x0 (= x0*v0)
{	return(x0*v0); }

template<typename Real>
inline VECTOR4<Real> get_vector4(VECTOR3<Real> vec)
{	VECTOR4<Real> temp(vec[0], vec[1], vec[2], 1.0); return temp;}

/// Vector 4 operators
template<typename Real> inline Real dot(const VECTOR4<Real> & v0, const VECTOR4<Real> & v1)
// return dot product of v0 and v1
{    return( v0[0]*v1[0] + v0[1]*v1[1] + v0[2]*v1[2] + v0[3]*v1[3]); }
//template<typename Real> inline VECTOR4<Real> cross(const VECTOR4<Real> & v0, const VECTOR4<Real> & v1)
template<typename Real> inline VECTOR4<Real>  operator +(const VECTOR4<Real>  & v0, const VECTOR4<Real>  & v1)
// return v0 + v1
{	return(VECTOR4<Real> (v0[0] + v1[0], v0[1] + v1[1], v0[2] + v1[2], v0[3] + v1[3])); }
template<typename Real> inline VECTOR4<Real>  operator -(const VECTOR4<Real>  & v0, const VECTOR4<Real> & v1)
// return v0 - v1
{	return(VECTOR4<Real>(v0[0] - v1[0], v0[1] - v1[1], v0[2] - v1[2], v0[3] - v1[3])); }
template<typename Real> inline Real operator *(VECTOR4<Real>& v0, VECTOR4<Real>& v1)
// return v0*v1T
{	return (v0[0] * v1[0] + v0[1] * v1[1] + v0[2] * v0[2] + v0[3] * v1[3]); }
template<typename Real> inline VECTOR4<Real> operator *(float x0, const VECTOR4<Real> & v0)
// return x0*v0
{	return(VECTOR4<Real>(x0*v0[0], x0*v0[1], x0*v0[2], x0*v0[3])); }
template<typename Real> inline VECTOR4<Real> operator *(const VECTOR4<Real> & v0, float x0)
// return v0*x0 (= x0*v0)
{	return(x0*v0); }

#if 0
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
#endif
