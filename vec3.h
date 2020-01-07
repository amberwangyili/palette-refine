#ifndef PALETTE_REFINE_VEC3_H
#define PALETTE_REFINE_VEC3_H


#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <float.h>
#include <numeric>      // std::iota
#include <algorithm>    // std::sort
#include <opencv2/opencv.hpp>
#include <fstream>
#include <sstream>

class vec3
{
public:

    vec3(float x_, float y_, float z_) { v[0] = x_; v[1] = y_; v[2] = z_; }

    vec3(vec3 const& p)
    {
        v[0] = (float)(p[0]);
        v[1] = (float)(p[1]);
        v[2] = (float)(p[2]);
    }

    vec3() { v[0] = 0; v[1] = 0; v[2] = 0; }

    inline  float x() const { return v[0]; }
    inline  float y() const { return v[1]; }
    inline  float z() const { return v[2]; }

    inline  float operator [] (unsigned int c) const
    {
        return v[c];
    }
    inline  float& operator [] (unsigned int c)
    {
        return v[c];
    }

    static vec3 Zero() { return vec3(0, 0, 0); }

    void setZero()
    {
        v[0] = 0;
        v[1] = 0;
        v[2] = 0;
    }

    void operator += (const vec3& other)
    {
        v[0] += other.x();
        v[1] += other.y();
        v[2] += other.z();
    }
    void operator -= (const vec3& other)
    {
        v[0] -= other.x();
        v[1] -= other.y();
        v[2] -= other.z();
    }

    // floathis is going to create problems if the compiler needs to resolve umbiguous casts, but it's the cleaner way to do it
    void operator *= (int s)
    {
        v[0] *= s;
        v[1] *= s;
        v[2] *= s;
    }
    void operator *= (unsigned int s)
    {
        v[0] *= s;
        v[1] *= s;
        v[2] *= s;
    }
    void operator *= (float s)
    {
        v[0] *= s;
        v[1] *= s;
        v[2] *= s;
    }
    void operator *= (double s)
    {
        v[0] *= s;
        v[1] *= s;
        v[2] *= s;
    }
    void operator /= (int s)
    {
        v[0] /= s;
        v[1] /= s;
        v[2] /= s;
    }
    void operator /= (unsigned int s)
    {
        v[0] /= s;
        v[1] /= s;
        v[2] /= s;
    }
    void operator /= (float s)
    {
        v[0] /= s;
        v[1] /= s;
        v[2] /= s;
    }
    void operator /= (double s)
    {
        v[0] /= s;
        v[1] /= s;
        v[2] /= s;
    }

    vec3 getOrthogonal() const
    {
        if (v[0] == 0)
        {
            return vec3(0, v[2], -v[1]);
        }
        else if (v[1] == 0)
        {
            return vec3(v[2], 0, -v[0]);
        }

        return vec3(v[1], -v[0], 0);
    }


    float sqrnorm() const
    {
        return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
    }
    float norm() const
    {
        return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    }

    void normalize()
    {
        float _n = norm();
        v[0] /= _n;
        v[1] /= _n;
        v[2] /= _n;
    }

    //float dot(const vec3& anotherVec) const {
    //	return ::dot(*this,anotherVec);
    //}
    const static vec3 null;
private:
    float v[3];
};




inline
float dot(const vec3 & p1, const vec3 & p2)
{
    return p1.x()* p2.x() + p1.y() * p2.y() + p1.z() * p2.z();
}
inline
vec3 cross(const vec3 & p1, const vec3 & p2) {
    return vec3(
            p1.y() * p2.z() - p1.z() * p2.y(),
            p1.z() * p2.x() - p1.x() * p2.z(),
            p1.x() * p2.y() - p1.y() * p2.x()
    );
}

inline
vec3 operator + (const vec3 & p1, const vec3 & p2)
{
    return vec3(p1.x() + p2.x(), p1.y() + p2.y(), p1.z() + p2.z());
}
inline
vec3 operator - (const vec3 & p1, const vec3 & p2)
{
    return vec3(p1.x() - p2.x(), p1.y() - p2.y(), p1.z() - p2.z());
}


inline
vec3 operator - (const vec3 & p2)
{
    return vec3(-p2.x(), -p2.y(), -p2.z());
}

inline
vec3 operator * (const vec3 & p, int s)
{
    return vec3(s * p.x(), s * p.y(), s * p.z());
}
inline
vec3 operator * (const vec3 & p, float s)
{
    return vec3(s * p.x(), s * p.y(), s * p.z());
}
inline
vec3 operator * (const vec3 & p, double s)
{
    return vec3(s * p.x(), s * p.y(), s * p.z());
}
inline
vec3 operator * (int s, const vec3 & p)
{
    return vec3(s * p.x(), s * p.y(), s * p.z());
}
inline
vec3 operator * (float s, const vec3 & p)
{
    return vec3(s * p.x(), s * p.y(), s * p.z());
}
inline
vec3 operator * (double s, const vec3 & p)
{
    return vec3(s * p.x(), s * p.y(), s * p.z());
}


inline
vec3 operator / (const vec3 & p, int s)
{
    return vec3(p.x() / s, p.y() / s, p.z() / s);
}
inline
vec3 operator / (const vec3 & p, float s)
{
    return vec3(p.x() / s, p.y() / s, p.z() / s);
}
inline
vec3 operator / (const vec3 & p, double s)
{
    return vec3(p.x() / s, p.y() / s, p.z() / s);
}


inline
float operator * (const vec3 & p1, const vec3 & p2)
{
    return p1.x()* p2.x() + p1.y() * p2.y() + p1.z() * p2.z();
}


inline
vec3 operator % (const vec3 & p1, const vec3 & p2)
{
    return cross(p1, p2);
}




inline std::ostream& operator << (std::ostream & s, vec3 const& p)
{
    s << p[0] << " \t" << p[1] << " \t" << p[2];
    return s;
}

using namespace std;
using namespace cv;

class Triangle;




class Ray {
public:
    vec3 origin;
    vec3 direction;
    Ray() { }
    Ray(vec3 s, vec3 d) {
        this->origin = s;
        this->direction = d;
        this->direction.normalize();
    }
    inline vec3 pointAtParameter(double t) const {
        return origin + direction * t;
    }
};

class Triangle {
public:
    vec3 a, b, c;
    double d;
    vec3 norm;
    Triangle() {}
    Triangle(vec3 a_, vec3 b_, vec3 c_) :a(a_), b(b_), c(c_) {
        norm = cross(b - a, c - a);
        norm.normalize();
        vec3 nnorm = -norm;
        vec3 center = (a + b + c) / 3;
        d = -dot(norm, a);
    }

public:
    inline static bool checkside(const vec3 & a, const vec3 & b, const vec3 & c, const vec3 & p) {
        vec3 ab = b - a, ac = c - a, ap = p - a;
        vec3 v1 = cross(ab, ac);
        vec3 v2 = cross(ab, ap);
        return dot(v1, v2) >= 0;
    }

    inline static bool in_tri(const vec3 & p, const Triangle & t) {
        return checkside(t.a, t.b, t.c, p) && checkside(t.b, t.c, t.a, p) && checkside(t.c, t.a, t.b, p);
    }

    //bool intersectionWithRayZ(const vec3& ray_origin) const {
    //	float eps = 1e-8;
    //	double t = -(d + dot(norm, ray_origin)) / (norm.z() + eps);
    //	if (t > eps) {
    //		vec3 point = ray_origin;
    //		point[2] += t;
    //		bool inter = in_tri(point, *this);
    //		return true;
    //	}
    //	return false;
    //}

    float intersectionParameter(const Ray& ray) const {
        float eps = 1e-8;

        float dot_product = dot(norm, ray.direction);
        if (fabs(dot_product) < eps)
            return NULL;

        double t = -(d + dot(norm, ray.origin)) / (dot_product + eps);
        if (t > eps) {
            vec3 point = ray.pointAtParameter(t);
            bool inter = in_tri(point, *this);
            if (inter) {
                return t;
            }
        }
        return NULL;
    }

};


struct Mesh {
public:
    vector<vec3> vertices;
    vector<vector<int>> faces;
public:
    void clear() { vertices.clear(); faces.clear();}
    int vertex_num() const { return int(vertices.size()); }
    int face_num() const { return int(faces.size()); }

    bool load_from_file(const string& filename);
    bool save_to_file(const string& filename);

    static bool isInside(const vector<Triangle>& triangles, const vec3& point);
    //bool isInsideFast(const vec3& point) const;
    void constructTriangles(vector<Triangle>& triangles) const;
};

bool Mesh::load_from_file(const string& filename) {
    this->clear();
    char s1;
    float f2, f3, f4;


    ifstream infile(filename.c_str());
    if (!infile.is_open()) {
        printf("fail loading");

        return false;
    }

    string sline;

    while (getline(infile, sline)) {
        istringstream sin(sline);
        sin >> s1 >> f2 >> f3 >> f4;
        if (s1 == 'v') {
            f2 = min(max(f2, 0.0f), 1.0f);
            f3 = min(max(f3, 0.0f), 1.0f);
            f4 = min(max(f4, 0.0f), 1.0f);
            vertices.emplace_back(vec3(f2, f3, f4));
        }
        else if (s1 == 'f') {
            vector<int> tri(3);
            tri[0] = int(f2);
            tri[1] = int(f3);
            tri[2] = int(f4);
            faces.emplace_back(tri);
        }
    }


    return true;
}


void Mesh::constructTriangles(vector<Triangle> & triangles) const {
    triangles.clear();
    for (int i = 0; i < (int)face_num(); ++i) {
        vec3 a = vertices[faces[i][0]];
        vec3 b = vertices[faces[i][1]];
        vec3 c = vertices[faces[i][2]];

        Triangle t = Triangle(a, b, c);
        triangles.push_back(t);
    }
}

bool Mesh::save_to_file(const string& filename) {
    ofstream ofs(filename.c_str());
    if (!ofs.is_open()) {
        return false;
    }
    for (int i = 0; i < (int)vertices.size(); ++i) {
        ofs << "v " << vertices[i].x() << " " << vertices[i].y() << " " << vertices[i].z() << "\n";
    }
    for (int j = 0; j < (int)faces.size(); ++j) {
        ofs << "f " << faces[j][0] << " " << faces[j][1] << " " << faces[j][2] << "\n";
    }
    ofs.close();
    return true;
}

bool Mesh::isInside(const vector<Triangle> & triangles,const vec3& pt) {

    Ray r = Ray(pt, vec3(1, 0, 1));

    int count = 0;

    for (int i = 0; i < (int)triangles.size(); ++i) {
        const Triangle& t = triangles[i];
        if (t.intersectionParameter(r) > 0){
            count += 1;
        }
    }
    if (count % 2 == 0) {
        return false;
    }
    else
        return true;
}



//bool Mesh::isInsideFast(const vec3& pt) const {
//	int count = 0;
//
//	for (int i = 0; i < (int)_triangles.size(); ++i) {
//		const Triangle& t = _triangles[i];
//		if (t.intersectionWithRayZ(pt)) {
//			count += 1;
//		}
//	}
//	if (count % 2 == 0) {
//		return false;
//	}
//	else
//		return true;
//}




#endif //PALETTE_REFINE_VEC3_H
