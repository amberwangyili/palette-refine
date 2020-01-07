
#ifndef PALETTE_REFINE_NEARESTPOINT_H
#define PALETTE_REFINE_NEARESTPOINT_H


#include "vec3.h"
using namespace std;

template <typename T> T clamp(const T& val, const T& val_min, const T& val_max)
{
    T tmp = std::min(val, val_max);
    tmp = std::max(val, val_min);
    return tmp;

}


inline vec3 closesPointOnTriangle(const vec3* triangle, const vec3& sourcePosition)
{
    vec3 edge0 = triangle[1] - triangle[0];
    vec3 edge1 = triangle[2] - triangle[0];
    vec3 v0 = triangle[0] - sourcePosition;

    float a = dot(edge0, edge0);
    float b = dot(edge0, edge1);
    float c = dot(edge1, edge1);
    float d = dot(edge0, v0);
    float e = dot(edge1, v0);

    float det = a * c - b * b;
    float s = b * e - c * d;
    float t = b * d - a * e;

    if (s + t < det)
    {
        if (s < 0.f)
        {
            if (t < 0.f)
            {
                if (d < 0.f)
                {
                    s = clamp(-d / a, 0.f, 1.f);
                    t = 0.f;
                }
                else
                {
                    s = 0.f;
                    t = clamp(-e / c, 0.f, 1.f);
                }
            }
            else
            {
                s = 0.f;
                t = clamp(-e / c, 0.f, 1.f);
            }
        }
        else if (t < 0.f)
        {
            s = clamp(-d / a, 0.f, 1.f);
            t = 0.f;
        }
        else
        {
            float invDet = 1.f / det;
            s *= invDet;
            t *= invDet;
        }
    }
    else
    {
        if (s < 0.f)
        {
            float tmp0 = b + d;
            float tmp1 = c + e;
            if (tmp1 > tmp0)
            {
                float numer = tmp1 - tmp0;
                float denom = a - 2 * b + c;
                s = clamp(numer / denom, 0.f, 1.f);
                t = 1 - s;
            }
            else
            {
                t = clamp(-e / c, 0.f, 1.f);
                s = 0.f;
            }
        }
        else if (t < 0.f)
        {
            if (a + d > b + e)
            {
                float numer = c + e - b - d;
                float denom = a - 2 * b + c;
                s = clamp(numer / denom, 0.f, 1.f);
                t = 1 - s;
            }
            else
            {
                s = clamp(-e / c, 0.f, 1.f);
                t = 0.f;
            }
        }
        else
        {
            float numer = c + e - b - d;
            float denom = a - 2 * b + c;
            s = clamp(numer / denom, 0.f, 1.f);
            t = 1.f - s;
        }
    }

    return triangle[0] + s * edge0 + t * edge1;
}
inline vec3 nearest_point(vec3 const& x, vector<vector<int> > const& triangles, vector<vec3> const& vertices)
{
    float min_distance = FLT_MAX;
    vec3 close_point;

    for (int i = 0; i < (int)triangles.size(); i++)
    {
        vec3 triangle_vertices[3];
        triangle_vertices[0] = vertices[triangles[i][0]];
        triangle_vertices[1] = vertices[triangles[i][1]];
        triangle_vertices[2] = vertices[triangles[i][2]];

        vec3 my_close_point = closesPointOnTriangle(triangle_vertices, x);
        float my_dis = (x - my_close_point).sqrnorm();
        if (my_dis < min_distance) {
            min_distance = my_dis;
            close_point = my_close_point;
        }
    }

    return close_point;
}


#endif //PALETTE_REFINE_NEARESTPOINT_H
