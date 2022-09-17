#ifndef _more_math_h_
#define _more_math_h_

namespace gs
{
    struct vec2 final
    {
        double x;
        double y;

        double length_squared() const;
        double length() const;

        double &operator[](int i);
        double operator[](int i) const;
    };

    bool operator==(const vec2 &u, const vec2 &v);
    bool operator!=(const vec2 &u, const vec2 &v);

    vec2 operator-(const vec2 &u);
    vec2 operator+(const vec2 &u, const vec2 &v);
    vec2 operator-(const vec2 &u, const vec2 &v);
    vec2 operator*(const vec2 &u, const vec2 &v);
    vec2 operator*(const vec2 &u, double t);
    vec2 operator*(double t, const vec2 &u);
    vec2 operator/(const vec2 &u, double t);

    vec2 &operator+=(vec2 &u, const vec2 &v);
    vec2 &operator-=(vec2 &u, const vec2 &v);
    vec2 &operator*=(vec2 &u, double t);
    vec2 &operator/=(vec2 &u, double t);

    vec2 normalize(const vec2 &u);
    double distance(const vec2 &u, const vec2 &v);
    double dot(const vec2 &u, const vec2 &v);

    struct vec3 final
    {
        double x;
        double y;
        double z;

        double length_squared() const;
        double length() const;

        double &operator[](int i);
        double operator[](int i) const;
    };

    bool operator==(const vec3 &u, const vec3 &v);
    bool operator!=(const vec3 &u, const vec3 &v);

    vec3 operator-(const vec3 &u);
    vec3 operator+(const vec3 &u, const vec3 &v);
    vec3 operator-(const vec3 &u, const vec3 &v);
    vec3 operator*(const vec3 &u, const vec3 &v);
    vec3 operator*(const vec3 &u, double t);
    vec3 operator*(double t, const vec3 &u);
    vec3 operator/(const vec3 &u, double t);

    vec3 &operator+=(vec3 &u, const vec3 &v);
    vec3 &operator-=(vec3 &u, const vec3 &v);
    vec3 &operator*=(vec3 &u, double t);
    vec3 &operator/=(vec3 &u, double t);

    vec3 normalize(const vec3 &u);
    double distance(const vec3 &u, const vec3 &v);
    double dot(const vec3 &u, const vec3 &v);
    vec3 cross(const vec3 &u, const vec3 &v);
    vec3 reflect(const vec3 &I, const vec3 &N);
    vec3 refract(const vec3 &I, const vec3 &N, double eta);

    struct ray final
    {
        // The origin of the ray, can be any vector.
        vec3 origin;

        // The direction of the ray, must be a unit vector.
        vec3 direction;
    };

    bool operator==(const ray &a, const ray &b);
    bool operator!=(const ray &a, const ray &b);

    // Uniform in [0, 1).
    double random();
    vec3 random_on_sphere();
    // Include equator.
    vec3 random_on_hemisphere(const vec3 &N);

    double sign(double x);

    struct aabb final
    {
        vec3 min;
        vec3 max;

        bool hit(const ray &r) const;
        bool hit(const ray &r, vec3 &p, vec3 &n, bool &in) const;
    };

    bool operator==(const aabb &u, const aabb &v);
    bool operator!=(const aabb &u, const aabb &v);

    aabb surrounding_box(const aabb &u, const aabb &v);

    struct sphere final
    {
        vec3 center;
        double radius;

        bool hit(const ray &r) const;
        bool hit(const ray &r, vec3 &p, vec3 &n, bool &in) const;
    };

    bool operator==(const sphere &u, const sphere &v);
    bool operator!=(const sphere &u, const sphere &v);
}

#endif /* _more_math_h_ */
