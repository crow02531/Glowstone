#include <more_math.hpp>

#include <random>
#include <time.h>
#include <bits/move.h>

using namespace gs;
using namespace std;

double vec3::length_squared() const
{
    return x * x + y * y + z * z;
}

double vec3::length() const
{
    return std::sqrt(length_squared());
}

double vec3::operator[](int i) const
{
    switch (i)
    {
    case 0:
        return x;
    case 1:
        return y;
    case 2:
        return z;
    default:
        throw 0;
    }
}

double &vec3::operator[](int i)
{
    switch (i)
    {
    case 0:
        return x;
    case 1:
        return y;
    case 2:
        return z;
    default:
        throw 0;
    }
}

bool gs::operator==(const vec3 &u, const vec3 &v)
{
    return u.x == v.x && u.y == v.y && u.z == v.z;
}

bool gs::operator!=(const vec3 &u, const vec3 &v)
{
    return !(u == v);
}

vec3 gs::operator-(const vec3 &u)
{
    return vec3{-u.x, -u.y, -u.z};
}

vec3 gs::operator+(const vec3 &u, const vec3 &v)
{
    return vec3{u.x + v.x, u.y + v.y, u.z + v.z};
}

vec3 gs::operator-(const vec3 &u, const vec3 &v)
{
    return vec3{u.x - v.x, u.y - v.y, u.z - v.z};
}

vec3 gs::operator*(const vec3 &u, const vec3 &v)
{
    return vec3{u.x * v.x, u.y * v.y, u.z * v.z};
}

vec3 gs::operator*(const vec3 &u, double t)
{
    return vec3{u.x * t, u.y * t, u.z * t};
}

vec3 gs::operator*(double t, const vec3 &u)
{
    return u * t;
}

vec3 gs::operator/(const vec3 &u, double t)
{
    return u * (1 / t);
}

vec3 &gs::operator+=(vec3 &u, const vec3 &v)
{
    u.x += v.x;
    u.y += v.y;
    u.z += v.z;

    return u;
}

vec3 &gs::operator-=(vec3 &u, const vec3 &v)
{
    u.x -= v.x;
    u.y -= v.y;
    u.z -= v.z;

    return u;
}

vec3 &gs::operator*=(vec3 &u, double t)
{
    u.x *= t;
    u.y *= t;
    u.z *= t;

    return u;
}

vec3 &gs::operator/=(vec3 &u, double t)
{
    return u *= (1 / t);
}

vec3 gs::normalize(const vec3 &u)
{
    return u / u.length();
}

double gs::distance(const vec3 &u, const vec3 &v)
{
    return (u - v).length();
}

double gs::dot(const vec3 &u, const vec3 &v)
{
    return u.x * v.x + u.y * v.y + u.z * v.z;
}

vec3 gs::cross(const vec3 &u, const vec3 &v)
{
    return vec3{u.y * v.z - u.z * v.y,
                u.z * v.x - u.x * v.z,
                u.x * v.y - u.y * v.x};
}

vec3 gs::reflect(const vec3 &I, const vec3 &N)
{
    return I - 2 * dot(N, I) * N;
}

vec3 gs::refract(const vec3 &I, const vec3 &N, double eta)
{
    double N_dot_I = dot(N, I);
    double k = 1 - eta * eta * (1 - N_dot_I * N_dot_I);

    if (k < 0)
        return vec3{};
    else
        return eta * I - (eta * N_dot_I + sqrt(k)) * N;
}

double gs::random()
{
    static uniform_real_distribution<double> dist(0, 1);
    static mt19937 gen(89000099);

    return dist(gen);
}

vec3 gs::random_on_sphere()
{
    double theta = random() * 2 * M_PI;
    double z = random() * 2 - 1;
    double s = sqrt(1 - z * z);

    return vec3{s * cos(theta), s * sin(theta), z};
}

vec3 gs::random_on_hemisphere(const vec3 &N)
{
    vec3 r = random_on_sphere();

    return dot(r, N) >= 0 ? r : -r;
}

double gs::sign(double x)
{
    return x > 0 ? 1 : (x == 0 ? 0 : -1);
}

bool aabb::hit(const ray &r) const
{
    auto dirfracX = 1 / r.direction.x;
    auto dirfracY = 1 / r.direction.y;
    auto dirfracZ = 1 / r.direction.z;

    auto t1 = (min.x - r.origin.x) * dirfracX;
    auto t2 = (max.x - r.origin.x) * dirfracX;
    auto t3 = (min.y - r.origin.y) * dirfracY;
    auto t4 = (max.y - r.origin.y) * dirfracY;
    auto t5 = (min.z - r.origin.z) * dirfracZ;
    auto t6 = (max.z - r.origin.z) * dirfracZ;

    auto tmin = fmax(fmax(fmin(t1, t2), fmin(t3, t4)), fmin(t5, t6));
    auto tmax = fmin(fmin(fmax(t1, t2), fmax(t3, t4)), fmax(t5, t6));

    return tmin < tmax && tmax > 0 && tmin != 0;
}

bool aabb::hit(const ray &r, vec3 &p, vec3 &n, bool &in) const
{
    auto dirfracX = 1 / r.direction.x;
    auto dirfracY = 1 / r.direction.y;
    auto dirfracZ = 1 / r.direction.z;

    auto t1 = (min.x - r.origin.x) * dirfracX;
    auto t2 = (max.x - r.origin.x) * dirfracX;
    auto t3 = (min.y - r.origin.y) * dirfracY;
    auto t4 = (max.y - r.origin.y) * dirfracY;
    auto t5 = (min.z - r.origin.z) * dirfracZ;
    auto t6 = (max.z - r.origin.z) * dirfracZ;

    auto tmin = fmax(fmax(fmin(t1, t2), fmin(t3, t4)), fmin(t5, t6));
    auto tmax = fmin(fmin(fmax(t1, t2), fmax(t3, t4)), fmax(t5, t6));

    if (tmin >= tmax || tmax <= 0 || tmin == 0)
        return false;

    if (in = tmin < 0)
    {
        p = r.origin + r.direction * tmax;

        if (t1 == tmax)
            n = vec3{1, 0, 0};
        else if (t2 == tmax)
            n = vec3{-1, 0, 0};
        else if (t3 == tmax)
            n = vec3{0, 1, 0};
        else if (t4 == tmax)
            n = vec3{0, -1, 0};
        else if (t5 == tmax)
            n = vec3{0, 0, 1};
        else if (t6 == tmax)
            n = vec3{0, 0, -1};
    }
    else
    {
        p = r.origin + r.direction * tmin;

        if (t1 == tmin)
            n = vec3{-1, 0, 0};
        else if (t2 == tmin)
            n = vec3{1, 0, 0};
        else if (t3 == tmin)
            n = vec3{0, -1, 0};
        else if (t4 == tmin)
            n = vec3{0, 1, 0};
        else if (t5 == tmin)
            n = vec3{0, 0, -1};
        else if (t6 == tmin)
            n = vec3{0, 0, 1};
    }

    return true;
}

aabb surrounding_box(const aabb &u, const aabb &v)
{
    vec3 min{fmin(u.min.x, v.min.x),
             fmin(u.min.y, v.min.y),
             fmin(u.min.z, v.min.z)};

    vec3 max{fmax(u.max.x, v.max.x),
             fmax(u.max.y, v.max.y),
             fmax(u.max.z, v.max.z)};

    return aabb{min, max};
}

bool sphere::hit(const ray &r, vec3 &p, vec3 &n, bool &in) const
{
    vec3 oc = r.origin - center;
    auto half_b = dot(oc, r.direction);
    auto c = oc.length_squared() - radius * radius;

    auto discriminant = half_b * half_b - c;
    if (discriminant <= 0)
        return false;
    auto sqrtd = sqrt(discriminant);

    auto root = (-half_b - sqrtd);
    if (root < 0)
    {
        root = (-half_b + sqrtd);
        if (root < 0)
            return false;
    }

    p = r.origin + r.direction * root;
    n = (p - center) / radius;

    // c = dot(r.direction, n);
    // in = c < 0 ? false : true;
    // n *= -sign(c);
    in = false;

    return true;
}
