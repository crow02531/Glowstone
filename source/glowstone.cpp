#include <glowstone.hpp>

#include <cmath>

using namespace gs;
using namespace std;

double scene::background(const ray &r) const
{
    return 0;
}

inline double fresnel_schlick(double eta_1, double eta_2, double cos_theta_i)
{
    double r0 = (eta_1 - eta_2) / (eta_1 + eta_2);
    r0 *= r0;

    return r0 + (1.f - r0) * pow(1.f - cos_theta_i, 5.f);
}

inline bool is_degenerated(const vec3 &v)
{
    return isnan(v.x) || isnan(v.y) || isnan(v.z);
}

double scene::trace(const ray &r) const
{
    ray current_ray = r;
    double path_factor = 1;
    double correct_factor = 1;
    double result = 0;

    while (true)
    {
        // get p_continue based on current path properties
        double p_continue = fmin(8 * path_factor, fmin(0.8 + 0.6 * path_factor, 1));

        if (random() < p_continue)
        {
            correct_factor /= p_continue;

            hit_record rec;
            hit(current_ray, rec);

            if (rec.success)
            {
                path_factor *= pow(rec.phase_1.transparency, distance(current_ray.origin, rec.point));
                result += path_factor * correct_factor * rec.phase_2.radiation;
                path_factor *= rec.phase_2.albedo;

                vec3 dir;

                if (rec.phase_2.specular)
                {
                    dir = refract(current_ray.direction, rec.normal, rec.phase_1.eta / rec.phase_2.eta);
                    double fs_r = fresnel_schlick(rec.phase_1.eta, rec.phase_2.eta, -dot(current_ray.direction, rec.normal));
                    double fs_t = sign(rec.phase_2.transparency) * dir.length_squared() * (1 - fs_r);

                    // get p_reflect based on current hit properties
                    double p_reflect = fs_r / (fs_r + fs_t);

                    if (random() < p_reflect)
                    {
                        dir = reflect(current_ray.direction, rec.normal);
                        path_factor *= fs_r;
                        correct_factor /= p_reflect;
                    }
                    else
                    {
                        path_factor *= fs_t;
                        correct_factor /= 1 - p_reflect;
                    }
                }
                else
                {
                    dir = normalize(rec.normal + random_on_sphere());
                    if (is_degenerated(dir))
                        dir = rec.normal;
                    path_factor *= dot(dir, rec.normal);
                    correct_factor /= (dot(dir, rec.normal) / M_PI) * 2 * M_PI;
                }

                current_ray = ray{rec.point + dir * 1e-5, dir};
            }
            else
            {
                result += path_factor * correct_factor * background(current_ray);

                break;
            }
        }
        else
            break;
    }

    return result;
}
