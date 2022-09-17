#ifndef _glowstone_h_
#define _glowstone_h_

#include <more_math.hpp>

namespace gs
{
    struct material final
    {
        // Either specular or scattering.
        bool specular;

        // Any value in (0, 1].
        double albedo;

        // Any value larger or equals to 0.
        double radiation;

        // Discribe the refractive index of the material.
        // Larger or equals to 1.
        double eta;

        // Any value in [0, 1].
        double transparency;
    } const MATERIAL_VOID = material{true, 1, 0, 1, 1};

    bool operator==(const material &a, const material &b);
    bool operator!=(const material &a, const material &b);

    struct hit_record final
    {
        // Hit or not.
        bool success;

        // The hit point, undefined when missed.
        vec3 point;

        // The normal of the hit surface, a unit vector point from
        // phase 2 to phase 1. Undefined when missed.
        vec3 normal;

        // The phase light comes from. Undefined when missed.
        material phase_1;

        // The opposite phase. Undefined when missed. Two phases
        // can't be the same.
        material phase_2;
    } const HIT_RECORD_MISS{};

    bool operator==(const hit_record &a, const hit_record &b);
    bool operator!=(const hit_record &a, const hit_record &b);

    class scene
    {
    public:
        // Dispatch a random experiment with the expectation of its result
        // equals to the spectral irradiance(in the unit W/(m*m*nm)) in the
        // space.
        double trace(const ray &r) const;

        virtual void hit(const ray &r, hit_record &rec) const = 0;
        virtual double background(const ray &r) const;
    };
}

#endif /* _glowstone_h_ */
