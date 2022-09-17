#include <glowstone.hpp>

#include <fstream>
#include <iostream>
#include <cmath>
#include <memory>
#include <vector>
#include <chrono>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

using namespace gs;

using std::make_shared;
using std::shared_ptr;

class hittable
{
public:
    virtual void hit(const ray &r, hit_record &rec, const material &outer) const = 0;
};

class obj_box final : public hittable
{
public:
    aabb geo;
    material mat;

public:
    obj_box(const aabb &a, const material &m) : geo(a), mat(m){};

    void hit(const ray &r, hit_record &rec, const material &outer) const
    {
        bool from_inside;

        if (geo.hit(r, rec.point, rec.normal, from_inside))
        {
            rec.success = true;

            if (from_inside)
            {
                rec.phase_1 = mat;
                rec.phase_2 = outer;
            }
            else
            {
                rec.phase_1 = outer;
                rec.phase_2 = mat;
            }
        }
        else
            rec = HIT_RECORD_MISS;
    }
};

class obj_sphere final : public hittable
{
public:
    sphere geo;
    material mat;

public:
    obj_sphere(const sphere &s, const material &m) : geo(s), mat(m){};

    void hit(const ray &r, hit_record &rec, const material &outer) const
    {
        bool from_inside;

        if (geo.hit(r, rec.point, rec.normal, from_inside))
        {
            rec.success = true;

            if (from_inside)
            {
                rec.phase_1 = mat;
                rec.phase_2 = outer;
            }
            else
            {
                rec.phase_1 = outer;
                rec.phase_2 = mat;
            }
        }
        else
            rec = HIT_RECORD_MISS;
    }
};

class simple_scene final : public scene
{
public:
    void hit(const ray &r, hit_record &rec) const
    {
        rec = HIT_RECORD_MISS;
        double t;

        for (const auto &obj : objects)
        {
            hit_record tmp;
            obj->hit(r, tmp, MATERIAL_VOID);

            if (tmp.success)
            {
                if (!rec.success)
                {
                    rec = tmp;
                    t = distance(tmp.point, r.origin);
                }
                else
                {
                    double t0 = distance(tmp.point, r.origin);
                    if (t0 < t)
                    {
                        t = t0;
                        rec = tmp;
                    }
                }
            }
        }
    }

    double background(const ray &r) const
    {
        return 0.1;
    }

    void add(shared_ptr<hittable> obj)
    {
        objects.push_back(obj);
    }

private:
    std::vector<shared_ptr<hittable>> objects;
};

void write_dat(const char *file_name, const double *src, int length)
{
    std::ofstream out(file_name, std::ios::binary);
    out.write((const char *)src, length * sizeof(double));
    out.close();
}

void read_dat(const char *file_name, double *dst, int length)
{
    std::ifstream in(file_name, std::ios::binary);
    in.read((char *)dst, length * sizeof(double));
    in.close();
}

#define WIDTH 800
#define HEIGHT 600
#define SIZE WIDTH *HEIGHT

double ray_color(const scene &s, const ray &r)
{
    static const int N = 10;

    double sum = 0;
    double max = -INFINITY;
    double min = INFINITY;

    for(int i = 0; i < N; ++i) {
        double result = s.trace(r);

        max = fmax(max, result);
        min = fmin(min, result);

        sum += result;
    }

    return (sum) / (N);
}

int main()
{
    simple_scene s;

    s.add(make_shared<obj_box>(aabb{{-0.5, -1.1, 7}, {0.5, -1, 7.5}}, material{true, 0.9, 0.7, 40, 0}));
    s.add(make_shared<obj_box>(aabb{{-0.5, -1.1, 4}, {0.5, -1, 4.5}}, material{true, 0.9, 0.7, 40, 0}));
    s.add(make_shared<obj_box>(aabb{{-0.5, -1.1, 3}, {0.5, -1, 3.5}}, material{true, 0.9, 0.7, 40, 0}));
    s.add(make_shared<obj_box>(aabb{{-0.5, -1.1, 2}, {0.5, -1, 2.5}}, material{true, 0.9, 0.7, 40, 0}));
    s.add(make_shared<obj_box>(aabb{{-2, -1.5, 0}, {2, -0.7, 9}}, material{true, 0.9, 0, 1.5, 0.9}));
    s.add(make_shared<obj_sphere>(sphere{{0, -0.19, 7}, 0.5}, material{true, 0.9, 0, 40, 0}));
    s.add(make_shared<obj_sphere>(sphere{{1, -0.19, 3}, 0.5}, material{false, 0.9, 0, 40, 0}));

    // s.add(make_shared<obj_box>(aabb{{-0.3, 0.8, 2}, {0.3, 0.81, 2.2}}, material{true, 0.01, 0.9, 40, 0}));
    // s.add(make_shared<obj_box>(aabb{{-1.5, -1.1, 0}, {1.5, -1, 4}}, material{true, 0.7, 0, 40, 0}));
    // s.add(make_shared<obj_box>(aabb{{-0.5, -1, 1.6}, {-0.1, -0.2, 2}}, material{true, 0.9, 0.1, 40, 0}));

    // s.add(make_shared<obj_rect>(rect{1, {-100, -100}, {100, 100}, -1}, material{0, false, const_0p9, const_0, 0}));
    // s.add(make_shared<obj_box>(aabb{{0.5, -0.5 + 1e-3, 3}, {1.5, 1.4, 4}}, material{0.9, false, const_0p9, const_0p3, const_1p5}));
    // s.add(make_shared<obj_sphere>(sphere{{0.5, 0.2, 6}, 1}, material{0, false, const_0p9, const_0, 0}));

    const double zNear = 0.01;
    const double U = zNear * WIDTH / HEIGHT;
    const double V = -zNear;

    std::cout << "Start." << std::endl;

    double *raw = new double[SIZE];

    auto begin = std::chrono::high_resolution_clock::now();
    for (int j = 0, pos = 0; j < HEIGHT; ++j)
    {
        for (int i = 0; i < WIDTH; ++i)
        {
            // std::cout << i << ' ' << j << std::endl << std::flush;
            double u = (double(i) + 0.5) / double(WIDTH) - 0.5;
            double v = (double(j) + 0.5) / double(HEIGHT) - 0.5;

            ray r{vec3{}, normalize(vec3{u * U, v * V, zNear})};

            raw[pos++] = ray_color(s, r);
        }
    }
    auto end = std::chrono::high_resolution_clock::now();

    char img_raw[SIZE];
    double max_raw = 0;
    for (int i = 0; i < SIZE; ++i)
    {
        img_raw[i] = 255 * fmin(sqrt(raw[i]), 1);

        max_raw = fmax(raw[i], max_raw);
    }

    /*double u = (double(423) + 0.5) / double(WIDTH) - 0.5;
    double v = (double(586) + 0.5) / double(HEIGHT) - 0.5;
    ray r{vec3{}, normalize(vec3{u * U, v * V, zNear})};

    double sum = 0;
    for (int i = 0; i < 1000000; ++i)
        sum += s.trace(r);
    double rv = sum / 1000000;

    std::ofstream out_("val_S_5000.csv");
    for (int j = 0; j < 10000; ++j)
    {
        sum = 0;
        for (int i = 0; i < 5000; ++i)
            sum += s.trace(r) - rv;
        out_ << sum << std::endl;
    }
    out_.close();*/

    stbi_write_png("image_raw.png", WIDTH, HEIGHT, 1, img_raw, WIDTH * 1);

    std::cout << "max_raw " << max_raw << std::endl;
    std::cout << "Done in " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "ms." << std::endl;
}
