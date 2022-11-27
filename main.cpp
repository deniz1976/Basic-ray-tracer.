/**
 * a basic ray tracer in C++17
 * Deniz Mutlu
 * Thanks for Marcus Mathiassen for the idea and source code of this algorithm
 * Its a simple algorithm and has no perspective
 */

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <cmath>

using namespace std;

struct Vector {
public:
    float x;
    float y;
    float z;

    Vector() {
        this->x = 0;
        this->y = 0;
        this->z = 0;
    }

    Vector(float x, float y, float z) {
        this->x = x;
        this->y = y;
        this->z = z;
    }

    Vector operator-(const Vector &vector) const { return {x - vector.x, y - vector.y, z - vector.z}; }

    Vector operator+(const Vector &vector) const { return {x + vector.x, y + vector.y, z + vector.z}; }

    Vector operator*(float d) const { return {x * d, y * d, z * d}; }

    Vector operator/(float d) const { return {x / d, y / d, z / d}; }

    [[std::nodiscard]] Vector normalize() const {
        float mg = sqrt(x * x + y * y + z * z);
        return {x / mg, y / mg, z / mg};
    }
};

inline float dot(const Vector &vector1, const Vector &vector2) {
    return (vector1.x * vector2.x + vector1.y * vector2.y + vector1.z * vector2.z);

}

struct Ray {
public:
    Vector origin;
    Vector direction;

    Ray(const Vector &origin, const Vector &direction) {
        this->origin = origin;
        this->direction = direction;
    }
};

struct Sphere {
public:
    Vector center;
    float radius;

    Sphere(const Vector &center, float radius) {
        if (radius > 0) { this->radius = radius; }
        else { this->radius = 0; }
        this->center = center;

    }

    [[nodiscard]] Vector getNormal(const Vector &pi) const {
        return (pi - center) / radius;
    }

    bool intersection(Ray ray, float &t) const {
        Vector origin = ray.origin;
        Vector direction = ray.direction;
        Vector originMinusCenter = origin - center;
        float b = 2 * dot(originMinusCenter, direction);
        float c = dot(originMinusCenter, originMinusCenter) - radius * radius;
        float delta = b * b - (4 * c);
        if (delta < 0) {
            return false;
        } else {
            delta = sqrt(delta);
            float t0 = (-b - delta) / 2;
            float t1 = (-b + delta) / 2;
            t = (t0 < t1) ? t0 : t1;
            return true;
        }
    }

    //added new operators for multiple light sources but just remember you cant create a sphere with negative radius
    //also there is a small problem with G mass of center,I will fix it.
    Sphere operator+(const Sphere &sphere) const {
        return {Vector(this->center.x * radius + sphere.center.x * radius, this->center.y * radius +
                                                                           sphere.center.y * radius,
                       this->center.z * radius + sphere.center.z * radius),
                (this->radius * this->radius + radius * radius) / 2};
    }

    Sphere operator-(const Sphere &sphere) const {
        return {Vector(this->center.x * radius - sphere.center.x * radius, this->center.y * radius -
                                                                           sphere.center.y * radius,
                       this->center.z * radius - sphere.center.z * radius),
                (this->radius * this->radius - radius * radius) / 2};
    }

};

struct Color {
public:
    float r;
    float g;
    float b;

    Color() {
        this->r = 0;
        this->g = 0;
        this->b = 0;
    }

    Color(float r, float g, float b) {
        this->r = r;
        this->g = g;
        this->b = b;
    }

    Color operator*(float d) const { return {r * d, g * d, b * d}; }

    Color operator+(Color color) const { return {(r + color.r) / 2, (g + color.g) / 2, (b + color.b) / 2}; }

};

int main() {

    const short width = 410;
    const short height = 410;
    std::ofstream out("out.ppm");
    out << "P3\n" << width << "\n" << height << "\n" << "255\n";
    Color pixel_color[410][410];
    Color white(255, 255, 255);

    /**
     * if you want to change the color of sphere just change r g b parameters of Color red.
     */

    Color sphereColor(255, 0, 0);
    Sphere sphere(Vector(width / 2.0, height / 2.0, 50), 20);
    Sphere light(Vector(0, 0, 50), 1);

    /**
     * multiple light sources option will be added in future but for now you can just create a new sphere object just
     * like below and add their coordinates with other spheres coordinates.
     * you have to remember that the radius of light source matters(So you have to use G(center of mass) while doing
     * +- operation.) and want to add two spheres you have to add x1 and x2 coordinates , y1 and y2 so on
     * there is no + or - operator for sphere struct but im going to add this operator in future.
     * and also remember image is 410x410 so these are your borders.
     */

//    Sphere light2(Vector(410,0,50),1);
//    Sphere light3((light.center+light2.center)/2,1);

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            Ray ray(Vector(float(i), float(j), -50), Vector(0, 0, 1));
            float t = 20000;
            if (sphere.intersection(ray, t)) {
                Vector P = ray.origin + ray.direction * t;
                Vector L = light.center - P;
                Vector N = sphere.getNormal(P);
                float dt = dot(L.normalize(), N.normalize());
                pixel_color[j][i] = sphereColor + (white * dt) * 1.0;
                if (pixel_color[j][i].r < 0) {
                    pixel_color[j][i].r = 0;
                }
                if (pixel_color[j][i].g < 0) {
                    pixel_color[j][i].g = 0;
                }
                if (pixel_color[j][i].b < 0) {
                    pixel_color[j][i].b = 0;
                }
                if (pixel_color[j][i].r > 255) {
                    pixel_color[j][i].r = 255;
                }
                if (pixel_color[j][i].g > 255) {
                    pixel_color[j][i].g = 255;
                }
                if (pixel_color[j][i].b > 255) {
                    pixel_color[j][i].b = 255;
                }
            }
            out << (int) pixel_color[j][i].r << std::endl;
            out << (int) pixel_color[j][i].g << std::endl;
            out << (int) pixel_color[j][i].b << std::endl;
        }
    }
    return 0;
}
