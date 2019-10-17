#include <aurora/Math.h>
#include <aurora/Utility.h>
#include <aurora/Vector.h>
#include <aurora/Color.h>
#include <aurora/Matrix.h>
#include <aurora/Image.h>

#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>

using namespace aurora;

float uniformRandom1D() {
    return uniformRandom();
}
Vector2 uniformRandom2D() {
    return Vector2(uniformRandom1D(), uniformRandom1D());
}

void stratifiedSample(int count, std::vector<Vector2> & samples) {
    samples.reserve(count);
    
    int size = std::sqrt(count);
    float inverseSize = 1.0 / size;
    
    for (int i = 0; i < count; i++) {
        Vector2 offset(i / size, i % size);
        Vector2 point = (offset + uniformRandom2D()) * inverseSize;
        
        samples.push_back(point);
    }
}

float gaussian1D(float sample, float width) {
    float radius = width * 0.5;
    return std::fmax(0, std::exp(-sample * sample) - std::exp(-radius * radius));
}
float gaussian2D(const Vector2 & sample, float width) {
    return gaussian1D(sample.x, width) * gaussian1D(sample.y, width);
}

struct Intersection {
    bool hit;
    float distance;
    int index;
    
    Intersection() {}
    Intersection(bool hit, float distance, int index) {
        this->hit = hit;
        this->distance = distance;
        this->index = index;
    }
};

struct Scene {
    std::vector<Shape *> shapes;
    
    Scene();
    Scene(const std::vector<Shape *> & shapes) {
        this->shapes = shapes;
    }
    
    bool intersects(const Ray & ray, Intersection & intersection) const {}
};

struct Ray {
    Vector3 origin;
    Vector3 direction;
    
    Ray() {}
    Ray(const Vector3 & origin, const Vector3 & direction) {
        this->origin = origin;
        this->direction = direction;
    }
    
    Vector3 point(float distance) const {
        return origin + direction * distance;
    }
};

enum BSDFType {
    Light = 0,
    Diffuse = 1,
    Specular = 2,
    None = 3
};

struct BSDF {
    BSDFType type;
    Color3 color;
    
    BSDF() {}
    BSDF(BSDFType type, const Color3 & color) {
        this->type = type;
        this->color = color;
    }
};

struct Shape {
    BSDF * bsdf;
    
    Shape() {}
    Shape(BSDF * bsdf) {
        this->bsdf = bsdf;
    }
    
    virtual bool intersects(const Ray & ray, Intersection & intersection) const = 0;
    virtual void calculateShaderGlobals(const Ray & ray, const Intersection & intersection,ShaderGlobals & shaderGlobals) const = 0;
    virtual float surfaceArea() const = 0;
};

struct Sphere : Shape {
    Vector3 position;
    float radius;
    
    Sphere() : Shape() {}
    Sphere(const Vector3 & position, float radius, BSDF * bsdf) : Shape(bsdf) {
        this->position = position;
        this->radius = radius;
    }
    
};

struct ShaderGlobals {
    Vector3 point;
    Vector3 normal;
    Vector2 uv;
    Vector3 tangentU;
    Vector3 tangentV;
    Vector3 viewDirection;
    Vector3 lightDirection;
    Vector3 lightPoint;
    Vector3 lightNormal;
    
    ShaderGlobals() {}
    ShaderGlobals(const Vector3 & point, const Vector3 & normal, const Vector2 & uv,const Vector3 & tangentU, const Vector3 & tangentV,const Vector3 & viewDirection, const Vector3 & lightDirection,const Vector3 & lightPoint, const Vector3 & lightNormal) {
        this->point = point;
        this->normal = normal;
        this->uv = uv;
        this->tangentU = tangentU;
        this->tangentV = tangentV;
        this->viewDirection = viewDirection;
        this->lightDirection = lightDirection;
        this->lightPoint = lightPoint;
        this->lightNormal = lightNormal;
    }
};

struct Triangle {
    Vertex vertices[3];

    Triangle() : Triangle() {}
    Triangle(Vertex vertices , BSDF* bsdf){
        this->vertices = vertices;
    }
}

struct Vertex {
    Vector3 position;
    Vector3 normal;
    Vector2 uv;

    Vertex() : Vertex(){}
    Vertex(const Vector3 & position , const Vector3 normal, const Vector2 uv){
        this->position = position;
        this->normal = normal;
        this->uv = uv;
    }
}

int main(int argc, char ** argv) {
    RenderOptions options(500, 500, 1, 4, 1, 1, 2.0, 2.2, 0);
    
    return 0;
}
