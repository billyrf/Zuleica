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

Color3 gamma(Color3 color, float value){
    float inverseGamma = (1 / value);
	return Color3(pow(color.r, inverseGamma), pow(color.g, inverseGamma), pow(color.b, inverseGamma));
}
	
Color3 exposure(Color3 color, float value){
	float power = pow(2, value);
	return Color3(color.r * power, color.g * power, color.b * power);
}
	
struct ShaderGlobals {
    Vector3 point;
    Vector3 normal;
    Vector3 tangentU;
    Vector3 tangentV;
    Vector3 viewDirection;
    Vector3 lightDirection;
    Vector3 lightPoint;
    Vector3 lightNormal;
    Vector2 textureCoordinate;
    Vector2 uv;
    
    ShaderGlobals() {}
    ShaderGlobals(const Vector3 & point, const Vector3 & normal, const Vector2 & uv,const Vector3 & tangentU, const Vector3 & tangentV,const Vector3 & viewDirection, const Vector3 & lightDirection,const Vector3 & lightPoint, const Vector3 & lightNormal, const Vector2 & textureCoordinate) {
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

struct Intersection {
    bool hit;
    float distance;
    int index;
    
    Intersection() {
        hit = false;
        distance = AURORA_INFINITY;
        index = -1;
    }
    Intersection(bool hit, float distance, int index) {
        this->hit = hit;
        this->distance = distance;
        this->index = index;
    }
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

struct Scene {
    std::vector<Shape *> shapes;
    
    Scene();
    Scene(const std::vector<Shape *> & shapes) {
        this->shapes = shapes;
    }
    
    bool intersects(const Ray & ray, Intersection & intersection) const {
        for (int i = 0; i < shapes.size(); i++) {
            Shape * shape = shapes[i];
            
            Intersection temp;
            shape->intersects(ray, temp);
            
            if (temp.hit && temp.distance < intersection.distance) {
                intersection.hit = temp.hit;
                intersection.distance = temp.distance;
                intersection.index = i;
            }
        }
        
        return intersection.hit;
    }
};

struct Sphere : Shape {
    Vector3 position;
    float radius;
    
    Sphere() : Shape() {}
    Sphere(const Vector3 & position, float radius, BSDF * bsdf) : Shape(bsdf) {
        this->position = position;
        this->radius = radius;
    }
    
    virtual bool intersects(const Ray & ray, Intersection & intersection) const {
        Vector3 l = position - ray.origin;
        float t = l.dot(ray.direction);
        
        if (t < 0)
            return false;
            
        float d2 = l.dot(l) - t * t;
        float r2 = radius * radius;
        
        if (d2 > r2)
            return false;
        
        float dt = std::sqrt(r2 - d2);
        
        float t0 = t - dt;
        float t1 = t + dt;
        
        if (t0 > t1)
            std::swap(t0, t1);
        
        if (t0 < 0) {
            t0 = t1;
            
            if (t0 < 0)
                return false;
        }
        
        intersection.hit = true;
        intersection.distance = t0;
        
        return true;
    }
    virtual void calculateShaderGlobals(const Ray & ray, const Intersection & intersection,ShaderGlobals & shaderGlobals) const {
        shaderGlobals.point = ray.point(intersection.distance);
        shaderGlobals.normal = (shaderGlobals.point - position).normalize();
        
        float theta = std::atan2(shaderGlobals.normal.x, shaderGlobals.normal.z);
        float phi = std::acos(shaderGlobals.normal.y);
        
        shaderGlobals.uv.x = theta * AURORA_INV_PI * 0.5;
        shaderGlobals.uv.y = phi * AURORA_INV_PI;
        
        shaderGlobals.tangentU.x = std::cos(theta);
        shaderGlobals.tangentU.y = 0;
        shaderGlobals.tangentU.z = -std::sin(theta);
        
        shaderGlobals.tangentV.x = std::sin(theta) * std::cos(phi);
        shaderGlobals.tangentV.y = -std::sin(phi);
        shaderGlobals.tangentV.z = std::cos(theta) * std::cos(phi);
        
        shaderGlobals.viewDirection = -ray.direction;
    }
    virtual float surfaceArea() const {
        return 4.0 * AURORA_PI * radius * radius;
    }
};

struct Vertex {
    Vector3 position;
    Vector3 normal;
    Vector2 textureCoordinate;

    Vertex(){}
    Vertex(const Vector3 & position , const Vector3 normal, const Vector2 textureCoordinate){
        this->position = position;
        this->normal = normal;
        this->textureCoordinate = textureCoordinate;
    }
};

struct Triangle : Shape {
    Vertex vertices[3];
    float area;

    Triangle() : Shape() {}
    Triangle(const Vertex & vertex0,const Vertex & vertex1, const Vertex & vertex2,BSDF * bsdf) : Shape(bsdf){
        this->vertices[0] = vertex0;
        this->vertices[1] = vertex1;
        this->vertices[2] = vertex2;

        updateSurfaceArea();
    }
    virtual bool intersects(const Ray & ray, Intersection & intersection) const {
        const Vector3 & v0 = vertices[0].position;
        const Vector3 & v1 = vertices[1].position;
        const Vector3 & v2 = vertices[2].position;
        
        Vector3 u = v1 - v0;
        Vector3 v = v2 - v0;
        
        Vector3 p = ray.direction.cross(v);
        float d = u.dot(p);
        
        if (std::abs(d) < AURORA_EPSILON)
            return false;
        
        Vector3 t = ray.origin - v0;
        float inverseD = 1.0f / d;
        
        float alpha = t.dot(p) * inverseD;
        
        if (alpha < 0.0f || alpha > 1.0f)
            return false;
            
        Vector3 q = t.cross(u);
        
        float beta = ray.direction.dot(q) * inverseD;
        
        if (beta < 0.0f || alpha + beta > 1.0f)
            return false;
        
        float t0 = v.dot(q) * inverseD;
        
        if (t0 < AURORA_EPSILON)
            return false;
        
        intersection.hit = true;
        intersection.distance = t0;
        
        return true;
    }
    virtual void calculateShaderGlobals( const Ray & ray, const Intersection & intersection, ShaderGlobals & shaderGlobals) const {
    	
        const Vector3 & p0 = vertices[0].position;
        const Vector3 & p1 = vertices[1].position;
        const Vector3 & p2 = vertices[2].position;
        
        const Vector3 & n0 = vertices[0].normal;
        const Vector3 & n1 = vertices[1].normal;
        const Vector3 & n2 = vertices[2].normal;
        
        const Vector2 & t0 = vertices[0].textureCoordinate;
        const Vector2 & t1 = vertices[1].textureCoordinate;
        const Vector2 & t2 = vertices[2].textureCoordinate;
        
        Vector3 b = barycentric(shaderGlobals.point, p0, p1, p2);
        
        shaderGlobals.point = ray.point(intersection.distance);
        
        shaderGlobals.normal = (n0 * b.x + n1 * b.y + n2 * b.z).normalize();
        shaderGlobals.textureCoordinate = t0 * b.x + t1 * b.y + t2 * b.z;
        
        shaderGlobals.uv = Vector2(b.x, b.y);
        
        calculateTangents(
            shaderGlobals.normal,
            shaderGlobals.tangentU,
            shaderGlobals.tangentV);
        
        shaderGlobals.viewDirection = -ray.direction;
    }
    virtual float surfaceArea() const {
        return area;
    }
    
    Triangle & updateSurfaceArea() {
        const Vector3 & v0 = vertices[0].position;
        const Vector3 & v1 = vertices[1].position;
        const Vector3 & v2 = vertices[2].position;
        
        Vector3 normal = calculateVectorArea(v0, v1, v2);
        area = 0.5f * normal.length();
        
        return *this;
    }

};

struct Film {
    float width;
    float height;
    
    Film() {}
    Film(float width, float height) {
        this->width = width;
        this->height = height;
    }
    
    float aspectRatio() const {
        return width / height;
    }
};

struct Camera {
    float fieldOfView;
    Film film;
    Matrix4 worldMatrix;
    
    Camera() {}
    Camera(float fieldOfView, const Film & film, const Matrix4 & worldMatrix) {
        this->fieldOfView =fieldOfView;
        this->film = film;
        this->worldMatrix = worldMatrix;
    }
    
    void lookAt(const Vector3 & position, const Vector3 & target, const Vector3 & up) {
        Vector3 w = (position - target).normalize();
        Vector3 u = up.cross(w).normalize();
        Vector3 v = w.cross(u);
        
        worldMatrix[0][0] = u.x;
        worldMatrix[0][1] = u.y;
        worldMatrix[0][2] = u.z;
        worldMatrix[0][3] = 0;
        
        worldMatrix[1][0] = v.x;
        worldMatrix[1][1] = v.y;
        worldMatrix[1][2] = v.z;
        worldMatrix[1][3] = 0;
        
        worldMatrix[2][0] = w.x;
        worldMatrix[2][1] = w.y;
        worldMatrix[2][2] = w.z;
        worldMatrix[2][3] = 0;
        
        worldMatrix[3][0] = position.x;
        worldMatrix[3][1] = position.y;
        worldMatrix[3][2] = position.z;
        worldMatrix[3][3] = 1.0;
    }
    Ray generateRay(float x, float y, const Vector2 & sample) const {
        float scale = std::tan(fieldOfView * 0.5);
        
        Vector3 pixel;
        
        pixel.x = (2.0 * (x + sample.x + 0.5) / film.width - 1.0) * scale * film.aspectRatio();
        pixel.y = (1.0 - 2.0 * (y + sample.y + 0.5) / film.height) * scale;
        pixel.z = -1.0;
        
        pixel *= worldMatrix;
        
        Vector3 position(worldMatrix[3][0], worldMatrix[3][1], worldMatrix[3][2]);
        Vector3 direction = (pixel - position).normalize();
        
        return Ray(position, direction);
    }
};

struct RenderOptions {
    int width;
    int height;
    int maximumDepth;
    int cameraSamples;
    int lightSamples;
    int diffuseSamples;
    float filterWidth;
    float gamma;
    float exposure;
    
    RenderOptions() {}
    RenderOptions(int width, int height, int maximumDepth,int cameraSamples, int lightSamples, int diffuseSamples, float filterWidth, float gamma, float exposure) {
        this->width = width;
        this->height = height;
        this->maximumDepth = maximumDepth;
        this->cameraSamples = cameraSamples;
        this->lightSamples = lightSamples;
        this->diffuseSamples = diffuseSamples;
        this->filterWidth = filterWidth;
        this->gamma = gamma;
        this->exposure = exposure;
    }
};

struct Renderer {
    RenderOptions * options;
    Camera * camera;
    Scene * scene;
    
    Renderer() {}
    Renderer(RenderOptions * options, Camera * camera, Scene * scene) {
        this->options = options;
        this->camera = camera;
        this->scene = scene;
    }
    
    Color3 trace(const Ray & ray, int depth) const {
        
        Intersection intersection;

        if (scene->intersects(ray, intersection))
            return Color3(1.0, 1.0, 1.0);

        return Color3(0, 0, 0);
        
    };
    
    Color3 render(Image3 * image) const {
		
		const Vector2 Point(0.5, 0.5);
		
        for (int i = 0; i < options->width; i++) {
            for (int j = 0; j < options->height; j++) {
				std::vector<Vector2> samples;
                stratifiedSample(options->cameraSamples,samples);
                
                Color3 color = Color3(0, 0, 0);
                float totalWeight = 0;
                
                for (int k = 0; k < options->cameraSamples; k++) {
                    Vector2 sample = (samples[k] - Point) * options->filterWidth;
                    Ray ray = camera->generateRay(i, j, sample);
                    
                    float weight = gaussian2D(sample, options->filterWidth);
                    
                    color += trace(ray, 0) * weight;
                    totalWeight += weight;
                }
                
                color /= totalWeight;
                
                image->setPixel(i, j, color);
                /*image[i][j] = saturate(gamma(exposure(color, exposureValue), gammaValue)) * 255;*/
            }
        }
       
    }
};

int main(int argc, char ** argv) {

	RenderOptions options(600, 600, 1, 4, 1, 1, 2.0, 2.0, 0);
	
    Film film(800, 600);
    
    Camera camera(radians(20.0), film, Matrix4());
    camera.lookAt(Vector3(0, 0, 35.0), Vector3(0, 0, 0), Vector3(0, 1.0, 0));
    
    BSDF * bsdf = new BSDF(BSDFType::Diffuse, Color3(1.0f, 1.0f, 1.0f));
    
    Shape * sphere = new Sphere(Vector3(0.0f, 0.0f, 0.0f), 1.0f, bsdf);
    
    Vertex vertex0(Vector3(-0.5f, -0.5f, 0.0f), Vector3(0.0f, 0.0f, 1.0f), Vector2(0.0f, 0.0f));
    Vertex vertex1(Vector3(0.5f, -0.5f, 0.0f), Vector3(0.0f, 0.0f, 1.0f), Vector2(1.0f, 0.0f));
    Vertex vertex2(Vector3(0.0f, 1.0f, 0.0f), Vector3(0.0f, 0.0f, 1.0f), Vector2(0.0f, 1.0f));
    
    Shape * triangle = new Triangle(vertex0, vertex1, vertex2, bsdf);
    
    std::vector<Shape *> shapes;
    
    shapes.push_back(sphere);
    shapes.push_back(triangle);
    
    Scene scene(shapes);
    
    Ray ray(Vector3(0.0f, 0.0f, 10.0f), Vector3(0.0f, 0.0f, -1.0f));
    Intersection intersection;
    
    scene.intersects(ray, intersection);
    
    std::cout << "Hit: " << intersection.hit << std::endl;
    std::cout << "Distance: " << intersection.distance << std::endl;
    std::cout << "Index: " << intersection.index << std::endl;
        
    if (intersection.hit) {
        Shape * shape = scene.shapes[intersection.index];
        
        ShaderGlobals shaderGlobals;
        shape->calculateShaderGlobals(ray, intersection, shaderGlobals);
        
        std::cout << "Point: " << shaderGlobals.point << std::endl;
        std::cout << "Normal: " << shaderGlobals.normal << std::endl;
        std::cout << "Texture coordinate: " << shaderGlobals.textureCoordinate << std::endl;
        std::cout << "UV: " << shaderGlobals.uv << std::endl;
        std::cout << "Tangent U: " << shaderGlobals.tangentU << std::endl;
        std::cout << "Tangent V: " << shaderGlobals.tangentV << std::endl;
        std::cout << "View direction: " << shaderGlobals.viewDirection << std::endl;
        std::cout << "Light direction: " << shaderGlobals.lightDirection << std::endl;
        std::cout << "Light point: " << shaderGlobals.lightPoint << std::endl;
        std::cout << "Light normal: " << shaderGlobals.lightNormal << std::endl;
    }
    
    Renderer renderer(&options, &camera, &scene);
    
    Image3 * image = new Image3(options.width, options.height);
   	renderer.render(image);
    
    if (writeImage("saida.ppm", image))
        std::cout << "Success." << std::endl;
    else
        std::cout << "Failure." << std::endl;
        
    delete bsdf;
    delete sphere;
    delete triangle;

    return 0;
}