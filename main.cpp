#include <iostream>
#include <math.h>
#include <GLUT/glut.h>
#include <vector>

#define WIDTH 256*2
#define HEIGHT 256*2
#define MAX_RAY_DEPTH 3
#define INFINITY 1e8

using namespace ::std;

double a = 0;
GLfloat image[WIDTH][HEIGHT][3];


class Vec3 {
public:
    float x,y,z;
    Vec3():x(0),y(0),z(0){}
    Vec3(float xx):x(xx),y(xx),z(xx){}
    Vec3(float xx, float yy, float zz) : x(xx), y(yy), z(zz) {}
    
    Vec3 operator+(const Vec3 &vector) {
		return Vec3(x + vector.x, y + vector.y, z + vector.z);
	}
	Vec3 operator-(const Vec3 &vector) {
		return Vec3(x - vector.x, y - vector.y, z - vector.z);
	}
    // this is used for calculate ambient diffuse specular
	Vec3 operator*(const Vec3 &vector) {
        return Vec3(x*vector.x,y*vector.y,z*vector.z);
	}
    
    Vec3 cross(const Vec3 &vector){
        return Vec3(y * vector.z - z * vector.y, z * vector.x - x * vector.z, x * vector.y - y*vector.x);
    }
	Vec3 operator*(const float val) {
		return Vec3(x*val, y *val, z *val);
	}
    
    bool operator==(const Vec3 &vector){
        return x==vector.x && y==vector.y && z == vector.z;
    }
    
    float dot(const Vec3 &vector){
        return x * vector.x + y * vector.y + z * vector.z;
    }
    
	Vec3 normalize() {
		float mod = sqrtf(x * x + y * y + z * z);
		return Vec3(x / mod, y / mod, z / mod);
	}
    
};

class Ray {
public:
    Vec3 origin,direction;
    Ray(){}
    Ray(Vec3 o,Vec3 d):origin(o),direction(d.normalize()){}
};

class Light {
public:
    Vec3 direction;
    //Ia : ambient light intesity ; I : intensity of light sourse
    Vec3 Ia , I;
    Light() {}
    Light(Vec3 d,Vec3 Ia, Vec3 I):direction(d.normalize()),Ia(Ia),I(I){}
};

class Material {
public:
    //Color ka: ambient , kd: diffuse , ks : specular
    Vec3 ka,kd,ks;
    // reflection:  for recursive ray tracing
    float reflection;
    // n_exp : power
    float n_exp;
    Material(){}
    Material(Vec3 ka,Vec3 kd, Vec3 ks ,float reflection, float n_exp):
            ka(ka),kd(kd),ks(ks),reflection(reflection),n_exp(n_exp){}
    
};

class Surface {
public:
    Material material;
    virtual bool intersect(Ray ray , float *t0 = NULL, float *t1 = NULL) = 0;
    virtual Vec3 getNormalAtIntersection(Ray ray,Vec3 phit) = 0;
};

class Sphere : public Surface{
public:
    Vec3 center;
    float radius;
    Sphere(Vec3 c,float r,Material m):center(c),radius(r){ material = m;}
    bool intersect(Ray ray , float *t0 = NULL, float *t1 = NULL) {
        Vec3 l = center - ray.origin;  // c -e
        float C = l.dot(l) - radius*radius; // A is equal to 1
        float B = l.dot(ray.direction);
        if(B < 0) return false;
        float D = B*B - C;
        if (D < 0.0f) return false;
        
        if (t0 != NULL && t1 != NULL) {
			*t0 = B - sqrtf(D);
			*t1 = B + sqrtf(D);
		}
        return true;
    }
    
    Vec3 getNormalAtIntersection(Ray ray,Vec3 phit) {
        return (phit - center).normalize();
    }
};

class Tetrahedron : public Surface {
public:
    Vec3 vertex[4];
	unsigned int hitTriangleStartVertex;
    Tetrahedron(Vec3 v1, Vec3 v2, Vec3 v3, Vec3 v4, Material m) {
		vertex[0] = v1; vertex[1] = v2; vertex[2] = v3; vertex[3] = v4;
		material = m;
		hitTriangleStartVertex = 0;
	}
    
    bool intersect(Ray ray, float *t0 = NULL, float *t1 = NULL) {
		bool result[4];
		for (int i = 0; i < 4; i++) {
			result[i] = intersectTriangle(i, ray, t0, t1);
		}
		return result[0] || result[1] || result[2] || result[3];
	}
    
    bool intersectTriangle(unsigned int startVertex,Ray ray , float *t0 = NULL, float *t1 = NULL) {
        Vec3& v0 = vertex[startVertex % 4];
		Vec3& v1 = vertex[(startVertex+1) % 4];
		Vec3& v2 = vertex[(startVertex+2) % 4];
        Vec3 v0v1 = v1 - v0;
        Vec3 v0v2 = v2 - v0;
        //no need to nomalize
        Vec3 N = v0v1.cross(v0v2).normalize();
     
        float nDotRay = N.dot(ray.direction);
        //ray is  parallel to triangle plane
        if (nDotRay == 0) {
            return false;
        }
        //
        float d = -v0.dot(N);
        float t = -(N.dot(ray.origin) +d) / nDotRay;
        // check whether the ray goes away from triangle
        if(t < 0.2) return false;
        //compute the intersection point
        Vec3 Phit = ray.origin + ray.direction *t;
        
        // edge0
        
        Vec3 v0p = Phit - v0;
        float v = N.dot(v0v1.cross(v0p));
        if (v<0) {
            return false;
        }
        // edge1
        Vec3 v1p = Phit - v1;
        Vec3 v1v2 = v2 - v1;
        float w = N.dot(v1v2.cross(v1p));
        if (w<0) {
            return false;
        }
        // edge2
        Vec3 v2p = Phit - v2;
        Vec3 v2v0 = v0 - v2;
        float u = N.dot(v2v0.cross(v2p));
        if (u<0) {
            return false;
        }
        *t0 = t;
        return true;
    }

    
    Vec3 getNormalAtIntersection(Ray ray,Vec3 phit) {
        Vec3 v0 = vertex[hitTriangleStartVertex % 4];
		Vec3 v1 = vertex[(hitTriangleStartVertex+1) % 4];
		Vec3 v2 = vertex[(hitTriangleStartVertex+2) % 4];
        Vec3 v0v1 = v1 - v0;
        Vec3 v0v2 = v2 - v0;
        Vec3 N = v0v1.cross(v0v2).normalize();
        
        if (N.dot(ray.direction) > 0) {
            return N * -1;
        }
        return N;
    }
    
};

class Scenario {
public:
    Vec3 lookAt;
    Light light;
    Vec3 Ia;
    Vec3 I;
    vector<Surface *> surfaces;
    GLfloat image[WIDTH][HEIGHT][3];
    
    Scenario(){}
    
    void setViewer(Vec3 &look){
        this->lookAt = look;
    }
    void addSurface(Surface *s){
        this->surfaces.push_back(s);
    }
    void addLight(Light &l){
        this->light = l;
    }
    void setIa(Vec3 &Ia){
        this->Ia = Ia;
    }
    void setI(Vec3 &I){
        this->I = I;
    }
    void render(){
        Vec3 planeBase(0.0);
        Vec3 horInc(1.0, 0.0, 0.0);
        Vec3 verInc(0.0, 1.0, 0.0);
        Vec3 depInc(0.0, 0.0, 1.0);

        //Light light(Vec3(-cos(a)*sqrt(2), 1.0, -sin(a)*sqrt(a)), Ia, I);
        Vec3 lookAt(-sin(a), 0.0, -cos(a));
        for (int i=0; i<HEIGHT; i++) {
            for (int j =0; j<WIDTH; j++) {
                GLfloat pos_x = (WIDTH / 2) * sin(a) - (j - WIDTH / 2) * cos(a);
                GLfloat pos_y = i - HEIGHT / 2;
                GLfloat pos_z = (WIDTH / 2) * cos(a) + (j - WIDTH / 2) * sin(a);
                
                Vec3 origin = planeBase + horInc * pos_x + verInc * pos_y + depInc * pos_z;
                
                Ray ray(origin, lookAt);
                Vec3 color = trace(ray, surfaces, 0,light,lookAt);
                image[i][j][0] = color.x;
                image[i][j][1] = color.y;
                image[i][j][2] = color.z;
            }
        }
    }
    Vec3 trace(Ray ray, vector<Surface *> &surfaces , const int &depth,Light light,Vec3 lookAt) {
        Vec3 color(0); // returned color
        
        float tnear = INFINITY;
        Surface *surface = NULL;
        // find neareast intersection
        for (unsigned int i =0; i< surfaces.size(); i++) {
            float t0 = INFINITY , t1 = INFINITY;
            if (surfaces[i]->intersect(ray,&t0,&t1)) {
                if(t0 < 0) t0 = t1;
                if (t0 < tnear) {
                    tnear = t0;
                    surface = surfaces[i];
                }
            }
        }
        
        if (!surface) return color;
        
        
        Vec3 phit = ray.origin + ray.direction * tnear;
        Vec3 nhit = surface->getNormalAtIntersection(ray,phit);
        float bias  = 1e-4; // add some bias to the point from which we will be tracing
        
        //if it is glazed surface
        if (surface->material.reflection >0 && depth < MAX_RAY_DEPTH) {
            Vec3 refldir = ray.direction - nhit *2 * ray.direction.dot(nhit);
            refldir.normalize();
            Vec3 reflection = trace(Ray(phit+nhit * bias,refldir), surfaces, depth + 1,light,lookAt);
            color =color + reflection;
            
        }
        //determine whether the object is in shadow
        Ray shadowRay(phit,light.direction * -1);
        bool isInShadow = false;
        for (int k =0; k<surfaces.size(); k++) {
            float t0,t1;
            if (surfaces[k]->intersect(shadowRay,&t0,&t1)) {
                isInShadow = true;
                break;
            }
        }
        // calculate ambient,diffuse,specular
        Vec3 ambient = surface->material.ka * light.Ia;
        if(isInShadow) return ambient;
        Vec3 N = surface->getNormalAtIntersection(ray,phit);
        Vec3 L = light.direction * -1;
        Vec3 V = lookAt * -1;
        Vec3 H = (L+V).normalize();
        
        Vec3 diffuse = surface->material.kd * light.I * max(0.0f,N.dot(L));
        Vec3 specular = surface->material.ks * light.I * powf(max(0.0f, N.dot(H)),surface->material.n_exp);
        color =color + ambient + diffuse + specular;
        return color;
    }

    
};

void init(Scenario &scenario) {
	glClearColor(0.0, 0.0, 0.0, 0.0);
    Material m1(Vec3(0.2,0.6,0.3),Vec3(0.2,0.6,0.3),Vec3(0.2,0.6,0.3),0,10);
    Sphere sphere1(Vec3(60.0, -20.0, -15.0), 60.0, m1);
    Material m2(Vec3(0.1, 0.4, 0.7),Vec3(0.1, 0.4, 0.7),Vec3(0.1, 0.4, 0.7),1,100);
    Sphere sphere2(Vec3(-40.0, 0.0, 15.0), 40.0, m2);
    Material m3(Vec3(0.4, 0.7, 0.1), Vec3(0.4, 0.7, 0.1), Vec3(0.5), 1,100.0);
    Tetrahedron tetra1(
                       Vec3(0.0, 40.0, 0.0),
                       Vec3(20.0, 100.0, -30.0),
                       Vec3(70.0, 40.0, -30.0),
                       Vec3(20.0, 40.0, -60.0),
                       m3);
    
    Vec3 lookAt=Vec3(0.0,0.0,-1.0).normalize();
    
    Vec3 Ia(0.4);
    Vec3 I(1.0);
    Light light(Vec3(-1.0, 1.0, -1.0), Ia, I);
    
    scenario.setViewer(lookAt);
    scenario.addLight(light);
    scenario.setIa(Ia);
    scenario.setI(I);
    scenario.addSurface(&sphere1);
    scenario.addSurface(&sphere2);
    scenario.addSurface(&tetra1);
	
    scenario.render();
}

void idle() {
     a = a + 0.1;
     glutPostRedisplay();
}

void display() {
	glClear(GL_COLOR_BUFFER_BIT);
    Scenario scenario;
    init(scenario);
	glDrawPixels(WIDTH, HEIGHT, GL_RGB, GL_FLOAT, scenario.image);
	glFlush();
    glutSwapBuffers();
}

int main(int argc, char** argv) {
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize(WIDTH, HEIGHT);
	glutInitWindowPosition(150, 120);
	glutCreateWindow("ray tracing");
	glutDisplayFunc(display);
    glutIdleFunc( idle );
	glutMainLoop();
	return 0;
}

