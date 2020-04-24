#include<bits/stdc++.h>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include"stb_image_write.h"
double pi = acos(-1);

struct Vector{
double x, y, z;
	Vector(){}
	Vector(double x, double y, double z):
		x(x),
		y(y),
		z(z){}
	Vector& operator+=(const Vector& b) {
		x += b.x;
		y += b.y;
		z += b.z;
		return * this ;
	}
};

Vector operator+(const Vector& a, const Vector& b) {
	return Vector(a.x+b.x, a.y+b.y, a.z+b.z);
}

Vector operator-(const Vector& a, const Vector& b) {
	return Vector(a.x-b.x, a.y-b.y, a.z-b.z);
}

Vector operator*(double t, const Vector& a) {
	return Vector(a.x*t, a.y*t, a.z*t);
}

Vector operator*(const Vector& a, double t) { //I want it to go both ways because sometimes I like to switch things around to add some spice. also I'm autistic :())
	return Vector(a.x*t, a.y*t, a.z*t);
}

double dot(const Vector& a, const Vector& b) {
	return double (a.x*b.x + a.y*b.y + a.z*b.z);
}

Vector pairwiseyoink(const Vector& a, const Vector& b) {
	return Vector (a.x*b.x, a.y*b.y, a.z*b.z);
}

struct Sphere{
	Vector C;
	Vector albedo;
	double R;
	Sphere(){}
	Sphere(const Vector& c, const Vector& albedo, double R):
		C(c),albedo(albedo), R(R)
	{}
};

struct Ray{
	Vector u;
	Vector O;
	Ray(){}
	Ray(const Vector& u, const Vector& O):
		u(u),
		O(O){}
};

struct Intersection{
	Vector point;
	Vector normal;
	double dist;
	bool hmm;
	Intersection(){}
	Intersection(const Vector& point, const Vector& normal, double t, bool b):
		point(point), normal(normal), dist(t), hmm(b)
	{}
};

struct Lightsource{
	Vector spot;
	double I;
	Vector hue;
	Lightsource(){}
	Lightsource(const Vector& spot, double i, const Vector& c):
		spot(spot), I(i), hue(c)
	{}
}


struct Scene{
	std::vector<Sphere> welt;
	std::vector<Lightsource> lights;
};

double norm(const Vector& V){
	return sqrt(dot(V, V));
}

Vector normalize(const Vector& V){
	return V*(1/norm(V));
}

Intersection intersect(const Sphere& sf, const Ray& ry){
	bool intersecB; //does it intersex this seks?
	double t;
	Vector OmC = ry.O-sf.C;
	double kk = dot(ry.u, OmC);
	double delta = kk*kk-(dot(OmC, OmC)-sf.R);
	if (delta >= 0){
		double t1 = -kk - sqrt(delta); //dot(u, C - O) = dot(u, -OmC) = -kk
		double t2 = -kk + sqrt(delta);
		if (t2 < 0){
			intersecB = false;
		}
		if (t1 >= 0){
			intersecB = true;
			t = t1;
		}
		if (t1 < 0){
			intersecB = true;
			t = t2;
		}
	}else{
		intersecB = false;
	}
	if(intersecB){
		Vector P = ry.O + t*ry.u;
		Vector N = (P-sf.C)*(1.0/norm(P-sf.C));
		return {P, N, t, intersecB};
	}
	return {{}, {}, 1e16, intersecB}; //{} works because we didn't use explicit
}


//Intersection: point, normal, dist

Intersection intersect(const Scene& sn, const Ray& ry){
	double curr_t = 1e16;
	bool intersex = 0;
	Sphere kept_ball;
	Vector bgcolor(0, 0, 0); 
	Intersection keep;


	for(auto& ball: sn.welt){
		Intersection Ic = intersect(ball, ry);
		double t = Ic.dist;
		if (t < curr_t){
			curr_t = t;
			keep = Ic;
			kept_ball = ball;
			intersex = 1;
			Vector surfacecolor = ball.albedo;
	}
	if(intersex){
		for(auto& ly: sn.lights){
			pairwiseyoink(surfacecolor, )
		}
	}else{
		return bgcolor;
	}
}



int main(){
	double fov = pi*60/180;
	Scene zouglou;
	zouglou.welt.push_back({Vector(0, 0, 1000), Vector(0xff, 0, 0xff), 940});
	zouglou.welt.push_back({Vector(0, 1000, 0), Vector(0xff, 0, 0), 940});
	zouglou.welt.push_back({Vector(0, 0,-1000), Vector(0, 0xff, 0), 940});
	zouglou.welt.push_back({Vector(0,-1000, 0), Vector(0, 0, 0xff), 940});
	zouglou.welt.push_back({Vector(0, 0, 0), Vector(0xff, 0xff, 0xff), 940});
	Vector Q(0, 0, 55);
	int w = 600;
	int h = 400; 
	uint8_t img[h][w][3]; //3 for r g b 
	for (int i=0; i<h; i++){
		for (int j=0; j<w; j++){
			Vector u(j+0.5-w/2, (h-i-1)+0.5-h/2, -w/(2*fov));
			Ray ry(normalize(u), Q);
			Intersection Ic = intersect(zouglou, ry);
			if (Ic.hmm){
				img[i][j][0] = 255;
				img[i][j][1] = 255;
				img[i][j][2] = 255;
			}
		}
	}
	stbi_write_png("zizirose.png", w, h, 3, img, w*3);
};