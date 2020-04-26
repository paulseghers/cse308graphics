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

std::ostream& operator<<(std::ostream& os, const Vector& V){
	os<<"("<<V.x<<", "<<V.y<<", "<<V.z<<")";
	return os;
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
	enum type {solid, reflec, trans} typ;
	double refractive_indeks;
	Sphere(){}
	Sphere(const Vector& c, const Vector& albedo, double R, type x, double n=1):
		C(c), albedo(albedo), R(R), typ(x), refractive_indeks(n)
	{}
};

struct Ray{
	Vector u;
	Vector O;
	Ray(){}
	Ray(const Vector& u, const Vector& O):
		u(u),
		O(O)
		{}
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
};

struct Camera{
	Vector Q;
	double fov;
	Camera(){}
	Camera(const Vector& q, double fov):
		Q(q), fov(fov)
	{}
};

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
	double pepsilon = 1e-4;
	bool intersecB; //does it intersex this seks?
	double t;
	Vector OmC = ry.O-sf.C;
	double kk = dot(ry.u, OmC);
	double delta = kk*kk-(dot(OmC, OmC)-sf.R*sf.R);
	bool flipneeded = false; //to make sure the normal is always pointing away from the sphere center
	if (delta >= 0){
		double t1 = -kk - sqrt(delta); //dot(u, C - O) = dot(u, -OmC) = -kk
		double t2 = -kk + sqrt(delta);
		if (t2 < pepsilon){
			intersecB = false;
		}
		if (t1 >= pepsilon){
			intersecB = true;
			t = t1;
		}
		if (t1 < pepsilon){
			if(t2 >= pepsilon){
				intersecB = true;
				t = t2;
				flipneeded = true;

			}
			else{
				intersecB = false;
			}
		}

	}else{
		intersecB = false;
	}
	if(intersecB){
		Vector P = ry.O + t*ry.u;
		Vector N = (P-sf.C)*(1.0/norm(P-sf.C));
		if(flipneeded){
			N = -1*N;
		}
		return {P, N, t, intersecB};
	}
	return {{}, {}, 1e16, intersecB}; //{} works because we didn't use explicit
}


//Intersection: point, normal, dist
Vector refLecc_this_pls(const Vector& incoming, const Vector& normal){
	return incoming - 2*dot(incoming, normal)*normal;
}

Vector refRacc_this_pls(const Vector& incoming, const Vector& normal, double n1surn2){
	double dotprod = dot(incoming, normal);
	Vector wTang = n1surn2 * (incoming - dotprod*normal);
	double chungusinsidethesqrt = 1 - n1surn2*n1surn2*(1-dotprod*dotprod);
	if (chungusinsidethesqrt<0){return refLecc_this_pls(incoming, normal);}
	std::cout<<"got in/out ;) \n";
	Vector wNorm = -1*normal*sqrt(chungusinsidethesqrt);
	return wTang+wNorm;
}


Vector intersect(const Scene& sn, const Ray& ry, int depth){
	if (depth<5){std::cout<<"ry: ("<<ry.u<<", "<<ry.O<<")"<<" dpth: "<<depth<<std::endl;}
	if (depth<0){return {0, 0, 0};} // check its ok
	double curr_t = 1e16;
	bool intersex = 0;
	Sphere kept_ball;
	Vector bgcolor(0, 0, 0);
	Vector ro = bgcolor; 
	Vector perceived=bgcolor;
	Intersection keep;
	double t;
	for(auto& ball: sn.welt){
		Intersection Ic = intersect(ball, ry);
		t = Ic.dist;
		if (t < curr_t){
			curr_t = t;
			keep = Ic;
			kept_ball = ball;
			//std::cout<<"kept_ball: center="<<ball.C<<" albedo="<<ball.albedo<<std::endl;
			intersex = 1;
			ro = ball.albedo;
			//std::cout<<ball.albedo<<std::endl;
		}
	}
	if(intersex){
		if(kept_ball.typ == Sphere::solid){
			for(auto& ly: sn.lights){
				double d = norm(ly.spot - keep.point);
				double dotprod_kanker = dot(keep.normal, (ly.spot - keep.point)*(1/d));
				double chungus1 = ly.I/(4*pi*d*d);
				Vector chungus2 = ro*(1/pi);
				//std::cout<<"dotprod= "<<dotprod_kanker<<"chungus1= "<<chungus1<<"chungus2= "<<chungus2<<std::endl;
				perceived += ly.I/(4*pi*d*d) * ro*(1/pi) * dotprod_kanker; //deze variable heest kanker dus het kanker is
			}
			return perceived;
			//return Vector(255, 255, 255);
		}
		else if(kept_ball.typ == Sphere::reflec){
			//std::cout<<"I GOT SO FAR AND TRIED SO HARD"<<std::endl;
			//return {0, 0, 0};
			return intersect(sn, Ray(refLecc_this_pls(ry.u, keep.normal), ry.O+t*ry.u), depth-1);
		}
		else if(kept_ball.typ == Sphere::trans){
			//keep.P
			double d1 = dot(keep.point, ry.u);
			double d2 = dot(kept_ball.C, ry.u);
			double n1surn2;
			if (d1 < d2){ //we're going in, boys
				n1surn2 = 1/kept_ball.refractive_indeks;
			}else{
				n1surn2 = kept_ball.refractive_indeks;
			}
			return intersect(sn, Ray(refRacc_this_pls(ry.u, keep.normal, n1surn2), ry.O+t*ry.u), depth-1);
		}
	}else{
		return bgcolor;
	}
}




int clamp(double c){
	if(c<0){
		return 0;
	}
	if(c>255){
		return 255;
	}
	return c;
}

int main(){
	double fov = pi*60/180;
	double intensity = 1e5;
	Scene zouglou;
	zouglou.welt.push_back({Vector(0, 0, 1000), Vector(0xff, 0, 0xff), 940, Sphere::solid});
	zouglou.welt.push_back({Vector(0, 1000, 0), Vector(0xff, 0, 0), 940, Sphere::solid});
	zouglou.welt.push_back({Vector(0, 0,-1000), Vector(0, 0xff, 0), 940, Sphere::solid});
	zouglou.welt.push_back({Vector(0,-1000, 0), Vector(0, 0, 0xff), 990, Sphere::solid});
	zouglou.welt.push_back({Vector(0, 0, 1), Vector(0xff, 0xff, 0xff), 9, Sphere::trans, 1.0/1.5});
	zouglou.welt.push_back({Vector(0, 0, 1), Vector(0xff, 0xff, 0xff), 10, Sphere::trans, 1.5});
	zouglou.lights.push_back({Vector(-10, 20, 40), intensity, Vector(0xff, 0xff, 0xff)});
	Vector Q(0, 0, 55);
	//std::cout<<zouglou.welt[0].C;
	int w = 600;
	int h = 400; 
	uint8_t img[h][w][3]; //3 for r g b 
	for (int i=0; i<h; i++){
		for (int j=0; j<w; j++){
			Vector u(j+0.5-w/2, (h-i-1)+0.5-h/2, -w/(2*fov));
			Ray ry(normalize(u), Q);
			Vector seen_color = intersect(zouglou, ry, 5);
			if (seen_color.x or seen_color.y or seen_color.z){
				//std::cout<<seen_color.x;
				img[i][j][0] = clamp(seen_color.x);
				img[i][j][1] = clamp(seen_color.y);
				img[i][j][2] = clamp(seen_color.z);
			}
		}
	}
	stbi_write_png("zizirose.png", w, h, 3, img, w*3);
};