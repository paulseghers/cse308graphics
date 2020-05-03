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

Vector crossprod(const Vector& a, const Vector& b){
	return Vector(
		a.y*b.z - a.z*b.y,
        a.z*b.x - a.x*b.z,
        a.x*b.y - a.y*b.x);
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
	double dist=1e32;
	bool hmm=0; //default values: we don't intersect, and we dot that very very far awae
	Intersection(){}
	Intersection(const Vector& point, const Vector& normal, double t, bool b):
		point(point), normal(normal), dist(t), hmm(b)
	{}
};

struct Lightsource{
	Vector spot;
	double I;
	Vector hue;
	double radius;
	Lightsource(){}
	Lightsource(const Vector& spot, double i, const Vector& c, double r):
		spot(spot), I(i), hue(c), radius(r)
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

struct triangle_indices{
	int vtx[3], norm[3], uv[3];
};

struct mesh{
	std::vector <triangle_indices> triangs;
	std::vector <Vector> vertices, norm, uv;
	uint8_t* image;
	int height;
	int width;
	//std::vector <Vector> colors;
};

struct Scene{
	std::vector<Sphere> welt;
	std::vector<Lightsource> lights;
	std::vector<mesh> meshWelt; //a world of meshes
};

double norm(const Vector& V){
	return sqrt(dot(V, V));
}

Vector normalize(const Vector& V){
	return V*(1/norm(V));
}

void mesh_scale(mesh& m, double lambda){
	for (auto& vtx : m.vertices){
		vtx = vtx*lambda;
	}
}

void mesh_translate(mesh& m, Vector s){
	for (auto& vtx : m.vertices){
		vtx = vtx+s;
	}
}


