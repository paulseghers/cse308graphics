struct Vector{
double x, y;
	Vector(){}
	Vector(double x, double y):
		x(x),
		y(y)
{}
	Vector& operator+=(const Vector& b) {
		x += b.x;
		y += b.y;
		return * this ;
	}
};

Vector operator+(const Vector& a, const Vector& b) {
	return Vector(a.x+b.x, a.y+b.y);
}

Vector operator-(const Vector& a, const Vector& b) {
	return Vector(a.x-b.x, a.y-b.y);
}

Vector operator*(double t, const Vector& a) {
	return Vector(a.x*t, a.y*t);
}

Vector operator*(const Vector& a, double t) { //I want it to go both ways because sometimes I like to switch things around to add some spice. also I'm autistic :())
	return Vector(a.x*t, a.y*t);
}

std::ostream& operator<<(std::ostream& os, const Vector& V){
	os<<"("<<V.x<<", "<<V.y<<")";
	return os;
}

double dot(const Vector& a, const Vector& b) {
	return double (a.x*b.x + a.y*b.y);
}

Vector pairwiseyoink(const Vector& a, const Vector& b) {
	return Vector (a.x*b.x, a.y*b.y);
}

double crossprod(const Vector& a, const Vector& b){
	return a.x*b.y - a.y*b.x;
}

Vector orthogonal(const Vector& v){
	return Vector(-v.y, v.x);
}

Vector intersexual(Vector& a, Vector& b, Vector& h, double th){
    double va = dot(a,h),
           vb = dot(b,h);
    double wa = (th-vb)/(va-vb),
           wb = 1.0-wa;
    return a*wa + b*wb;
}

double norm(Vector v){
	return sqrt(dot(v, v));
}

struct polygre{
	std::vector<Vector> corners;
	Vector centerpunkt;
	polygre(Vector p){
		corners = std::vector<Vector>{ Vector{0, 0}, {1, 0}, {1, 1}, {0, 1}};
		centerpunkt = p;
		for_each(corners.begin(), corners.end(), [&p](Vector& v){v = v-p;});
	} /*										          	   |   */
	void intersex_with_halfplane(Vector M, Vector h){ // --h-> M plane line
		int n = corners.size();					          /*   |   */
		double th = dot(M, h);
		auto indapussy = [&](Vector & v){return dot(v, h) < th; };
		bool zoukneeded = false, veryzouky = false;
		for (auto& c : corners){
			if (dot(c, h) > th) zoukneeded = true;
            else veryzouky = true;
		}
		if (!zoukneeded) return;
        if (!veryzouky) { corners = {}; return; }

		// std::cout<<"we did an intersex"<<std::endl;
		int i, j;
		for(i = 0; indapussy(corners[i]); i++); //i is now the first guy out of da pusi
		for(j = i; !indapussy(corners[j%n]); j++); //first guy back in da pusi :D
		for(i = j; indapussy(corners[i%n]); i++);  //i is now the first guy out of da pusi but FOR REAL THIS TIME >:)		

		Vector jntersext = intersexual(corners[(j-1)%n], corners[j%n], h, th);
		Vector intersext = intersexual(corners[(i-1)%n], corners[i%n], h, th);//we modulo n every. fucking. where. just to be sure :') no alarms and no surprises
		std::vector<Vector> newcorns;
		newcorns.push_back(intersext); 
		newcorns.push_back(jntersext);
		for (int k = j; k < i; k++){
			newcorns.push_back(corners[k%n]); //don't forget the mod n (this was a surprise)
		}
		std::swap(corners, newcorns);
	}
    double territory() {
        double ter = 0;
        for (int i=1; i<(int)corners.size()-1; i++) {
            ter += crossprod(corners[i]-corners[0], corners[i+1]-corners[0])/2;
        }
        return ter;
    }
    double secondmom() {
        int N = corners.size();
        double mom = 0;
        for (int i=0; i<N; i++) {
            Vector a = corners[i], b = corners[(i+1) % N];
            mom += crossprod(a,b) * ( dot(a,a)+dot(a,b)+dot(b,b)) / 12;
        }
        return mom;
    }
};

std::ostream& operator<<(std::ostream& os, const polygre& P){
	os<<"(";
	for (auto & p : P.corners){
		os<<p<<", ";
	}
	os<<")";
	return os;
}
