#include<bits/stdc++.h>

static std::default_random_engine sheitan ( 10 ) ; // random s e e d = 10
static std::uniform_real_distribution<double> uniform ( 0 , 1 ) ;
Vector random_cos(){
	double r1 = uniform(sheitan);
	double r2 = uniform(sheitan);
	double x = cos(2*pi*r1)*sqrt(1-r2);
	double y = sin(2*pi*r1)*sqrt(1-r2);
	double z = sqrt(r2);
	return {x, y, z};
}

Vector eel_slap(const Vector& Normal, const Vector& sampled){
	double x = Normal.x;
	double y = Normal.y;
	double z = Normal.z; //c'est tous des doubles HEIN
	Vector T2;
	if ((std::abs(x) <= std::abs(y)) && (std::abs(x) <= std::abs(z))){T2 = {0, z, -y}; }
	if ((std::abs(y) <= std::abs(x)) && (std::abs(y) <= std::abs(z))){T2 = {z, 0, -x}; }
	if ((std::abs(z) <= std::abs(x)) && (std::abs(z) <= std::abs(y))){T2 = {y, -x, 0}; }
	//std::cout << Normal << T2 << std::endl;
	T2 = normalize(normalize(T2));
	Vector T1 = crossprod(T2,Normal);
	Vector res = sampled.z*Normal+sampled.x*T1+sampled.y*T2;
	/* std::cout<<"in:"<<std::endl;
	std::cout<<"we were fed random: "<<sampled<<std::endl;
	std::cout<<"T1: "<<T1<<"and T2: "<<T2<<std::endl; 
	std::cout<<"Normal: "<<Normal<<", res: "<<res<<std::endl; */
	return res;
} 

void BoxMuller(double stdev, double &x , double &y){ //Müller-Milch ...mmmmh!
	double r1 = uniform (sheitan);
	double r2 = uniform (sheitan);
	x = sqrt (-2 * log ( r1 ) ) * cos ( 2 * pi * r2 ) * stdev;
	y = sqrt (-2 * log ( r1 ) ) * sin ( 2 * pi * r2 ) * stdev;
}

