#include<bits/stdc++.h>
#include"voronoi_header.h"
#include"svg_put.cpp"

static std::default_random_engine sheitan ( 10 ) ; // random s e e d = 10
static std::uniform_real_distribution<double> uniform ( 0 , 1 ) ;


Vector randomkaka(){
	double x = uniform(sheitan);
	double y = uniform(sheitan);
	return Vector(x, y);
}

Vector dohalfplaenpls(Vector orig, Vector point2){
	return (1/2)*(point2-orig);
}

bool checkhalfplane( Vector v, Vector checkpunkt){ //checkpunkt: the punkt for which we go "hmmmm"
	return (dot(v, checkpunkt) <= dot(v, v));
}

int main(){
	int numpoints = 10;
	std::vector<Vector> zouglou;
	for (int i = 0; i < numpoints; i++){
		zouglou.push_back(randomkaka());
	}
	std::vector<polygre> polyglongs;
	for (int i = 0; i < numpoints; i++){
		polygre Pet(zouglou[i]);
		for (int j = 0; j < numpoints; j++){
			if (j != i){
				Pet.intersex_with_halfplane((1.0/2.0)*(zouglou[j]-zouglou[i]));
			}
		}
		polyglongs.push_back(Pet);
		std::cout<<"for "<<zouglou[i]<<" the polygon is "<<Pet<<std::endl;
	}
	save_svg(polyglongs, "snoglou.svg","white");

	std::cout<<"nice cock bro"<<std::endl;
	return 0;
}

