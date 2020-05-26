#include<bits/stdc++.h>
#include"types_and_helperfunctions.h"
#define STB_IMAGE_IMPLEMENTATION
#include"stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include"stb_image_write.h"

int width, height;
static std::default_random_engine sheitan (std::time(0)) ; // random s e e d = 10
static std::uniform_real_distribution<double> uniform ( 0 , 1 ) ;

Vector vectors_on_unitshpere(){
	double r1 = uniform (sheitan);
	double r2 = uniform (sheitan);
	double x = cos(2*pi*r1)*sqrt(r2*(1-r2));
	double y = sin(2*pi*r1)*sqrt(r2*(1-r2));
	double z = 1- 2*r2;
	return Vector(x, y, z);
}

std::vector<Vector> reeed(char* fname){//een std vektor met mijne vektoren binnen
	int ncomps;
	uint8_t* image = stbi_load(fname, &width, &height, &ncomps, 3);
	std::vector <Vector> res(width*height); //my text editor is autistic so it makes my a vector a function :D since I'm drunk this suprised me !
	for (int i = 0 ; i < height ; i++){
		for (int j = 0 ; j <                       width ; j++){
			int k = i*width+j;
			res[k] = Vector(image[3*k+0], image[3*k+1], image[3*k+2]);
		}
	}
	return res;
}

void wriii(char* where_to_poo, const std::vector<Vector>& imvec){
	uint8_t* image = new uint8_t[3*width*height];
	for (int i = 0 ; i < height ; i++){
		for (int j = 0 ; j < width ; j++){
			int k = i*width+j;
			image[3*k+0]=imvec[k].x, 
			image[3*k+1]=imvec[k].y,
			image[3*k+2]=imvec[k].z;
		}
	}
	stbi_write_png(where_to_poo, width, height, 3, image, width*3); //this function is very autistic but that's ok
}

void do_slicy_stuff(const std::vector<Vector>& paint, std::vector<Vector>& face){// SLICED OPTIMAL TRANSPORT MATCHING :DDD
	for (int k = 69; k < 420-390+69; k++){
		Vector zouk = vectors_on_unitshpere();
		int n = width*height;
		std::vector <std::pair<double, int>> facesque(n), paintesque(n);
		for (int i = 0; i < n; i++){
			facesque[i]  ={dot(zouk, face[i]), i};
			paintesque[i]={dot(zouk, paint[i]), i};
		}
		std::sort(facesque.begin()  , facesque.end());//in C++ 20+ this can be done spicily but where not there yet
		std::sort(paintesque.begin(), paintesque.end());

		/*light-colored thigs always arrive at the begginng of the array*/
		for (int j = 0; j < n; j++){
			face[facesque[j].second] += zouk * (paintesque[j].first-facesque[j].first);
		}
	}
}

int main(){
	std::vector<Vector> patrick = reeed("pattern.png");
	std::vector<Vector> to = reeed("target.png");
	do_slicy_stuff(patrick, to);
	wriii("new_image.png", to);
	std::cout<<"done :)"<<std::endl;
	return 0;
}