#include<bits/stdc++.h>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include"stb_image_write.h"
#include<omp.h>
//#pragma once
double pi = acos(-1);

/* ##### FUN FUNCTIONS */
#include"types_and_helperfunctions.h"
#include"spicy_lighting_stuff.h"
#include"haslers_cat_eater.h"
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

Intersection intersect(const mesh& m, const Ray& ry, Vector& color){
	Intersection Ic;
	for (auto& tri : m.triangs){
		Vector a = m.vertices[tri.vtx[0]];
		Vector b = m.vertices[tri.vtx[1]];
		Vector c = m.vertices[tri.vtx[2]];
		Vector e1 = b-a;
		Vector e2 = c-a;
		Vector  N = crossprod(e1, e2);
		Vector AO = a-ry.O; 
		double Beta = dot(e2, crossprod(AO, ry.u))/dot(ry.u, N);
		double gamma= -dot(e1, crossprod(AO, ry.u))/dot(ry.u, N);
		double alpha= 1-Beta-gamma;
		double t = dot(AO, N)/dot(ry.u, N);
		if(0<=alpha && alpha<=1 && 0<=Beta && Beta<=1 && 0<=gamma && gamma<=1 && t < Ic.dist){
			Ic.dist = t;
			Ic.hmm  = true;
			//Ic.normal= normalize(N);
			Ic.point = ry.O+ry.u*t;
			Vector a = m.norm[tri.norm[0]];
			Vector b = m.norm[tri.norm[1]];
			Vector c = m.norm[tri.norm[2]];
			Ic.normal = normalize(alpha*a+Beta*b+gamma*c);

			Vector aA = m.uv[tri.uv[0]];
			Vector bB = m.uv[tri.uv[1]];
			Vector cC = m.uv[tri.uv[2]];
			Vector pos = alpha*aA+Beta*bB+gamma*cC;
			double x = pos.x, y = pos.y;
			x -= std::floor(x);
			y -= std::floor(y);
			int bigX = x*m.width;
			int bigY = (1-y)*m.height;
			color.x = m.image[3*(bigY*m.width+bigX)+0];
			color.y = m.image[3*(bigY*m.width+bigX)+1];
			color.z = m.image[3*(bigY*m.width+bigX)+2];
			color = color*(1.0/255);

			/*std::cout<<"punkt: "<<Ic.point<<std::endl;
			std::cout<<"a: "<<a<<std::endl;
			std::cout<<"b: "<<b<<std::endl;
			std::cout<<"c: "<<c<<std::endl;
			std::cout<<"ry: "<<ry.u<<", "<<ry.O<<std::endl;
			std::cout<<"alpha: "<<alpha<<std::endl;
			std::cout<<"beta: "<<Beta<<std::endl;
			std::cout<<"gamma: "<<gamma<<std::endl<<std::endl;*/
		}

	}
return Ic;
}

//Intersection: point, normal, dist
Vector refLecc_this_pls(const Vector& incoming, const Vector& normal){
	Vector kungus = incoming - 2*dot(incoming, normal)*normal;
	//std::cout<<"reFLecc gives: "<<kungus<<std::endl;
	return kungus;
}

Vector refRacc_this_pls(const Vector& incoming, const Vector& normal, double n1surn2){
	double dotprod = dot(incoming, normal);
	Vector wTang = n1surn2 * (incoming - dotprod*normal);
	double chungusinsidethesqrt = 1 - n1surn2*n1surn2*(1-dotprod*dotprod);
	if (chungusinsidethesqrt<0){return refLecc_this_pls(incoming, normal);}
	//std::cout<<"got in/out ;) \n";
	Vector wNorm = -1*normal*sqrt(chungusinsidethesqrt);
	//std::cout<<"refRacc gives: "<<wTang+wNorm<<std::endl;
	return wTang+wNorm;
}


Vector intersect(const Scene& sn, const Ray& ry, int depth){
	//if (depth<5){std::cout<<"ry: ("<<ry.u<<", "<<ry.O<<")"<<" dpth: "<<depth<<std::endl;}
	//std::cout<<"the ry (u, O): "<<ry.u<<", "<<ry.O<<std::endl;
	if (norm(ry.O) > 1e10){
		exit(0);
	}
	if (depth<0){return {0, 0, 0};} // check its ok
	
// std::cout << ry.u;
	double curr_t = 1e16;
	bool intersex = 0;
	Sphere kept_ball;
	Vector bgcolor(0, 0, 0);
	Vector ro = bgcolor; 
	Vector perceived=bgcolor;
	Intersection keep;
	for (auto& m : sn.meshWelt){
		//Vector color = {0, 0, 0};
		Intersection Ic = intersect(m, ry, ro);
			if(Ic.hmm){
				//return color*1;
				//return Vector(1, 1, 1)*1e32;
				for(auto& ly: sn.lights){
					double d = norm(ly.spot - Ic.point);
					double dotprod_kanker = dot(Ic.normal, (ly.spot - Ic.point)*(1/d));
					double chungus1 = ly.I/(4*pi*d*d);
					Vector chungus2 = ro*(1/pi);
					//std::cout<<"dotprod= "<<dotprod_kanker<<"chungus1= "<<chungus1<<"chungus2= "<<chungus2<<std::endl;
					Vector checking_vec = Ic.point-ly.spot; //we need a vec to do a check
					
					Vector checkvec_normalized = normalize(checking_vec);
					Vector new_origin =  ly.spot+eel_slap(checkvec_normalized, random_cos())*ly.radius;
					Ray checking_ry = Ray(normalize(Ic.point-new_origin), new_origin);
					double george = norm(Ic.point-new_origin);
					bool theressmthindawae = 0;
					for(auto& sphe: sn.welt){
						Intersection checking_intersec = intersect(sphe, checking_ry);
						double d = checking_intersec.dist;
						//std::cout<<"big G: "<<george<<std::endl;
						//std::cout<<"big d: "<<d<<std::endl;
						if (d < george-0.0001){
							theressmthindawae = 1;
						}
					}
					if (!theressmthindawae)
					for (auto& ding: sn.meshWelt){
						Vector c = Vector(360, 420, 69);
						Intersection checking_intersec = intersect(ding, checking_ry, c);
						double d = checking_intersec.dist;
						if (d < george-0.0001){
							theressmthindawae = 1;
						}

					}
					if (!theressmthindawae)
						perceived += ly.I/(4*pi*d*d) * ro*(1/pi) * dotprod_kanker;
					//perceived += ly.I/(4*pi*d*d) * ro*(1/pi) * dotprod_kanker; //deze variable heest kanker dus het kanker is
				}
				return perceived;
			}
			


	}
	for(auto& ball: sn.welt){
		Intersection Ic = intersect(ball, ry);
		double t = Ic.dist;
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
				Vector checking_vec = keep.point-ly.spot; //we need a vec to do a check
				
				Vector checkvec_normalized = normalize(checking_vec);
				Vector new_origin =  ly.spot+eel_slap(checkvec_normalized, random_cos())*ly.radius;
				Ray checking_ry = Ray(normalize(keep.point-new_origin), new_origin);
				double george = norm(keep.point-new_origin);
				bool theressmthindawae = 0;
				for(auto& sphe: sn.welt){
					Intersection checking_intersec = intersect(sphe, checking_ry);
					double d = checking_intersec.dist;
					//std::cout<<"big G: "<<george<<std::endl;
					//std::cout<<"big d: "<<d<<std::endl;
					if (d < george-0.0001){
						theressmthindawae = 1;
					}
				}
				if (!theressmthindawae)
				for (auto& ding: sn.meshWelt){
					Vector c = Vector(360, 420, 69);
					Intersection checking_intersec = intersect(ding, checking_ry, c);
					double d = checking_intersec.dist;
					if (d < george-0.0001){
						theressmthindawae = 1;
					}

				}
				if (!theressmthindawae)
				perceived += ly.I/(4*pi*d*d) * ro*(1/pi) * dotprod_kanker; //deze variable heest kanker dus het kanker is
			}
			//zici
			Vector sampled = random_cos();
			//assert(dot(keep.normal, ry.u) < 0);
			Ray diffuse_ry = Ray(eel_slap(keep.normal, sampled), keep.point);
			return perceived+pairwiseyoink(intersect(sn, diffuse_ry, depth-1),kept_ball.albedo);
			//return Vector(255, 255, 255);
		}
		else if(kept_ball.typ == Sphere::reflec){
			//std::cout<<"we do reflec"<<std::endl;
			//return {0, 0, 0};
			return intersect(sn, Ray(refLecc_this_pls(ry.u, keep.normal), ry.O+curr_t*ry.u), depth-1);
		}
		else if(kept_ball.typ == Sphere::trans){
			//std::cout<<"we do trains"<<std::endl;
			//keep.P
			double d1 = dot(keep.point, ry.u);
			double d2 = dot(kept_ball.C, ry.u);
			double n1surn2;
			if (d1 < d2){ //we're going in, boys
				n1surn2 = 1/kept_ball.refractive_indeks;
			}else{
				n1surn2 = kept_ball.refractive_indeks;
			}
			//std::cout<<"t: "<<t<<std::endl;
			return intersect(sn, Ray(refRacc_this_pls(ry.u, keep.normal, n1surn2), ry.O+curr_t*ry.u), depth-1);
		}
	}
		//std::cout<<"ray (u, O): "<<ry.u<<", "<<ry.O<<std::endl;
		//std::cout<<"we got here"<<std::endl;
		return bgcolor;
	
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

int gamma_correct(double c){
	if(c<0)return 0;
	if(c>1)return 255;
	return 255*pow(c, 1/2.2);
}

int main(){
	//double fov = pi*60/180;
	double fov = pi*60/180;
	double intensity = 1e5;
	Scene zouglou;
	zouglou.welt.push_back({Vector(0, 0, 1000), Vector(0xff, 0, 0xff)*(1.0/255.0), 940, Sphere::solid}); //if vectors are too big they make light
	zouglou.welt.push_back({Vector(0, 1000, 0), Vector(0xff, 0, 0)*(1.0/255.0), 940, Sphere::solid});
	zouglou.welt.push_back({Vector(0, 0,-1000), Vector(0, 0xff, 0)*(1.0/255.0), 940, Sphere::solid});
	zouglou.welt.push_back({Vector(0,-1000, 0), Vector(0, 0, 0xff)*(1.0/255.0), 990, Sphere::solid});
	zouglou.welt.push_back({Vector(-1000, 0, 0), Vector(0, 0xff, 0xff)*(1.0/255.0), 940, Sphere::solid});
	zouglou.welt.push_back({Vector( 1000, 0, 0), Vector(0xff, 0xff, 0)*(1.0/255.0), 940, Sphere::solid});
	//zouglou.meshWelt.push_back(eat_da_mesh());
	//zouglou.welt.push_back({Vector(25, 0, 0), Vector(0xff, 0xff, 0xff), 10, Sphere::trans, 1.5});
	zouglou.welt.push_back({Vector(0, 0, 0), Vector(0xff, 0xff, 0xff)*(1.0/255.0), 10, Sphere::solid});
	zouglou.welt.push_back({Vector(25, 0, 0), Vector(0xff, 0xff, 0xff)*(1.0/255.0), 10, Sphere::reflec});
	zouglou.welt.push_back({Vector(-20, 1, 0), Vector(0xff, 0xff, 0xff), 9, Sphere::trans, 1.0/1.5});
	zouglou.welt.push_back({Vector(-20, 1, 0), Vector(0xff, 0xff, 0xff), 10, Sphere::trans, 1.5});

	zouglou.lights.push_back({Vector(-20, 20, 25), intensity, Vector(0xff, 0xff, 0xff)*(1.0/255.0), 10});
	
/*	for (int v = 6; v<10; v++){
		Sphere s = zouglou.welt[v];
		std::cout<<"Sphere: "<<s.C<<", "<<s.albedo<<", "<<s.R<<", "<<s.typ<<", "<<s.refractive_indeks<<std::endl;
	}*/

	Vector Q(0, 0, 55);
	//std::cout<<zouglou.welt[0].C;
	int w = 1920;
	int h = 1080; 
	uint8_t img[h][w][3]; //3 for r g b 
#pragma omp parallel for
	for (int i=0; i<h; i++){
		for (int j=0; j<w; j++){
/*			i = h/2.0;
			j = 3*w/4.0;*/
			Vector u(j+0.5-w/2, (h-i-1)+0.5-h/2, -w/(2*fov));
			Vector seen_color = {0, 0, 0};
			for (int k = 0; k <24; k++){
				double dx, dy;
				BoxMuller(0.9, dx, dy);
				Vector pepsilon2d = {dx, dy, 0};
				Ray ry(normalize(u+pepsilon2d), Q);
				seen_color += intersect(zouglou, ry, 5);
			}
			seen_color = seen_color*(1.0/24);
			if (seen_color.x or seen_color.y or seen_color.z){
				//std::cout<<seen_color.x;
				img[i][j][0] = gamma_correct(seen_color.x);
				img[i][j][1] = gamma_correct(seen_color.y);
				img[i][j][2] = gamma_correct(seen_color.z);
			}
		}
	}
	stbi_write_png("zizirose.png", w, h, 3, img, w*3);
};