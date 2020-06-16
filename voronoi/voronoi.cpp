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



std::vector<polygre> powerdiagre(std::vector<Vector> pts, const double* w){
	int numpoints = pts.size();
	std::vector<polygre> polyglongs;
	for (int i = 0; i < numpoints; i++){
		polygre Pet(pts[i]);
		for (int j = 0; j < numpoints; j++){
			if (j != i){
                Vector h = pts[j] - pts[i];
				Vector M = .5 * (1. + (w[i] - w[j])/dot(h,h)) * h;
				Pet.intersex_with_halfplane(M,h);
			}
		} 
		polyglongs.push_back(Pet);
		// std::cout<<"for "<<pts[i]<<" the polygon is "<<Pet<<std::endl;
	}
	return polyglongs;
}

#include <lbfgs.c>

struct lbfgs_kek {
    std::vector<double> lamb;
	std::vector<Vector> zouglou;
};

static lbfgsfloatval_t evaluate(
    void *instance,
    const lbfgsfloatval_t *x,
    lbfgsfloatval_t *g,
    const int n,
    const lbfgsfloatval_t step
    )
{

    std::vector<double> lamb = ((lbfgs_kek*)instance)->lamb;
    std::vector<Vector> zouglou = ((lbfgs_kek*)instance)->zouglou;
    auto polyglongs = powerdiagre(zouglou, x);
    lbfgsfloatval_t fx = 0.0;
    for (int i = 0; i < n; i ++) {
        double T = polyglongs[i].territory(),
               S = polyglongs[i].secondmom();
        g[i] = T - lamb[i];
        fx += x[i] * g[i] - S;
    }
    return fx;
}
static int progress(
    void *instance,
    const lbfgsfloatval_t *x,
    const lbfgsfloatval_t *g,
    const lbfgsfloatval_t fx,
    const lbfgsfloatval_t xnorm,
    const lbfgsfloatval_t gnorm,
    const lbfgsfloatval_t step,
    int n,
    int k,
    int ls
    )
{
    printf("Iteration %d:\n", k);
    printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);
    printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
    printf("\n");
    return 0;
}

int main(){
	int numpoints = 100;
    std::vector<double> lamb(numpoints);
	std::vector<Vector> zouglou(numpoints);
    Vector Panopticon = {0.5,0.5};
    double food = 0;
	for (int i = 0; i < numpoints; i++) {
		zouglou[i]= randomkaka();
		food += lamb[i] = exp(-5*norm(zouglou[i] - Panopticon));
    }
	for (int i = 0; i < numpoints; i++) {
		lamb[i] /= food;
    }

    double fx;
    lbfgs_parameter_t param;
    lbfgs_parameter_init(&param);
    lbfgs_kek args = { lamb, zouglou };
    double x[numpoints];
	for (int i = 0; i < numpoints; i++) {
        x[i] = 0;
    }
    int ret = lbfgs(numpoints, x, &fx, evaluate, progress, &args, &param);

    auto polyglongs = powerdiagre(zouglou, x);
	save_svg(polyglongs, "snoglou.svg","white");

	std::cout<<"nice cock bro"<<std::endl;
	return 0;
}

