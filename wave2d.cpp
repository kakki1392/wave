#include "wave2d.h"
#include <cmath>
using namespace arma;
using namespace std;

Wave2d::Wave2d(){
	isAtInitial = true;
	L = 1.0;
	c_0 = 1.0;
	tau = L/c_0;
	t = 0.0;
	N = 100;
	//dt = 0.0001;
	h = 1.0/((double) (N-1));
	//r = (dt*dt)/(h*h);
	r = 0.1;
	dt = sqrt(r*h*h);

	x = zeros<vec>(N);
	y = zeros<vec>(N);
	u = zeros<mat>(N,N);
	u_previous = zeros<mat>(N,N);
	alpha = zeros<mat>(N,N);
	beta_x = zeros<mat>(N,N);
	beta_y = zeros<mat>(N,N);
	fill_xy();
	createDefaultInitial();
}

Wave2d::~Wave2d(){

}


void Wave2d::createCustomInitial(){
	for(size_t i = 0; i<N; i++){
		for(size_t j = 0; j<N; j++){
			u(i,j) = u_initial(x(i), y(j));
		}
	}
}

void Wave2d::createDefaultInitial(){
	for(size_t i = 1; i<(N-1); i++){
		for(size_t j = 1; j<(N-1); j++){
			u(i,j) = exp(-1000.0*((x(i)-0.5)*(x(i)-0.5) + (y(j)-0.5)*(y(j)-0.5)));
		}
	}
}


double Wave2d::calc_alpha(size_t i, size_t j){
	double f_value = f(x(i), y(j));
	return r*f_value*f_value;
}

double Wave2d::calc_beta_x(size_t i, size_t j){
	return r*h*f(x(i), y(i))*f_x(x(i),y(j));
}

double Wave2d::calc_beta_y(size_t i, size_t j){
	return r*h*f(x(i), y(j))*f_y(x(i),y(j));
}

void Wave2d::fill_xy(){
	for(size_t i = 0; i < N; i++){
		x(i) = ((double) i) * h;
	}
	y = x;
}

void Wave2d::fill_alpha_beta(){
	for(size_t i = 0; i<N; i++){
		for(size_t j=0; j<N; j++){
			alpha(i,j) = calc_alpha(i,j);
			beta_x(i,j) = calc_beta_x(i,j);
			beta_y(i,j) = calc_beta_y(i,j);
		}
	}
}

//Printing

void Wave2d::print_alpha(){
	alpha.print("alpha");
}

void Wave2d::print_beta_x(){
	beta_x.print("beta_x");
}

void Wave2d::print_beta_y(){
	beta_y.print("beta_y");
}

void Wave2d::print_x(){
	x.print("x");
}

void Wave2d::print_y(){
	y.print("y");
}

void Wave2d::print_u(){
	u.raw_print("u");
}

void Wave2d::print_u_previous(){
	u_previous.raw_print("u_previous");
}

//PLOTTING
void Wave2d::plot_2d(){
	gplt.heatmap(N,N,u);
}

void Wave2d::plot_3d(){
	gplt.xyzstream(N,x,y,u);
}

void Wave2d::initialize(){
	t = 0.0;
	isAtInitial = true;
	fill_alpha_beta();
}

void Wave2d::iterate(size_t it){
	for(size_t i=0; i<it; i++){
		iterate_single();
		t = t + dt;
	}
}

void Wave2d::iterate_single(){
	mat unew = zeros<mat>(N,N);
	if(isAtInitial){
		isAtInitial = false;
		//Do corners
		double a_LL = alpha(1,1);
		double a_LR = alpha(N-2,1);
		double a_UL = alpha(1,N-2);
		double a_UR = alpha(N-2,N-2);
		
		unew(1,1) = 0.5*((a_LL+beta_x(1,1))*u(2,1) + (a_LL+beta_y(1,1))*u(1,2) + (2.0-4.0*a_LL)*u(1,1));
		unew(N-2,1) = 0.5*((a_LR-beta_x(N-2,1))*u(N-3,1) + (a_LR+beta_y(N-2,1))*u(N-2,2) + (2.0-4.0*a_LR)*u(N-2,1));
		unew(1,N-2) = 0.5*((a_UL+beta_x(1,N-2))*u(2,N-2) + (a_UL-beta_y(1,N-2))*u(1,N-3) + (2.0-4.0*a_UL)*u(1,N-2));
		unew(N-2,N-2) = 0.5*((a_UR-beta_x(N-2,N-2))*u(N-3,N-2) + (a_UR-beta_y(N-2,N-2))*u(N-2,N-3) + (2.0-4.0*a_UR)*u(N-2,N-2));

		//Do boundaries
		for(size_t k = 2; k<(N-2); k++){
			double a_k_1 = alpha(k,1);
			double b_k_1x = beta_x(k,1);
			unew(k,1) =  0.5*((a_k_1+b_k_1x)*u(k+1,1) + (a_k_1-b_k_1x)*u(k-1,1) + (a_k_1+beta_y(k,1))*u(k,2) + (2.0-4.0*a_k_1)*u(k,1));
			double a_k_N2 = alpha(k,N-2);
			double b_k_N2x = beta_x(k,N-2);
			unew(k,N-2) = 0.5*((a_k_N2+b_k_N2x)*u(k+1,N-2) + (a_k_N2-b_k_N2x)*u(k-1,N-2) + (a_k_N2-beta_y(k,N-2))*u(k,N-3) + (2.0-4.0*a_k_N2)*u(k,N-2));
			double a_1_k = alpha(1,k);
			double b_1_ky = beta_y(1,k);
			unew(1,k) = 0.5*((a_1_k+beta_x(1,k))*u(2,k) + (a_1_k+b_1_ky)*u(1,k+1) + (a_1_k-b_1_ky)*u(1,k-1) + (2.0-4.0*a_1_k)*u(1,k));
			double a_N2_k = alpha(N-2,k);
			double b_N2_ky = beta_y(N-2,k);
			unew(N-2,k) = 0.5*((a_N2_k-beta_x(N-2,k))*u(N-3,k) + (a_N2_k+b_N2_ky)*u(N-2,k+1) + (a_N2_k-b_N2_ky)*u(N-2,k-1) + (2.0-4.0*a_N2_k)*u(N-2,k));
		}

		//Do central
		for(size_t i = 2; i < (N-2); i++){
			for(size_t j = 2; j < (N-2); j++){
				double a = alpha(i,j);
				double b_x = beta_x(i,j);
				double b_y = beta_x(i,j);
				unew(i,j) = 0.5*((a+b_x)*u(i+1,j)+(a-b_x)*u(i-1,j)+(a+b_y)*unew(i,j+1)+(a-b_y)*u(i,j-1)+(2.0-4.0*a)*u(i,j));
			}
		}
		u_previous = u;
		u = unew;
	}else{
		//Do corners
		double a_LL = alpha(1,1);
		double a_LR = alpha(N-2,1);
		double a_UL = alpha(1,N-2);
		double a_UR = alpha(N-2,N-2);
		
		unew(1,1) = (a_LL+beta_x(1,1))*u(2,1) + (a_LL+beta_y(1,1))*u(1,2) + (2.0-4.0*a_LL)*u(1,1);
		unew(N-2,1) = (a_LR-beta_x(N-2,1))*u(N-3,1) + (a_LR+beta_y(N-2,1))*u(N-2,2) + (2.0-4.0*a_LR)*u(N-2,1);
		unew(1,N-2) = (a_UL+beta_x(1,N-2))*u(2,N-2) + (a_UL-beta_y(1,N-2))*u(1,N-3) + (2.0-4.0*a_UL)*u(1,N-2);
		unew(N-2,N-2) = (a_UR-beta_x(N-2,N-2))*u(N-3,N-2) + (a_UR-beta_y(N-2,N-2))*u(N-2,N-3) + (2.0-4.0*a_UR)*u(N-2,N-2);

		//Do boundaries
		for(size_t k = 2; k<(N-2); k++){
			double a_k_1 = alpha(k,1);
			double b_k_1x = beta_x(k,1);
			unew(k,1) =  (a_k_1+b_k_1x)*u(k+1,1) + (a_k_1-b_k_1x)*u(k-1,1) + (a_k_1+beta_y(k,1))*u(k,2) + (2.0-4.0*a_k_1)*u(k,1);
			double a_k_N2 = alpha(k,N-2);
			double b_k_N2x = beta_x(k,N-2);
			unew(k,N-2) = (a_k_N2+b_k_N2x)*u(k+1,N-2) + (a_k_N2-b_k_N2x)*u(k-1,N-2) + (a_k_N2-beta_y(k,N-2))*u(k,N-3) + (2.0-4.0*a_k_N2)*u(k,N-2);
			double a_1_k = alpha(1,k);
			double b_1_ky = beta_y(1,k);
			unew(1,k) = (a_1_k+beta_x(1,k))*u(2,k) + (a_1_k+b_1_ky)*u(1,k+1) + (a_1_k-b_1_ky)*u(1,k-1) + (2.0-4.0*a_1_k)*u(1,k);
			double a_N2_k = alpha(N-2,k);
			double b_N2_ky = beta_y(N-2,k);
			unew(N-2,k) = (a_N2_k-beta_x(N-2,k))*u(N-3,k) + (a_N2_k+b_N2_ky)*u(N-2,k+1) + (a_N2_k-b_N2_ky)*u(N-2,k-1) + (2.0-4.0*a_N2_k)*u(N-2,k);
		}

		//Do central
		for(size_t i = 2; i < (N-2); i++){
			for(size_t j = 2; j < (N-2); j++){
				double a = alpha(i,j);
				double b_x = beta_x(i,j);
				double b_y = beta_x(i,j);
				unew(i,j) = (a+b_x)*u(i+1,j)+(a-b_x)*u(i-1,j)+(a+b_y)*unew(i,j+1)+(a-b_y)*u(i,j-1)+(2.0-4.0*a)*u(i,j);
			}
		}
		unew = unew - u_previous;
		u_previous = u;
		u = unew;
	}
}

void Wave2d::ezIterate_single(){
	mat unew = zeros<mat>(N,N);
	if(isAtInitial){
		isAtInitial = false;
		for(size_t i=1; i<(N-1); i++){
			for(size_t j=1; j<(N-1); j++){
				unew(i,j) = 0.5*r*(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1));
			}
		}
		unew = unew + (1.0-2.0*r)*u;
		u_previous = u;
		u = unew;
	}else{
		for(size_t i=1; i<(N-1); i++){
			for(size_t j=1; j<(N-1); j++){
				unew(i,j) = r*(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1));
			}
		}
		unew = unew - u_previous + (2.0-4.0*r)*u;
		u_previous = u;
		u = unew;
	}
}

void Wave2d::ezIterate(size_t it){
	for(size_t i=0; i<it; i++){
		ezIterate_single();
		t = t + dt;
	}
}

void Wave2d::iterate_test_single(){
	mat unew = zeros<mat>(N,N);
	if(isAtInitial){
		isAtInitial = false;
		for(size_t i=1; i<(N-1); i++){
			for(size_t j=1; j<(N-1); j++){
				double a = alpha(i,j);
				double bx = beta_x(i,j);
				double by = beta_y(i,j);
				unew(i,j) = 0.5*((a+bx)*u(i+1,j) + (a-bx)*u(i-1,j) + (a+by)*u(i,j+1) + (a-by)*u(i,j-1) - 4.0*a*u(i,j));
			}
		}
		unew = unew + u;
		u_previous = u;
		u = unew;
	}else{
		for(size_t i=1; i<(N-1); i++){
			for(size_t j=1; j<(N-1); j++){
				double a = alpha(i,j);
				double bx = beta_x(i,j);
				double by = beta_y(i,j);
				unew(i,j) = (a+bx)*u(i+1,j) + (a-bx)*u(i-1,j) + (a+by)*u(i,j+1) + (a-by)*u(i,j-1) - 4.0*a*u(i,j);
			}
		}
		unew = unew + 2.0*u - u_previous;
		u_previous = u;
		u = unew;
	}
}

void Wave2d::iterate_test(size_t it){
	for(size_t i=0; i<it; i++){
		iterate_test_single();
		t = t + dt;
	}
}

//SUBCLASSES

ConstantWave::ConstantWave(): Wave2d(){
	createDefaultInitial();
}

ConstantWave::~ConstantWave(){

}

double ConstantWave::f(double & x, double & y){
	return 1.0;
}

double ConstantWave::f_x(double & x, double & y){
	return 0.0;
}

double ConstantWave::f_y(double & x, double & y){
	return 0.0;
}

double ConstantWave::u_initial(double & x, double & y){
	return 1.0;
}
