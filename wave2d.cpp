#include "wave2d.h"
#include <cmath>
using namespace arma;
using namespace std;

Wave2d::Wave2d(){
	isAtInitial = true;
	L = 1.0;
	c_0 = 1.0;
	tau = L/C;
	t = 0.0;
	N = 200;
	dt = 0.0001;
	h = 1.0/((double) (N-1));
	r = (dt*dt)/(h*h);

	x = zeros<vec>(N);
	y = zeros<vec>(N);
	fill_xy();
}

Wave2d::~Wave2d(){

}


void Wave2d::createCustomInitial(){
	for(size_t i = 0; i<N; i++){
		for(size_t j = 0; j<N; j++){
			u(i,j) = u_initial(x(i), y(i));
		}
	}
}

void Wave2d::print_dt_SI(){
	cout << dt*tau << endl;
}

void Wave2d::print_dx_SI(){
	cout << dx*L << endl;
}


double Wave2d::calc_alpha(size_t i, size_t j){
	double f_value = f(x(i), y(i));
	return r*f_value*f_value;
}

double Wave2d::calc_beta_x(size_t i, size_t j){
	return r*h*f(x(i), y(i))*f_x(x(i),y(j));
}

double Wave2d::calc_beta_y(size_t i, size_t j){
	return r*h*f(x(i), y(i))*f_y(x(i),y(j));
}

void Wave2d::fill_xy(){
	for(size_t i = 0; i < N; i++){
		x(i) = ((double) i) * dx;
	}
	y = x;
}

//PLOTTING
void Wave2d::plot(){
	gplt.xystream(N,x,u);
}

void Wave2d::initialize(){
	t = 0.0;
	isAtInitial = true;
	u = zeros<mat>(N,N);
	u_previous = zeros<mat>(N,N);
	fill_alpha_beta();
}

void Wave2d::iterate(size_t it){
	mat unew = zeros<mat>(N,N);
	if(isAtInitial){
		isAtInitial = false;
		//Do corners
		double a_LL = alpha(1,1);
		double a_LR = alpha(N-1,1);
		double a_UL = alpha(1,N-1);
		double a_UR = alpha(N-1,N-1);
		
		unew(1,1) = 0.5*((a_LL+beta_x(1,1))*u(2,1) + (a_LL+beta_y(1,1))*u(1,2) + (2.0+4.0*a_LL)*u(1,1));
		unew(N-2,1) = 0.5*((a_LR-beta_x(N-2,1))*u(N-3,1) + (a_LR+beta_y(N-2,1))*u(N-2,2) + (2.0+4.0*a_LR)*u(N-2,1));
		unew(1,N-2) = 0.5*((a_UL+beta_x(1,N-2))*u(2,N-2) + (a_UL-beta_y(1,N-2))*u(1,N-3) + (2.0+4.0*a_UL)*u(1,N-2));
		unew(N-2,N-2) = 0.5*((a_UR-beta_x(N-2,N-2))*u(N-3,N-2) + (a_UR-beta_y(N-2,N-2))*u(N-2,N-3) + (2.0+4.0*a_UR)*u(N-2,N-2));

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
				unew(i,j) = 0.5*((a+b_x)*u(i+1,j)+(a-b_x)*u(i-1,j)+(a+b_y)*unew(i,j+1)+(a-b_y)*u(i,j-1)+(2.0+4.0*a)*u(i,j));
			}
		}
		u_previous = u;
		u = unew;
	}else{
		//Do corners
		double a_LL = alpha(1,1);
		double a_LR = alpha(N-1,1);
		double a_UL = alpha(1,N-1);
		double a_UR = alpha(N-1,N-1);
		
		unew(1,1) = (a_LL+beta_x(1,1))*u(2,1) + (a_LL+beta_y(1,1))*u(1,2) + (2.0+4.0*a_LL)*u(1,1);
		unew(N-2,1) = (a_LR-beta_x(N-2,1))*u(N-3,1) + (a_LR+beta_y(N-2,1))*u(N-2,2) + (2.0+4.0*a_LR)*u(N-2,1);
		unew(1,N-2) = (a_UL+beta_x(1,N-2))*u(2,N-2) + (a_UL-beta_y(1,N-2))*u(1,N-3) + (2.0+4.0*a_UL)*u(1,N-2);
		unew(N-2,N-2) = (a_UR-beta_x(N-2,N-2))*u(N-3,N-2) + (a_UR-beta_y(N-2,N-2))*u(N-2,N-3) + (2.0+4.0*a_UR)*u(N-2,N-2);

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
				unew(i,j) = (a+b_x)*u(i+1,j)+(a-b_x)*u(i-1,j)+(a+b_y)*unew(i,j+1)+(a-b_y)*u(i,j-1)+(2.0+4.0*a)*u(i,j);
			}
		}
		unew = unew - u_previous;
		u_previous = u;
		u = unew;
	}
}

