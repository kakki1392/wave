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
	fill_XY();
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


double Wave2d::alpha(size_t i, size_t j){
	double f_value = f(x(i), y(i));
	return r*f_value*f_value;
}

double Wave2d::beta_x(size_t i, size_t j){
	return r*h*f(x(i), y(i))*f_x(x(i),y(j));
}

double Wave2d::beta_y(size_t i, size_t j){
	return r*h*f(x(i), y(i))*f_y(x(i),y(j));
}

void Wave2d::createXY(){
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
	createDiagonals();
	createB();
	createDeltaFunction();
}

void Wave2d::iterate(size_t it){
	if(isAtInitial){
		//DO STUFF
		isAtInitial = false;
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
		//Down
		for(size_t kx = 2; k<(N-2); k++){
			//unew(kx,1) =  ;
		}

		for(size_t kx = 2; k<(N-2); k++){
			//unew(kx,1) =  ;
		}

		for(size_t kx = 2; k<(N-2); k++){
			//unew(kx,1) =  ;
		}

		for(size_t kx = 2; k<(N-2); k++){
			//unew(kx,1) =  ;
		}

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


