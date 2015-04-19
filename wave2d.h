#ifndef WAVE2D_H
#define WAVE2D_H
#include <armadillo>
#include "gnuplotting.h"

using namespace arma;
using namespace std;

class Wave2d{
	public:
		Wave2d();
		~Wave2d();
		virtual double f(double & x, double & y) = 0;
		virtual double f_x(double & x, double & y) = 0;
		virtual double f_y(double & x, double & y) = 0;
		virtual double u_initial(double & x, double & y) = 0;

		void createCustomInitial();
		void createDefaultInitial();

		void print_alpha();
		void print_beta_x();
		void print_beta_y();
		void print_x();
		void print_y();
		void print_u();
		void print_u_previous();

		Gnuplotting gplt;

		void plot_2d();
		void plot_3d();

		void initialize();    //Initializes matrices and vectors 
		void iterate(size_t it);       //Iterates the system through time
		void ezIterate_single();
		void ezIterate(size_t it);
		void iterate_single();

		void iterate_test_single();
		void iterate_test(size_t it);

	private:
		void fill_xy();	
		void fill_alpha_beta();
		double calc_alpha(size_t i, size_t j);
		double calc_beta_x(size_t i, size_t j);
		double calc_beta_y(size_t i, size_t j);

		//SYSTEM PARAMETERS
		double L;   //size of system
		double c_0; //wave velocity
		double tau; //characteristic time scale, tau=L/c_0


		bool isAtInitial;

		double t;   //current time of system
		size_t N;   //Number of spatial points in [0,1]
		double h;   //x=ih, i=0,....N-1, y=jh, j=0;....N-1
		double dt;  //timestep
		double r;   //r=(dt/h)^2

		mat alpha;   //arma::matrix containing alpha(x,y) values
		mat beta_x;  //arma::matrix containing beta_y(x,y) values
		mat beta_y;  //arma::matrix containing beta_x(x,y) values

		vec x;  //arma::vector containing all x positions
		vec y;  //arma::vector containing all y positions

		mat u;  //arma::matrix containing wave height u(x,y,t)
		mat u_previous; //arma::matrix containing wave height u(x,y,t-dt)

};



class ConstantWave: public Wave2d{
	public:
		ConstantWave();
		~ConstantWave();
		virtual double f(double & x, double & y);
		virtual double f_x(double & x, double & y);
		virtual double f_y(double & x, double & y);
		virtual double u_initial(double & x, double & y);
};

#endif
