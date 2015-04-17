#include "wave2d.h"
#include <armadillo>
#include <iostream>
#include "gnuplotting.h"
#include <unistd.h>

using namespace arma;
using namespace std;

int  main(){
	cout.precision(11);
	cout.setf(ios::fixed);
	ConstantWave A,B;
	A.initialize();
	B.initialize();
	
	wall_clock timer;
	timer.tic();
	A.ezIterate(10000000);
	cout << "ez: " << timer.toc() << endl;
	timer.tic();
	B.iterate_test(10000000);
	cout << "test: " << timer.toc() << endl;

	A.print_alpha();
	A.print_beta_x();
	A.print_beta_y();
	A.print_u();
	A.print_u_previous();
	B.print_u();
	B.print_u_previous();





	return 0;
}
