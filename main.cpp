#include "wave2d.h"
#include "gnuplotting.h"
#include <unistd.h>

using namespace arma;
using namespace std;

int  main(){
	ConstantWave test;
	test.initialize();
	test.plot();
	/*
	test.print_x();
	test.print_y();
	test.print_alpha();
	test.print_beta_x();
	test.print_beta_y();
	*/
	
	for(size_t i=0; i<200; i++){
		test.ezIterate(2);
		test.plot();
		usleep(1e5);
	}















	return 0;
}
