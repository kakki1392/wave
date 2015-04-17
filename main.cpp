#include "wave2d.h"
#include "gnuplotting.h"
#include <unistd.h>

using namespace arma;
using namespace std;

int  main(){
	ConstantWave test;
	test.initialize();
	test.gplt.cmd("reset");
	test.gplt.cmd("set pm3d map");
	test.gplt.cmd("set cbrange [-1:1]");
	test.gplt.cmd("set palette gray");
	for(int i=0; i<1000; i++){
		test.ezIterate(1);
		test.plot();
		usleep(1e5);
	}

	/*
	test.print_x();
	test.print_y();
	test.print_alpha();
	test.print_beta_x();
	test.print_beta_y();
	*/















	return 0;
}
