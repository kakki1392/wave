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
	
	/*
	wall_clock timer;
	timer.tic();
	A.ezIterate(100);
	cout << "ez: " << timer.toc() << endl;
	timer.tic();
	B.iterate_test(100);
	cout << "test: " << timer.toc() << endl;

	A.print_alpha();
	A.print_beta_x();
	A.print_beta_y();
	A.print_u();
	A.print_u_previous();
	B.print_u();
	B.print_u_previous();
	*/

	A.gplt.cmd("set term gif animate delay 5 1");
	A.gplt.cmd("set output 'gauss.gif'");
	A.gplt.cmd("set pm3d map");
	//A.gplt.cmd("set pm3d at s hidden3d 100");
	//A.gplt.cmd("set style line 100 lt 5 lw 0.5");
	//A.gplt.cmd("unset hidden3d");
	//A.gplt.cmd("unset surf");
	A.gplt.cmd("set palette gray");
	A.gplt.cmd("set cbrange [-0.2:0.2]");
	//A.gplt.cmd("set hidden3d");
//	A.gplt.cmd("set zrange [-0.7:0.7]");
	for(int i=0; i<1000; i++){
		A.iterate_test(1);
		A.plot_2d();
	//	usleep(1e4);
	}
	A.gplt.cmd("unset output");






	return 0;
}
