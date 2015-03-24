// test.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "head.h"
#include "fluid.h"
#include "OPENGL_DRIVER.h"
#include <iostream>
using namespace std;


void fun() {
	body.solve();
	cout << body.particles[0].x[0] <<' '<< body.particles[0].x[1]<<' ' <<body.particles[0].x[2]<<endl;
}

int _tmain(int argc, char** argv)
{
	OPENGL_DRIVER::target[0] = 0;
	OPENGL_DRIVER::target[1] = 0;
	OPENGL_DRIVER::target[2] = 0;
	OPENGL_DRIVER(&argc, argv);

/*	for(int i = 0; i < 10; i++) {
		fun();
	}
	int a;
	cin>>a;
	*/
	return 0;
}

