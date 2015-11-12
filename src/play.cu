#include <iostream>
#include <dlib/timing.h>
#include "matrix.h"

int main() {

	Matrix<host<float>> hm1(500,2000);
	Matrix<host<float>> hm2(500,2000);
	dlib::timing::start(1, "copy");	
	Matrix<device<float>> m1 = hm1;
	Matrix<device<float>> m2 = hm2; 
	dlib::timing::stop(1);
	dlib::timing::start(2, "DEV mult");	
	Matrix<device<float>> m3 = m1 * m2;
	dlib::timing::stop(2);
	dlib::timing::start(3, "HOST mult");
	auto hm3 = hm1 * hm2;
	dlib::timing::stop(3);

	std::cout << m3(0, 0) << std::endl;
	std::cout << hm3(0, 0) << std::endl;
	dlib::timing::print();

	return 0;
}
