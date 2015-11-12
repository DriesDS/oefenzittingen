#include "cg.hpp"
#include "parse.hpp"

#include <algorithm>
#include <functional>
#include <iostream>
#include <iomanip>
#include <limits>
#include <vector>
#include <cmath>

void matvec1( std::vector<double> const& x, std::vector<double>& y ) {
	//int n = x.size();
	for (std::size_t i=0; i<x.size(); ++i) {
		y[i] = x[i]*(i+1);
	}
}

struct matvec2 {
	double m;
	matvec2(const double m = 2) : m(m) { }
	
	template <typename X, typename Y>
	void operator() ( X const& x, Y& y ) const {
		for (int i=0; i<x.size(); ++i) {
			y[i] = x[i]/m ;
		}
	}
};

auto matvec3 = []( std::vector<double> const& x, std::vector<double>& y ) {
	for (std::size_t i=0; i<x.size(); ++i) {
		y[i] = x[i]*(i+1);	
	}	
} ;

template <class Vector>
void matvec4( Vector const& x, Vector& y ) {
	for (std::size_t i=0; i<x.size(); ++i) {
		y[i] = x[i]*(i+1);
	}
}

template <class Vector>
void matvec( Vector const& x, Vector& y ) {
	int n = x.size();
	if (n > 1) {
		y[0] = (-2*x[0]+x[1])*pow(n+1,2);
		y[n-1] = (-2*x[n-1] + x[n-2])*pow(n+1,2);
		for (std::size_t i=1; i<x.size()-1; ++i) {
			y[i] = (x[i-1] -2*x[i] + x[i+1])*pow(n+1,2);
		}
	} else {
		// scalar version of matrix A
		y[0] = -2*x[0]*pow(n+1,2);
	}
}

int main(int argc, char* argv[]) {
	std::size_t n = 128;
	long double ild;
	if(argc > 1) parse(argv[1],n);
	
	std::vector<float> xf(n-1) ;
	std::vector<double> xd(n-1) ;
	std::vector<long double> xld(n-1) ;
	std::vector<long double> bld(n-1) ;
	std::vector<long double> s(n-1) ;
	std::vector<long double> errf(n-1);
	std::vector<long double> errd(n-1);
	std::vector<long double> errld(n-1);

	std::cout << std::setprecision(std::numeric_limits<long double>::digits10+1) 
		<< std::scientific;

	//initialize x and b appropriately
	for (std::size_t i=0; i<xf.size(); ++i) {
		ild = static_cast<long double>(i+1);
		xld[i] = 0;
		xf[i] = 0;
		xd[i] = 0;
		bld[i] = -(3*ild/n+pow(ild/n,2))*exp(ild/n);
		s[i] = (ild/n)*exp(ild/n)-pow(ild/n,2)*exp(ild/n);
	}

	tws::cg(matvec<std::vector<float>>,xf,bld);
	tws::cg(matvec<std::vector<double>>,xd,bld);
	tws::cg(matvec<std::vector<long double>>,xld,bld);
	
	std::transform(xf.begin(),xf.end(),s.begin(),errf.begin(),std::minus<long double>());
	std::transform(xd.begin(),xd.end(),s.begin(),errd.begin(),std::minus<long double>());
	std::transform(xld.begin(),xld.end(),s.begin(),errld.begin(),std::minus<long double>());
	for(std::size_t i=0; i<errf.size(); ++i) {
		if(errf[i] < 0) errf[i] *= -1;
		if(errd[i] < 0) errd[i] *= -1;
	        if(errld[i] < 0) errld[i] *= -1;
	}

	long double maxerrf = *max_element(errf.begin(),errf.end());
	long double maxerrd = *max_element(errd.begin(),errd.end());
	long double maxerrld = *max_element(errld.begin(),errld.end());
	std::cout << n << "\t" << maxerrf << "\t" << maxerrd << "\t" << maxerrld << std::endl;

}
