#ifndef _COLOR_CONVERSIONS_H_
#define _COLOR_CONVERSIONS_H_
#include <algorithm>
#include <cmath>

static const double LAB_NORM_RADIUS = 60; //radius of th esphere in Lab space to normalize colors to.

struct color{
    double r,g,b;
    color(){};
    color( double r_, double g_, double b_ )
	:r(r_), g(g_), b(b_){};
    color operator *( double k )const{
	return color(r*k, g*k,b*k);
    };
    color operator + (const color & c)const{
	return color(r+c.r, g+c.g, b+c.b);
    };
    color operator + (double d )const{
	return color( r+d, g+d, b+d);
    };
    color normalize()const{
	using namespace std;
        double m = max( max( fabs(r-0.5), fabs(g-0.5) ),
			max( fabs(b-0.5), 1e-5) ); // 0 < m < 0.5
        double k = 0.5 / m;
	return (*this + -0.5)*k+0.5;
    };
    static int bounded( int rgb ){
	if ( rgb < 0 ) return 0;
	if ( rgb > 255 ) return 255;
	return rgb;
    }
    int int_r()const{
	return bounded( int(floor(r*255)) );
    };
    int int_g()const{
	return bounded( int(floor(g*255)) );
    };
    int int_b()const{
	return bounded( int(floor(b*255)) );
    };
};


struct lab_color{
    double L,a,b;

    lab_color( double L_, double a_,  double b_ ):
	L(L_), a(a_), b(b_){};
    lab_color(){};

    lab_color operator *( double k )const{
	return lab_color(L*k, a*k, b*k);
    };
    lab_color operator + (const lab_color & c)const{
	return lab_color( L+c.L, a+c.a, b+c.b);
    };
    lab_color normalize()const;
};

lab_color rgb2lab( const color & rgb );
color lab2rgb( const lab_color & lab );

#endif /* _COLOR_CONVERSIONS_H_ */
