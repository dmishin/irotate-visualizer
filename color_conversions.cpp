#include "color_conversions.hpp"
#include <cmath>
using namespace std;

double ref_X =  95.047;
double ref_Y = 100.000;
double ref_Z = 108.883;

double rgb2xyz_f( double rgb ){
    if ( rgb > 0.04045 ) 
	return pow ( ( rgb + 0.055 ) / 1.055, 2.4 );
    else
	return rgb / 12.92;
}

double xyz2lab_f( double xyz ){
    if ( xyz > 0.008856 ) 
	return pow( xyz ,  1.0/3 );
    else
	return ( 7.787 * xyz ) + ( 16 / 116.0 );
}

lab_color rgb2lab( const color & rgb )
{
    ///First ocnvert to XYZ
    double var_R = rgb.r;
    double var_G = rgb.g;
    double var_B = rgb.b;

    var_R = rgb2xyz_f( var_R);
    var_G = rgb2xyz_f( var_G);
    var_B = rgb2xyz_f( var_B);

    var_R = var_R * 100;
    var_G = var_G * 100;
    var_B = var_B * 100;

    //Observer. = 2°, Illuminant = D65
    double X,Y,Z;
    X = var_R * 0.4124 + var_G * 0.3576 + var_B * 0.1805;
    Y = var_R * 0.2126 + var_G * 0.7152 + var_B * 0.0722;
    Z = var_R * 0.0193 + var_G * 0.1192 + var_B * 0.9505    ;

    //Now convert to lab
    double var_X = X / ref_X;          //ref_X =  95.047   Observer= 2°, Illuminant= D65
    double var_Y = Y / ref_Y;          //
    double var_Z = Z / ref_Z;          //

    var_X = xyz2lab_f( var_X );
    var_Y = xyz2lab_f( var_Y );
    var_Z = xyz2lab_f( var_Z );

    return lab_color( 
	( 116 * var_Y ) - 16,
	500 * ( var_X - var_Y ),
	200 * ( var_Y - var_Z ) );
}


double lab2xyz_f( double lab ){
    if ( pow(lab, 3) > 0.008856 ) 
	return pow( lab, 3 );
    else
	return ( lab - 16 / 116.0 ) / 7.787;
}

double xyz2rgb_f( double xyz )
{
    if ( xyz > 0.0031308 ) 
	return 1.055 * ( pow( xyz , 1 / 2.4 ) ) - 0.055;
    else                     
	return 12.92 * xyz;
}

color lab2rgb( const lab_color & lab )
{
    double var_Y = ( lab.L + 16 ) / 116;
    double var_X = (lab.a) / 500 + var_Y;
    double var_Z = var_Y - lab.b / 200;

    var_X = lab2xyz_f( var_X );
    var_Y = lab2xyz_f( var_Y );
    var_Z = lab2xyz_f( var_Z );

    double X = ref_X * var_X;     //ref_X =  95.047     Observer= 2°, Illuminant= D65
    double Y = ref_Y * var_Y;     //ref_Y = 100.000
    double Z = ref_Z * var_Z;     //ref_Z = 108.883
    
    //now convert bact to the rgb
    var_X = X / 100;        //X from 0 to  95.047      (Observer = 2°, Illuminant = D65)
    var_Y = Y / 100;        //Y from 0 to 100.000
    var_Z = Z / 100;        //Z from 0 to 108.883

    double var_R = var_X *  3.2406 + var_Y * -1.5372 + var_Z * -0.4986;
    double var_G = var_X * -0.9689 + var_Y *  1.8758 + var_Z *  0.0415;
    double var_B = var_X *  0.0557 + var_Y * -0.2040 + var_Z *  1.0570;

    var_R = xyz2rgb_f( var_R );
    var_G = xyz2rgb_f( var_G );
    var_B = xyz2rgb_f( var_B );
    
    return color( var_R, var_G, var_B );
}


double sqr( double x ){ return x*x; };
lab_color lab_color::normalize()const
{
    double n = sqrt( sqr( L-50) + sqr(a)+ sqr(b) );
    if (n < 1e-6 ) return *this;
    double in = LAB_NORM_RADIUS/n;
    return lab_color( (L-50)*in+50, a*in, b*in );
}
