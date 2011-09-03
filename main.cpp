#include <png++/png.hpp>
#include <iostream>
#include <math.h>
#include <vector>
#include <set>
#include <algorithm>
#include <cassert>
#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <locale>
#include "color_conversions.hpp"

using namespace std;

//Integer rotation
struct irot{
    double k1, k2;
    irot( double angle ){
	k1 = tan(angle*0.5);
	k2 = -sin(angle);
    };
    int iround( double x ){
	return int(floor(x+0.5));
    };
    void rotate( int& x, int & y ){
        x += iround( y*k1 );
        y += iround( x*k2 );
        x += iround( y*k1 );
    };
};

struct decreaser{ size_t operator()(size_t x ){ return x - 1; }; };

struct diagram{
    size_t size;
    size_t max_iters;
    
    typedef std::vector<size_t> index_vec_t;
    typedef std::vector<std::pair<int,int> > point_vector_t;
    index_vec_t data;

    diagram( size_t size_ ){
	size = size_;
	max_iters = 0; 
    }

    bool in_range( int x, int y ) const {
	return ( x >= 0  && x < (int)size  && 
		 y >= 0  && y < (int)size );
    };

    size_t & at( int x, int y){
	return data[x+y*size];
    };
    size_t at( int x, int y)const{
	return data[x+y*size];
    };
    size_t calculate( double angle, int dx, int dy, size_t repeats=1 );
};

size_t diagram::calculate( double angle, int dx, int dy, size_t repeats ){
    size_t next_idx = 0; //starting numeration from 1.    
    point_vector_t points;
    size_t NO_ORBIT = (size_t)(-1);

    data.resize( size*size );
    std::fill( data.begin(), data.end(), NO_ORBIT );
	
    int xc, yc; 
    yc = size / 2;
    xc = size / 2;

    max_iters = size/2+max(abs(dx),(dy)) * 100; //guard to avoid hanging, if algorithm goes into infinite cycle.

    irot rotate( angle );

    for ( int yy = 0; yy < (int)size; ++ yy ){
	cout<< yy << " of "<< size << "    \r"; cout.flush();
	for ( int xx = 0; xx < (int)size; ++xx ){ 
	    int p0x=xx-xc+dx, p0y=yy-yc+dy;
	    if ( at( xx, yy ) != NO_ORBIT ) continue; //this point is already processed.
            
	    points.push_back( make_pair(xx, yy));
	    size_t orbit = NO_ORBIT; //index of the current orbit. 0 - not found.
	    int px = p0x, py = p0y; //current point
	    size_t iters = 0;
	    while ( ++iters < max_iters ){
		for( size_t i_rep=0; i_rep < repeats; ++ i_rep) //rotate several times
		    rotate.rotate( px, py );
		int px_s = px + xc - dx, py_s = py + yc - dy; //screen coordinates
		if ( orbit == NO_ORBIT && in_range(px_s, py_s) ){
		    if ( px == p0x && py==p0y ) 
			break; //cycle closed
		    orbit = at( px_s, py_s );
		    points.push_back( make_pair(px_s,py_s) );
		}
	    }
                
	    if ( orbit == NO_ORBIT){ //all cycle is not indexed.
		orbit = next_idx; //assign new orbit
		next_idx += 1;
	    }
	    point_vector_t::iterator ip, ep = points.end();
	    for ( ip = points.begin(); ip != ep; ++ip){
		at(ip->first, ip->second) = orbit;
	    }
	    points.clear();
	}
    }
    return next_idx;
};

struct build_graph{
    const diagram & diag;
    size_t num_groups;
    typedef std::set< size_t > intermediate_graph_record_t;
    typedef std::vector< size_t > graph_record_t;

    intermediate_graph_record_t * temp_graph;
    graph_record_t * graph;

    build_graph( const diagram& diagram_, size_t num_groups_ )
	:diag(diagram_), num_groups( num_groups_ ){	
	temp_graph = new intermediate_graph_record_t[ num_groups ];
	graph = NULL;

	int w,h;
	w=h= diag.size;
	for( int y=0; y < h; ++ y ){
	    for( int x=0; x < w; ++ x){
		size_t v = diag.at(x,y);
		if (x > 0) doadd( v, diag.at(x-1,y) );
		if (y > 0) doadd( v, diag.at(x,y-1) );
	    }
	}
	compile_graph();
	delete[] temp_graph;
	temp_graph = NULL;
	cout<<"Graph compiled"<<endl;
    }
    ~build_graph(){
	delete[] graph;
    };
    void compile_graph();
    graph_record_t & neighbores( size_t orbit ){
	assert( orbit < num_groups );
	return graph[ orbit ];
    };
    static void insert( intermediate_graph_record_t & rec, size_t value ){
	rec.insert( value );
    };

    void doadd( size_t v0, size_t v1){
        if (v1 != v0 ){
	    insert( temp_graph[v1], v0 );
	    insert( temp_graph[v0], v1 );
	}
    };
};
void build_graph::compile_graph()
{
    graph = new graph_record_t[ num_groups ];
    for( size_t i = 0; i < num_groups; ++i){
	graph[i].resize( temp_graph[i].size() );
	copy( temp_graph[i].begin(), temp_graph[i].end(), 
	      graph[i].begin() );
	sort(graph[i].begin(), graph[i].end());
    }
}

////////////////////////////////////////////////////////////////////////////////
// RGB colorization
////////////////////////////////////////////////////////////////////////////////
struct colorize{
    build_graph & graph;
    typedef std::vector< color > colormap_t;
    colormap_t colormap;
    colorize( build_graph & graph_ )
	:graph( graph_ ),
	 colormap( graph_.num_groups ){
    };
    color meancolor( size_t orbit );
    void iterate();
    void operator()(size_t num_iters);
    void write_image( const diagram & dia, const char * ofile );
    void random_init();
    void auto_levels();
};

void colorize::auto_levels()
{
    double r0=1,g0=1,b0=1;
    double r1=0,g1=0,b1=0;
    for( colormap_t::const_iterator i = colormap.begin(); i != colormap.end(); ++ i){
	r0 = min( r0, i->r );
	g0 = min( g0, i->g );
	b0 = min( b0, i->b );
	r1 = max( r1, i->r );
	g1 = max( g1, i->g );
	b1 = max( b1, i->b );
    }
    double kr = 1 / max( 1e-5, r1-r0 );
    double kg = 1 / max( 1e-5, g1-g0 );
    double kb = 1 / max( 1e-5, b1-b0 );
    for( colormap_t::iterator i = colormap.begin(); i != colormap.end(); ++i) {
	i->r = (i->r - r0)*kr;
	i->g = (i->g - g0)*kg;
	i->b = (i->b - b0)*kb;
    }
}
	
void colorize::random_init()
{
    for( size_t i = 0; i < colormap.size(); ++i){
	color randc( rand()/(double)RAND_MAX,
		     rand()/(double)RAND_MAX,
		     rand()/(double)RAND_MAX );
	colormap[i]=randc.normalize();
    }
}

void colorize::write_image( const diagram & dia, const char * ofile )
{
    png::image< png::rgb_pixel > image( dia.size, dia.size );
    for( size_t y = 0; y<dia.size; ++y){
	for( size_t x = 0; x<dia.size; ++x){
	    color & c = colormap[ dia.at(x,y) ];
	    image[y][x] = png::rgb_pixel( 
		c.int_r(), c.int_g(), c.int_b() );
	}
    }
    image.write( ofile );
}
color colorize::meancolor( size_t orbit )
{
    color rval(0,0,0);
    build_graph::graph_record_t::const_iterator i, e;
    build_graph::graph_record_t & neighbores( graph.neighbores( orbit));
    e = neighbores.end();
    for( i = neighbores.begin(); i != e; ++i ){
	rval = rval + colormap[ *i ];
    }
    if (neighbores.size() == 0 ){
	return color(0.5, 0.5, 0.5);
    }else{
	return rval * (1.0/neighbores.size());
    };
}
    
void colorize::iterate()
{
    colormap_t new_colormap( colormap.size() );
    for( size_t i = 0; i < colormap.size(); ++ i){
	new_colormap[i] = meancolor( i ).normalize();
    };
    colormap.swap(new_colormap);
    auto_levels();
}

void colorize::operator()(size_t num_iters)
{
    random_init();
    for( size_t i =0; i < num_iters; ++i){
	cout<<"Step "<<i<<" of "<<num_iters<<"    \r";
	cout.flush();
	iterate();
    };
}

////////////////////////////////////////////////////////////////////////////////
// LAB colorization
////////////////////////////////////////////////////////////////////////////////
struct colorize_lab{
    build_graph & graph;
    typedef std::vector< lab_color > colormap_t;
    colormap_t colormap;
    colorize_lab( build_graph & graph_ )
	:graph( graph_ ),
	 colormap( graph_.num_groups ){
    };
    lab_color meancolor( size_t orbit );
    void iterate();
    void operator()(size_t num_iters);
    void write_image( const diagram & dia, const char * ofile );
    void random_init();
    void auto_levels();
};

void colorize_lab::auto_levels()
{
    double r0=1,g0=1,b0=1;
    double r1=0,g1=0,b1=0;
    for( colormap_t::const_iterator i = colormap.begin(); i != colormap.end(); ++ i){
	color c( lab2rgb( *i ) );
	r0 = min( r0, c.r );
	g0 = min( g0, c.g );
	b0 = min( b0, c.b );
	r1 = max( r1, c.r );
	g1 = max( g1, c.g );
	b1 = max( b1, c.b );
    }
    double kr = 1 / max( 1e-5, r1-r0 );
    double kg = 1 / max( 1e-5, g1-g0 );
    double kb = 1 / max( 1e-5, b1-b0 );
    for( colormap_t::iterator i = colormap.begin(); i != colormap.end(); ++i) {
	color c = lab2rgb( *i );
	color c1( (c.r - r0)*kr,
		  (c.g - g0)*kg,
		  (c.b - b0)*kb );
	*i = rgb2lab( c1 );
    }
}
	

	
void colorize_lab::random_init()
{
    for( size_t i = 0; i < colormap.size(); ++i){
	lab_color randc( rand()/(double)RAND_MAX * 100,
			 rand()/(double)RAND_MAX * 100 - 50,
			 rand()/(double)RAND_MAX * 100 - 50 );
	colormap[i]= randc.normalize();
    }
}

void colorize_lab::write_image( const diagram & dia, const char * ofile )
{
    png::image< png::rgb_pixel > image( dia.size, dia.size );
    for( size_t y = 0; y<dia.size; ++y){
	for( size_t x = 0; x<dia.size; ++x){
	    color c = lab2rgb( colormap[ dia.at(x,y) ] ).normalize();
	    image[y][x] = png::rgb_pixel( 
		c.int_r(), c.int_g(), c.int_b() );
	}
    }
    image.write( ofile );
}
lab_color colorize_lab::meancolor( size_t orbit )
{
    lab_color rval(0,0,0);
    build_graph::graph_record_t::const_iterator i, e;
    build_graph::graph_record_t & neighbores( graph.neighbores( orbit));
    e = neighbores.end();
    for( i = neighbores.begin(); i != e; ++i ){
	rval = rval + colormap[ *i ];
    }
    if (neighbores.size() == 0 ){
	return lab_color(0.5, 0, 0);
    }else{
	return rval * (1.0/neighbores.size());
    };
}
    
void colorize_lab::iterate()
{
    colormap_t new_colormap( colormap.size() );
    for( size_t i = 0; i < colormap.size(); ++ i){
	new_colormap[i] = meancolor( i ).normalize();
    };
    colormap.swap(new_colormap);
}

void colorize_lab::operator()(size_t num_iters)
{
    random_init();
    for( size_t i =0; i < num_iters; ++i){
	cout<<"Step "<<i<<" of "<<num_iters<<"    \r";
	cout.flush();
	iterate();
    };
}

namespace po = boost::program_options;

enum ColorModelT{
    COLOR_MODEL_RGB,
    COLOR_MODEL_LAB,
    COLOR_MODEL_UNKNOWN
};
ColorModelT convert_cm( const string & cm_name ){
    string sname;
    transform( cm_name.begin(), cm_name.end(), back_inserter( sname ), ::tolower );
    if ( sname == "rgb" )
	return COLOR_MODEL_RGB;
    if ( sname == "lab" )
	return COLOR_MODEL_LAB;
    return COLOR_MODEL_UNKNOWN;
}

int main(int argc, char *argv[])
{
    using boost::lexical_cast;
    using boost::bad_lexical_cast;

    size_t diagram_size = 2048;
    size_t smoothing_steps = 20;
    double angle = M_PI/7*2;
    string image_name = "output.png";
    string s_angle;
    string s_color_model = "rgb";
    ColorModelT color_model;
    size_t repeats = 1;
    int dx = 0, dy = 0;

    //parsing and setting program options
    po::options_description desc("Allowed options");
    desc.add_options()
	("help", "Show help message")
	("size,s", po::value<size_t>( &diagram_size )->default_value(512), 
	 "Size of the image to generate. Default is 512.")
	("smoothing,S", po::value<size_t>( &smoothing_steps )->default_value(10), 
	 "How many smoothing steps to perform. Default is 10.")
	("output,o", po::value<string>( &image_name )->default_value("output.png"), 
	 "Output file name, PNG format is used.")
	("angle,a", po::value<string>( &s_angle ), "Angle, given as multiplier to PI. for example, '2/5' gives 2/5PI")
	("repeat,r", po::value<size_t>( &repeats )->default_value(1), "How many times to repeat rotation (default is 1)" )
	("dx,x", po::value<int>( & dx )->default_value(0), "Horizontal offset, default is 0: center of the image is rendered" )
	("dy,y", po::value<int>( & dy )->default_value(0), "Vertical offset, default is 0: center of the image is rendered" )
	("color,C", po::value<string>( &s_color_model )->default_value("rgb"), 
	 "Color model to use for smoothing. Possbile values are RGB (default) and LAB" );
    
    po::variables_map vm;
    try{
	po::store( po::parse_command_line( argc, argv, desc ), vm );
	po::notify( vm );
    }catch( exception e ){
	cerr<<"Failed to parse command-line options:"<<e.what()<<endl;
	return 1;
    }
    if( vm.count("help") ){
	cout << desc <<endl;
	return 0;
    }
    if ( diagram_size > 10000 ){
	cerr << "Size is too big"<<endl;
	return 1;
    }
    color_model = convert_cm( s_color_model );
    if (color_model == COLOR_MODEL_UNKNOWN ){
	cerr << "Unknown name of the oclor model: "<<s_color_model<<endl;
	return 1;
    }

    
    if (s_angle.empty()){
	angle = M_PI/5*2;
    }else{
	try{
	    size_t slash_pos = s_angle.find( "/" );
	    if (slash_pos == s_angle.npos ){
		//no slash
		double angle_k = lexical_cast<double>(s_angle);
		angle = angle_k * M_PI;
		cout << "Using angle "<<angle_k<<" PI"<<endl;
	    }else{
		double num = lexical_cast<double>( s_angle.substr( 0, slash_pos ) );
		double den = lexical_cast<double>( s_angle.substr( slash_pos+1 ) );
		angle = M_PI * num / den;
		cout << "Using angle "<<num<<"/"<<den<<" PI"<<endl;
	    }
	}catch( bad_lexical_cast e ){
	    cerr<<"Failed to parse angle:"<<e.what()<<endl;
	    return 1;
	}
    }
    

    srand( time(NULL) );

    cout<<"Searching for all orbits..."<<endl;
    diagram dia( diagram_size );
    size_t orbits = dia.calculate( angle, dx, dy, repeats );
    std::cout<<"Found "<<orbits<<" orbits"<<std::endl;

    cout<<"Now building orbit adjacency graph"<<endl;
    build_graph graph( dia, orbits );
    cout<<"Done"<<endl;

    cout<<"Now colorizing "<<smoothing_steps<<" steps"<<endl;
    if (color_model == COLOR_MODEL_RGB ){
	cout << "Using RGB color model for smoothing"<<endl;
	colorize colorizer( graph );
	colorizer( smoothing_steps );
	cout <<"Done colorizing"<<endl;
	cout<<"Writing image "<<image_name<<endl;
	colorizer.write_image( dia, image_name.c_str() );
    }else if (color_model == COLOR_MODEL_LAB ){
	cout << "Using LAB color model for smoothing"<<endl;
	colorize_lab colorizer( graph );
	colorizer( smoothing_steps );
	cout <<"Done colorizing"<<endl;
	cout<<"Writing image "<<image_name<<endl;
	colorizer.write_image( dia, image_name.c_str() );
    }
    
    return 0;
}
