#include <algorithm>
#include <chrono>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <math.h>

const std::string USAGE =
    "./surface_path pre.data post.data -d x y [options]\n"
    " [options]\n"
    " -ppath print indivual points\n"
    " -h  Print this help \n";

const uint32_t h_unit = 30;
const uint32_t v_unit = 11;

class Point{
public:
    double x;
    double y;
};

/// file I/O

void readFile(std::vector<uint8_t> &data,
	      std::string file_name,
	      uint32_t dx,
	      uint32_t dy)
{
    std::cout << "Reading "<<file_name <<" dims=[" << dx <<"x"<<dy <<"]\n";
    //std::vector<uint8_t> data;
    uint8_t max_v, min_v;
    std::fstream file;
    file.open(file_name, std::fstream::in | std::fstream::binary);
    for (uint32_t y = 0; y < dy; ++y)
	for (uint32_t x = 0; x < dx; ++x)
	    {
		uint8_t buff;
		file.read((char*)(&buff), sizeof(buff));
		data.push_back(buff);
		if (x+y == 0){max_v = buff; min_v = buff;}
		else {
		    max_v = std::max(max_v, buff);
		    min_v = std::min(min_v, buff);
		}
	    }
    printf("max: %hu min: %hu\n", max_v, min_v);
}

/// mesh operations

void getGridPoints(std::vector<Point> &out,
		   Point p1,
		   Point p2)
{
    Point p_start = {0, 0};
    Point p = {p2.x - p1.x, p2.y -p1.y};
    Point p_offset = p1;
    out.resize(0);
    out.push_back(p_offset);
    
    if (p.x == 0){
	// traverse along y
	for (uint32_t i=1; i<p.y; i++){
	    Point a = {p_offset.x, i+p_offset.y};
	    out.push_back(a);
	}
    }else if(p.y == 0){
	for (uint32_t i=1; i<p.x; i++){
	    Point a = {p_offset.x+i, p_offset.y};
	    out.push_back(a);
	}
    }else{
	double k = double(p.y)/double(p.x);
	while ((p_start.x < p.x) || (p_start.y < p.y)){
	    double x_next = p_start.x > int(p_start.x)? std::floor(p_start.x+1) : p_start.x+1;
	    double y_next = p_start.y > int(p_start.y)? std::floor(p_start.y+1) : p_start.y+1;
	    if (abs((x_next-p_start.x)*k) < abs(y_next-p_start.y)){
		// x_next is nearest
		p_start.y += (x_next-p_start.x)*k;
		p_start.x = x_next;
		//std::cout << "x="<<x_next<<"\n";
	    }else{
		// y_next is nearest
		p_start.x += (y_next-p_start.y)/k;
		p_start.y = y_next;
		//std::cout << "y="<<y_next<<"\n";
	    }
	    //std::cout <<"p "<< p_start.x <<" "<<p_start.y <<"\n";
	    out.push_back({p_start.x+p_offset.x, p_start.y+p_offset.y});
	}
    }
}

// utils

long toGlobal(uint32_t x, uint32_t y, uint32_t dx){
    return dx*y + x;
}

double toRealWorldSurfaceDist(double p1x,
			      double p1y,
			      double p2x,
			      double p2y,
			      uint8_t p1_val,
			      uint8_t p2_val)
{
    double base_x = (p2x - p1x) * h_unit;
    double base_y = (p2y - p1y) * h_unit;
    double base_z = (p2_val - p1_val) * v_unit;
    return sqrt(base_x*base_x + base_y*base_y + base_z*base_z);
}


double interpolate(uint8_t v1,
		   uint8_t v2,
		   double dist1,
		   double dist2,
		   double dist_src){
    double scale = (dist_src - dist1) / (dist2 - dist1);
    return v1+(v2-v1)*scale;
}

inline double Det(double a, double b, double c, double d)
{
        return a*d - b*c;
}


////  cell operations

void getCellVal(uint8_t *pts,
		uint32_t x_start,
		uint32_t y_start,
		uint32_t dx,
		std::vector<uint8_t> &data)
{
    pts[0] = data[toGlobal(x_start, y_start, dx)];
    pts[1] = data[toGlobal(x_start+1, y_start, dx)];
    pts[2] = data[toGlobal(x_start+1, y_start+1, dx)];
    pts[3] = data[toGlobal(x_start, y_start+1, dx)];
}

double getPointVal(Point p, 
		uint32_t x_start,
		uint32_t y_start,
		uint8_t *cell)
{
    // 0 1
    // 3 2
    uint8_t val = 0;
    if (p.x == x_start){
	// interpolate pts 0 and 3
	val = interpolate(cell[0], cell[3], y_start, y_start+1, p.y);
    }else if (p.x == x_start+1){
	// interpolate pts 1 and 2
	val = interpolate(cell[1], cell[2], y_start, y_start+1, p.y);
    }else if (p.y == y_start){
	// interpolate pts 0 and 1
	val = interpolate(cell[0], cell[1], x_start, x_start+1, p.x);
    }else if (p.y == y_start+1){
	// interpolate pts 2 and 3
	val = interpolate(cell[2], cell[3], x_start+1, x_start, p.x);
    }else{
	// diagonal
	val = interpolate(cell[1], cell[3], x_start+1, x_start, p.x);
    }
    return val;
}

bool getMiddlePoint(Point p1, Point p2, Point &ret){
    // 0 1
    // 3 2
    // assume triangulation [0 3 1] [1 3 2]
    uint32_t x_start = std::min(p1.x, p2.x);
    uint32_t y_start = std::min(p1.y, p2.y);
    
    // return false when overlaps at end points
    if ((p1.x == x_start+1) && (p1.y == y_start)) return false;
    if ((p2.x == x_start+1) && (p2.y == y_start)) return false;
    if ((p1.x == x_start) && (p1.y == y_start+1)) return false;
    if ((p2.x == x_start) && (p2.y == y_start+1)) return false;

    // https://stackoverflow.com/questions/16524096/how-to-calculate-the-point-of-intersection-between-two-lines
    // http://mathworld.wolfram.com/Line-LineIntersection.html
    Point p3 = {x_start+1, y_start};
    Point p4 = {x_start, y_start+1};
    double detL1 = Det(p1.x, p1.y, p2.x, p2.y);
    double detL2 = Det(p3.x, p3.y, p4.x, p4.y);
    double x1mx2 = p1.x - p2.x;
    double x3mx4 = p3.x - p4.x;
    double y1my2 = p1.y - p2.y;
    double y3my4 = p3.y - p4.y;

    double xnom = Det(detL1, x1mx2, detL2, x3mx4);
    double ynom = Det(detL1, y1my2, detL2, y3my4);
    double denom = Det(x1mx2, y1my2, x3mx4, y3my4);
    if(denom == 0.0)//Lines don't seem to cross
    {
        return false;
    }

    double ixOut = xnom / denom;   
    double iyOut = ynom / denom; 
    if(!isfinite(ixOut) || !isfinite(iyOut)) //Probably a numerical issue
        return false;
    ret.x = ixOut;
    ret.y = iyOut;
    return true;
}

double getCellDiff(Point p1,
		   Point p2,
		   std::vector<uint8_t> &pre,
		   std::vector<uint8_t> &post,
		   uint32_t dx,
		   uint32_t dy)
{
    if ((abs(p1.x - p2.x) > 1) || (abs(p1.y - p2.y) > 1))
	std::cerr << " 2 points must be in the same cell! \n";
    // read cell ccw
    uint8_t pts[4];
    uint32_t x_start = std::min(p1.x, p2.x);
    uint32_t y_start = std::min(p1.y, p2.y);
    double ret, ret2;
    // get diagonal point
    Point p_middle;
    bool hasMiddle = getMiddlePoint(p1, p2, p_middle);
    if (hasMiddle)
	std::cout << "tri point ["<< p_middle.x <<" "<<p_middle.y << "] \n";
    
    // pre
    getCellVal(pts, x_start, y_start, dx, pre);
    double p1_val = getPointVal(p1, x_start, y_start, pts);
    double p2_val = getPointVal(p2, x_start, y_start, pts);
    double pm_val = getPointVal(p_middle, x_start, y_start, pts);
    if (!hasMiddle)
	ret = toRealWorldSurfaceDist(p1.x, p1.y, p2.x, p2.y, p1_val, p2_val);
    else ret = toRealWorldSurfaceDist(p1.x, p1.y, p_middle.x, p_middle.y, p1_val, pm_val)
	     + toRealWorldSurfaceDist(p2.x, p2.y, p_middle.x, p_middle.y, p2_val, pm_val);

    // post
    getCellVal(pts, x_start, y_start, dx, post);
    p1_val = getPointVal(p1, x_start, y_start, pts);
    p2_val = getPointVal(p2, x_start, y_start, pts);

    // get cell difference
    if (!hasMiddle)
	ret2 = toRealWorldSurfaceDist(p1.x, p1.y, p2.x, p2.y, p1_val, p2_val);
    else {
	pm_val = getPointVal(p_middle, x_start, y_start, pts);
	ret2 = toRealWorldSurfaceDist(p1.x, p1.y, p_middle.x, p_middle.y, p1_val, pm_val)
	     + toRealWorldSurfaceDist(p2.x, p2.y, p_middle.x, p_middle.y, p2_val, pm_val);
    }
    ret = ret2 - ret;

    //std::cout << "cell pts values for ";
    //std::cout <<"["<< p1.x <<" "<<p1.y <<" "<<p2.x <<" "<<p2.y<<"] \n"
    //for (int i=0; i<4; i++)
    //	std::cout << unsigned(pts[i]) <<" ";
    //std::cout <<std::endl;

    return ret;
}



int main(int argc, char **argv)
{
    if (argc < 2) std::cout << USAGE<< std::endl;

    std::vector<std::string> args(argv, argv + argc);
    std::string file_name, file_name_post;
    uint32_t dx, dy;
    uint32_t p1x=0, p1y=0, p2x=0, p2y=0;
    bool print_path = false;
    
    // cmd line parse
    for (int i=1; i<argc; i++){
	if (args[i] == "-h") {
	    std::cout << USAGE<<std::endl;;
	    return 0;
	}else if(args[i] == "-ppath") print_path = true;
	    
	if (i == 1) {
	    file_name = args[i];
	    file_name_post = args[++i];
	}
	else if (args[i] == "-d"){
	    dx = stoi(args[++i]);
	    dy = stoi(args[++i]);
	}
    }

    // read from disk
    std::vector<uint8_t> pre_data, post_data;
    readFile(pre_data, file_name, dx, dy);
    readFile(post_data, file_name_post, dx, dy);

    std::vector<Point> pts;
    bool quit = false;

    std::cout << "Start program, type quit to stop\n";

    while (!quit){		
	std::string in;
	std::cout <<"\n[user input 4 nums]: ";
	std::getline(std::cin, in);

	if (in == "quit") {quit = true; break;}
	    
	std::string temp;
	std::vector <int> nums;
	std::stringstream ss;
	ss << in;
	while(ss >> temp)
	    nums.push_back(stoi(temp));
	if (nums.size() < 4){ std::cout << "input 4 ints, p1x p1y p2x p2y\n"; continue; }
	p1x = nums[0];
	p1y = nums[1];
	p2x = nums[2];
	p2y = nums[3];

	Point s = {p1x,p1y};
	Point e = {p2x,p2y};
	getGridPoints(pts, s, e);

	std::cout << "\n";
	double diff = 0;
	for (uint32_t i =0; i<pts.size()-1; i++){
	    auto v = getCellDiff(pts[i], pts[i+1], pre_data, post_data, dx, dy);
	    diff += v;
	    if (print_path)
		std::cout << "from ["<<pts[i].x <<", "<<pts[i].y<<"] \n" 
			  << "to   ["<<pts[i+1].x <<", "<<pts[i+1].y<<"] \n"
			  << "cell diff="<< v  <<" \n\n" ;
	}
	
	std::cout << "diff=" << diff <<"\n";
    }
    
    return 0;
}
