#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"

#include <string>
#include <iostream>
#include <algorithm>
#include <boost/filesystem.hpp>
#include "math.h"
#include <utility>
using namespace utility;
using namespace boost::filesystem;
using namespace std;

//TODO: ADD 1/8" margin around all outer radii!

double INCHES_TO_MM = 25.4;
double MATERIAL_WIDTH = 96.0 * INCHES_TO_MM;
double MATERIAL_LENGTH = 48.0 * INCHES_TO_MM;
double MATERIAL_HEIGHT = (23.0 / 32.0) * INCHES_TO_MM; 


double sq(double x) {
    return x * x;
}

class Vec2d { 
    public:
    double x, y;
    Vec2d(double x, double y) {
        this->x = x;
        this->y = y;
    }
    double dot(Vec2d &other) {
        return x * other.x + y * other.y;
    }
    double dist(Vec2d &pos) {
        return sqrt(sq(x - pos.x) + sq(y - pos.y));
    }
    double norm() {
        return dist(Vec2d(0,0));
    Vec2d normalize() {
        return Vec2d(x, y) * (1.0 / norm());
    }
    //Returns CCW rotated vector by 90 degrees
    Vec2d perp() {
        return Vec2d(y, -x);
    }
    Vec2d antiperp() {
        return Vec2d(-y, x);
    }
    //Returns the rotation direction (+1 = CCW) to another vector
    int rot_direction(Vec2d &other) {
        Mat2d matr(x, y, other.x, other.y);
        return sgn(matr.det());
    }
    bool operator==(Vec2d other) {
        return this.x == other.x && this.y == other.y;
    }
};

Vec2d operator+(const Vec2d &one, const Vec2d &two) {
    return Vec2d(one.x + two.x, one.y + two.y);
}
Vec2d operator*(const Vec2d &one, double scaling) {
    return Vec2d(one.x * scaling, one.y * scaling);
}
Vec2d operator*(double scaling, const Vec2d &one) {
    return one * scaling;
}
Vec2d operator-(const Vec2d &one, const Vec2d &two) {
    return one + (-1.0 * two);
}

class Shape {
    public:
    virtual Vec2d closest_point_to(Vec2d p);
    //Computes the minimal distance from the shape to a given point
    double distance_to(Vec2d p) {
        return closest_point_to(p).dist(p);
    }
};

class Ring : public Shape { 
    public:
    Vec2d pos;
    double inner_radius;
    double outer_radius;
    Ring(Vec2d pos_init, double inner_radius, double outer_radius) : 
        pos(pos_init.x, pos_init.y) { 
        this->inner_radius = inner_radius;
        this->outer_radius = outer_radius;
    }
    Ring move_to(Vec2d p) {
        return Ring(p, inner_radius, outer_radius);
    }
    //NOTE: Only tests for bounding circle intersections
    bool intersects(Ring other) {
        return pos.dist(other.pos) < (outer_radius + other.outer_radius);
    }
    double distance_to(Ring other) {
        return fmax(0, pos.dist(other.pos) - outer_radius - other.outer_radius);
    }
    Vec2d closest_point_to(Vec2d p) {
        if (strict_in(p)) {
            return p;
        }
        return (p - pos).normalize() * outer_radius;
    }
    bool strict_in(Vec2d test) {
        double dist = pos.dist(test);
        return dist > inner_radius && dist < outer_radius;
    }
    bool operator==(Ring other) {
        return (pos == other.pos) && (inner_radius == other.inner_radius) && 
            (outer_radius == other.outer_radius);
    }
};

class Slice { 
    public:
    //Bounding ring of mesh in mesh coordinates
    Ring bounding_ring;
    //Position of mesh center relative to material
    Vec2d pos;
    //Mesh represented (there is a unique slice per mesh)
    tinyobj::mesh_t &mesh;
    Slice(tinyobj::mesh_t &src_mesh) : 
        bounding_ring(Vec2d(0,0),0,0),
        pos(0,0), mesh(src_mesh) {
    }
};

class Mat2d {
    double a;
    double b;
    double c;
    double d;
    public:
    Mat2d(double a, double b, double c, double d) {
        this->a = a; this->b = b; this->c = c; this->d = d;
    }
    double det() {
        return a*c-b*d;
    }
    Mat2d inverse() {
        return Mat2d(d, -b, -c, a) * (1.0 / det());
    }
    Vec2d operator*(const Vec2d& in) {
        return Vec2d(a*in.x + b*in.y, c*in.x, d*in.y);
    }
    Mat2d operator*(double in) {
        return Mat2d(a * in, b * in, c * in, d * in);
    }
};
    

class LineSeg2d : public Shape {
    public:
    Vec2d pt1;
    Vec2d pt2;
    LineSeg2d(Vec2d pt1_init, Vec2d pt2_init) : pt1(pt1_init.x, pt1_init.y), 
        pt2(pt2_init.x, pt2_init.y) {}
    Vec2d perp() {
        return (pt1 - pt2).perp();
    }
    Line2d bisector() {
        return line_through((pt1 - pt2) * 0.5, perp());
    }
    Line2d toLine() {
        return line_through(pt1, tangent());
    }
    Vec2d tangent() {
        return pt2 - pt1;
    }
    Vec2d closest_point_to(Vec2d p) {
        Line2d l = line_through(p, perp());
        Vec2d closest = l.intersect(toLine());
        //Now, clamp to the line
        double line_seg_len = pt1.dist(pt2);
        double dist_1 = closest.dist(p1);
        double dist_2 = closest.dist(p2);
        if (dist_1 + dist_2 > line_seg_len) {
            if (dist_1 > dist_2) {
                return p2;
            }
            return p1;
        }
        return closest;
    }
    bool operator==(LineSeg2d other) {
        return pt1 == other.pt1 && pt2 == other.pt2;
    }
};

Line2d line_through(Vec2d p, Vec2d tangent) {
    return Line2d(tangent.x, tangent.y, tangent.dot(p));
}

//Line of the form ax+by=c
class Line2d {
    public:
    double a;
    double b;
    double c;
    Line2d(double a, double b, double c) {
        this->a = a; this->b = b; this->c = c;
    }
    Vec2d intersect(Line2d other) {
        Mat2d m(a, b, other.a, other.b);
        m = m.inverse();
        Vec2d k(c, other.c);
        return m*k;
    }
};

void swap(vector<Vec2d> &v, int i, int j) {
    Vec2d tmp = v[i];
    v[i] = v[j];
    v[j] = tmp;
}

void permute(vector<Vec2d> &in) {
    for (int i = 0; i < in.size(); i++) {
        int j = (rand() % (in.size() - i)) + i;
        swap(in, i, j);
    }
}
   
Ring calc_circle_directly(vector<Vec2d> boundaryPoints) {
    if (boundaryPoints.size() == 2) {
        Vec2d center = 0.5 * (boundaryPoints[0] + boundaryPoints[1]);
        return Ring(center, 0, center.dist(boundaryPoints[0]));
    }
    else if (boundaryPoints.size() == 3) {
        LineSeg2d l1(boundaryPoints[0], boundaryPoints[1]);
        LineSeg2d l2(boundaryPoints[1], boundaryPoints[2]);
        Vec2d center = l1.bisector().intersect(l2.bisector());
        return Ring(center, 0, center.dist(boundaryPoints[0]));
    }
    else {
        cerr << "Internal Error: Boundary points vector of wrong size\n";
        return Ring(Vec2d(0,0),0,0);
    }
}

//Citation: Emo Welzl
//and http://www.sunshine2k.de/coding/java/Welzl/Welzl.html
Ring bounding_circle_helper(vector<Vec2d> planePoints, vector<Vec2d> boundaryPoints) {
    if (planePoints.size() + boundaryPoints.size() <= 3) {
        vector<Vec2d> combined = planePoints.copy();
        combined.insert(combined.end(), boundaryPoints.begin(), boundaryPoints.end());
        return calc_circle_directly(combined);
    }
    else {
        permute(planePoints);
        Vec2d p = planePoints.pop_back();
        Ring r = bounding_circle_helper(planePoints, boundaryPoints);
        if (!p.strict_in(r)) {
            vector<Vec2d> enlargedBoundary = boundaryPoints.copy();
            enlargedBoundary.push_back(p);
            r = bounding_circle_helper(planePoints, enlargedBoundary);
            enlargedBoundary.clear();
        }
        return r;
    }
}

//Gets the bounding ring of a slice that already has been populated
//with bounding circle information
void get_bounding_ring(Slice &slice, vector<LineSeg2d> &lines) {
    Vec2d center = slice.bounding_ring.pos;
    if (intersects(center, lines)) {
        //Must NOT be a ring, do nothing
        return;
    }
    double inner = 0.0;
    for (auto l : lines) {
        inner = fmax(inner, l.closest_point_to(center));
    }
    slice.bounding_ring.inner_radius = inner;
}

//Computes whether or not a given point intersects any 2d cross-section of the slice
//Assumes that the vector of lines passed is in triples of triangle lines
bool intersects(Vec2d p, vector<LineSeg2d> &lines) {
    auto getDir = [](Vec2d p, LineSeg2d l) -> int {
        return l.tangent().rot_direction(p - l.pt1);
    };

    for (int i = 0; i < lines.size() / 3; i++) {
        LineSeg2d l1 = lines[3*i];
        LineSeg2d l2 = lines[3*i+1];
        LineSeg2d l3 = lines[3*i+2];
        //Technique drawn from http://www.blackpawn.com/texts/pointinpoly/
        int rot_direction = getDir(p, l1);
        if (rot_direction == getDir(p, l2) && rot_direction = getDir(p, l3)) {
            //If the vector from a triangle vertex to a point lies on the same
            //side of every line, it must be in the triangle.
            return true;
        }
    }
    return false;
}

//Gets all 2d line segments of triangles in an order that respects triangles.
vector<LineSeg2d> get_line_segs(Slice &slice) {
    vector<LineSeg2d> ret;
    for (int i = 0; i < slice.mesh.positions.size() / 9; i++) {
        Vec2d p1 = Vec2d(slice.mesh.positions[9*i], slice.mesh.positions[9*i+1]);
        Vec2d p2 = Vec2d(slice.mesh.positions[3+9*i], slice.mesh.position[3+9*i+1]);
        Vec2d p3 = Vec2d(slice.mesh.positions[6+9*i], slice.mesh.position[6+9*i+1]);
        ret.push_back(LineSeg2d(p1, p2));
        ret.push_back(LineSeg2d(p2, p3));
        ret.push_back(LineSeg2d(p3, p1));
    }
    return ret;
}

void compute_bounds(Slice &slice) {
    vector<Vec2d> points;
    //TODO: avoid this by overriding Vec2d with projection
    for (int i = 0; i < slice.mesh.positions.size / 3; i++) {
        points.push_back(Vec2d(slice.mesh.positions[3*i+0], slice.mesh.positions[3*i+1]));
    }
    slice.bounding_ring = bounding_circle_helper(points, vector<Vec2d>());

    vector<LineSeg2d> lineSegs = get_line_segs(slice);
    
    slice.bounding_ring = get_bounding_ring(slice, lineSegs);

}

//TODO: Ensure d_min kept with the list of corner placements to avoid re-calculation.
double min_dist(Ring c, Shape u, Shape v, &vector<pair<&Slice, Ring>> cfg, &vector<LineSeg2d> edges) {
    //d_min will hold the distance from the ring's center to the closest shape in the
    //configuration (including sides), but exclude u and v
    d_min = nan("");
    for (auto l : edges) {
        //TODO: How to test for equality here?
        if (true /*(not (l == u)) && (not (l == v)) */) {
            d_min = fmin(d_min, l.distance_to(c.pos));
        }
    }
    for (pair<&Slice, Ring> p : cfg) {
        Ring r = p.second;
        //TODO: How to test for equality here?
        if (true /*(not (r == u)) && (not (r == v)) */) {
            d_min = fmin(d_min, r.distance_to(c.pos));
        }
    }
    d_min = fmax(0, d_min - c.outer_radius);

    return d_min;
}

double hole_degree(Ring c, double d_min) {
    return (1 - (d_min / c.outer_radius));
}

struct Placement {
    int index;
    Vec2d pos;
    double d_min;
};



//Citation: This is a modified version of http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.99.5620&rep=rep1&type=pdf
vector<pair<&Slice, Ring>> ring_pack(vector<&Slice> in) {
    vector<pair<&Slice, Ring>> config;
    vector<Placement> placements;
    vector<LineSeg2d> edges;
    

    Vec2d v1 = Vec2d(0,0);
    Vec2d v2 = Vec2d(MATERIAL_WIDTH, 0);
    Vec2d v3 = Vec2d(MATERIAL_WIDTH, MATERIAL_LENGTH);
    Vec2d v4 = Vec2d(0, MATERIAL_LENGTH);

    edges.emplace_back(v1, v2);
    edges.emplace_back(v2, v3);
    edges.emplace_back(v3, v4);
    edges.emplace_back(v4, v1);

    auto radiusCompare = [](&Slice a, &Slice b) {
        return a.bounding_ring.outer_radius >= b.bounding_ring.outer_radius;
    }

    //Sorted in order of decreasing radius
    std::sort(in.begin(), in.end(), radiusCompare);

    //Instead of what's done in the paper, we arbitrarily pick an initial
    //configuration with the two smallest rings at opposite ends of the rectangle.
    //This reduces the time complexity from O(n^10) to O(n^8).
    Slice& s1 = in.pop_back();
    Slice& s2 = in.pop_back();
    double r1 = s1.bounding_ring.outer_radius;
    double r2 = s2.bounding_ring.outer_radius;
    config.emplace_back(s1, s1.bounding_ring.move_to(Vec2d(r1, r1)));
    config.emplace_back(s2, s2.bounding_ring.move_to(Vec2d(MATERIAL_WIDTH - r2, MATERIAL_LENGTH - r2));
    
    //Generate initial placements
    for (int i = 0; i < in.size(); i++) {
        vector temp = gen_placements(i, in[i].bounding_ring.outer_radius, config, edges);
        placements.insert(placements.end(), temp.begin(), temp.end());
    }

    auto addPlacement = [](Placement accepted) {
        Slice &s = in[accepted.index];
        Ring ci = s.bounding_ring.move_to(accepted.pos);

        config.emplace_back(s, ci);

        //TODO: Make this cache-friendlier?

        //Great, we modified "config", but now, we also need to modify both "placements"
        //and "in". We'll swap in[p.index] to the back of in, then pop it. 
        //Since we need to loop over "placements" anyway, we'll modify the indices as needed
        vector<Placement> new_placements;
        for (Placement p : placements) {
            //Ignore placements involving ci
            if (p.index == accepted.index) {
                continue;
            }
            Ring p_bounds = in[p.index].bounding_ring.move_to(p.pos);
            //Ignore placements that would overlap ci
            if (ci.intersects(p_bounds)) {
                continue;
            }
            //Update d_min for all remaining placements
            p.d_min = fmin(d_min, ci.distance_to(p_bounds));

            if (p.index == in.size() - 1) {
                //In this case, change the index to the index of the accepted.
                p.index = accepted.index; 
            }
            new_placements.push_back(p);
        }
        placements.clear();
        placements = new_placements;
        //Now, perform the swap and pop on 'in'
        in[accepted.index] = in[in.size() - 1];
        in.pop_back();

        //Okay, but we're still not done -- now we must add new placements involving ci.
        
    };

    while (placements.size() > 0) {
        //Stores the index into the corner placements with the maximum lambda
        int max_index = 0;
        int max_lambda = 0;
        //Find the best corner placement
        for (i = 0; i < placements.size(); i++) {
            double lambda = hole_degree(in[placements[i].index].bounding_ring, placements[i].d_min);
            if (lambda > max_lambda) {
                max_index = i;
                max_lambda = lambda;
            }
        }
        addPlacement(placements[max_index]);
    }
}

bool eps_equals(double a, double b) {
    return (a - b) < 0.001;
}

vector<double> quad_solve(double a, double b, double c) {
    vector<double> result;
    discrim = b*b-4*a*c;
    //TODO: use eps_equals!
    if (eps_equals(discrim, 0)) {
        result.push_back(-b / (2*a));
    }
    else if (discrim > 0) {
        result.push_back((-b + sqrt(discrim)) / 2*a);
        result.push_back((-b - sqrt(discrim)) / 2*a);
    }
    return result;
}

//Gets the center of a circle tangent to a circle and a horizontal line, with radius r
vector<Vec2d> get_tangent_circle_center_horiz(Ring one, double y_h, double r) {
    delta_x = sqrt(sq(one.outer_radius + r) - sq(y_h = r - one.pos.y));
    vector<Vec2d> result;
    result.emplace_back(one.pos.x + delta_x, y_h + r);
    result.emplace_back(one.pos.x - delta_x, y_h + r);
}
//Same as before. TODO: Make this not as dumb
vector<Vec2d> get_tangent_circle_center_vert(Ring one, double x_h, double r) {
    return perp(
        get_tangent_circle_center_horiz(Ring(one.pos.antiperp(), one.inner_radius, one.outer_radius),
                                           x_h, r));
}

    

vector<Vec2d> perp(vector<Vec2d> in) {
    vector<Vec2d> result(in.size());
    for (int i = 0; i < in.size(); i++) { result[i] = in.perp(); }
    return result;
}

//Gets the center of a circle tangent to the two given ones, with radius r.
//To do this, it uses a (horribly ugly) series of formulas.
vector<Vec2d> get_tangent_circle_center(Ring one, Ring two, double r) {
    double r1 = one.outer_radius;
    double r2 = two.outer_radius;
    double x1 = one.pos.x;
    double x2 = two.pos.x;
    double y1 = one.pos.y;
    double y2 = two.pos.y;
    double ydiff = y2 - y1;
    if (eps_equals(ydiff, 0)) {
        //Rotate the whole plane 90 degrees, and try again.
        return perp(
        get_tangent_circle_center(Ring(one.pos.antiperp(), one.inner_radius, one.outer_radius),
                                  Ring(two.pos.antiperp(), two.inner_radius, two.outer_radius),
                                  r));
    }
    //From subtracting the two "tangency" equations, we get a linear
    //equation of the form y = mx + k
    double m = (x1 - x2) / ydiff;
    double k = ((2*r+r1+r2)*(r1-r2) + x2*x2 + y2*y2 - x1*x1 - y1*y1) / (2.0 * ydiff);
    //Now, we substitute into (x-x1)^2+(y-y1)^2=(r+r1)^2 to get a quadratic equation in x
    double a = 1 + m*m;
    double b = 2*(m*(k-y1)-x1);
    double c = (x1*x1 + sq(k-y1) - sq(r+r1));
    vector<double> x = quad_solve(a, b, c);
    vector<Vec2d> result(x.size());
    std::transform(x.begin(), x.end(), result.begin(), [](double x) { return Vec2d(x, m*x + k) });
    return result;
}

bool placement_valid(Vec2d pos, double r, &vector<pair<&Slice, Ring>> cfg) {
    if (pos.x + r > MATERIAL_WIDTH || pos.x - r < 0 ||
        pos.y + r > MATERIAL_LENGTH || pos.y - r < 0) {
        return false;
    }
    Ring test = Ring(pos, 0, r);
    for (int i = 0; i < cfg.size(); i++) {
        if (cfg[i].second.intersects(test)) {
            return false;
        }
    }
    return true;
}

Placement pos_to_placement(Vec2d p, Shape u, Shape v, vector<pair<Slice&, Ring>> &cfg, vector<LineSeg2d> &edges) {
    Placement result;
    result.index = index;
    result.pos = p;
    result.d_min = min_dist(Ring(p, 0, r), u, v, cfg, edges);
    return result;
}

void add_centers(vector<Vec2d> centers, double r, vector<pair<&Slice, Ring>> &cfg, vector<LineSeg2d> &edges)

vector<Placement> gen_placements_involving(Ring one, int index, double r, &vector<pair<&Slice, Ring>> cfg,
                                           &vector<LineSeg2d> edges) {
}

//Generate all possible placements for a given index
vector<Placement> gen_placements(int index, double r, &vector<pair<&Slice, Ring>> cfg, &vector<LineSeg2d> edges) {
    auto pos_to_placement = [](Vec2d p, Shape u, Shape v) -> Placement { 
        Placement result;
        result.index = index;
        result.pos = p;
        result.d_min = min_dist(Ring(p, 0, r), u, v, cfg, edges);
        return result;
    }

    vector<Placement> result;

    auto add_center = [](Vec2d center, Shape one, Shape two) {
        if (placement_valid(center, r, cfg)) {
            result.push_back(pos_to_placement(center, one, two));
        }
    }
    
    auto add_centers = [](vector<Vec2d> centers, Shape one, Shape two) {
        for (Vec2d center : centers) { add_center(center, one, two) }
    }

    //Compute circle-circle and circle-line corner placements  
    for (int i = 0; i < cfg.size(); i++) {
        Ring one = cfg[i].second;
        for (int j = 0; j < cfg.size(); j++) {
            //TODO: Can this check be optimized by storing pairwise distances with the config?
            //Make sure that a corner placement is remotely possible
            Ring two = cfg[j].second;
            if (one.distance_to(second) > 2*r) {
                continue;
            }
            else {
                add_centers(get_tangent_circle_center(one, two, r));
            }
        }
        
        //Now, circle-line case
        add_centers(get_tangent_circle_horiz(one, 0, r));
        add_centers(get_tangent_circle_vert(one, 0, r));
        add_centers(get_tangent_circle_horiz(one, MATERIAL_LENGTH, r));
        add_centers(get_tangent_circle_vert(one, MATERIAL_WIDTH, r));
    }
    //Finally, test the line-line placements
    vector<Vec2d> corner;
    corner.emplace_back(r, r);
    corner.emplace_back(MATERIAL_WIDTH - r, r);
    corner.emplace_back(MATERIAL_WIDTH - r, MATERIAL_LENGTH - r);
    corner.emplace_back(r, MATERIAL_LENGTH - r);
    add_centers(corner);

    return result;
}


int main(int argc, char** argv) {
    if (argc < 2) {
        cout << "Usage: 3d-slice-layout path\n";
        return 1;
    }
    path source_dir(argv[1]);
    if (!exists(source_dir) || !is_directory(source_dir)) {
        cout << "The path input must be a directory!\n";
        return 1;
    }

    //For now, don't even bother to load materials
    //TODO: should we do this?
    vector<tinyobj::shape_t> shapes;
    vector<tinyobj::material_t> materials;

    for (directory_entry& entry : directory_iterator(source_dir)) {
        if (entry.path().extension() == ".obj") {
            vector<tinyobj::shape_t> temp_shapes;
            vector<tinyobj::material_t> temp_matls;
            string err;
            bool ret = tinyobj::LoadObj(temp_shapes, temp_matls, err, entry.path().c_str());
            if (!err.empty()) {
                std::cerr << err << std::endl;
            }
            if (!ret) {
                return 1;
            }

            shapes.insert(shapes.end(), temp_shapes.begin(), temp_shapes.end());
        }
    }

    return 0;
}



