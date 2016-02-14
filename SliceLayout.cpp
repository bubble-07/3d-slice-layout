#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"

#include <string>
#include <iostream>
#include <unordered_map>
#include <algorithm>
#include <boost/filesystem.hpp>
#include "math.h"
#include <utility>
using namespace boost::filesystem;
using namespace std;

//TODO: ADD 1/8" margin around all outer radii!

double INCHES_TO_MM = 25.4;
double MATERIAL_WIDTH = 96.0 * INCHES_TO_MM;
double MATERIAL_LENGTH = 48.0 * INCHES_TO_MM;
double MATERIAL_HEIGHT = (23.0 / 32.0) * INCHES_TO_MM; 

//To be set later
int BOTTOM = 0;
int TOP = 0;
int LEFT = 0;
int RIGHT = 0;

bool eps_equals(double a, double b) {
    return (fabs(a - b)) < 0.00001;
}

double sq(double x) {
    return x * x;
}
int sgn(double x) {
    return (x > 0) - (x < 0);
}

//I don't like that std::vector doesn't return the popped element.
//Still has UB if in is empty, though
template<typename T> T pop_back(vector<T> &in) {
    T ret = in[in.size() - 1];
    in.pop_back();
    return ret;
}

template<typename T> void swap(vector<T> &v, int i, int j) {
    T tmp = v[i];
    v[i] = v[j];
    v[j] = tmp;
}

class Mat2d;

class Vec2d { 
    public:
    double x, y;
    Vec2d(double x, double y) {
        this->x = x;
        this->y = y;
    }
    double dot(const Vec2d &other) {
        return x * other.x + y * other.y;
    }
    double dist(const Vec2d &pos) {
        return sqrt(sq(x - pos.x) + sq(y - pos.y));
    }
    double norm() {
        return dist(Vec2d(0,0));
    }
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
    int rot_direction(const Vec2d &other); 

    bool operator==(Vec2d other) {
        return x == other.x && y == other.y;
    }
    const Vec2d operator*(double scaling) {
        return Vec2d(x * scaling, y * scaling);
    }
};

Vec2d operator+(const Vec2d &one, const Vec2d &two) {
    return Vec2d(one.x + two.x, one.y + two.y);
}
Vec2d operator*(double scaling, const Vec2d &one) {
    return Vec2d(one.x * scaling, one.y * scaling);
}
Vec2d operator-(const Vec2d &one, const Vec2d &two) {
    return one + (-1.0 * two);
}

class Mat2d {
    double a;
    double b;
    double c;
    double d;
    public:
    Mat2d(double a_init, double b_init, double c_init, double d_init) {
        a = a_init; b = b_init; c = c_init; d = d_init;
    }
    bool singular() {
        return eps_equals(det(), 0);
    }
    double det() {
        return a*c-b*d;
    }
    Mat2d inverse() {
        return Mat2d(d, -b, -c, a) * (1.0 / det());
    }
    private:
    Vec2d linsolvehelper(Vec2d expected, bool swapped) {
        //Technical specification: inverse() * expected, but numerically stable.
        if (swapped || fabs(a) >= fabs(c)) {
            double y = (expected.y - ((c*expected.x)/a)) / (d - (c*b)/a);
            double x = ((expected.x - b*y) / a);
            return Vec2d(x, y);
        }
        else {
            Vec2d temp = Mat2d(c, d, a, b).linsolvehelper(Vec2d(expected.y, expected.x), true);
            return Vec2d(temp.y, temp.x);
        }
    }
    public:
    Vec2d linsolve(Vec2d expected) {
        return linsolvehelper(expected, false);
    }
    Vec2d linsolve2(Vec2d expected) {
        return linsolven(expected, 2);
    }
    Vec2d linsolven(Vec2d expected, int n) {
        Vec2d x = linsolve(expected);
        for (int i = 0; i < n - 1; i++) {
            Vec2d actual = operator*(x);
            Vec2d dx = linsolve(expected - actual);
            x = x + dx;
        }
        return x;
    }
    Vec2d operator*(const Vec2d& in) {
        return Vec2d(a*in.x + b*in.y, c*in.x + d*in.y);
    }
    Mat2d operator*(double in) {
        return Mat2d(a * in, b * in, c * in, d * in);
    }
};

int Vec2d::rot_direction(const Vec2d &other) {
    Mat2d matr(x, y, other.x, other.y);
    return sgn(matr.det());
}   

class Vec3d {
    public:
    double x, y, z;
    Vec3d(double x, double y, double z) {
        this->x = x;
        this->y = y;
        this->z = z;
    }
    double dot(const Vec3d &other) {
        return x * other.x + y * other.y + z * other.z;
    }
    double norm() {
        return sqrt(x * x + y * y + z * z);
    }
    Vec3d normalize() {
        return Vec3d(x, y, z) * (1.0 / norm());
    }
    const Vec3d operator*(double scaling) {
        return Vec3d(x * scaling, y * scaling, z * scaling);
    }
    Vec3d cross(const Vec3d& o) {
        return Vec3d(y*o.z - z * o.y, z * o.x - x * o.z, x * o.y - y * o.x);
    }
};
Vec3d operator+(const Vec3d &one, const Vec3d &two) {
    return Vec3d(one.x + two.x, one.y + two.y, one.z + two.z);
}
Vec3d operator*(double scaling, const Vec3d &one) {
    return Vec3d(one.x * scaling, one.y * scaling, one.z * scaling);
}
Vec3d operator-(const Vec3d &one, const Vec3d &two) {
    return one + (-1.0 * two);
}

class Shape {
    public:
    enum Type { RING, LINESEG };
    Type type;
    virtual Vec2d closest_point_to(Vec2d p) = 0;
    //Computes the minimal distance from the shape to a given point
    double distance_to(Vec2d p) {
        return closest_point_to(p).dist(p);
    }
    virtual bool operator==(Shape& other) {
        return false;
    }
    bool operator!=(Shape& other) {
        return !(operator==(other));
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
        this->type = Shape::Type::RING;
    }
    string to_string() {
        return "(" + std::to_string(pos.x) + "," + std::to_string(pos.y) + "," + std::to_string(outer_radius) + ")";
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
    using Shape::distance_to;

    bool in(Vec2d test) {
        double dist = pos.dist(test);
        return dist >= inner_radius && dist <= outer_radius;
    }
    bool fuzz_in(Vec2d test) {
        double dist = pos.dist(test);
        return dist >= inner_radius - 0.001 && dist <= outer_radius + 0.001;
    }

    bool strict_in(Vec2d test) {
        double dist = pos.dist(test);
        return dist > inner_radius && dist < outer_radius;
    } 
    Vec2d closest_point_to(Vec2d p) {
        if (strict_in(p)) {
            return p;
        }
        return (p - pos).normalize() * outer_radius;
    }
    
    bool operator==(Shape& o) {
        if (o.type == Shape::Type::RING) {
            Ring& other = (Ring&) o;
            return (pos == other.pos) && (inner_radius == other.inner_radius) && 
                (outer_radius == other.outer_radius);
        }
        return false;
    }
};
//Line of the form ax+by=c
class Line2d {
    public:
    double a;
    double b;
    double c;
    Line2d(double a, double b, double c) {
        this->a = a;
        this->b = b;
        this->c = c;
    }
    //Construct a line using <X, c> = b
    Line2d(Vec2d c, double b) : Line2d(c.x, c.y, b) {}

    Vec2d intersect(Line2d other) {
        Mat2d m(a, b, other.a, other.b);
        //m = m.inverse();
        Vec2d k(c, other.c);
        //Vec2d result = m*k;
        Vec2d result = m.linsolven(k, 4);
        return result;
    }
    bool parallel(Line2d other) {
        return Mat2d(a, b, other.a, other.b).singular();
    }
    Vec2d tangent() {
        return Vec2d(b, -a);
    }
    Vec2d normal() {
        return Vec2d(a, b);
    }
    Vec2d closest_point_to(Vec2d p) {
        return intersect(Line2d(tangent(), p.dot(tangent())));
    }
};

Line2d line_through(Vec2d p, Vec2d tangent) {
    Vec2d perp = tangent.perp();
    return Line2d(perp, perp.dot(p));
}

class LineSeg2d : public Shape {
    public:
    Vec2d pt1;
    Vec2d pt2;
    LineSeg2d(Vec2d pt1_init, Vec2d pt2_init) : pt1(pt1_init.x, pt1_init.y), 
        pt2(pt2_init.x, pt2_init.y) {
        this->type = Shape::Type::LINESEG;     
    }
    bool parallel(LineSeg2d other) {
        Vec2d t = other.tangent();
        return Mat2d(tangent().x, tangent().y, t.x, t.y).singular();
    }
    Vec2d perp() {
        return tangent().perp();
    }
    Line2d bisector() {
        return line_through((pt1 + pt2) * 0.5, perp());
    }
    Line2d toLine() {
        return line_through(pt1, tangent());
    }
    Vec2d tangent() {
        return (pt2 - pt1);
    }
    Vec2d closest_point_to(Vec2d p) {
        Line2d l = line_through(p, perp());
        Vec2d closest = l.intersect(toLine());
        //Now, clamp to the line
        double line_seg_len = pt1.dist(pt2);
        double dist_1 = closest.dist(pt1);
        double dist_2 = closest.dist(pt2);
        if (dist_1 + dist_2 > line_seg_len) {
            if (dist_1 > dist_2) {
                return pt2;
            }
            return pt1;
        }
        return closest;
    }
    bool operator==(Shape& o) {
        if (o.type == Shape::Type::LINESEG) {
            LineSeg2d& other = (LineSeg2d&) o;
            return pt1 == other.pt1 && pt2 == other.pt2;
        }
        return false;
    }
};

//Input: file, mesh, prior vertex index. Output: new vertex index.
int export_mesh(FILE* out, tinyobj::mesh_t m, int index) {

    auto comma_sep_triple = [](vector<float> &in, int base) -> string {
        return to_string(in[base]) + " " +
               to_string(in[base + 1]) + " " +
               to_string(in[base + 2]);
    };
    //The Obj format is STUPID and has indices starting from 1.
    auto s = [](int in) -> string {
        return to_string(in + 1) + "//" + to_string(in + 1);
    };

    //First, export all vertices
    for (int i = 0; i < m.positions.size(); i += 3) {
        fputs(("v " + comma_sep_triple(m.positions, i) + "\n").c_str(), out);
    }
    //Then, export all vertex normals
    for (int i = 0; i < m.normals.size(); i += 3) {
        fputs(("vn " + comma_sep_triple(m.normals, i) + "\n").c_str(), out);
    }
    //Export all faces (trivial)
    for (int i = index; i < (m.positions.size() / 3) + index; i += 3) {
        fputs(("f " + s(i) + " " + s(i + 1) + " " + s(i + 2) + "\n").c_str(), out);
    }
    return index + (m.positions.size() / 3);
}

template<typename T> void append(vector<T> &one, vector<T> &two) {
    one.insert(one.end(), two.begin(), two.end());
}

//Converts the silly index form to a raw sequential form
tinyobj::mesh_t& standardize_mesh(tinyobj::mesh_t& in) {
    vector<float> positions;
    vector<float> normals;
    vector<unsigned int> indices;
    for (int i = 0; i < in.indices.size(); i++) {
        positions.push_back(in.positions[3 * in.indices[i]]);
        positions.push_back(in.positions[3 * in.indices[i] + 1]);
        positions.push_back(in.positions[3 * in.indices[i] + 2]);
        
        normals.push_back(in.normals[3 * in.indices[i]]);
        normals.push_back(in.normals[3 * in.indices[i] + 1]);
        normals.push_back(in.normals[3 * in.indices[i] + 2]);

        indices.push_back(i);
    }
    in.positions.clear();
    in.normals.clear();
    in.indices.clear();

    in.positions = positions;
    in.normals = normals;
    in.indices = indices;

    return in;
}

//Merges the second mesh's data onto the first's
tinyobj::mesh_t& merge_meshes(tinyobj::mesh_t& first, tinyobj::mesh_t& second) {
    append(first.positions, second.positions);
    append(first.normals, second.normals);
    //For now, we don't care about material ids, texcoords, indices
    //TODO: Should we care?
    return first;
}

class Slice { 
    public:
    //Bounding ring of mesh in mesh coordinates
    Ring bounding_ring;
    //Mesh represented (there is a unique slice per mesh)
    tinyobj::mesh_t mesh;
    //Name of the slice (as a composite object)
    string name;

    Slice(string src_name, tinyobj::mesh_t src_mesh) : 
        bounding_ring(Vec2d(0,0),0,0),
        mesh(src_mesh), name(src_name){}

    //Exports an object with the prior vertex index
    int export_obj(FILE* dest, int index) {
        fputs(("s 1\n\no " + name + "\n").c_str(), dest);
        return export_mesh(dest, mesh, index);
    }

    //Modifies the bounds of the slice,
    //translating the mesh if necessary
    void modifyBounds(Ring r) {
        //In the mesh, the normals, texcoords, and indices remain the same
        //The only thing that changes is the list of positions.
        double dx = r.pos.x - bounding_ring.pos.x;
        double dy = r.pos.y - bounding_ring.pos.y;

        for (int i = 0; i < mesh.positions.size() / 3; i++) {
            mesh.positions[3*i] += dx;
            mesh.positions[3*i + 1] += dy;
        }
    }

    //Determines if the given slice needs to be flipped (based on having all parts 2.5d-machined
    //from the positive z axis.
    bool needsFlip() {
        double cumulative_normal_dir = 0;
        //This method is based on computing an cumulative normal (weighted by surface area of triangles)
        for (int i = 0; i < mesh.normals.size(); i += 3) {
            Vec3d normal = Vec3d(mesh.normals[i], mesh.normals[i+1], mesh.normals[i+2]);
            //Assumes the normals are paired with the triangles
            //Get triangle vertices
            int t1 = i;
            int t2 = t1 + 3;
            int t3 = t2 + 3;
            Vec3d v1 = Vec3d(mesh.positions[t1], mesh.positions[t1 + 1], mesh.positions[t1 + 2]);
            Vec3d v2 = Vec3d(mesh.positions[t2], mesh.positions[t2 + 1], mesh.positions[t2 + 2]);
            Vec3d v3 = Vec3d(mesh.positions[t3], mesh.positions[t3 + 1], mesh.positions[t3 + 2]);
            double area = (v1 - v2).cross(v3 - v2).norm() * 0.5;
            cumulative_normal_dir += sgn(normal.z) * fabs(area);
        }
        return cumulative_normal_dir < 0;
    }

    //Fits the mesh within the material bounds
    void fitZ() {
    //In the mesh, the normals + positions must flip.
        //First, compute the maximum and minimum z values over mesh positions
        double min_z = nan("");
        double max_z = nan("");
        for (int i = 0; i < mesh.positions.size() / 3; i++) {
            double z = mesh.positions[3*i+2];
            min_z = fmin(min_z, z);
            max_z = fmax(max_z, z);
        }
        //TODO: throw an exception/return special value if it cannot be done
        double dz = 0 - min_z;
        for (int i = 0; i < mesh.positions.size() / 3; i++) {
            mesh.positions[3*i + 2] += dz;
        }
    }

    //Flips the slice (over z) , maintaining positioning
    //in the material bounding box
    void flip() {
        double centerz = MATERIAL_HEIGHT / 2.0;
        //Flip positions of points
        for (int i = 0; i < mesh.positions.size() / 3; i++) {
            mesh.positions[3*i + 2] = centerz - mesh.positions[3*i + 2];
        }
        //Now, flip triangles (change orientation)
        for (int i = 0; i < mesh.positions.size(); i += 9) {
            swap(mesh.positions, i, i + 3);
            swap(mesh.positions, i + 1, i + 4);
            swap(mesh.positions, i + 2, i + 5);
        }
        //Flip normal z components
        for (int i = 0; i < mesh.normals.size() / 3; i++) {
            mesh.normals[3*i + 2] *= -1.0;
        }
    }
    
    //Flips and fits within the material as necesary
    void flipAndFit() {
        if (needsFlip()) {
            cerr << "FLIP!" << endl;
            flip();
        }
        fitZ();
    }
};


void permute(vector<Vec2d> &in) {
    for (int i = 0; i < in.size(); i++) {
        int j = (rand() % (in.size() - i)) + i;
        swap(in, i, j);
    }
}

//deletes (if you don't care about order)
void del(vector<Vec2d> &in, int index) {
    swap(in, index, in.size() - 1);
    in.pop_back();
}
void del_two(vector<Vec2d> &in, int ind1, int ind2) {
    swap(in, ind1, in.size() - 1);
    swap(in, ind2, in.size() - 2);
    in.pop_back();
    in.pop_back();
}

Ring calc_circle_two(Vec2d one, Vec2d two) {
    Vec2d center = 0.5 * (one + two);
    return Ring(center, 0, fmax(center.dist(one), center.dist(two)));
}

Ring calc_circle_three(Vec2d one, Vec2d two, Vec2d three) {
    LineSeg2d l1(one, two);
    LineSeg2d l2(two, three);
    Vec2d center = l1.bisector().intersect(l2.bisector());
    return Ring(center, 0, center.dist(one));
}

//Indices indicate the indices of non-boundary points in the algorithm below
class Candidate {
    public:
    Ring ring;
    int index1 = -1;
    int index2 = -1;
    bool success;
    Candidate() : ring(Vec2d(0,0),0,0) {
        success = false;
    }
    Candidate(Ring r, int index) : ring(r.pos, r.inner_radius, r.outer_radius) {
        success = false;
        index1 = index;
    }
    Candidate(Ring r, int index1, int index2) : ring(r.pos, r.inner_radius, r.outer_radius) {
        success = false;
        this->index1 = index1;
        this->index2 = index2;
    }
};

//Citation: Emo Welzl
Ring bounding_circle(vector<Vec2d> planePoints) {
    //Special case: only one point
    if (planePoints.size() == 1) {
        return Ring(planePoints[0], 0, 0);
    }
    permute(planePoints);

    vector<Vec2d> support;
    Vec2d v1 = planePoints[0];
    Vec2d v2 = planePoints[1];
    support.push_back(v1);
    support.push_back(v2);
    Ring result = calc_circle_two(v1, v2);

    //Forces a given candidate to be "successful" -- that is, contain all of the supporting points
    //This approach seems more numerically stable than the usual approach.
    auto force_success = [&support](Candidate& in) {
        double r = 0;
        for (Vec2d &p : support) {
            r = fmax(r, p.dist(in.ring.pos));
        }
        in.ring.outer_radius = r;
        in.success = true;
    };

    auto pick_better = [&support, &force_success](Candidate current, Candidate in) -> Candidate {
        force_success(in);
        if (!current.success || in.ring.outer_radius < current.ring.outer_radius) {
            return in;
        }
        return current;
    };
    auto update = [&support, &result](Candidate in) {
        if (in.index2 == -1) {
            del(support, in.index1);
        }
        else {
            del_two(support, in.index1, in.index2);
        }
        result = in.ring;
    };
    auto circle_two = [&support](int i) -> Ring {
        return calc_circle_two(support[i], support[support.size() - 1]);
    };
    auto circle_three = [&support](int i, int j) -> Ring {
        return calc_circle_three(support[i], support[j], support[support.size() - 1]);
    };
    bool finished = false;
    
    while (!finished) {
        bool restart = false;
        for (int i = 0; i < planePoints.size() && !restart; i++) {
            if (!result.fuzz_in(planePoints[i])) {
                //Uh-oh. New support points!
                support.push_back(planePoints[i]);
                Candidate test = Candidate();
                switch (support.size()) {
                    case 2:
                        result = calc_circle_two(support[0], support[1]);
                        break;
                    case 3:
                        //First, try all options with 2 supporting points
                        test = pick_better(test, Candidate(circle_two(0), 1));
                        test = pick_better(test, Candidate(circle_two(1), 0));
                        if (test.success) {
                            update(test);
                            break;
                        }
                        //Didn't work? Must have 3 supporting points.
                        result = calc_circle_three(support[0], support[1], support[2]);
                        break;
                    case 4:
                        //Try all 2-permutations.
                        test = pick_better(test, Candidate(circle_two(0), 1, 2));
                        test = pick_better(test, Candidate(circle_two(1), 0, 2));
                        test = pick_better(test, Candidate(circle_two(2), 0, 1));
                        //Try all 3-permutations.
                        test = pick_better(test, Candidate(circle_three(0, 1), 2));
                        test = pick_better(test, Candidate(circle_three(1, 2), 0));
                        test = pick_better(test, Candidate(circle_three(0, 2), 1));
                        update(test);
                        break;
                }
                restart = true;
            }
        }
        if (!restart) {
            finished = true;
        }
    }
    return result;
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
        if (rot_direction == getDir(p, l2) && rot_direction == getDir(p, l3)) {
            //If the vector from a triangle vertex to a point lies on the same
            //side of every line, it must be in the triangle.
            return true;
        }
    }
    return false;
}

//Gets the bounding ring of a slice that already has been populated
//with bounding circle information
void get_bounding_ring(Slice* slice, vector<LineSeg2d> &lines) {
    Vec2d center = slice->bounding_ring.pos;
    if (intersects(center, lines)) {
        //Must NOT be a ring, do nothing
        return;
    }
    double inner = 0.0;
    for (auto l : lines) {
        inner = fmax(inner, l.distance_to(center));
    }
    slice->bounding_ring.inner_radius = inner;
}

//Gets all 2d line segments of triangles in an order that respects triangles.
vector<LineSeg2d> get_line_segs(Slice *slice) {
    vector<LineSeg2d> ret;
    for (int i = 0; i < slice->mesh.positions.size() / 9; i++) {
        Vec2d p1 = Vec2d(slice->mesh.positions[9*i], slice->mesh.positions[9*i+1]);
        Vec2d p2 = Vec2d(slice->mesh.positions[3+9*i], slice->mesh.positions[3+9*i+1]);
        Vec2d p3 = Vec2d(slice->mesh.positions[6+9*i], slice->mesh.positions[6+9*i+1]);
        ret.push_back(LineSeg2d(p1, p2));
        ret.push_back(LineSeg2d(p2, p3));
        ret.push_back(LineSeg2d(p3, p1));
    }
    return ret;
}

void compute_bounds(Slice *slice) {
    vector<Vec2d> points;
    //TODO: avoid this by overriding Vec2d with projection
    for (int i = 0; i < slice->mesh.positions.size() / 3; i++) {
        points.push_back(Vec2d(slice->mesh.positions[3*i+0], slice->mesh.positions[3*i+1]));
    }
    slice->bounding_ring = bounding_circle(points);

    vector<LineSeg2d> lineSegs = get_line_segs(slice);

    get_bounding_ring(slice, lineSegs);
}

//TODO: Ensure d_min kept with the list of corner placements to avoid re-calculation.
double min_dist(Ring c, Shape &u, Shape &v, vector<pair<Slice*, Ring>> &cfg, vector<LineSeg2d> &edges) {
    //d_min will hold the distance from the ring's center to the closest shape in the
    //configuration (including sides), but exclude u and v
    double d_min = nan("");
    for (auto l : edges) {
        if (l != u && l != v) {
            d_min = fmin(d_min, l.distance_to(c.pos));
        }
    }
    for (pair<Slice*, Ring> p : cfg) {
        Ring r = p.second;
        if (r != u && r != v) {
            d_min = fmin(d_min, r.distance_to(c.pos));
        }
    }
    d_min = fmax(0, d_min - c.outer_radius);

    return d_min;
}

double hole_degree(Ring c, double d_min) {
    return (1 - (d_min / c.outer_radius));
}

class Placement {
    public:
    Placement(int index, Vec2d pos_init, double d_min) : pos(pos_init.x, pos_init.y) {
        this->index = index;
        this->d_min = d_min;
    }
    int index;
    Vec2d pos;
    double d_min;
};

vector<double> quad_solve(double a, double b, double c) {
    vector<double> result;
    double discrim = b*b-4*a*c;
    if (eps_equals(discrim, 0)) {
        result.push_back(-b / (2*a));
    }
    else if (discrim > 0) {
        result.push_back((-b + sqrt(discrim)) / 2*a);
        result.push_back((-b - sqrt(discrim)) / 2*a);
    }
    return result;
}

vector<Vec2d> perp(vector<Vec2d> in) {
    vector<Vec2d> result;
    for (int i = 0; i < in.size(); i++) { 
        result.emplace_back(in[i].perp());
    }
    return result;
}

vector<Vec2d> get_tangent_circle_center(Ring one, Ring two, double r) {
    cerr << one.to_string() + two.to_string() + "," + to_string(r) << endl;
    Vec2d c = two.pos - one.pos;
    if (eps_equals(c.x, 0) && eps_equals(c.y, 0)) {
        return vector<Vec2d>();
    }
    double d = c.norm();
    c = c.normalize();

    double r1 = one.outer_radius;
    double r2 = two.outer_radius;
    double cos = (sq(r + r1) - sq(r + r2) + sq(d)) / (2.0 * d * (r1 + r));
    double sin = sqrt(1 - sq(cos));
    cerr << "cos = " << cos << " sin = " << sin << " d = " << d << endl;
    Mat2d rotMat1 = Mat2d(cos, sin, -sin, cos);
    Mat2d rotMat2 = Mat2d(cos, -sin, sin, cos);
    c = c * (r1 + r);
    vector<Vec2d> result;
    result.push_back(one.pos + rotMat1 * c);
    result.push_back(one.pos + rotMat2 * c);
    return result;
}

/*
//Gets the center (x) of a circle tangent to the two given ones, with radius r.
vector<Vec2d> get_tangent_circle_center(Ring one, Ring two, double r) {
    cerr << one.to_string() + two.to_string() + "," + to_string(r) << endl;
    Vec2d x1 = one.pos;
    Vec2d x2 = two.pos;
    double r1 = one.outer_radius;
    double r2 = two.outer_radius;
    //First, compute the line the tangent circle must lie on, of the form <X, c> = b
    Vec2d c = x1 - x2;
    if (eps_equals(c.x, 0) && eps_equals(c.y, 0)) {
        //Special case: dealing with the same ring
        return vector<Vec2d>();
    }

    double b = (x1.dot(x1) - x2.dot(x2) + sq(r + r2) - sq(r + r1)) / 2.0;
    cerr << to_string(b) + ",";
    Line2d l(c, b);
    //Now, compute the closest point to x1 on the line and the distance of x1 to the line
    Vec2d p = l.closest_point_to(x1);
    cerr << to_string(p.x) + "," + to_string(p.y) << endl;
    double d = x1.dist(p);
    //Then p, x1, x defines a right triangle with side lengths d, h, r + r1
    double h = sqrt(sq(r + r1) - sq(d));
    Vec2d t = l.tangent().normalize();
    vector<Vec2d> result;
    result.push_back(p + t*h);
    result.push_back(p - t*h);
    return result;
}
*/

/*
//Gets the center of a circle tangent to the two given ones, with radius r.
//To do this, it uses a (horribly ugly) series of formulas.
vector<Vec2d> get_tangent_circle_center(Ring one, Ring two, double r) {
    cerr << one.to_string() + two.to_string() + "," + to_string(r) << endl;
    double r1 = one.outer_radius;
    double r2 = two.outer_radius;
    double x1 = one.pos.x;
    double x2 = two.pos.x;
    double y1 = one.pos.y;
    double y2 = two.pos.y;
    double ydiff = y2 - y1;
    
    if (eps_equals(ydiff, 0)) {
        if (eps_equals(x1 - x2, 0)) {
            return vector<Vec2d>();
        }
    }
   if (fabs(ydiff) <= fabs(x1 - x2)) {
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
    vector<Vec2d> result(x.size(), Vec2d(0,0));
    std::transform(x.begin(), x.end(), result.begin(), [&m, &k](double x) { return Vec2d(x, m*x + k); });
    for (Vec2d &v : result) {
        cerr << "Points: " + to_string(v.x) + "," + to_string(v.y) << endl;
    }
    return result;
}
*/

bool placement_valid(Vec2d pos, double r, vector<pair<Slice*, Ring>> &cfg) {
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

Placement pos_to_placement(int index, Vec2d p, double r, Shape &u, Shape &v, 
                           vector<pair<Slice*, Ring>> &cfg, vector<LineSeg2d> &edges) {
    return Placement(index, p, min_dist(Ring(p, 0, r), u, v, cfg, edges));
}

void add_center(int index, Vec2d center, double r, Shape &one, Shape &two, 
        vector<pair<Slice*, Ring>> &cfg, vector<LineSeg2d> &edges, vector<Placement> &result) {
    if (placement_valid(center, r, cfg)) {
        result.push_back(pos_to_placement(index, center, r, one, two, cfg, edges));
    }
}


//Gets the center of a circle tangent to a circle and a horizontal line, with radius r
vector<Vec2d> get_tangent_circle_center_horiz(Ring one, double y_h, double r) {
    double delta_x = 2 * sqrt(one.outer_radius * r);
    vector<Vec2d> result;
    result.emplace_back(one.pos.x + delta_x, y_h + r);
    result.emplace_back(one.pos.x - delta_x, y_h + r);
    result.emplace_back(one.pos.x + delta_x, y_h - r);
    result.emplace_back(one.pos.x - delta_x, y_h - r);
    return result;
}
//Same as before. TODO: Make this not as dumb
vector<Vec2d> get_tangent_circle_center_vert(Ring one, double x_h, double r) {
    return perp(
        get_tangent_circle_center_horiz(Ring(one.pos.antiperp(), one.inner_radius, one.outer_radius),
                                           x_h, r));
}


vector<Placement> gen_placements_involving(Ring one, int index, double r, vector<pair<Slice*, Ring>> &cfg,
                                           vector<LineSeg2d> &edges) {
    vector<Placement> result;
    
    auto add_centers = [&index, &cfg, &r, &edges, &result](vector<Vec2d> centers, Shape &one, Shape &two) {
        for (Vec2d center : centers) { 
            add_center(index, center, r, one, two, cfg, edges, result); 
            //cerr << "Placement at " << center.x << " " <<  center.y << "\n";
        }
    };

    //Handle circle-circle case
    for (int j = 0; j < cfg.size(); j++) {
        //TODO: Can this check be optimized by storing pairwise distances with the config?
        //Make sure that a corner placement is remotely possible
        Ring two = cfg[j].second;
        if (one.distance_to(two) > 2*r) {
            continue;
        }
        else {
            //cout << "CIRCLES" << "\n";
            add_centers(get_tangent_circle_center(one, two, r), one, two);
        }
    }
    
    //cout << "CIRCLE-LINES" << "\n";
    //Now, circle-line case
    add_centers(get_tangent_circle_center_horiz(one, 0, r), one, edges[BOTTOM]);
    add_centers(get_tangent_circle_center_vert(one, 0, r), one, edges[LEFT]);
    add_centers(get_tangent_circle_center_horiz(one, MATERIAL_LENGTH, r), one, edges[TOP]);
    add_centers(get_tangent_circle_center_vert(one, MATERIAL_WIDTH, r), one, edges[RIGHT]);
    return result;

}

//Generate all possible placements for a given index
vector<Placement> gen_placements(int index, double r, vector<pair<Slice*, Ring>> &cfg, vector<LineSeg2d> &edges) {
    vector<Placement> result;

    auto add = [&index, &r, &cfg, &edges, &result](Vec2d center, Shape &one, Shape &two) {
        add_center(index, center, r, one, two, cfg, edges, result);
    };

    //Compute circle-circle and circle-line corner placements  
    for (int i = 0; i < cfg.size(); i++) {
        Ring one = cfg[i].second;
        vector<Placement> temp = gen_placements_involving(one, index, r, cfg, edges);
        result.insert(result.end(), temp.begin(), temp.end());
        temp.clear();
    }

    //Finally, test the line-line placements
    add(Vec2d(r, r), edges[BOTTOM], edges[LEFT]);
    add(Vec2d(MATERIAL_WIDTH - r, r), edges[BOTTOM], edges[RIGHT]);
    add(Vec2d(MATERIAL_WIDTH - r, MATERIAL_LENGTH - r), edges[RIGHT], edges[TOP]);
    add(Vec2d(r, MATERIAL_LENGTH - r), edges[TOP], edges[LEFT]);

    return result;
}


//Citation: This is a modified version of http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.99.5620&rep=rep1&type=pdf
vector<pair<Slice*, Ring>> ring_pack(vector<Slice*> in) {
    vector<pair<Slice*, Ring>> config;
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

    //TODO: I don't like any of the above, or the below.
    //Here, we set the (constant) indices for edges
    BOTTOM = 0;
    RIGHT = 1;
    TOP = 2;
    LEFT = 3;
    

    auto radiusCompare = [](Slice *a, Slice *b) {
        return a->bounding_ring.outer_radius >= b->bounding_ring.outer_radius;
    };

    //Sorted in order of decreasing radius
    std::sort(in.begin(), in.end(), radiusCompare);

    //Instead of what's done in the paper, we arbitrarily pick an initial
    //configuration with the two smallest rings at opposite ends of the rectangle.
    //This reduces the time complexity from O(n^10) to O(n^8), at the possible loss of some packing efficiency.
    Slice* s1 = pop_back(in);
    Slice* s2 = pop_back(in);
    double r1 = s1->bounding_ring.outer_radius;
    double r2 = s2->bounding_ring.outer_radius;
    config.emplace_back(s1, s1->bounding_ring.move_to(Vec2d(r1, r1)));
    config.emplace_back(s2, s2->bounding_ring.move_to(Vec2d(MATERIAL_WIDTH - r2, MATERIAL_LENGTH - r2)));
    
    //Generate initial placements
    for (int i = 0; i < in.size(); i++) {
        vector<Placement> temp = gen_placements(i, in[i]->bounding_ring.outer_radius, config, edges);
        placements.insert(placements.end(), temp.begin(), temp.end());
    }

    auto addPlacement = [&in, &config, &placements, &edges](Placement accepted) {
        Slice *s = in[accepted.index];
        Ring ci = s->bounding_ring.move_to(accepted.pos);

        cerr << "Circle center: " + to_string(ci.pos.x) + "," + to_string(ci.pos.y) << endl;
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
            Ring p_bounds = in[p.index]->bounding_ring.move_to(p.pos);
            //Ignore placements that would overlap ci
            if (ci.intersects(p_bounds)) {
                continue;
            }
            //Update d_min for all remaining placements
            p.d_min = fmin(p.d_min, ci.distance_to(p_bounds));

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
        for (int i = 0; i < in.size(); i++) {
            Slice *s = in[i];
            vector<Placement> temp = gen_placements_involving(ci, i, s->bounding_ring.outer_radius, 
                                                              config, edges);
            placements.insert(placements.end(), temp.begin(), temp.end());
            temp.clear();
        }
    };

    while (placements.size() > 0) {
        //Stores the index into the corner placements with the maximum lambda
        int max_index = 0;
        int max_lambda = 0;
        //Find the best corner placement
        for (int i = 0; i < placements.size(); i++) {
            double lambda = hole_degree(in[placements[i].index]->bounding_ring, placements[i].d_min);
            if (lambda > max_lambda) {
                max_index = i;
                max_lambda = lambda;
            }
        }
        addPlacement(placements[max_index]);
    }
    return config;
}

//Gets a slice number from a name of the form "Slice_slice-1_Part_0000_"... 
//(As generated by 123D Make)
int getSliceNum(string name) {
    string numStart = name.substr(12);
    int i;
    std::istringstream(numStart) >> i;
    return i;
}

string getSliceName(string filename, string nameinfile) {
    return filename + "_" + to_string(getSliceNum(nameinfile));
}
    

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cout << "Usage: 3d-slice-layout sourcepath (optional: destpath) \n";
        return 1;
    }
    path source_dir(argv[1]);
    if (!exists(source_dir) || !is_directory(source_dir)) {
        cerr << "The path input must be a directory!\n";
        return 1;
    }
    FILE* deststream = stdout;
    if (argc == 3) {
        deststream = fopen(argv[2], "w");
        if (deststream == NULL) {
            cerr << "Cannot open output file!" << endl;
        }
    }

    //We need to bunch together meshes belonging to the same slice in the same file
    //This is a dumb way to do it.
    unordered_map<string, vector<tinyobj::shape_t>> shapes;

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
            string filename = entry.path().filename().string();
            for (tinyobj::shape_t &shape : temp_shapes) {
                string newName = getSliceName(filename, shape.name);
                shape.name = newName;
                if (shapes.count(newName) == 0) {
                    shapes[newName] = vector<tinyobj::shape_t>();
                }
                shapes[newName].push_back(shape);
            }
        }
    }

    //Now, iterate over all of the buckets to get a final list of Slices
    vector<Slice*> slices = vector<Slice*>();
    for (auto kv : shapes) {
        string name = kv.first;
        vector<tinyobj::shape_t> slice_shapes = kv.second;
        tinyobj::mesh_t combined_mesh = standardize_mesh(slice_shapes[0].mesh);
        for (int i = 1; i < slice_shapes.size(); i++) {
            merge_meshes(combined_mesh, standardize_mesh(slice_shapes[i].mesh)); 
        }
        if (combined_mesh.positions.size() == 0) {
            continue;
        }
        slices.push_back(new Slice(name, combined_mesh));
    }

    //Compute the bounds of each material, and flip and fit each one within the material
    for (Slice* slice : slices) {
        slice->flipAndFit();
        compute_bounds(slice);
        cerr << slice->bounding_ring.to_string() << endl;
    }

    vector<pair<Slice*, Ring>> packed = ring_pack(slices);
    
    //Generate the circle chart and translate the underlying mesh.
    cout << "<svg xmlns = \"http://www.w3.org/2000/svg\" verstion=\"1.1\">" << endl;
    for (auto sr : packed) {
        Vec2d c = sr.second.pos;
        //Write out the circle to SVG
        cout << "<circle cx=\"" << to_string(c.x) << "\" cy=\"" << to_string(c.y) << "\" r=\"";
        cout << to_string(sr.second.outer_radius) << "\" fill=\"orange\" />" << endl;
        
        //Write out labels to SVG
        cout << "<text x=\"" << to_string(c.x - sr.second.outer_radius);
        cout << "\" y=\"" << to_string(c.y) << "\" font-size=\"15\">";
        cout << sr.first->name << endl;
        cout << "</text>" << endl;

        //Translate the mesh
        sr.first->modifyBounds(sr.second);
        cerr << "Stupid point: " << sr.first->mesh.positions[0] << ", " <<
            sr.first->mesh.positions[1] << endl;
    }
    cout << "</svg>" << endl;
    
    //Great. Now everything should be fitted, and in the right position. Write .obj file format to stdout.
    int index = 0;
    for (Slice* slice : slices) {
        index = slice->export_obj(deststream, index);
        delete slice;
    }

    if (argc == 3) {
        fclose(deststream);
    }

    return 0;
}



