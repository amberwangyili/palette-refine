#include <iostream>
#include <nlopt.h>
#include <nlopt.hpp>
#include <opencv2/opencv.hpp>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <random>
#include <stdio.h>
#include <time.h>
#include "vec3.h"
#include "nlopt.h"
#include "nearestPoint.h"
#include "cxxopt.h"

using namespace std;
using namespace cv;
using namespace nlopt;



float outsidehull_points_distance(const Mesh &mesh, const vector<vec3> &outsidepoints)
{
    float total_distance = 0;
    for (int i = 0; i < outsidepoints.size(); ++i) {
        vec3 nearestPoint = nearest_point(outsidepoints[i], mesh.faces, mesh.vertices);
        total_distance += (nearestPoint - outsidepoints[i]).norm();
    }
    return total_distance;

}

void compute_outside_points(const Mesh& mesh, const vector<vec3>& points, vector<vec3> &outside_points, vector<vec3> &inside_points)
{
    vector<Triangle> triangles;
    mesh.constructTriangles(triangles);
    outside_points.reserve(points.size());
    inside_points.reserve(points.size());

    vector<bool> isInside;
    isInside.resize(points.size());

    for (int l = 0; l < points.size(); ++l) {
        const vec3 & pt = points[l];
        isInside[l] = Mesh::isInside(triangles, pt);
    }

    for (int l = 0; l < points.size(); ++l) {
        if (isInside[l]) {
            inside_points.push_back(points[l]);
        } else
            outside_points.push_back(points[l]);
    }
}

vector<int> sort_indexes(const vector<float> & v) {

    // initialize original index locations
    vector<int> idx(v.size());
    iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    sort(idx.begin(), idx.end(),
         [&v](int i1, int i2) {return v[i1] < v[i2]; });

    return idx;
}

vector<int> sort_indexes_n(const vector<float> & v, int n_elements) {

    // initialize original index locations
    vector<int> idx(v.size());
    iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    nth_element(idx.begin(), idx.begin() + n_elements, idx.end(),
                [&v](int i1, int i2) {return v[i1] < v[i2]; });

    return idx;
}

vec3 compute_specific_center_point(const Mesh& mesh, const vector<vec3>& inside_points, int point_num, int vertex_index) {
    vector<float> dist(inside_points.size());
    vec3 vertex = mesh.vertices[vertex_index];

    for (int i = 0; i < inside_points.size(); ++i) {
        dist[i] = (inside_points[i] - vertex).norm();
    }

    vector<int> indx = sort_indexes_n(dist, point_num);

    vec3 center;
    for (int k = 0; k < point_num; ++k) {
        center += inside_points[indx[k]];
    }
    center = vec3(center.x() / (float)point_num, center.y() / (float)point_num, center.z() / (float)point_num);
    return center;
}



vector<vec3> compute_center_points(const Mesh& mesh, const vector<vec3>& inside_points, int point_num) {

    vector<vec3> center_points(mesh.vertex_num());

    for (int l = 0; l < mesh.vertex_num(); ++l) {
        center_points[l] = compute_specific_center_point(mesh, inside_points, point_num, l);
    }
    return center_points;
}


pair<int,double> find_nearest(const Mesh& mesh, const vec3& point){
    double dist = MAXFLOAT;
    int parent;

    for(int i = 0 ; i< mesh.vertex_num(); i++){
        if( dist >(mesh.vertices[i] - point ).norm()){
            parent = i;
            dist = (mesh.vertices[i] - point ).norm();
        }
    }
    return make_pair(parent,dist);
}

struct point_dist{
    vec3 point;
    double dist;
    point_dist ( vec3 point_ , double dist_): point(point_) , dist(dist_) {}
    bool operator <(const point_dist &x) const{
        return dist < x.dist;
    }
};

vector<vec3> compute_center_points_unique(const Mesh& mesh, const vector<vec3>& inside_points, int point_num) {


    vector<vec3> center_points(mesh.vertex_num());
    vector<vector<point_dist>> included_points(mesh.vertex_num(),vector<point_dist>());

    int parent;
    double dist;

    for(auto &it:inside_points){

        pair<int,double> res = find_nearest(mesh,it);
        parent = res.first;
        dist = res.second;
        included_points[parent].push_back(point_dist(it,dist));

    }


    for(int i = 0; i<mesh.vertex_num(); i++){
        vec3 center;
        int counter = 0;
        nth_element(included_points[i].begin(),included_points[i].begin()+point_num,included_points[i].end());
        for(auto &it : included_points[i] ){
            center += (it).point;
            counter ++;
            if(counter >= point_num)break;
        }
        center = vec3(center.x() / (float)point_num, center.y() / (float)point_num, center.z() / (float)point_num);
        center_points[i] = center;
    }
    return center_points;
}



struct Data {
    vector<vec3> center_point;
    int indx;
    Mesh mesh;
    vector<vec3> points;
    float lambda;
    int point_num;
    int option;
    int unique;
};



Mat load_image(const string & img_filepath) {
    Mat img = imread(img_filepath);
    if (img.rows == 0 || img.cols == 0)
    {
        cout << "load image from" << img_filepath << "failes (file not exist).\n";
    }
    return img;
}


struct Result {
    float lambda;
    float total_error;
    float reconstruct_error;
    float representative_error;
};



void compute_loss_function_base(const Mesh& mesh, const vector<vec3>& outside_points, int n_total_point,
                                const vector<vec3>& centerNNPoints, float lambda, Result& result) {

    float reconstruct_error = outsidehull_points_distance(mesh, outside_points) / n_total_point;
    float represent_error = 0;
    for (int i = 0; i < mesh.vertex_num(); i++) {
        represent_error += (mesh.vertices[i] - centerNNPoints[i]).norm();
    }
    represent_error /= mesh.vertex_num();

    float total_error = reconstruct_error * lambda + represent_error;

    result.total_error = total_error;
    result.lambda = lambda;
    result.reconstruct_error = reconstruct_error;
    result.representative_error = represent_error;
}


void compute_loss_function(const Mesh& mesh, const vector<vec3>& points, float lambda, int nn_point_num, int option,int unique, Result& result)
{
    vector<vec3> insidepoint, outsidepoint;
    compute_outside_points(mesh, points, outsidepoint, insidepoint);
    vector<vec3> centerpoints;

    if(unique == 1){
        if (option == 1) {
            int center_point_num = nn_point_num < (int)insidepoint.size()  ? nn_point_num : (int)insidepoint.size();
            centerpoints = compute_center_points_unique(mesh, insidepoint, center_point_num);

        }
        else {
            centerpoints = compute_center_points_unique(mesh, points, nn_point_num);
        }

    }else{
        if (option == 1) {
            int center_point_num = nn_point_num < (int)insidepoint.size()  ? nn_point_num : (int)insidepoint.size();
            centerpoints = compute_center_points(mesh, insidepoint, center_point_num);
        }
        else {
            centerpoints = compute_center_points(mesh, points, nn_point_num);
        }
    }


    compute_loss_function_base(mesh, outsidepoint, points.size(), centerpoints, lambda, result);
}

double cost_func(const std::vector<double> & x, std::vector<double> & grad, void* data)
{

    const Data* dt = reinterpret_cast<Data*>(data);
    double target_dist = x[0];
    vec3 target_point = (1 - target_dist) * dt->center_point[dt->indx] + target_dist * dt->mesh.vertices[dt->indx];

    Mesh newmesh = dt->mesh;
    newmesh.vertices[dt->indx] = target_point;

    Result result;
    {
        vector<vec3> insidepoint, outsidepoint;

        compute_outside_points(newmesh, dt->points, outsidepoint, insidepoint);
        int nn_point_num = dt->point_num;
        vector<vec3> centerpoints = dt->center_point;

        if(dt->unique == 1){
            if (dt->option == 1) {
                int center_point_num = nn_point_num < (int)insidepoint.size() ? nn_point_num : (int)insidepoint.size();

                centerpoints[dt->indx] = compute_center_points_unique(newmesh, insidepoint, center_point_num)[dt->indx];
            }
            else {
                centerpoints[dt->indx] = compute_center_points_unique(newmesh, dt->points, nn_point_num)[dt->indx];
            }
        }
        else{
            if (dt->option == 1) {
                int center_point_num = nn_point_num < (int)insidepoint.size() ? nn_point_num : (int)insidepoint.size();
                centerpoints[dt->indx] = compute_specific_center_point(newmesh, insidepoint, center_point_num,dt->indx);
            }
            else {
                centerpoints[dt->indx] = compute_specific_center_point(newmesh, dt->points, nn_point_num, dt->indx);
            }
        }

        compute_loss_function_base(newmesh, outsidepoint, dt->points.size(), centerpoints, dt->lambda, result);

     }


    return result.total_error;

}


int main(int argc, char* argv[])
{


    int random_sample_num_per_side = 200;
    float neighbor_point_ratio = 0.001;
    float lambda = 20;
    int iteration = 2;
    int option = 1;
    int unique = 0;



    std::string pic_path = "./image/02.png",obj_path = "./obj_original/02-mesh-obj-files.obj",prefix = "02";
    try
    {
        cxxopts::Options options("Palette Refinement", "Refine color palette of an image");
        options
                .positional_help("[pic] [obj] [prefix]")
                .show_positional_help();
        options.add_options()
                ("p, pic", "original picture", cxxopts::value<std::string>(pic_path))
                ("o, obj", "convex hull .obj", cxxopts::value<std::string>(obj_path))
                ("f, prefix", "prefix", cxxopts::value<std::string>(prefix))
                ("s, sample", "number of samples", cxxopts::value<int>(random_sample_num_per_side))

                ("c, option", "center point selection method", cxxopts::value<int>(option))

                ("r, ratio", "ratio between neighbor_point and random points", cxxopts::value<float>(neighbor_point_ratio))
                ("l, lambda", "lambda", cxxopts::value<float>(lambda))
                ("i, iter", "number of iteration", cxxopts::value<int>(iteration))
                ("u, unique", "unique parent of a point", cxxopts::value<int>(unique))
                ("h, help", "Print help");

        options.parse_positional({ "pic", "obj", "prefix" });
        auto parseResult = options.parse(argc, argv);
        if (parseResult.count("h"))
        {
            std::cout << options.help() << std::endl;
            exit(0);
        }
    }
    catch (const cxxopts::OptionException& e)
    {
        std::cout << "error parsing options: " << e.what() << std::endl;
        exit(1);
    }

    std::cout << "image:" << pic_path << "\n";
    std::cout << "ori palette:" << obj_path << "\n";

    Mesh mesh;
    bool bLoadedMesh = mesh.load_from_file(obj_path);
    if (!bLoadedMesh) {
        std::cout << "Load from " << obj_path << " failed.\n";
        return 0;
    }
    Mat img = load_image(pic_path);
    if (img.cols == 0) {
        std::cout << "Load from " << pic_path << " failed.\n";
        return 0;
    }

    vector<vec3> points;
    std::random_device rd;
    auto seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 g(seed);
    std::uniform_real_distribution<float> dis(0.0f, 0.99f);
    // stratified sample image
    int sample_width = std::min<int>(random_sample_num_per_side, img.cols);
    int sample_height = std::min<int>(random_sample_num_per_side, img.rows);
    float ratio_width = static_cast<float>(img.cols) / static_cast<float>(sample_width);
    float ratio_hegiht = static_cast<float>(img.rows) / static_cast<float>(sample_height);
    //std::cout << "ratio_width = " << ratio_width << ", ratio_height = " << ratio_hegiht << std::endl;
    for (int i = 0; i < sample_height; ++i)
        for (int j = 0; j < sample_width; ++j)
        {
            int x = std::min<int>(std::floor((i + dis(g)) * ratio_hegiht), img.rows - 1);
            int y = std::min<int>(std::floor((j + dis(g)) * ratio_width), img.cols - 1);
            //std::cout << "(x, y) = " << "(" << x << ", " << y << ")" << std::endl;
            const uchar* pointer = img.ptr<uchar>(x);
            vec3 pt(pointer[3 * y + 2] / 255.0f, pointer[3 * y + 1] / 255.0f, pointer[3 * y] / 255.0f);
            points.push_back(std::move(pt));
        }


    int point_num = std::max<int>(neighbor_point_ratio * sample_width * sample_height, 2);
    Mesh refind_palette = mesh;

    for (int j = 0; j < iteration;j++)
    {

        for (int i = 0; i < mesh.vertex_num(); i++)
        {
            vector<vec3> centerpoints;
            if(unique ==1){
                centerpoints = compute_center_points_unique(refind_palette, points, point_num);
            }else{
                centerpoints = compute_center_points(refind_palette, points, point_num);
            }

            int indx = i;
            Data dt = { centerpoints,indx,refind_palette,points,lambda,point_num,option,unique };
            float initialguess = 0;
            double minval = 999999;

            double initial_x[11];
            double cost[11];

            for (int k = 0; k <= 10; k++)
            {
                initial_x[k] = 0.1f * k;
                vector<double> x = { initial_x[k] };
                cost[k] = cost_func(x, x, &dt);
                //cout << cost_func(x, x, &dt) << endl;
            }

            for (int k = 0; k <= 10; k++) {
                if (minval > cost[k]) {
                    minval = cost[k];
                    initialguess = initial_x[k];
                }
            }

            std::vector<double> x = { initialguess };

            nlopt::opt opt(nlopt::LN_COBYLA, 1);   // algorithm, n
            opt.set_lower_bounds(0);
            opt.set_upper_bounds(1);
            opt.set_min_objective(cost_func, &dt);
            opt.set_xtol_rel(1e-4);

            double minf;
            nlopt::result result = opt.optimize(x, minf);
            refind_palette.vertices[dt.indx] = (1 - x[0]) * dt.center_point[dt.indx] + x[0] * dt.mesh.vertices[dt.indx];
        }
    }

    if(unique == 1){
        refind_palette.save_to_file("./obj_refine_unique/" + prefix + "-2-unique-mesh-obj-files.obj");
    }else{
        refind_palette.save_to_file("./obj_refine/" + prefix + "-2-mesh-obj-files.obj");
    }


    Result result;
    compute_loss_function(refind_palette, points, lambda, point_num,option,unique,result);

    std::cout << "Total Loss = " << result.total_error << "\n";
    std::cout << "lambda = " << result.lambda << "\n";
    std::cout << "Reconstrcut Loss = " << result.reconstruct_error << "\n";
    std::cout << "Sparse Loss = " << result.representative_error<< "\n";





    return 0;
}


