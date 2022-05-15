#include <algorithm>
#include <vector>
#include "rasterizer.hpp"
#include <opencv2/opencv.hpp>
#include <math.h>


rst::pos_buf_id rst::rasterizer::load_positions(const std::vector<Eigen::Vector3f>& positions)
{
    auto id = get_next_id();
    pos_buf.emplace(id, positions);

    return { id };
}

rst::ind_buf_id rst::rasterizer::load_indices(const std::vector<Eigen::Vector3i>& indices)
{
    auto id = get_next_id();
    ind_buf.emplace(id, indices);

    return { id };
}

rst::col_buf_id rst::rasterizer::load_colors(const std::vector<Eigen::Vector3f>& cols)
{
    auto id = get_next_id();
    col_buf.emplace(id, cols);

    return { id };
}

auto to_vec4(const Eigen::Vector3f& v3, float w = 1.0f)
{
    return Vector4f(v3.x(), v3.y(), v3.z(), w);
}


static bool insideTriangle(float x, float y, const Vector3f* _v)
{
    // TODO : Implement this function to check if the point (x, y) is inside the triangle represented by _v[0], _v[1], _v[2]

    // check if each pixel is inside the triangle.
    Eigen::Vector3f Point = Vector3f(x, y, 1);// each point

    //Three sides of a triangle
    Eigen::Vector3f tri_side1 = _v[1] - _v[0];
    Eigen::Vector3f tri_side2 = _v[2] - _v[1];
    Eigen::Vector3f tri_side3 = _v[0] - _v[2];

    //The result of the cross product
    Eigen::Vector3f ans1 = tri_side1.cross(Point - _v[0]);
    Eigen::Vector3f ans2 = tri_side2.cross(Point - _v[1]);
    Eigen::Vector3f ans3 = tri_side3.cross(Point - _v[2]);

    //whether the results are in the same direction
    bool flag = (ans1[2] >= 0 && ans2[2] >= 0 && ans3[2] >= 0)
        || (ans1[2] < 0 && ans2[2] < 0 && -ans3[2] < 0);
    return flag;
}

static std::tuple<float, float, float> computeBarycentric2D(float x, float y, const Vector3f* v)
{
    float c1 = (x * (v[1].y() - v[2].y()) + (v[2].x() - v[1].x()) * y + v[1].x() * v[2].y() - v[2].x() * v[1].y()) / (v[0].x() * (v[1].y() - v[2].y()) + (v[2].x() - v[1].x()) * v[0].y() + v[1].x() * v[2].y() - v[2].x() * v[1].y());
    float c2 = (x * (v[2].y() - v[0].y()) + (v[0].x() - v[2].x()) * y + v[2].x() * v[0].y() - v[0].x() * v[2].y()) / (v[1].x() * (v[2].y() - v[0].y()) + (v[0].x() - v[2].x()) * v[1].y() + v[2].x() * v[0].y() - v[0].x() * v[2].y());
    float c3 = (x * (v[0].y() - v[1].y()) + (v[1].x() - v[0].x()) * y + v[0].x() * v[1].y() - v[1].x() * v[0].y()) / (v[2].x() * (v[0].y() - v[1].y()) + (v[1].x() - v[0].x()) * v[2].y() + v[0].x() * v[1].y() - v[1].x() * v[0].y());
    return { c1,c2,c3 };
}

void rst::rasterizer::draw(pos_buf_id pos_buffer, ind_buf_id ind_buffer, col_buf_id col_buffer, Primitive type)
{
    auto& buf = pos_buf[pos_buffer.pos_id];
    auto& ind = ind_buf[ind_buffer.ind_id];
    auto& col = col_buf[col_buffer.col_id];

    float f1 = (50 - 0.1) / 2.0;
    float f2 = (50 + 0.1) / 2.0;

    Eigen::Matrix4f mvp = projection * view * model;
    for (auto& i : ind)
    {
        Triangle t;
        Eigen::Vector4f v[] = {
                mvp * to_vec4(buf[i[0]], 1.0f),
                mvp * to_vec4(buf[i[1]], 1.0f),
                mvp * to_vec4(buf[i[2]], 1.0f)
        };
        //Homogeneous division
        for (auto& vec : v) {
            vec /= vec.w();
        }
        //Viewport transformation
        for (auto& vert : v)
        {
            vert.x() = 0.5 * width * (vert.x() + 1.0);
            vert.y() = 0.5 * height * (vert.y() + 1.0);
            vert.z() = vert.z() * f1 + f2;
        }

        for (int i = 0; i < 3; ++i)
        {
            t.setVertex(i, v[i].head<3>());
            t.setVertex(i, v[i].head<3>());
            t.setVertex(i, v[i].head<3>());
        }

        auto col_x = col[i[0]];
        auto col_y = col[i[1]];
        auto col_z = col[i[2]];

        t.setColor(0, col_x[0], col_x[1], col_x[2]);
        t.setColor(1, col_y[0], col_y[1], col_y[2]);
        t.setColor(2, col_z[0], col_z[1], col_z[2]);

        rasterize_triangle_MSAA(t);
    }
}

//Screen space rasterization
void rst::rasterizer::rasterize_triangle(const Triangle& t) {
    auto v = t.toVector4();

    // TODO : Find out the bounding box of current triangle.
    // iterate through the pixel and find if the current pixel is inside the triangle
    //create the Bounding Box
    float x_min = std::min(t.v[0].x(), std::min(t.v[1].x(), t.v[2].x()));
    float x_max = std::max(t.v[0].x(), std::max(t.v[1].x(), t.v[2].x()));
    float y_min = std::min(t.v[0].y(), std::min(t.v[1].y(), t.v[2].y()));
    float y_max = std::max(t.v[0].y(), std::max(t.v[1].y(), t.v[2].y()));
    int xmin = std::floor(x_min);
    int xmax = std::ceil(x_max);
    int ymin = std::floor(y_min);
    int ymax = std::ceil(y_max);

    for (int i = xmin; i < xmax; i++) {
        for (int j = ymin; j < ymax; j++) {
            if (insideTriangle(i + 0.5, j + 0.5, t.v)) {
                //Computational depth interpolation
                auto [alpha, beta, gamma] = computeBarycentric2D(i + 0.5, j + 0.5, t.v); //Compute interpolated centroids
                float w_reciprocal = 1.0 / (alpha / v[0].w() + beta / v[1].w() + gamma / v[2].w());
                float z_interpolated = alpha * v[0].z() / v[0].w() + beta * v[1].z() / v[1].w() + gamma * v[2].z() / v[2].w();
                z_interpolated *= w_reciprocal; //z-buffer

                if (z_interpolated < depth_buf[get_index(i, j)]) {
                    //Update depth buffer
                    depth_buf[get_index(i, j)] = z_interpolated;
                    //set color
                    set_pixel(Vector3f(i, j, 1), t.getColor());
                }
            }
        }
    }
}

//use MSAA to rasterization
void rst::rasterizer::rasterize_triangle_MSAA(const Triangle& t) {
    auto v = t.toVector4();

    //create the Bounding Box
    float x_min = std::min(t.v[0].x(), std::min(t.v[1].x(), t.v[2].x()));
    float x_max = std::max(t.v[0].x(), std::max(t.v[1].x(), t.v[2].x()));
    float y_min = std::min(t.v[0].y(), std::min(t.v[1].y(), t.v[2].y()));
    float y_max = std::max(t.v[0].y(), std::max(t.v[1].y(), t.v[2].y()));
    int xmin = std::floor(x_min);
    int xmax = std::ceil(x_max);
    int ymin = std::floor(y_min);
    int ymax = std::ceil(y_max);

    //MSAA list
    std::vector<float> dis{ 0.25,0.25,0.75,0.75,0.25 };

    float cnt;//Record the num in pixel
    float z_buffer = INT_MAX;

    for (int i = xmin; i < xmax; i++) {
        for (int j = ymin; j < ymax; j++) {
            cnt = 0;
            //Four-part of one pixel
            for (int d = 0; d < 4; d++) {
                if (insideTriangle(i + dis[d], j + dis[d + 1], t.v)) {
                    //Computational depth interpolation
                    auto [alpha, beta, gamma] = computeBarycentric2D(i + dis[d], j + dis[d + 1], t.v); //Compute interpolated centroids
                    float w_reciprocal = 1.0 / (alpha / v[0].w() + beta / v[1].w() + gamma / v[2].w());
                    float z_interpolated = alpha * v[0].z() / v[0].w() + beta * v[1].z() / v[1].w() + gamma * v[2].z() / v[2].w();
                    z_interpolated *= w_reciprocal; //z-buffer

                    z_buffer = std::min(z_buffer, z_interpolated);
                    cnt += 0.25;
                }
            }
            //if this pixel cnt >0 ,set color
            if (cnt > 0 && z_buffer < depth_buf[get_index(i, j)]) {
                //Update depth buffer
                depth_buf[get_index(i, j)] = z_buffer;
                //set color,According to the number of cnt
                set_pixel(Vector3f(i, j, 0), t.getColor() * cnt);
            }
        }
    }
}

void rst::rasterizer::set_model(const Eigen::Matrix4f& m)
{
    model = m;
}

void rst::rasterizer::set_view(const Eigen::Matrix4f& v)
{
    view = v;
}

void rst::rasterizer::set_projection(const Eigen::Matrix4f& p)
{
    projection = p;
}

void rst::rasterizer::clear(rst::Buffers buff)
{
    if ((buff & rst::Buffers::Color) == rst::Buffers::Color)
    {
        std::fill(frame_buf.begin(), frame_buf.end(), Eigen::Vector3f{ 0, 0, 0 });
    }
    if ((buff & rst::Buffers::Depth) == rst::Buffers::Depth)
    {
        std::fill(depth_buf.begin(), depth_buf.end(), std::numeric_limits<float>::infinity());
    }
}

rst::rasterizer::rasterizer(int w, int h) : width(w), height(h)
{
    frame_buf.resize(w * h);
    depth_buf.resize(w * h);
}

int rst::rasterizer::get_index(int x, int y)
{
    return (height - 1 - y) * width + x;
}

void rst::rasterizer::set_pixel(const Eigen::Vector3f& point, const Eigen::Vector3f& color)
{
    //old index: auto ind = point.y() + point.x() * width;
    auto ind = (height - 1 - point.y()) * width + point.x();
    frame_buf[ind] = color;

}