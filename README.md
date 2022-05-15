# GAMES101-Assignment02 

Assignment 2 requires rasterizing triangles. 

In the process of completion, I need to create a bounding box, sample the internal pixels, and finally determine the pixel colors through Z-buffer and frame-buffer. Of course, I also completed the improvement task, using the super-sampling algorithm to anti-aliasing on triangles.

I'll start with a review of what I've learned in class.

**Pixel: Picture of element**

​		A pixel is a little square with color(red, green, blue);

​		When we use the viewport matrix to transform the object in the XY plane, the pixel's indices are from (0,0) to (width-1, height -1). And the Pixel(x,y) is located at (x+0.5, y+0.5); 


**Sampling: Determine whether the pixel is in the triangle**

   Check if each pixel center is inside the triangle.

​	Using the triangle vertices, use the right-hand rule to determine if the pixel is in the triangle.  

**Bounding Box: reduce the pixels that need to be checked**


​		The bounding box is determined according to triangle vertices.



**Rasterization: Drawing to Raster Displays**

​	We need to reduce the effect of Aliasing caused by sampling.

​	**Antialiased Sampling: undersampling will introduce aliasing**

- Filter before sampling (Pre-filtering -> Sampling)
- Filtering: Get rid of the frequency content we don't need.
- Sampling: Repeating frequency Contents.

​	**Antialiasing By Averaging Values in Pixel Area:** 

- Convolve f(x,y) by a 1-pixel box-blur. 
- Then sample at every pixel.



**Supersampling(MSAA?):** 

​		Approximate the effect of the 1-pixel box filter by sampling multiple locations within a pixel and averaging their values.

1. Take N*N samples in each pixel. (Increase sampling rate)
2. Average the N*N samples inside each pixel.



**Expand Contents: **

1. The cost of MSAA: 
2. FXAA(Fast Approximate AA): 
3. TAA(Temporal AA):  



Back to the Assignment2.

First of all, we need to do the MVP transformation just like assignment 1. As requested, I filled the function **get_projection_matrix ()**.

```CPP
Eigen::Matrix4f get_projection_matrix(float eye_fov, float aspect_ratio,
    float zNear, float zFar)
{
    // Define Matrix
    Eigen::Matrix4f projection = Eigen::Matrix4f::Identity();
    Eigen::Matrix4f M_persp2ortho = Eigen::Matrix4f::Identity();
    Eigen::Matrix4f M_ortho = Eigen::Matrix4f::Identity();//orthographic Matrix

    float angle = eye_fov * MY_PI / 180.0; // half fov angle

    auto t = -zNear * tan(angle / 2);//use fov calculate h
    auto r = t * aspect_ratio;// use aspect calculate r
    auto l = -r;
    auto b = -t;

    M_persp2ortho << zNear, 0, 0, 0,
        0, zNear, 0, 0,
        0, 0, zNear + zFar, -zNear * zFar,
        0, 0, 1, 0;

    M_ortho << 2 / (r - l), 0, 0, -(r + l) / 2,
        0, 2 / (t - b), 0, -(t + b) / 2,
        0, 0, 2 / (zNear - zFar), -(zNear + zFar) / 2,
        0, 0, 0, 1;

    projection = M_ortho * M_persp2ortho * projection;
    return projection;
}
```

Next, I need to create a Bounding Box in the function called **rasterize_triangle(const Triangle& t)**:

The input parameter t is an instance of the triangle class.

```CPP
class Triangle{
public:
    Vector3f v[3]; /*the original coordinates of the triangle, v0, v1, v2 in counter clockwise order*/
    /*Per vertex values*/
    Vector3f color[3]; //color at each vertex;
    Vector2f tex_coords[3]; //texture u,v
    Vector3f normal[3]; //normal vector for each vertex
}
```

We just need the vertex data there. (v[3])，the upper and lower bounds of triangle vertices are the Bounding Box scopes.

```CPP
 //create the Bounding Box
    float x_min = std::min(t.v[0].x(), std::min(t.v[1].x(), t.v[2].x()));
    float x_max = std::max(t.v[0].x(), std::max(t.v[1].x(), t.v[2].x()));
    float y_min = std::min(t.v[0].y(), std::min(t.v[1].y(), t.v[2].y()));
    float y_max = std::max(t.v[0].y(), std::max(t.v[1].y(), t.v[2].y()));
    int xmin = std::floor(x_min);
    int xmax = std::ceil(x_max);
    int ymin = std::floor(y_min);
    int ymax = std::ceil(y_max);
```

Next,   I should fill the function **"insideTriangle()"**to check if each pixel is inside the triangle.

```CPP
static bool insideTriangle(int x, int y, const Vector3f* _v)
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
    bool flag = (ans1[2] >= 0 && ans2[2] >= 0 && ans3[2] >= 0) || (ans1[2] < 0 && ans2[2] < 0 && -ans3[2] < 0);
    return flag;
}
```

Next, I need to compare the depth value and depth buffer of each pixel to set the color.

```CPP
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
                    set_pixel(Vector3f(i, j, 0), t.getColor());
                }
 }
```

Two triangles should display correctly now.

![QQ图片20220515153350](https://user-images.githubusercontent.com/68177870/168466393-1970cf08-d27c-40af-8535-acdb4c153dd6.png)


We need to use MSAA for anti-aliasing. Create a new function :

```CPP
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
    std::vector<float> dis{0.25,0.25,0.75,0.75,0.25};

    float cnt;//Record the num in pixel
    float z_buffer=INT_MAX;

    for (int i = xmin; i < xmax; i++) {
        for (int j = ymin; j < ymax; j++) {
            cnt = 0;
            //Four-part of one pixel
            for (int d = 0; d < 4; d++) {
                if (insideTriangle(i + dis[d], j + dis[d+1], t.v)) {
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
            if (cnt > 0 && z_buffer< depth_buf[get_index(i, j)]) {
                //Update depth buffer
                depth_buf[get_index(i, j)] = z_buffer;
                //set color,According to the number of cnt
                set_pixel(Vector3f(i, j, 0), t.getColor()*cnt);
            }
        }
    }
}

```

Divide a pixel into four pixels for interpolation. Sets the color based on the number of pixels in the triangle.

Let's compare the effects: 

Base:
![image](https://user-images.githubusercontent.com/68177870/168466448-5a209696-f0c3-4e49-9579-f5011d6950cf.png)



MSAA:

![image](https://user-images.githubusercontent.com/68177870/168466556-f59e9528-6470-441a-b1df-215256a4f1ed.png)


Obviously, MSAA is better, but it adds a lot of running time. We can see a black line in the MSAA result ,but I haven't found a solution yet.

