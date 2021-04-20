# Geometry Processing - Poisson Surface Reconstruction

This is my implementation of Poisson Surface Reconstruction in [CSC419/CSC2520 Geometry Processing](https://github.com/alecjacobson/geometry-processing-csc2520/).

## Build & Execution
```
git submodule update --init --recursive
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make 
./mesh-reconstruction [path to point cloud]
```

Press `p` to show the point cloud
![](images/elephant-points-normals.jpg)
Press `m` to show the point cloud
![](images/elephant-mesh.jpg)

