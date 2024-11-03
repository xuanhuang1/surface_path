# surface_path


this program takes two 2D raw elevation value files with matching dimensions. The program runs interactively, taking 2 end points (4 integers for grid indices) and computes the surface distance

to compile:
```
mkdir build
cd build
cmake .. 
make
```

to run:
```
./surface_path pre.data post.data -d x y [options]
 -h  Print this help 
 [options]
 -ppath print indivual points
```