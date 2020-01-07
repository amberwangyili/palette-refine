# refined geometric approach for palette-based image decomposition

Code for the paper "An improved Geometric Approach for Palette-based Image Decomposition and Recoloring"

## Requirements

------

The following setup steps are developed and tested on MacOS High Sierra 10.13.6 and with minimum requirements of [CMake](<https://cmake.org/download/>) 3.12 and C++ 14 standard.

**Used Libs:** 

- [`opencv`](<https://github.com/opencv/opencv>)
- [`nlopt`](<https://nlopt.readthedocs.io/en/latest/#download-and-installation>)

## Build with CMake

------

```bash
$ cd <palette-refine>
$ cmake .
$ make
```



## Usage

directly use a toy example

```bash
$ ./palette_refine
```

or

```bash
$./palette_refine [OPTION...] [pic] [obj] [prefix]

  -p, --pic arg     original picture
  -o, --obj arg     convex hull .obj
  -f, --prefix arg  prefix
  -s, --sample arg  number of samples
  -c, --option arg  center point selection method
  -r, --ratio arg   ratio between neighbor_point and random points
  -l, --lambda arg  lambda
  -i, --iter arg    number of iteration
  -u, --unique arg  unique parent of a point
  -h, --help        Print help
```

