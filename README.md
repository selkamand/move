
<!-- README.md is generated from README.Rmd. Please edit that file -->

# move

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/move)](https://CRAN.R-project.org/package=move)
<!-- badges: end -->

> \[!WARNING\]  
> This package is in early development and not yet ready for use

The move R package provides functions for basic 3d & multidimensional
transformation of points, vectors, and planes.

## Installation

You can install the development version of move like so:

``` r
if (!require("remotes"))
    install.packages("remotes")

remotes::install_github("selkamand/move")
```

## Quick Start

### Rotate a vector around an axis (Rodrigues)

``` r
library(move)

v <- c(1, 0, 0)
axis <- c(0, 0, 1)
angle <- pi/2  # 90 degrees in radians

rotate_vector_around_axis(v, axis, angle)
#> [1] 0 1 0
```

### Rotate a vector to align with a target direction

``` r
v <- c(2, -3, 4)    # any vector
target <- c(1,  1, 1)    # direction we want

v_rot <- rotate_vector_to_align_with_target(v, target)
sqrt(sum(v^2))     # original length
#> [1] 5.385165
sqrt(sum(v_rot^2)) # preserved length
#> [1] 5.385165
```

### Rotate a vector into a plane (not projection)

``` r
v <- c(1, 2, 3)
n <- c(0, 1, -1)     # plane normal

v_in_plane <- rotate_vector_into_a_plane(v, n)
# v_in_plane lies in the plane; ||v|| is preserved (unlike a projection)
```

### Project Vectors

``` r
# Vector projection of a onto b
project_vector_into_vector(c(2, 3, 4), c(1, 0, 0))
#> [1] 2 0 0

# Scalar projection (signed length of a along b)
compute_scalar_projection(c(2, 3, 4), c(1, 0, 0))
#> [1] 2
```

### Translate positions

``` r
translate_position_by_vector(c(1, 2, 3), c(0.5, 0, -1))
#> [1] 1.5 2.0 2.0
#> c(1.5, 2, 2)

translate_position_in_direction(c(0, 0, 0), c(1, 1, 0), 10)
#> [1] 7.071068 7.071068 0.000000
#> approximately c(7.071068, 7.071068, 0)

compute_translation_vector(c(1, 2, 3), c(4, 6, 3))
#> [1] 3 4 0
#> c(3, 4, 0)
```

### Measuring

``` r
a <- c(1, 0, 0)
b <- c(0,1,0)

measure_angle_between_vectors(a, b)
#> [1] 1.570796

a %angle% b
#> [1] 1.570796
```

### Conversion

``` r
radians_to_degrees(pi)
#> [1] 180
degrees_to_radians(180)
#> [1] 3.141593
```

### Plane representation conversions

Planes can be represented in two equivalent forms:

1.  **Point–normal form**  
    Defined by a 3D point `point` that lies on the plane and a normal
    vector `normal`.  
    Equation:  
    $$$
    \[
    \hat{n} \cdot (x - p) = 0
    \]
    $$\$

2.  **Normal–offset (Hesse) form**  
    Defined by a **unit** normal vector `normal` and a **signed offset**
    (distance) `offset` from the origin along that normal.  
    Equation:  
    $$$
    \[
    \hat{n} \cdot x = s
    \]
    $$\$ where `s` = `offset`.

These forms can be converted with the following functions:

``` r
# Point–normal -> Normal–offset
x <- convert_plane_point_normal_to_normal_offset(normal = c(0,0,2), point = c(2,3,5))
x$normal  # c(0,0,1)
#> [1] 0 0 1
x$offset  # 5
#> [1] 5

# Normal–offset -> Point–normal (closest point to origin)
convert_plane_normal_offset_to_point(normal = c(0,0,1), offset = 5)
#> $point
#> [1] 0 0 5
#> 
#> $normal
#> [1] 0 0 1
```
