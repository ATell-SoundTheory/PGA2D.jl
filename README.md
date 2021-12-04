# PGA2D.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ATell-SoundTheory.github.io/PGA2D.jl/dev)
[![Build Status](https://github.com/ATell-SoundTheory/PGA2D.jl/workflows/CI/badge.svg)](https://github.com/ATell-SoundTheory/PGA2D.jl/actions)

This Julia package implements a two-dimensional plane-based geometric algebra on top of [CliffordAlgebras.jl](https://github.com/ATell-SoundTheory/CliffordAlgebras.jl). The multivectors of the algebra can encode points, lines, directions and motors. This package provides a number of named operations on these objects to make using the algebra easier and more intuitive. In addition, we also offer a [Plots.jl](http://docs.juliaplots.org/latest/) recipe for visualising the supported geometric primitives.

Create points, lines and directions.
```
julia> using PGA2D

julia> point(0,0)
point(0.0, 0.0) ≜ +1×e1e2 ∈ Cl(2, 0, 1)

julia> line(1,1,0)
line(0.7071067811865475, 0.7071067811865475, 0.0) ≜ +1×e1+1×e2 ∈ Cl(2, 0, 1)

julia> direction(2,1)
direction(2, 1) ≜ +1×e0e1+2×e2e0 ∈ Cl(2, 0, 1)
```

Line coordinates are given as the coefficients $(a,b,c)$ in the line equation $a x + b y + c = 0$. Coordinates can be extracted from MultiVectors:

```
julia> point_coordinates(point(-2,1))
(-2.0, 1.0)

julia> line_coordinates(line(-1,1,-3))
(-1, 1, -3)

julia> direction_coordinates(direction(1,0))
(1, 0)
```

However, if the MultiVector does not encode the object we are asking for, we will receive a DomainError:

```
julia> direction_coordinates(line(-1,1,-3))
ERROR: DomainError with line(-0.7071067811865475, 0.7071067811865475, -2.1213203435596424) ≜ -1×e1+1×e2-3×e0 ∈ Cl(2, 0, 1):
 is not a direction.
```

We can make sure we have the right kind of object by using the following methods:

```
julia> p = point(2,1)
point(2.0, 1.0) ≜ +1×e1e2+1×e0e1+2×e2e0 ∈ Cl(2, 0, 1)

julia> is_point(p)
true

julia> is_line(p)
false

julia> is_direction(p)
false
```

If you would like to work with the MultiVectors directly, you can use the exported symbols `pga2d` for the global algebra instance and `PGA2DMV` for the MultiVector type. The basis vectors for Clifford(2,0,1) are also exported as `e0`, `e1`, `e2`, `e01`, `e12`, `e20`, `e012` or alternatively `𝐈`. For full support of the underlying Clifford Algebra, you can import `CliffordAlgebras.jl`.

The recommended use of this package is to rely on the available higher-level functions that act on points, directions and lines. Because the type of all objects is that of the MultiVector `PGA2DMV`, the geometric functions can not check for the proper type. You can provide a line in place of a point and get a MultiVector as a result. To make sure you are aware of the expected object, the nomenclature of the functions is designed to indicate both the kind of object returned and the expected argument types.

The function `l_ortho_lp` returns a line, as indicated by the leading `l`. The line is orthogonal to the line and passes through the point provided as arguments and described by the `lp` suffix:

```
julia> p = point(1,1)
point(1.0, 1.0) ≜ +1×e1e2+1×e0e1+1×e2e0 ∈ Cl(2, 0, 1)

julia> l = line(1,1,0)
line(0.7071067811865475, 0.7071067811865475, 0.0) ≜ +1×e1+1×e2 ∈ Cl(2, 0, 1)

julia> l_ortho_lp(l,p)
line(-0.7071067811865475, 0.7071067811865475, 0.0) ≜ -1×e1+1×e2 ∈ Cl(2, 0, 1)
```

There are many more functions that allow us to intersect two lines, find the line passing through to two points, find distances and angles, bisect angles and distances, calculate areas and transform objects using motors. Functions that are missing should be easily composable from the functions in this package.


