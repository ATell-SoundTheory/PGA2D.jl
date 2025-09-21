```@meta
CurrentModule = PGA2D
```

# PGA2D

PGA2D implements the 2D plane-based geometric algebra Cl(2,0,1) on top of CliffordAlgebras.jl.
It provides high-level constructors and operations for points, lines, directions, and motors.

## Quick start

```julia
using PGA2D

p = point(0, 0)
q = point(1, 0)
r = point(0, 1)

l = join_pp(p, q)         # line through p and q
m = meet_ll(join_pp(q,r), join_pp(p,r))  # their intersection (point r)

dist = dist_pp(p, q)      # 1.0
θ = angle_ll(join_pp(p,q), join_pp(p,r)) # π/2
```

Coordinates helpers:

```julia
point_coordinates(point(-2, 1))         # (-2.0, 1.0)
direction_coordinates(direction(1, 0))  # (1, 0)
line_coordinates(line(1, 1, 0))         # (≈0.707, ≈0.707, 0.0)

# Non-throwing variants return `nothing` if the type doesn’t match
try_point_coordinates(direction(1,0))   # nothing
```

Distance and angle helpers:

```julia
dist_lp(join_pp(p,q), r)  # distance from point r to line pq
angle_dd(direction(1,0), direction(0,1))  # π/2
```

### Triangle centers

Incenter (intersection of internal angle bisectors) and centroid (intersection of medians):

```julia
P1, P2, P3 = point(0,0), point(1,0), point(0,1)
I = incenter_ppp(P1, P2, P3)
Cg = meet_ll(join_pp(P3, p_bisect_pp(P1,P2)), join_pp(P1, p_bisect_pp(P2,P3)))
```

## Plotting (optional)

PGA2D ships a Plots.jl recipe as an optional extension.
If Plots is present in the environment, the recipe is loaded automatically on Julia ≥ 1.9.

```julia
using Plots
plot([p, q, r, l]; aspect_ratio=:equal)
```

## Index and API

```@index
```

```@autodocs
Modules = [PGA2D]
```
