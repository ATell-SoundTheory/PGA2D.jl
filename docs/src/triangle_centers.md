# Triangle centers

This page collects constructions for classic triangle centers using PGA2D.

```@meta
CurrentModule = PGA2D
```

## Setup

```julia
P1, P2, P3 = point(0,0), point(1,0), point(0,1)
```

## Incenter

Intersection of internal angle bisectors:

```julia
I = incenter_ppp(P1, P2, P3)
point_coordinates(I)
```

Expected for the right unit triangle: `(1 - √2/2, 1 - √2/2)`.

## Circumcenter

Intersection of perpendicular bisectors:

```julia
Cc = circumcenter_ppp(P1, P2, P3)
point_coordinates(Cc)  # (0.5, 0.5)
```

## Circumcircle

```julia
Cc2, R = circumcircle_ppp(P1, P2, P3)
point_coordinates(Cc2), R
```

## Incircle

```julia
Ic, r = incircle_ppp(P1, P2, P3)
point_coordinates(Ic), r
```

## Orthocenter

Intersection of altitudes:

```julia
H = orthocenter_ppp(P1, P2, P3)
point_coordinates(H)   # (0.0, 0.0)
```

## Centroid (for comparison)

Intersection of medians:

```julia
Cg = meet_ll(join_pp(P3, p_bisect_pp(P1,P2)), join_pp(P1, p_bisect_pp(P2,P3)))
point_coordinates(Cg)   # (1/3, 1/3)
```

## Plotting demo (optional)

The plotting recipe is provided via an optional extension and loads automatically when Plots is available on Julia ≥ 1.9.
This example is not doctested or executed during CI, keeping builds deterministic.

```julia
using PGA2D, Plots

# Triangle vertices and centers
P1, P2, P3 = point(0,0), point(1,0), point(0,1)
I  = incenter_ppp(P1, P2, P3)
Cc = circumcenter_ppp(P1, P2, P3)
H  = orthocenter_ppp(P1, P2, P3)
Cg = meet_ll(join_pp(P3, p_bisect_pp(P1,P2)), join_pp(P1, p_bisect_pp(P2,P3)))

# Triangle sides
L12, L23, L31 = join_pp(P1,P2), join_pp(P2,P3), join_pp(P3,P1)

plot(legend=:outerright, aspect_ratio=:equal, size=(500,400))
scatter!([P1,P2,P3], label="Vertices")
plot!([L12,L23,L31], label=["L12" "L23" "L31"])  # recipe handles arrays

scatter!([I], label="Incenter")
scatter!([Cc], label="Circumcenter")
scatter!([H], label="Orthocenter")
scatter!([Cg], label="Centroid")
# Optionally save the figure locally:
# savefig("triangle_centers.png")
```

