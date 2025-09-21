module PGA2D

using CliffordAlgebras

import Base.join, Base.show

export pga2d, PGA2DT, PGA2DMV
export e0,e1,e2, e01, e12, e20, e012, 𝐈
export direction, point, line
export point_coordinates
export direction_coordinates
export line_coordinates
export try_point_coordinates
export try_direction_coordinates
export try_line_coordinates
export is_point, is_direction
export is_line, is_motor
export normalize
export join_pp
export meet_ll
export l_ortho_lp
export p_proj_lp
export l_para_lp
export d_ortho_l
export dist_pp
export dist_ll
export dist_lp
export angle_ll
export angle_dl
export angle_dd
export angle_ll_signed
export angle_dl_signed
export angle_dd_signed
export dist_orient_lp
export l_bisect_ll
export l_bisect_pp
export p_bisect_pp
export reflect_l
export area_ppp
export area_loop
export length_loop
export motor_pa
export motor_ll
export move_pa
export move_ll
export move_m
export line_pd
export incenter_ppp
export circumcenter_ppp
export orthocenter_ppp
export try_incenter_ppp
export try_circumcenter_ppp
export try_orthocenter_ppp
export incircle_ppp
export circumcircle_ppp

const pga2d = CliffordAlgebra(:PGA2D)
const PGA2DT = typeof(pga2d)
const PGA2DMV = MultiVector{PGA2DT}

const e0 = pga2d.e0
const e1 = pga2d.e1
const e2 = pga2d.e2

const e01 = pga2d.e0e1
const e12 = pga2d.e1e2
const e20 = pga2d.e2e0

const e012 = pga2d.e1e2e0
const 𝐈 = e012



"""
    point(x::Real, y::Real)

Constructs a MultiVector that encodes a point (x,y).
"""
point(x::Real, y::Real) = direction(x, y) + e12

"""
    direction(x::Real, y::Real)

Constructs a MultiVector that encodes a direction (x,y).
"""
direction(x::Real, y::Real) = x * e20 + y * e01

"""
    line(a::Real, b::Real, c::Real)

Constructs a MultiVector that encodes a line a·x + b·y + c = 0.
"""
line(a::Real, b::Real, c::Real) = a * e1 + b * e2 + c * e0

# Convenience constructors
point(t::Tuple{<:Real,<:Real}) = point(t[1], t[2])
direction(t::Tuple{<:Real,<:Real}) = direction(t[1], t[2])
line(t::Tuple{<:Real,<:Real,<:Real}) = line(t[1], t[2], t[3])

"""
    point_coordinates(P::PGA2DMV)

Extracts the point coordinates (x,y) from a MultiVector P. If the MultiVector does not encode a point a DomainError is thrown.
"""
function point_coordinates(P::PGA2DMV)
    w = P.e1e2
    if isgrade(P,2) && !iszero(w)
        wi = inv(w)
        return wi.*(P.e2e0, P.e0e1)
    else
        throw(DomainError(P, " is not a point."))
    end
end

"""
    direction_coordinates(d::PGA2DMV)

Extracts the direction coordinates (x,y) from a MultiVector P. If the MultiVector does not encode a direction a DomainError is thrown.
"""
function direction_coordinates(d::PGA2DMV)
    w = d.e1e2
    if isgrade(d,2) && iszero(w)
        return (d.e2e0, d.e0e1)
    else
        throw(DomainError(d, " is not a direction."))
    end
end

"""
    line_coordinates(l::PGA2DMV)

Extracts the line coordinates (a,b,c) from a MultiVector l. If the MultiVector does not encode a line a DomainError is thrown.
"""
function line_coordinates(l::PGA2DMV)
    if isgrade(l,1) && (!iszero(l.e1) || !iszero(l.e2))
        return (l.e1, l.e2, l.e0)
    else
        throw(DomainError(l, " is not a line."))
    end
end

"""
    is_point(P::PGA2DMV)

Returns true iff the MultiVector P encodes a point.
"""
function is_point(P::PGA2DMV)
    isgrade(P,2) && !iszero(P.e1e2)
end

"""
    is_direction(d::PGA2DMV)

Returns true iff the MultiVector d encodes a direction.
"""
function is_direction(d::PGA2DMV)
    isgrade(d,2) && iszero(d.e1e2) && !iszero(d)
end

"""
    is_line(l::PGA2DMV)

Returns true iff the MultiVector l encodes a line.
"""
function is_line(l::PGA2DMV)
    isgrade(l,1) && (!iszero(l.e1) || !iszero(l.e2))
end

"""
    is_motor(m::PGA2DMV)

Returns true iff the MultiVector m encodes a motor.
"""
function is_motor(m::PGA2DMV)
    iszero(m.e1) && iszero(m.e2) && iszero(m.e0) && iszero(m.e1e2e0) && !iszero(m)
end


function show(io::IO, x::PGA2DMV)
    if is_point(x)
        print(io, "point", point_coordinates(x), " ≜ ")
    elseif is_direction(x)
        print(io, "direction", direction_coordinates(x), " ≜ ")
    elseif is_line(x)
        print(io, "line", line_coordinates(normalize(x)), " ≜ ")
    elseif is_motor(x)
        m = normalize(x)
        P = grade(m,2)
        if is_point(P)
            orientation = sign(P.e1e2)
            α = 2 * orientation * atan(norm(P), scalar(m))
            print(io, "motor: rotation $(α/π)π around $P", " ≜ ")
        elseif is_direction(P)
            d = 2 * sqrt(P.e0e1^2 + P.e2e0^2)
            d_ortho = (P.e2e0,-P.e0e1).*(2/d)
            print(io, "motor: translation $d direction $(d_ortho)", " ≜ ")
        end
    end
    show_multivector(io, x)
end

"""
    normalize(x::PGA2DMV)

Scales the MultiVector x so that it has unit norm. Throws a DomainError for zero-norm inputs.
"""
function normalize(x::PGA2DMV)
    n = norm(x)
    if n == 0
        throw(DomainError(x, " cannot be normalized (zero norm)."))
    end
    x / n
end

"""
    try_point_coordinates(P::PGA2DMV) -> Union{Tuple{<:Real,<:Real},Nothing}

Non-throwing variant of `point_coordinates`. Returns `(x,y)` or `nothing` if `P` is not a point.
"""
function try_point_coordinates(P::PGA2DMV)
    w = P.e1e2
    if isgrade(P,2) && !iszero(w)
        wi = inv(w)
        return wi.*(P.e2e0, P.e0e1)
    end
    nothing
end

"""
    try_direction_coordinates(d::PGA2DMV) -> Union{Tuple{<:Real,<:Real},Nothing}

Non-throwing variant of `direction_coordinates`. Returns `(x,y)` or `nothing` if `d` is not a direction.
"""
function try_direction_coordinates(d::PGA2DMV)
    w = d.e1e2
    if isgrade(d,2) && iszero(w) && !iszero(d)
        return (d.e2e0, d.e0e1)
    end
    nothing
end

"""
    try_line_coordinates(l::PGA2DMV) -> Union{Tuple{<:Real,<:Real,<:Real},Nothing}

Non-throwing variant of `line_coordinates`. Returns `(a,b,c)` or `nothing` if `l` is not a line.
"""
function try_line_coordinates(l::PGA2DMV)
    if isgrade(l,1) && (!iszero(l.e1) || !iszero(l.e2))
        return (l.e1, l.e2, l.e0)
    end
    nothing
end

"""
    meet_ll(l1::PGA2DMV, l2::PGA2DMV)

Calculates the point of intersection of the two line-encoding MultiVectors l1 and l2.
If the lines are parallel the result is a direction parallel to the lines. If the lines coincides, the result is 0.
"""
meet_ll(l1::PGA2DMV, l2::PGA2DMV) = l1 ∧ l2

"""
    join_pp(P1::PGA2DMV, P2::PGA2DMV)

Calculates the line through the points P1 and P2. If P1 and P2 coincide, the result is 0.
"""
join_pp(P1::PGA2DMV, P2::PGA2DMV) = P1 ∨ P2

"""
    l_ortho_lp(l::PGA2DMV, P::PGA2DMV)

Calculates the line that is orthogonal to l and passes through P.
"""
l_ortho_lp(l::PGA2DMV, P::PGA2DMV) = l ⋅ P

"""
    p_proj_lp(l::PGA2DMV, P::PGA2DMV)

Orthogonally projects the point P onto the line l.
"""
p_proj_lp(l::PGA2DMV, P::PGA2DMV) = (l ⋅ P) * l

"""
    l_para_lp(l::PGA2DMV, P::PGA2DMV)

Constructs a parallel line to the line l in the point P.
"""
l_para_lp(l::PGA2DMV, P::PGA2DMV) = (l ⋅ P ) * P

"""
    d_ortho_l(l::PGA2DMV)

Calculates the direction orthogonal to l.
"""
d_ortho_l(l::PGA2DMV) = l * 𝐈

"""
    dist_pp(P1::PGA2DMV, P2::PGA2DMV)

Calculates the Euclidean distance between P1 and P2.
"""
dist_pp(P1::PGA2DMV, P2::PGA2DMV) = norm(P1 ∨ P2) / (norm(P1) * norm(P2))

"""
    dist_ll(l1::PGA2DMV, l2::PGA2DMV)

Calculates the orthogonal Euclidean distance between parallel lines `l1` and `l2`.

Notes
- This quantity is well-defined as the perpendicular distance only when the two lines are parallel.
- For almost-parallel lines with a small angle φ between them, `dist_ll` equals the parallel offset multiplied by `cos(φ)`. Hence, as the lines deviate from perfect parallelity, the returned value decreases by a factor `≈ 1 - φ^2/2` (second-order in φ).
- For non-parallel lines (intersecting), this is not the minimal point-to-line distance (which is 0); it returns the component of the offset along the common normal as if the lines were treated as nearly parallel. No attempt is made to detect parallelism.

Example
Let `l1: y = 0` and `l2: y = d + m x`. Then `l1 = line(0, 1, 0)` and `l2 = line(m, -1, d)`. One finds
`dist_ll(l1, l2) = |d| / √(1 + m^2) = |d| cos(φ)`, where `φ = atan(m)` is the angle between the two lines.
Thus for `m → 0` (nearly parallel) the value approaches `|d|` with a relative error `≈ m^2/2`.
"""
dist_ll(l1::PGA2DMV, l2::PGA2DMV) = norm(dual(l1 ∧ l2)) / (norm(l1) * norm(l2))

"""
    angle_ll(l1::PGA2DMV, l2::PGA2DMV)

Calculates the angle between two lines l1 and l2.
"""
angle_ll(l1::PGA2DMV, l2::PGA2DMV) = begin
    s = norm(l1 ∧ l2) / (norm(l1) * norm(l2))
    asin(clamp(s, -1.0, 1.0))
end

"""
    angle_dl(d::PGA2DMV, l2::PGA2DMV)

Calculates the angle between the direction d and the line l.
"""
angle_dl(d::PGA2DMV, l::PGA2DMV) = begin
    s = norm(dual(d ∧ l)) / (norm(d) * norm(l))
    asin(clamp(s, -1.0, 1.0))
end

"""
    angle_dl_signed(d::PGA2DMV, l::PGA2DMV)

Signed angle from direction `d` to the tangent of line `l` in (-π, π]. The line's tangent is chosen as (-b, a) for `l = a e1 + b e2 + c e0`.
"""
function angle_dl_signed(d::PGA2DMV, l::PGA2DMV)
    (dx, dy) = direction_coordinates(d)
    (a, b, _) = line_coordinates(l)
    # Tangent direction of the line: base (b, -a) so that for y=0 (a=0,b=1) tangent is +x.
    tx, ty = b, -a
    # Choose orientation to align tangent with the reference direction d
    if dx * tx + dy * ty < 0
        tx = -tx; ty = -ty
    end
    dot = dx * tx + dy * ty
    cross = dx * ty - dy * tx
    atan(cross, dot)
end

"""
    dist_lp(l::PGA2DMV, P::PGA2DMV)

Unsigned Euclidean distance between line `l` and point `P`.
"""
dist_lp(l::PGA2DMV, P::PGA2DMV) = abs(dist_orient_lp(l, P))

"""
    angle_dd(d1::PGA2DMV, d2::PGA2DMV)

Angle between directions `d1` and `d2`.
Computes from Euclidean coordinates to avoid null-norm issues.
"""
function angle_dd(d1::PGA2DMV, d2::PGA2DMV)
    (x1, y1) = direction_coordinates(d1)
    (x2, y2) = direction_coordinates(d2)
    dot = x1 * x2 + y1 * y2
    cross = x1 * y2 - y1 * x2
    # atan2-style: robust and returns principal angle in [0, π]
    atan(abs(cross), dot)
end

"""
    angle_dd_signed(d1::PGA2DMV, d2::PGA2DMV)

Signed angle from direction `d1` to `d2` in (-π, π], using `atan2` of the 2D cross/dot of Euclidean coordinates.
"""
function angle_dd_signed(d1::PGA2DMV, d2::PGA2DMV)
    (x1, y1) = direction_coordinates(d1)
    (x2, y2) = direction_coordinates(d2)
    dot = x1 * x2 + y1 * y2
    cross = x1 * y2 - y1 * x2
    atan(cross, dot)
end

"""
    dist_orient_lp(l::PGA2DMV, P::PGA2DMV)

Calculates the oriented Euclidean distance between the line l and the point P.
"""
dist_orient_lp(l::PGA2DMV, P::PGA2DMV) = scalar(P ∨ l) / (norm(P) * norm(l))

"""
    l_bisect_ll(l1::PGA2DMV, l2::PGA2DMV)

Calculates both angle-bisecting lines between l1 and l2 and return them in a tuple of MultiVectors.
"""
function l_bisect_ll(l1::PGA2DMV, l2::PGA2DMV)
    l1n = normalize(l1)
    l2n = normalize(l2)
    (l1n + l2n, l1n - l2n)
end

"""
    angle_ll_signed(l1::PGA2DMV, l2::PGA2DMV)

Signed angle from line `l1` to `l2` in (-π, π], defined via their outward normals.
"""
function angle_ll_signed(l1::PGA2DMV, l2::PGA2DMV)
    (a1, b1, _) = line_coordinates(l1)
    (a2, b2, _) = line_coordinates(l2)
    dot = a1 * a2 + b1 * b2
    cross = a1 * b2 - b1 * a2
    atan(cross, dot)
end

"""
    l_bisect_pp(P1::PGA2DMV, P2::PGA2DMV)

Calculates the orthogonal bisecting line between the points P1 and P2.
"""
function l_bisect_pp(P1::PGA2DMV, P2::PGA2DMV) 
    P1n = normalize(P1)
    P2n = normalize(P2)
    (P1n + P2n) ×₋ (P1n ∨ P2n)
end


"""
    p_bisect_pp(P1::PGA2DMV, P2::PGA2DMV)

Calculates the midpoint between P1 and P2.
"""
p_bisect_pp(P1::PGA2DMV, P2::PGA2DMV) = normalize(P1) + normalize(P2)


"""
    reflect_l(x::PGA2DMV, l::PGA2DMV)

Reflects any MultiVector encoded object x on the line l.
"""
reflect_l(x::PGA2DMV, l::PGA2DMV) = l * x * l

"""
    area_ppp(P1::PGA2DMV, P2::PGA2DMV, P3::PGA2DMV)

Calculates the oriented area of the triangle P1,P2,P3.
"""
area_ppp(P1::PGA2DMV, P2::PGA2DMV, P3::PGA2DMV) = scalar((P1 ∨ P2 ∨ P3) / (2 * norm(P1) * norm(P2) * norm(P3)))

"""
    area_loop(points)

Calculates the oriented area of the loop defined by the iterable points. 
"""
function area_loop(points)
    area = 0
    for k = 2:(length(points)-1)
        area += area_ppp(points[1], points[k], points[k+1])
    end
    area
end

"""
    length_loop(points)

Calculates the length of the loop defined by the iterable points.
"""
function length_loop(points)
    Pp = last(points)
    Ppn = normalize(Pp)
    l = zero(eltype(Ppn))
    for P in points
        Pn = P / norm(P)
        l += norm( Ppn ∨ Pn )
        Ppn = Pn
    end
    l
end

"""
    motor_pa(P::PGA2DMV, α::Real)

Calculates the motor for a rotation by the angle α around P. If P is a direction, then the motor describes a translation by α orthogonal to the direction.
"""
function motor_pa(P::PGA2DMV, α::Real)
    if is_direction(P)
        1 + α/(2*sqrt(P.e0e1^2+P.e2e0^2)) * P
    else
        exp(α/2*normalize(P))
    end
end

"""
    motor_ll(l1::PGA2DMV, l2::PGA2DMV)

Calculates the motor that maps the lines l1 to l2.
"""
function motor_ll(l1::PGA2DMV, l2::PGA2DMV)
    l1n = normalize(l1)
    l2n = normalize(l2)
    l2l1 = l2n * l1n
    l2l1n = normalize(l2l1)
    1 + l2l1n
end

"""
    move_pa(x::PGA2DMV, P::PGA2DMV, α::Real)

Moves x with motor_pa(x,P,α).
"""
move_pa(x::PGA2DMV, P::PGA2DMV, α::Real) = motor_pa(P,α) ≀ x

"""
    move_ll(x::PGA2DMV, l1::PGA2DMV, l2::PGA2DMV)

Moves x with motor_ll(l1,l2).
"""
move_ll(x::PGA2DMV, l1::PGA2DMV, l2::PGA2DMV) = motor_ll(l1,l2) ≀ x

"""
    move_m(x::PGA2DMV, m::PGA2DMV)

Moves x with the motor m.
"""
move_m(x::PGA2DMV, m::PGA2DMV) = m ≀ x

"""
    line_pd(d::PGA2DMV, P::PGA2DMV) -> PGA2DMV

Construct the line passing through point `P` with tangent direction `d`.
"""
function line_pd(d::PGA2DMV, P::PGA2DMV)
    (dx, dy) = direction_coordinates(d)
    (x0, y0) = point_coordinates(P)
    # Line normal is perpendicular to direction: (a, b) = (dy, -dx)
    a = dy
    b = -dx
    c = -(a * x0 + b * y0)
    line(a, b, c)
end


"""
    incenter_ppp(P1::PGA2DMV, P2::PGA2DMV, P3::PGA2DMV) -> PGA2DMV

Returns the incenter of triangle (P1, P2, P3) as the intersection of internal angle bisectors.
"""
function incenter_ppp(P1::PGA2DMV, P2::PGA2DMV, P3::PGA2DMV)
    # Triangle sides as lines
    L12 = join_pp(P1, P2)
    L23 = join_pp(P2, P3)
    L31 = join_pp(P3, P1)
    # Internal angle bisectors at vertices P1 and P2
    B1a, B1b = l_bisect_ll(L31, L12)
    B2a, B2b = l_bisect_ll(L12, L23)
    candidates = (
        meet_ll(B1a, B2a),
        meet_ll(B1a, B2b),
        meet_ll(B1b, B2a),
        meet_ll(B1b, B2b),
    )
    # Oriented side signs for the opposite vertex of each side
    s12 = sign(dist_orient_lp(L12, P3))
    s23 = sign(dist_orient_lp(L23, P1))
    s31 = sign(dist_orient_lp(L31, P2))
    # Select the candidate that lies inside the triangle (same side as opposite vertices)
    for X in candidates
        is_point(X) || continue
        sx12 = sign(dist_orient_lp(L12, X))
        sx23 = sign(dist_orient_lp(L23, X))
        sx31 = sign(dist_orient_lp(L31, X))
        if sx12 == s12 && sx23 == s23 && sx31 == s31
            return X
        end
    end
    # Fallback (should not happen for non-degenerate triangles)
    first(candidates)
end

"""
    circumcenter_ppp(P1::PGA2DMV, P2::PGA2DMV, P3::PGA2DMV) -> PGA2DMV

Returns the circumcenter of triangle (P1, P2, P3) as the intersection of two perpendicular bisectors.
"""
function circumcenter_ppp(P1::PGA2DMV, P2::PGA2DMV, P3::PGA2DMV)
    L12 = join_pp(P1, P2); M12 = p_bisect_pp(P1, P2); B12 = l_ortho_lp(L12, M12)
    L23 = join_pp(P2, P3); M23 = p_bisect_pp(P2, P3); B23 = l_ortho_lp(L23, M23)
    meet_ll(B12, B23)
end

"""
    orthocenter_ppp(P1::PGA2DMV, P2::PGA2DMV, P3::PGA2DMV) -> PGA2DMV

Returns the orthocenter of triangle (P1, P2, P3) as the intersection of two altitudes.
"""
function orthocenter_ppp(P1::PGA2DMV, P2::PGA2DMV, P3::PGA2DMV)
    L23 = join_pp(P2, P3); A1 = l_ortho_lp(L23, P1)
    L31 = join_pp(P3, P1); A2 = l_ortho_lp(L31, P2)
    meet_ll(A1, A2)
end

"""
    try_incenter_ppp(P1::PGA2DMV, P2::PGA2DMV, P3::PGA2DMV) -> Union{PGA2DMV,Nothing}

Non-throwing/degeneracy-friendly incenter; returns `nothing` if it cannot be determined.
"""
function try_incenter_ppp(P1::PGA2DMV, P2::PGA2DMV, P3::PGA2DMV)
    L12 = join_pp(P1, P2)
    L23 = join_pp(P2, P3)
    L31 = join_pp(P3, P1)
    B1a, B1b = l_bisect_ll(L31, L12)
    B2a, B2b = l_bisect_ll(L12, L23)
    candidates = (
        meet_ll(B1a, B2a),
        meet_ll(B1a, B2b),
        meet_ll(B1b, B2a),
        meet_ll(B1b, B2b),
    )
    s12 = sign(dist_orient_lp(L12, P3))
    s23 = sign(dist_orient_lp(L23, P1))
    s31 = sign(dist_orient_lp(L31, P2))
    for X in candidates
        is_point(X) || continue
        sx12 = sign(dist_orient_lp(L12, X))
        sx23 = sign(dist_orient_lp(L23, X))
        sx31 = sign(dist_orient_lp(L31, X))
        if sx12 == s12 && sx23 == s23 && sx31 == s31
            return X
        end
    end
    nothing
end

"""
    try_circumcenter_ppp(P1::PGA2DMV, P2::PGA2DMV, P3::PGA2DMV) -> Union{PGA2DMV,Nothing}

Non-throwing/degeneracy-friendly circumcenter; returns `nothing` if it cannot be determined.
"""
function try_circumcenter_ppp(P1::PGA2DMV, P2::PGA2DMV, P3::PGA2DMV)
    L12 = join_pp(P1, P2); M12 = p_bisect_pp(P1, P2); B12 = l_ortho_lp(L12, M12)
    L23 = join_pp(P2, P3); M23 = p_bisect_pp(P2, P3); B23 = l_ortho_lp(L23, M23)
    X = meet_ll(B12, B23)
    is_point(X) ? X : nothing
end

"""
    try_orthocenter_ppp(P1::PGA2DMV, P2::PGA2DMV, P3::PGA2DMV) -> Union{PGA2DMV,Nothing}

Non-throwing/degeneracy-friendly orthocenter; returns `nothing` if it cannot be determined.
"""
function try_orthocenter_ppp(P1::PGA2DMV, P2::PGA2DMV, P3::PGA2DMV)
    L23 = join_pp(P2, P3); A1 = l_ortho_lp(L23, P1)
    L31 = join_pp(P3, P1); A2 = l_ortho_lp(L31, P2)
    X = meet_ll(A1, A2)
    is_point(X) ? X : nothing
end

"""
    incircle_ppp(P1::PGA2DMV, P2::PGA2DMV, P3::PGA2DMV) -> (center::PGA2DMV, radius::Real)

Returns the incenter and inradius of triangle (P1, P2, P3).
"""
function incircle_ppp(P1::PGA2DMV, P2::PGA2DMV, P3::PGA2DMV)
    C = incenter_ppp(P1, P2, P3)
    L12 = join_pp(P1, P2)
    r = dist_lp(L12, C)
    return (center = C, radius = r)
end

"""
    circumcircle_ppp(P1::PGA2DMV, P2::PGA2DMV, P3::PGA2DMV) -> (center::PGA2DMV, radius::Real)

Returns the circumcenter and circumradius of triangle (P1, P2, P3).
"""
function circumcircle_ppp(P1::PGA2DMV, P2::PGA2DMV, P3::PGA2DMV)
    C = circumcenter_ppp(P1, P2, P3)
    r = dist_pp(C, P1)
    return (center = C, radius = r)
end


end
