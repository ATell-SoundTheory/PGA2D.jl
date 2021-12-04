module PGA2D

using CliffordAlgebras

import Base.join, Base.show

export pga2d, PGA2DT, PGA2DMV
export e0,e1,e2, e01, e12, e20, e012, 𝐈
export direction, point, line
export point_coordinates
export direction_coordinates
export line_coordinates
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
export angle_ll
export angle_dl
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

Constructs a MultiVector that encodes a line ax + bx + c = 0.
"""
line(a::Real, b::Real, c::Real) = a * e1 + b * e2 + c * e0

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
        print(io, "motor", " ≜ ")
    end
    show_multivector(io, x)
end

"""
    normalize(x::PGA2DMV)

Scales the MultiVector x so that it has unit norm.
"""
normalize(x::PGA2DMV) = x / norm(x)

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

Calculates the orthogonal Euclidean distance between parallel lines l1 and l2.
"""
dist_ll(l1::PGA2DMV, l2::PGA2DMV) = norm(dual(l1 ∧ l2)) / (norm(l1) * norm(l2))

"""
    angle_ll(l1::PGA2DMV, l2::PGA2DMV)

Calculates the angle between two lines l1 and l2.
"""
angle_ll(l1::PGA2DMV, l2::PGA2DMV) = asin( norm(l1 ∧ l2) / (norm(l1) * norm(l2)))

"""
    angle_dl(d::PGA2DMV, l2::PGA2DMV)

Calculates the angle between the direction d and the line l.
"""
angle_dl(d::PGA2DMV, l::PGA2DMV) = asin( norm(dual(d ∧ l)) / (norm(d) * norm(l)))

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

Calculates the oriented area of the loop defined by the iteratable points. 
"""
function area_loop(points)
    area = 0
    for k = 2:(length(points)-1)
        area += triangle_area(points[1],points[k],points[k+1])
    end
    area
end

"""
    length_loop(points)

Calculates the length of the loop defined by the iteratable points.
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
motor_pa(P::PGA2DMV, α::Real) = exp(α/2*P)

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

using Plots

@recipe function multivector_plot_recipe(mv::PGA2DMV)
    if is_point(mv)
        seriestype := :path
        markershape --> :circle
        markersize --> 5
        (x,y) = point_coordinates(mv)
        [x],[y]
    elseif is_line(mv)
        seriestype := :straightline
        (a,b,c) = line_coordinates(mv)
        if abs(a) > abs(b)
            [-c/a, -c/a+b],[0,-a]
        else
            [0,-b],[-c/b, -c/b+a]
        end
    else
        nothing
    end
end


end
