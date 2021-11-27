module PGA2D

using CliffordAlgebras

import Base.join, Base.show

export pga2d
export e0,e1,e2, e01, e12, e20, e012
export direction, point, line
export point_coordinates
export direction_coordinates
export line_coordinates
export is_point
export is_direction
export is_line
export join, meet
export line_orthogonal
export point_project
export line_point
export direction_orthogonal
export distance_points
export distance_lines
export angle_lines
export angle_direction_line
export oriented_distance_line_point
export angle_bisectors
export perpendicular_bisector
export triangle_area
export loop_area
export loop_length
export motor, move
export motor_from_lines
export reflect_line
export normalize

const pga2d = CliffordAlgebra(:PGA2D)
const PGA2DT = typeof(pga2d)
const PGA2DMV = MultiVector{PGA2DT}

const e0 = pga2d.e0
const e1 = pga2d.e1
const e2 = pga2d.e2

const e01 = pga2d.e0e1
const e20 = pga2d.e2e0
const e12 = pga2d.e1e2

const e012 = pga2d.e1e2e0
const 𝐈 = e012

direction(x::Real, y::Real) = x * e20 + y * e01
point(x::Real, y::Real) = direction(x, y) + e12
line(a::Real, b::Real, c::Real) = a * e1 + b * e2 + c * e0

function point_coordinates(P::PGA2DMV)
    w = P.e1e2
    if isgrade(P,2) && !iszero(w)
        wi = inv(w)
        return wi.*(P.e2e0, P.e0e1)
    else
        throw(DomainError(P, " is not a point."))
    end
end

function direction_coordinates(d::PGA2DMV)
    w = d.e1e2
    if isgrade(d,2) && iszero(w)
        return (d.e2e0, d.e0e1)
    else
        throw(DomainError(d, " is not a direction."))
    end
end

function line_coordinates(l::PGA2DMV)
    if isgrade(l,1) && (!iszero(l.e1) || !iszero(l.e2))
        return (l.e1, l.e2, l.e0)
    else
        throw(DomainError(l, " is not a line."))
    end
end

function is_point(P::PGA2DMV)
    isgrade(P,2) && !iszero(P.e1e2)
end

function is_direction(d::PGA2DMV)
    isgrade(d,2) && iszero(d.e1e2)
end

function is_line(l::PGA2DMV)
    isgrade(l,1) && (!iszero(l.e1) || !iszero(l.e2))
end

function is_motor(m::PGA2DMV)
    # TODO    
    false
end

function show(io::IO, x::PGA2DMV)
    if is_point(x)
        print(io, "point", point_coordinates(x), " ≜ ")
    elseif is_direction(x)
        print(io, "direction", direction_coordinates(x), " ≜ ")
    elseif is_line(x)
        print(io, "line", line_coordinates(x), " ≜ ")
    elseif is_motor(x)
        # TODO
        print(io, "motor", " ≜ ")
    end
    show_multivector(io, x)
end


meet(l1::PGA2DMV, l2::PGA2DMV) = l1 ∧ l2
join(P1::PGA2DMV, P2::PGA2DMV) = P1 ∨ P2

line_orthogonal(l::PGA2DMV, P::PGA2DMV) = l ⋅ P
point_project(l::PGA2DMV, P::PGA2DMV) = (l ⋅ P) * l
line_point(l::PGA2DMV, P::PGA2DMV) = (l ⋅ P ) * P

direction_orthogonal(l::PGA2DMV) = l * 𝐈

distance_points(P1::PGA2DMV, P2::PGA2DMV) = norm(P1 ∨ P2) / (norm(P1) * norm(P2))
distance_lines(l1::PGA2DMV, l2::PGA2DMV) = norm(dual(l1 ∧ l2)) / (norm(l1) * norm(l2))

angle_lines(l1::PGA2DMV, l2::PGA2DMV) = asin( norm(l1 ∧ l2) / (norm(l1) * norm(l2)))
angle_direction_line(d::PGA2DMV, l::PGA2DMV) = asin( norm(dual(d ∧ l)) / (norm(d) * norm(l)))

oriented_distance_line_point(l::PGA2DMV, P::PGA2DMV) = scalar(P ∨ l) / (norm(P) * norm(l))

angle_bisectors(l1::PGA2DMV, l2::PGA2DMV) = (l1 + l2, l1 - l2) .* inv(norm(l1) * norm(l2))

function perpendicular_bisector(P1::PGA2DMV, P2::PGA2DMV) 
    P1n = normalize(P1)
    P2n = normalize(P2)
    (P1n + P2n) * (P1n ∨ P2n)
end

triangle_area(P1::PGA2DMV, P2::PGA2DMV, P3::PGA2DMV) = scalar((P1 ∨ P2 ∨ P3) / (2 * norm(P1) * norm(P2) * norm(P3)))

function loop_area(points)
    area = 0
    for k = 2:(length(points)-1)
        area += triangle_area(points[1],points[k],points[k+1])
    end
    area
end


function loop_length(points)
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

motor(P::PGA2DMV, α::Real) = exp(α/2*P)

move(x::PGA2DMV, P::PGA2DMV, α::Real) = motor(P,α) ≀ x

function motor_from_lines(l1::PGA2DMV, l2::PGA2DMV)
    l1n = normalize(l1)
    l2n = normalize(l2)
    l2l1 = l2n * l1n
    l2l1n = normalize(l2l1)
    1 + l2l1n
end

reflect_line(x::PGA2DMV, l::PGA2DMV) = l * x * l

normalize(x::PGA2DMV) = x / norm(x)

end
