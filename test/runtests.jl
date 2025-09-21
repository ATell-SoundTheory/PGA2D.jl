using PGA2D
using Test

@testset "PGA2D.jl" begin

    @testset "Coordinates" begin
        @test point_coordinates(point(1.0,2.0)) == (1.0,2.0)
        @test direction_coordinates(direction(1.0,2.0)) == (1.0,2.0)        
        @test line_coordinates(line(1.0,2.0,3.0)) == (1.0,2.0,3.0)
        @test try_point_coordinates(point(1.0,2.0)) == (1.0,2.0)
        @test try_direction_coordinates(direction(1.0,2.0)) == (1.0,2.0)
        @test try_line_coordinates(line(1.0,2.0,3.0)) == (1.0,2.0,3.0)
        @test try_point_coordinates(direction(1,2)) === nothing
        @test try_direction_coordinates(point(1,2)) === nothing
        @test try_line_coordinates(direction(1,2)) === nothing
        
        @test_throws DomainError point_coordinates(direction(1.0,2.0))
        @test_throws DomainError point_coordinates(line(1.0,2.0,3.0))
        @test_throws DomainError direction_coordinates(point(1.0,2.0))
        @test_throws DomainError direction_coordinates(line(1.0,2.0,3.0))
        @test_throws DomainError line_coordinates(point(1.0,2.0))
        @test_throws DomainError line_coordinates(direction(1.0,2.0))

        @test is_point(point(1,2))
        @test is_line(line(1,2,3))
        @test is_direction(direction(1,2,))

        @test !is_point(direction(1,2))
        @test !is_point(line(1,2,3))
        @test !is_direction(point(1,2))
        @test !is_direction(line(1,2,3))
        @test !is_line(point(1,2))
        @test !is_line(direction(1,2))
    end

    @testset "Geometric operations" begin
        P1 = point(0,0)
        P2 = point(2,0)
        P3 = point(3,1)
        P4 = point(1,2)

        l1 = join_pp(P1,P2)
        l2 = join_pp(P2,P3)
        l3 = join_pp(P3,P4)
        l4 = join_pp(P4,P1)

        @test point_coordinates(meet_ll(l1,l2)) == point_coordinates(P2)
        @test point_coordinates(meet_ll(l2,l3)) == point_coordinates(P3)
        @test point_coordinates(meet_ll(l3,l4)) == point_coordinates(P4)
        @test point_coordinates(meet_ll(l4,l1)) == point_coordinates(P1)

        @test point_coordinates(meet_ll(join_pp(P1,P2),join_pp(P3,P4))) == (5.0,0.0)

       # Distance and angles
       l = join_pp(P1,P2)
       @test dist_lp(l, point(0,1)) ≈ 1.0
       @test angle_ll(join_pp(P1,P2), join_pp(P1,P4)) > 0
       d1 = direction(1,0)
       d2 = direction(0,1)
       @test angle_dd(d1, d2) ≈ π/2 atol=1e-12
    @test angle_dd_signed(d1, d2) ≈ π/2 atol=1e-12
    @test angle_dd_signed(d2, d1) ≈ -π/2 atol=1e-12
       # Clamping robustness: nearly parallel lines
       lε = join_pp(point(0,0), point(1,1e-14))
       θ = angle_ll(l, lε)
       @test 0 <= θ <= π/2
    # Signed line-line and dir-line
    @test angle_ll_signed(l, join_pp(P1, P4)) > 0
    @test angle_dl_signed(d1, l) ≈ 0 atol=1e-12
  @test angle_dl_signed(d2, l) ≈ -π/2 atol=1e-12
    # line through point and direction
  lpd = line_pd(direction(1,0), point(0,1))
  @test dist_lp(lpd, point(0,1)) ≈ 0 atol=1e-12
  @test abs(angle_dl_signed(direction(1,0), lpd)) ≈ 0 atol=1e-12
   end  

   @testset "Normalize guard" begin
       # Construct a concrete zero multivector (scalar 0 in the algebra)
       z = zero(point(0,0))
       @test_throws DomainError normalize(z)
   end  

   @testset "Triangle centers" begin
       P1, P2, P3 = point(0,0), point(1,0), point(0,1)
       I = incenter_ppp(P1,P2,P3)
       @test is_point(I)
       (x, y) = point_coordinates(I)
       r = 1 - sqrt(2)/2
       @test x ≈ r atol=1e-12
       @test y ≈ r atol=1e-12

       Cc = circumcenter_ppp(P1,P2,P3)
       @test is_point(Cc)
       (xc, yc) = point_coordinates(Cc)
        @test xc ≈ 0.5 atol=1e-12
        @test yc ≈ 0.5 atol=1e-12

       H = orthocenter_ppp(P1,P2,P3)
       @test is_point(H)
       (xh, yh) = point_coordinates(H)
        @test xh ≈ 0.0 atol=1e-12
        @test yh ≈ 0.0 atol=1e-12

      # Circles
    inc = incircle_ppp(P1,P2,P3)
    @test is_point(inc.center)
    (xic, yic) = point_coordinates(inc.center)
      @test xic ≈ r atol=1e-12
      @test yic ≈ r atol=1e-12
    @test inc.radius ≈ r atol=1e-12

    circ = circumcircle_ppp(P1,P2,P3)
    @test is_point(circ.center)
    (xcc, ycc) = point_coordinates(circ.center)
      @test xcc ≈ 0.5 atol=1e-12
      @test ycc ≈ 0.5 atol=1e-12
    @test circ.radius ≈ sqrt(0.5) atol=1e-12

    # try_* centers should succeed here
    @test is_point(try_incenter_ppp(P1,P2,P3))
    @test is_point(try_circumcenter_ppp(P1,P2,P3))
    @test is_point(try_orthocenter_ppp(P1,P2,P3))
   end

end
