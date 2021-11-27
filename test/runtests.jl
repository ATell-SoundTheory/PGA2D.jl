using PGA2D
using Test

@testset "PGA2D.jl" begin

    @testset "Coordinates" begin
        @test point_coordinates(point(1.0,2.0)) == (1.0,2.0)
        @test direction_coordinates(direction(1.0,2.0)) == (1.0,2.0)        
        @test line_coordinates(line(1.0,2.0,3.0)) == (1.0,2.0,3.0)
        
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

        l1 = join(P1,P2)
        l2 = join(P2,P3)
        l3 = join(P3,P4)
        l4 = join(P4,P1)

        @test point_coordinates(meet(l1,l2)) == point_coordinates(P2)
        @test point_coordinates(meet(l2,l3)) == point_coordinates(P3)
        @test point_coordinates(meet(l3,l4)) == point_coordinates(P4)
        @test point_coordinates(meet(l4,l1)) == point_coordinates(P1)

        @test point_coordinates(meet(join(P1,P2),join(P3,P4))) == (5.0,0.0)

   end  

end
