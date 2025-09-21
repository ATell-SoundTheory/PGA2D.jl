module PGA2DPlotsExt

using PGA2D
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
