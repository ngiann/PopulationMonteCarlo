##############################################
function plot_ellipse(μ, Σ, clr="b", label="")
##############################################

    # http://www.visiondummy.com/2014/04/draw-error-ellipse-representing-covariance-matrix/

    # calculate eigevanlues
    E = eigvals(Σ)
    @assert(E[1]<=E[2])

    # calculate eigenvectors, columns are the eigenvectors
    V = eigvecs(Σ)

    ϕ = atan(V[2,2], V[1,2]) # was atan2

    if(ϕ < 0.0)
        ϕ = ϕ + 2*pi;
    end

    chisquare_val = sqrt(2.41) # 70% confidence interval
    theta_grid    = LinRange(0.0, 2*pi, 500)
    a = chisquare_val*sqrt(maximum(E))
    b = chisquare_val*sqrt(minimum(E))

    # the ellipse in x and y coordinates
    ellipse_x_r  = a*cos.( theta_grid )
    ellipse_y_r  = b*sin.( theta_grid )

    # Define a rotation matrix
    R = [ cos.(ϕ) sin.(ϕ); -sin.(ϕ) cos.(ϕ) ];


    # let's rotate the ellipse to angle ϕ
    r_ellipse = [ellipse_x_r ellipse_y_r]  * R

    plot(μ[1], μ[2], "x", markersize=12, markeredgewidth=1, color=clr)
    plot(r_ellipse[:,1] .+ μ[1], r_ellipse[:,2] .+ μ[2], clr, linewidth=1, label=label)

    # plot ellipse axis (eigenvectors)
    a_pnts = [vec(LinRange(0.0, a, 500))  vec(LinRange(0.0, a, 500))*0.0    ]*R
    b_pnts = [vec(LinRange(0.0, b, 500))*0.0      vec(LinRange(0.0, b, 500))]*R
    plot(a_pnts[:,1] .+ μ[1], a_pnts[:,2] .+ μ[2], @sprintf("--%s",clr), linewidth=1)
    plot(b_pnts[:,1] .+ μ[1], b_pnts[:,2] .+ μ[2], @sprintf("--%s",clr), linewidth=1)
    axis("equal")
end
