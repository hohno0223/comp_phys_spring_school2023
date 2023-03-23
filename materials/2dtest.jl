function udot(x,y)
    r = sqrt(x^2+y^2)
    return -x/r^3
end

function vdot(x,y)
    r = sqrt(x^2+y^2)
    return -y/r^3
end

function calc_E(x,y,u,v)
    return (u^2+v^2)/2 - 1/sqrt(x^2+y^2)
end

using Plots
ENV["GKSwstype"] = "nul"



function main()
    x = 1.0
    y = 0.0
    u = 0.0
    v = 1.0
    n = 10000
    tmax = 20
    Δt = tmax/n
    t = 0.0
    filename = "anim_v$v.gif"
    nprint = 50

    u_m = u - udot(x,y)*Δt/2
    v_m = v - vdot(x,y)*Δt/2
    E0 = calc_E(x,y,u,v)
    anim = Animation()

    for i=2:n
        u_p = u_m + udot(x,y)*Δt
        v_p = v_m + vdot(x,y)*Δt

        x += u_p*Δt
        y += v_p*Δt
        t += Δt

        E = calc_E(x,y,u,v)
        dE = abs(E-E0)/abs(E)

        u_m = u_p
        v_m = v_p

        u = (u_p + u_m)/2
        v = (v_p + v_m)/2


        if i % nprint == 0
            println("$t $x $y $E $dE")
            plt = scatter([0,x],[0,y],label="t = $t",xlims=(-2,2),ylims=(-2,2),aspect_ratio = 1)
            frame(anim,plt)
        end
    end
    gif(anim, filename, fps = 30)



end 
main()