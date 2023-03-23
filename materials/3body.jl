struct Body
    r::Vector{Float64}
    v::Vector{Float64}
    f::Vector{Float64}
    mass::Float64
    function Body(r,v,mass)
        f = zero(r)
        return new(r,v,f,mass)
    end
end

mutable struct SolarSystem
    bodies::Vector{Body}
    nbody::Int64
    function SolarSystem()
        bodies = Body[]
        nbody = 0
        return new(bodies,nbody)
    end
end

function Base.push!(s::SolarSystem,body::Body)
    push!(s.bodies,body)
    s.nbody = length(s.bodies)
end

function Base.getindex(s::SolarSystem,i)
    return s.bodies[i]
end

function calc_kinetic_e(s::SolarSystem)
    e = 0.0
    nbody = s.nbody
    for i=1:nbody
        mass = s[i].mass
        v = s[i].v
        e += (mass/2)*dot(v,v)
    end
    return e
end

function calc_potential_e(s::SolarSystem)
    e = 0.0
    nbody = s.nbody
    for i=1:nbody
        r_i = s[i].r
        mass_i = s[i].mass
        for j=i+1:nbody
            r_j = s[j].r
            mass_j = s[j].mass
            rij = r_i - r_j
            #println(G*mass_i*mass_j/norm(rij))
            e += G*mass_i*mass_j/norm(rij)
        end
    end
    return e
end

function calc_force!(s::SolarSystem)
    nbody = s.nbody
    for i=1:nbody
        s[i].f .= 0
    end
    for i=1:nbody
        r_i = s[i].r
        mass_i = s[i].mass
        for j=i+1:nbody
            r_j = s[j].r
            mass_j = s[j].mass
            rij = r_i - r_j
            s[i].f .+= -G*mass_i*mass_j*rij/norm(rij)^3
            s[j].f .+= -s[i].f
        end
    end
end

function update!(s::SolarSystem,Δt)
    nbody = s.nbody
    for i=1:nbody
        mass = s[i].mass
        s[i].v .+= s[i].f*Δt/(2*mass)
        s[i].r .+= s[i].v*Δt
    end

    calc_force!(s)

    for i=1:nbody
        mass = s[i].mass
        s[i].v .+= s[i].f*Δt/(2*mass)
    end

end


const G = 6.67384e-11
const mass_Sun = 1.989e30 #[kg]
const mass_Earth = 5.972e24 #[kg]
const mass_Moon = 7.348e22 #[kg]

const radius_Earth = 1.496e11 #[m]
const velocity_Earth = 2.978e4 #[m/s]
const radius_Moon = 3.844e8 #[m]
const velocity_Moon = 1.022e3 #[m/s]

using LinearAlgebra


using Plots

function main()
    solarsystem = SolarSystem()
    sun = Body([0,0],[0,0],mass_Sun )
    push!(solarsystem,sun)

    r_Earth = [radius_Earth,0]
    v_Earth = [0,velocity_Earth]
    earth = Body(r_Earth,v_Earth ,mass_Earth )
    push!(solarsystem,earth)

    r_Moon = [radius_Earth+radius_Moon ,0]
    v_Moon = [0,velocity_Earth+velocity_Moon]
    moon = Body(r_Moon,v_Moon ,mass_Moon )
    push!(solarsystem,moon)


    t = 0.0
    tmax = 3.2e7
    n = 10^6
    Δt = tmax/n

    E0 = calc_kinetic_e(solarsystem ) + calc_potential_e(solarsystem )
    anim1 = Animation()
    anim2 = Animation()
    filename1 = "3body_SE.gif"
    filename2 = "3body_EM.gif"
    #fp = open("3body_data.txt","w")
    for i = 2:n
        
        if i % div(n,100) == 0
            E = calc_kinetic_e(solarsystem ) + calc_potential_e(solarsystem )
            dE = abs(E-E0)/abs(E0)
            println("$t $E $dE")
            r_E = solarsystem[2].r
            r_M = solarsystem[3].r
            plt1 = scatter([0,r_E[1],r_M[1]],[0,r_E[2],r_M[2]],xlims=(-2*radius_Earth,2*radius_Earth),ylims=(-2*radius_Earth,2*radius_Earth),label="t = $t Earth",aspect_ratio = 1, markersize=1)
            plt1 = scatter!([r_M[1]],[r_M[2]],xlims=(-2*radius_Earth,2*radius_Earth),ylims=(-2*radius_Earth,2*radius_Earth),label="Moon",aspect_ratio = 1,markersize=1)
            frame(anim1,plt1)
            plt2 = scatter([0,r_M[1]-r_E[1]],[0,r_M[2]-r_E[2]],xlims=(-2*radius_Moon,2*radius_Moon),ylims=(-2*radius_Moon,2*radius_Moon),label="t = $t Moon",aspect_ratio = 1)
            frame(anim2,plt2)
        end
        update!(solarsystem,Δt)
        t += Δt
    end
    gif(anim1, filename1, fps = 60)
    gif(anim2, filename2, fps = 60)




    
end

main()