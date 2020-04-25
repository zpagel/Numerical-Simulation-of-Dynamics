
# Let's feed the integration function an initial wavefunction ψ0(x), 
# a potential V(x,t), a position vector xs for gridding position, 
# a time vector ts for gridding time. 
using LinearAlgebra

function integrate_SE(ψ0,V,xs,ts)
    nTs = length(ts);
    nXs = length(xs);
    dx = xs[2] - xs[1];
    dt = ts[2] - ts[1];
    
    #initialize state gridded in position and time
    ψf = (1+1im).*zeros((nXs, nTs));
    
    #initial conditions
    ψf[:,1] = (1+0im).*ψ0.(xs);
    
    #α is a parameter in the state evolution matrices
    α = 0.5im*dt/dx^2;
    
    for i=2:nTs
        #potential vector for this timestep
        V_i = V.(xs,ts[i]);
        #diagonal entries for U_1 and U_2
        ξ_i = 1 .+ (0.5im*dt).*(2/dx^2 .+ V_i);
        γ_i = 1 .- (0.5im*dt).*(2/dx^2 .+ V_i);
        #U_1 and U_2 matrices
        U_1 = SymTridiagonal(ξ_i,-α.*ones(nXs-1));
        U_2 = SymTridiagonal(γ_i,α.*ones(nXs-1));
        #Solve for state vector at next timestep. Note that I don't explicitly calculate inv(U_1), which 
        #would take much longer. This line of code let's Julia use algorithms that benefite from a) not
        #having to explicitly calvulate inv(U_1), and also fast algorithms for linear equations using
        #symmetric tridiagonal matrices. In total, you get a massive speedup.
        ψf[:,i] = U_1\(U_2*ψf[:,i-1]);
    end
    
    return ψf
end

# Now let's test this out with a Gaussian wavefunction in free space!

function ψ0(x)
    σ = 2;
    μ = 0;
    a = 1/sqrt(2pi*σ^2)*exp(-(x-μ)^2/(2σ^2));
    return a
end

#free space wavefunction --> no potential.
function V(x,t)
    return 0
end

xs = -10:.1:10;
ts = 0:0.05:5;
@time ψ = integrate_SE(ψ0,V,xs,ts);

# pcolormesh is the best plotting function I could find to make a 2D density plot. 
# It's based on Python's matplotlib library, so we have to use PyPlot package for things to go smoothly. 
using PyPlot
#density plot of dynamics
pcolormesh(abs.(ψ))

# Cold wavefunction means it's velocity spread is much less than 1, in this case I use 0.05*v_r.
function ψ0_cold(x)
    dv = 0.05; #units of recoil velocity
    σ = 1/(2*dv); #Heisenberg limited wavefunction
    μ = 0;
    a = 1/sqrt(2pi*σ^2)*exp(-(x-μ)^2/(2σ^2));
    return a
end

# Gives the potential as a function of x and t. This is a traveling cosine wave with a Gaussian intensity profile in time.
function V_Bragg(x,t)
    tau = .2596;
    n = 4;
    int = 3.0871;
    t_pulse = 10tau;
    f_Bragg = 4*n;
    # do one Bragg pulse
    🐢 = 2pi*int*exp(-(t-t_pulse)^2/(tau^2))*cos(f_Bragg*t - 2*x);
    return 🐢
end

# test out
xs = -50:.1:150;
ts = 0:0.005:10;
@time ψ = integrate_SE(ψ0_cold,V_Bragg,xs,ts);

using PyPlot
#density plot of dynamics
pcolormesh(abs.(ψ'))
title("Fourth order Bragg diffraction")

using LinearAlgebra

#define a new integrate_SE_BBS function specific to the Bloch beamsplitter. In addition to the 
# initial condition ψ, potential V, and the position and time grids xs and ts, it also takes in
# an array params. See the definition of params below to see which variables go where. This was 
# the easiest way I could figure out to get all of the necessary parameters into each function 
# without too many variables flying around.
function integrate_SE_BBS(ψ0,V,xs,ts,params)
    nTs = length(ts);
    nXs = length(xs);
    dx = xs[2] - xs[1];
    dt = ts[2] - ts[1];
    ψf = (1+1im).*zeros((nXs, nTs));
    ψf[:,1] = (1+0im).*[ψ0(xs[j],params) for j=1:nXs];
    α = 0.5im*dt/dx^2;
    
    for i=2:nTs
        V_i = [V(xs[j],ts[i],params) for j=1:nXs];
        ξ_i = 1 .+ (0.5im*dt).*(2/dx^2 .+ V_i);
        γ_i = 1 .- (0.5im*dt).*(2/dx^2 .+ V_i);
        U_1 = SymTridiagonal(ξ_i,-α.*ones(nXs-1));
        U_2 = SymTridiagonal(γ_i,α.*ones(nXs-1));
        ψf[:,i] = U_1\(U_2*ψf[:,i-1]);
    end
    
    return ψf
end

#unit step/heaviside function, because they don't have a built-in Julia on for some reason.
function unit_step(x)
    if x<0
        return 0
    else
        return 1
    end
end

function ψ0_cold(x,params)
    σ = params[12];    #Heisenberg limited wavefunction, velocity spread σ
    μ = 0;
    a = 1/sqrt(2pi*σ^2)*exp(-(x-μ)^2/(2σ^2));
    return a
end


#BBS potential is made up of an adiabatic lattice load, then a frequency ramp, then an adiabatic
# unload, and finally a free propogation where the potential is zero.

function V_BBS(x,t,params)
    bdry = params[1];
    t_propogate = params[2];
    t_load = params[3];
    n_bloch = params[2];
    n_ramp = params[5];
    n_lat_dep = params[6];
    grav_ramp = params[7];
    ramp_rate = params[8];
    bloch_period = params[9];
    t_bloch = params[10];
    dv = params[11];
    dx = params[12];
    
    #BBS Potential. unit_step(x) defined above. 🐢s all the way down. I'll switch to 🍍 at some point.
    
    🐢 = n_lat_dep * cos(2*x)*(t/t_load*unit_step(t)*unit_step(t_load-t) +
        cos(ramp_rate*(t-t_load)^2/2)*unit_step(t-t_load)*unit_step(t_load + t_bloch - t) + 
        (t_load - t + t_load + t_bloch)/t_load*cos(ramp_rate*t_bloch^2/2 + ramp_rate*t_bloch*
            (t - t_load - t_bloch))*unit_step(t - t_load - t_bloch)*unit_step(2*t_load + t_bloch - t));
    return 🐢
end

#reduce the size of the output state to reduce plotting time. This function keeps every $nx^{th}$ point
# in space, and every $nt^{th}$ point in time, so it reduces the wavefuction by a factor of nx*nt. 
# Note that nx and nt should probably be integers to avoid problems.
function reduce_size(ψ, nx, nt)
    old_x = length(ψ[:,1]);
    old_t = length(ψ[1,:]);
    new_x = Int(floor(old_x/nx));
    new_t = Int(floor(old_t/nt));
    
    ψ_new = ψ[nx:nx:new_x*nx, nt:nt:new_t*nt];
    println("Size of new wavefunction: ", size(ψ_new))
    return ψ_new
end


#parameters for simulation
bdry = 75;
t_propogate = 65;
t_load = pi;
n_bloch = 1;
n_ramp = 0.1;
n_lat_dep = 3.5;
grav_ramp = 0.8576;
ramp_rate = grav_ramp * n_ramp;
bloch_period = 8/ramp_rate;
t_bloch = n_bloch*bloch_period;
dv = 0.05;
dx = 1/(2*dv);

#put parameters into params vector to pass into other functions
params = [bdry,t_propogate,t_load,n_bloch,n_ramp,n_lat_dep,grav_ramp,ramp_rate,bloch_period,t_bloch,dv,dx];

#Define grid
xs = -bdry*dx:0.04:bdry*dx;
ts = 0:0.04:(2*t_load + t_bloch + t_propogate);

#Solve for wavefunction
println("Metrics on solving for wavefunction, including time to solve: ")
@time ψ = integrate_SE_BBS(ψ0_cold,V_BBS,xs,ts,params);
println("And you just used up ", sizeof(ψ)/10^9, " Gb to store ψ.")
println("")

#Reduce the wavefunction size for plotting. 
ψ_new = reduce_size(ψ, 15, 2);
println("On the other hand, ψ_new only uses ", sizeof(ψ_new)/10^9, " Gb to store.")


using PyPlot, Colors  # I don't think I ever use Colors...?
#Make my own colormap 🍎!
#I'm trying to replicate the CMYK colormap in Mathematica. This colormap also exists in Julia, but
#I don't think it's compatible with the matplotlib functions. Anyways, it's cooler to make your own.

#I haven't made this very general yet, so changing n_colors doesn't automatically update everything else. 
n_colors = 4;
c = zeros((n_colors, 4));

#format is [r,g,b,α], where α is the transparency.
c[1,:] = [.4, .75, 1, .7];
c[2,:] = [1, .4, 0.45, .7];
c[3,:] = [1, 1, .55, .6];
c[4,:] = [0, 0, 0, 1.0];

#Set the limits for when the colors change. For example, if lim1 = .15, then the color will reach that of c[2,:]
# by the time the function reaches 0.15*(f_max - f_min). 
lim1 = 0.15;
lim2 = 0.5;

#You need a weird data structure for making a colormap, see the PyPlot matplotlib documentation.
ctable = [[(0.0, c[1,i], c[1,i]), (lim1, c[2,i], c[2,i]), (lim2, c[3,i], c[3,i]), (1.0, c[4,i], c[4,i])] for i=1:4];

🍎 = ColorMap("🍎", ctable[1], ctable[2], ctable[3], ctable[4])

A = [i for i=1:1000, j=1:1]

pcolormesh(A, cmap = 🍎)

#Plot the BBS dynamics

using Pkg
#Pkg.add("Colors")
using PyPlot
#import Colors.delta


#density plot of dynamics
fig = figure("Plotting state evolution and timing",figsize=(20,20))
subplot(221)
main = pcolormesh(abs.(ψ_new'),cmap = 🍎)
title("Bloch Beamsplitter evolution")

#plot the timing sequence on the RHS of the graph.
subplot(222)
timing = PyPlot.plot([0;1],[0;0],[0,1],[t_load,t_load],[0,1],[t_load+t_bloch,t_load+t_bloch],
    [0,1],[2*t_load+t_bloch,2*t_load+t_bloch],[0,1],
    [2*t_load+t_bloch+t_propogate,2*t_load+t_bloch+t_propogate])
xlim(0,1)
ylim(0,2*t_load+t_bloch+t_propogate)
fig[:canvas][:draw]() # Update the figure
suptitle("Timing sequence")

using LinearAlgebra

function integrate_SE_BBS_lite(ψ0,V,xs,ts,params)
    #dumb variables
    nTs = length(ts);
    nXs = length(xs);
    reduced_nXs = Int(floor(nXs/params[13]));
    reduced_nTs = Int(floor(nTs/params[14]));
    
    #current and previous state vectors.
    ψ_current = (1+1im).*zeros(nXs);
    ψ_previous = (1+1im).*zeros(nXs);
    ψf = (1+0.0im).*zeros((reduced_nXs, reduced_nTs));
    
    #Initial condition
    ψ_previous= (1+0im).*[ψ0(xs[j],params) for j=1:nXs];
    
    #more dumb varriables
    dx = xs[2] - xs[1];
    dt = ts[2] - ts[1];
    α = 0.5im*dt/dx^2;
    
    for i=2:nTs
        V_i = [V(xs[j],ts[i],params) for j=1:nXs];
        ξ_i = 1 .+ (0.5im*dt).*(2/dx^2 .+ V_i);
        γ_i = 1 .- (0.5im*dt).*(2/dx^2 .+ V_i);
        U_1 = SymTridiagonal(ξ_i,-α.*ones(nXs-1));
        U_2 = SymTridiagonal(γ_i,α.*ones(nXs-1));
        
        #rewrite state vectors
        ψ_current = U_1\(U_2*ψ_previous);
        ψ_previous = ψ_current;
        
        #only store in output vector if timestep is a multiple of n_reduce_t, and only store spatial values
        #that are a multiple of n_reduce_x
        if i%params[14] == 0 && Int(i/params[14])<=reduced_nTs
            ψf[:,Int(i/params[14])] = [ψ_current[Int(i)] for i=params[13]:params[13]:params[13]*reduced_nXs];
        end
    end
    
    return ψf
end


bdry = 55;
t_propogate = 35;
t_load = pi;
n_bloch = 1;
n_ramp = 0.1;
n_lat_dep = 2.5;
grav_ramp = 0.8576;
ramp_rate = grav_ramp * n_ramp;
bloch_period = 8/ramp_rate;
t_bloch = n_bloch*bloch_period;
dv = 0.05;
dx = 1/(2*dv);
n_reduce_x = 15;
n_reduce_t = 2;

params = [bdry,t_propogate,t_load,n_bloch,n_ramp,n_lat_dep,grav_ramp,ramp_rate,bloch_period,
    t_bloch,dv,dx,n_reduce_x,n_reduce_t];

xs = -bdry*dx:0.04:bdry*dx;
ts = 0:0.04:(2*t_load + t_bloch + t_propogate);

println("Metrics on solving for wavefunction, including time to solve: ")
@time ψ = integrate_SE_BBS_lite(ψ0_cold,V_BBS,xs,ts,params);
println("You used ", sizeof(ψ)/10^9, " Gb memory to store ψ.")

using PyPlot
#import Colors.delta

#define weird x & y vectors for making the timing plot on the right. 
x = 1:length(ψ[:,1])
y = 1:length(ψ[1,:])
dy = .1;
dx = .1;

#density plot of dynamics
fig= figure("Plotting state evolution and timing",figsize=(20,20))

subplot(221)
#main plot
main = pcolormesh(abs.(ψ'),cmap = 🍎)
axis("off")
title("Very High Lattice Depth", fontsize = 24.0)
xlabel("Position")
ylabel("Time")

#subplot on the side for showing timing of sequence
subplot(264)
timing1 = PyPlot.plot([0,1],[0,0],color = :black)
timing2 = PyPlot.plot([0,1],[t_load,t_load],color = :black)
timing3 = PyPlot.plot([0,1],[t_load+t_bloch,t_load+t_bloch],color = :black)
timing4 = PyPlot.plot([0,1],[2*t_load+t_bloch,2*t_load+t_bloch],color = :black)
timing5 = PyPlot.plot([0,1],[2*t_load+t_bloch+t_propogate,2*t_load+t_bloch+t_propogate],color = :black)
#annotate("Look, data!",
#    xy=[x[convert(Int64,floor(length(x)/4.1))];y[convert(Int64,floor(length(y)/4.1))]],
#    xytext=[x[convert(Int64,floor(length(x)/4.1))]+0.1dx;y[convert(Int64,floor(length(y)/4.1))]+0.1dy],
#    xycoords="data")
xlim(0,1)
ylim(0 - 0.1,2*t_load+t_bloch+t_propogate + 0.1)
annotate("Frequency ramping",
    fontsize = 18.0,
    xy=[0,0],
    xytext=[-0.05,t_load + t_propogate/2 + 10],
    xycoords="data")
annotate("Lattice Load",
    fontsize = 18.0,
    xy=[0,0],
    xytext=[0.1,t_load + 15],
    xycoords="data")
annotate("Lattice Unload",
    fontsize = 18.0,
    xy=[0,0],
    xytext=[0.1,t_load + 15 + t_bloch + t_load],
    xycoords="data")
annotate("Free Propogation",
    fontsize = 18.0,
    xy=[0,0],
    xytext=[0,t_load + 15 + t_bloch + t_load + t_propogate/3],
    xycoords="data")
PyPlot.arrow(0.5,
    16,
    -0.05,
    -10,
    head_width=0.1,
    width=0.03,
    head_length=3,
    overhang=0,
    head_starts_at_zero="true",
    facecolor="black")
PyPlot.arrow(0.5,
    16 + t_bloch + t_load,
    -0.05,
    -10,
    head_width=0.1,
    width=0.03,
    head_length=3,
    overhang=0,
    head_starts_at_zero="true",
    facecolor="black")
axis("off")
fig[:canvas][:draw]() # Update the figure

#function to export the figure. Can put in whatever file extension you want for the image. Note that it's
# very tine consuming to export large image files. 
#@time PyPlot.savefig("BBS_very_high_lat_dep.pdf")

# Doesn't really work as planned... On second thought, this was a terrible idea
function V_BBS_four_port(x,t,params)
    bdry = params[1];
    t_propogate = params[2];
    t_load = params[3];
    n_bloch = params[2];
    n_ramp = params[5];
    n_lat_dep = params[6];
    grav_ramp = params[7];
    ramp_rate = params[8];
    bloch_period = params[9];
    t_bloch = params[10];
    dv = params[11];
    dx = params[12];
    
    #BBS Potential. unit_step(x) defined above. 🐢s all the way down
    
    🐢 = n_lat_dep * cos(2*x)*(t/t_load*unit_step(t)*unit_step(t_load-t) +
        0.5*(cos(ramp_rate*(t-t_load)^2/2) + 
            cos(2ramp_rate*(t-t_load)^2/2))*unit_step(t-t_load)*unit_step(t_load + t_bloch - t) + 
        0.5*(t_load - t + t_load + t_bloch)/t_load*(cos(ramp_rate*t_bloch^2/2 + ramp_rate*t_bloch*
            (t - t_load - t_bloch)) + cos(2ramp_rate*t_bloch^2/2 + 2ramp_rate*t_bloch*
            (t - t_load - t_bloch)))*unit_step(t - t_load - t_bloch)*unit_step(2*t_load + t_bloch - t));
    return 🐢
end

#I want to try seeing multiple Bloch oscillations/ BO in a higher Bloch band

function V_BO(x,t,params)
    bdry = params[1];
    t_propogate = params[2];
    t_load = params[3];
    n_bloch = params[2];
    n_ramp = params[5];
    n_lat_dep = params[6];
    grav_ramp = params[7];
    ramp_rate = params[8];
    bloch_period = params[9];
    t_bloch = params[10];
    dv = params[11];
    dx = params[12];
    ω_0 = params[15];
    
    #Bloch oscillation potential
    
    🐢 = n_lat_dep * (cos(2*x + ω_0*t)*t/t_load*unit_step(t)*unit_step(t_load-t) +
        cos(2*x + ω_0*t_load + (ramp_rate*(t-t_load) + ω_0)*(t-t_load)/2)*unit_step(t-t_load)*unit_step(t_load + t_bloch - t) + 
        (t_load - t + t_load + t_bloch)/t_load*cos(2*x + ω_0*t_load + (ramp_rate*(t_bloch) + ω_0)*(t_bloch)/2 + (ramp_rate*t_bloch + ω_0)*
            (t - t_load - t_bloch))*unit_step(t - t_load - t_bloch)*unit_step(2*t_load + t_bloch - t));
    return 🐢
end

bdry = 55;
t_propogate = 10;
t_load = 5pi;
n_bloch = 3;
n_ramp = 0.2;
n_lat_dep = 5.0;
grav_ramp = 0.8576;
ramp_rate = grav_ramp * n_ramp;
bloch_period = 8/ramp_rate;
t_bloch = n_bloch*bloch_period;
dv = 0.05;
dx = 1/(2*dv);
n_reduce_x = 100;
n_reduce_t = 20;
ω_0 = 6;

params = [bdry,t_propogate,t_load,n_bloch,n_ramp,n_lat_dep,grav_ramp,ramp_rate,bloch_period,
    t_bloch,dv,dx,n_reduce_x,n_reduce_t,ω_0];

#xs = -bdry*dx:0.04:bdry*dx;
xs = -160*dx:0.04:15*dx;
ts = 0:0.005:(2*t_load + t_bloch + t_propogate);

println("Metrics on solving for wavefunction, including time to solve: ")
@time ψ = integrate_SE_BBS_lite(ψ0_cold,V_BO,xs,ts,params);
println("You used ", sizeof(ψ)/10^9, " Gb memory to store ψ.")

using PyPlot
#import Colors.delta

#ψ_new = reduce_size(ψ, 10,10)

#define weird x & y vectors for making the timing plot on the right. 
x = 1:length(ψ[:,1])
y = 1:length(ψ[1,:])
dy = .1;
dx = .1;

#density plot of dynamics
fig= figure("Plotting state evolution and timing",figsize=(20,20))

subplot(221)
#main plot
main = pcolormesh(abs.(ψ'),cmap = 🍎)
axis("off")
title("Bloch Oscillations!", fontsize = 24.0)
xlabel("Position")
ylabel("Time")

#subplot on the side for showing timing of sequence
subplot(264)
timing1 = PyPlot.plot([0,1],[0,0],color = :black)
timing2 = PyPlot.plot([0,1],[t_load,t_load],color = :black)
timing3 = PyPlot.plot([0,1],[t_load+t_bloch,t_load+t_bloch],color = :black)
timing4 = PyPlot.plot([0,1],[2*t_load+t_bloch,2*t_load+t_bloch],color = :black)
timing5 = PyPlot.plot([0,1],[2*t_load+t_bloch+t_propogate,2*t_load+t_bloch+t_propogate],color = :black)
#annotate("Look, data!",
#    xy=[x[convert(Int64,floor(length(x)/4.1))];y[convert(Int64,floor(length(y)/4.1))]],
#    xytext=[x[convert(Int64,floor(length(x)/4.1))]+0.1dx;y[convert(Int64,floor(length(y)/4.1))]+0.1dy],
#    xycoords="data")
xlim(0,1)
ylim(0 - 0.1,2*t_load+t_bloch+t_propogate + 0.1)
annotate("Frequency ramping",
    fontsize = 18.0,
    xy=[0,0],
    xytext=[-0.05,t_load + t_propogate/2 + 10],
    xycoords="data")
annotate("Lattice Load",
    fontsize = 18.0,
    xy=[0,0],
    xytext=[0.1,t_load + 15],
    xycoords="data")
annotate("Lattice Unload",
    fontsize = 18.0,
    xy=[0,0],
    xytext=[0.1,t_load + 15 + t_bloch + t_load],
    xycoords="data")
annotate("Free Propogation",
    fontsize = 18.0,
    xy=[0,0],
    xytext=[0,t_load + 15 + t_bloch + t_load + t_propogate/3],
    xycoords="data")
PyPlot.arrow(0.5,
    16,
    -0.05,
    -10,
    head_width=0.1,
    width=0.03,
    head_length=3,
    overhang=0,
    head_starts_at_zero="true",
    facecolor="black")
PyPlot.arrow(0.5,
    16 + t_bloch + t_load,
    -0.05,
    -10,
    head_width=0.1,
    width=0.03,
    head_length=3,
    overhang=0,
    head_starts_at_zero="true",
    facecolor="black")
axis("off")
fig[:canvas][:draw]() # Update the figure

#function to export the figure. Can put in whatever file extension you want for the image. Note that it's
# very tine consuming to export large image files. 
#@time PyPlot.savefig("BBS_very_high_lat_dep.pdf")

using LinearAlgebra

function integrate_SE_lite(ψ0,V,xs,ts,params)
    #dumb variables
    nTs = length(ts);
    nXs = length(xs);
    reduced_nXs = Int(floor(nXs/params[13]));
    reduced_nTs = Int(floor(nTs/params[14]));
    
    #current and previous state vectors.
    ψ_current = (1+1im).*zeros(nXs);
    ψ_previous = (1+1im).*zeros(nXs);
    ψf = (1+0.0im).*zeros((reduced_nXs, reduced_nTs));
    
    #Initial condition
    ψ_previous= (1+0im).*[ψ0(xs[j],params) for j=1:nXs];
    
    #more dumb varriables
    dx = xs[2] - xs[1];
    dt = ts[2] - ts[1];
    α = 0.5im*dt/dx^2;
    
    #This is where the 💰 is made
    for i=2:nTs
        V_i = [V(xs[j],ts[i],params) for j=1:nXs];
        ξ_i = 1 .+ (0.5im*dt).*(2/dx^2 .+ V_i);
        γ_i = 1 .- (0.5im*dt).*(2/dx^2 .+ V_i);
        U_1 = SymTridiagonal(ξ_i,-α.*ones(nXs-1));
        U_2 = SymTridiagonal(γ_i,α.*ones(nXs-1));
        
        #rewrite state vectors
        ψ_current = U_1\(U_2*ψ_previous);
        ψ_previous = ψ_current;
        
        #only store in output vector if timestep is a multiple of n_reduce_t, and only store spatial values
        #that are a multiple of n_reduce_x
        if i%params[14] == 0 && Int(i/params[14])<=reduced_nTs
            ψf[:,Int(i/params[14])] = [ψ_current[Int(i)] for i=params[13]:params[13]:params[13]*reduced_nXs];
        end
    end
    
    return ψf
end

function ψ0_hot(x,params)
    σ = params[12]    #Heisenberg limited wavefunction, velocity spread σ
    μ = 0;
    a = 1/sqrt(2pi*σ^2)*exp(-(x-μ)^2/(2σ^2));
    return a
end

#Let's do an adiabatic load, then hold for time t_bloch, then adiabatic unload, then free propogation
#Note that there is no force/acceleration in the potential: atoms are just sitting in a still lattice.
function V_lattice(x,t,params)
    bdry = params[1];
    t_propogate = params[2];
    t_load = params[3];
    n_bloch = params[2];
    n_ramp = params[5];
    n_lat_dep = params[6];
    grav_ramp = params[7];
    ramp_rate = params[8];
    bloch_period = params[9];
    t_bloch = params[10];
    dv = params[11];
    dx = params[12];
    
    #🐒 ftw
    
    🐒 = n_lat_dep * cos(2*x)*(t/t_load*unit_step(t)*unit_step(t_load-t) +
        1*unit_step(t-t_load)*unit_step(t_load + t_bloch - t) + 
        (t_load - t + t_load + t_bloch)/t_load*unit_step(t - t_load - t_bloch)*unit_step(2*t_load + t_bloch - t));
    return 🐒
end

bdry = 1000;
t_propogate = 15;
t_load = pi;
n_bloch = 2;
n_ramp = 1.0;
n_lat_dep = 5.0;
grav_ramp = 0.8576;      #for ramping due to gravity, bloch period is ~7 in these units
ramp_rate = grav_ramp * n_ramp;
bloch_period = 8/ramp_rate;
t_bloch = n_bloch*bloch_period;
dv = 1.0;
dx = 1/(2*dv);
n_reduce_x = 30;
n_reduce_t = 4;

params = [bdry,t_propogate,t_load,n_bloch,n_ramp,n_lat_dep,grav_ramp,ramp_rate,bloch_period,
    t_bloch,dv,dx,n_reduce_x,n_reduce_t];

xs = -bdry*dx:0.005:bdry*dx;
ts = 0:0.04:(2*t_load + t_bloch + t_propogate);

println("Metrics on solving for wavefunction, including time to solve: ")
@time ψ = integrate_SE_lite(ψ0_hot,V_lattice,xs,ts,params);
println("You just used up ", sizeof(ψ)/10^9, " Gb to store ψ.")

using PyPlot
#import Colors.delta

#define weird x & y vectors for making the timing plot on the right. 
x = 1:length(ψ[:,1])
y = 1:length(ψ[1,:])
dy = .1;
dx = .1;

#density plot of dynamics
fig= figure("Plotting state evolution and timing",figsize=(20,20))

subplot(221)
#main plot
main = pcolormesh(abs.(ψ'),cmap = 🍎)
axis("off")
title("Hot Atom in Lattice", fontsize = 24.0)
xlabel("Position")
ylabel("Time")

#subplot on the side for showing timing of sequence
subplot(264)
timing1 = PyPlot.plot([0,1],[0,0],color = :black)
timing2 = PyPlot.plot([0,1],[t_load,t_load],color = :black)
timing3 = PyPlot.plot([0,1],[t_load+t_bloch,t_load+t_bloch],color = :black)
timing4 = PyPlot.plot([0,1],[2*t_load+t_bloch,2*t_load+t_bloch],color = :black)
timing5 = PyPlot.plot([0,1],[2*t_load+t_bloch+t_propogate,2*t_load+t_bloch+t_propogate],color = :black)
#annotate("Look, data!",
#    xy=[x[convert(Int64,floor(length(x)/4.1))];y[convert(Int64,floor(length(y)/4.1))]],
#    xytext=[x[convert(Int64,floor(length(x)/4.1))]+0.1dx;y[convert(Int64,floor(length(y)/4.1))]+0.1dy],
#    xycoords="data")
xlim(0,1)
ylim(0 - 0.1,2*t_load+t_bloch+t_propogate + 0.1)
annotate("Frequency ramping",
    fontsize = 18.0,
    xy=[0,0],
    xytext=[-0.05,t_load + t_propogate/2 + 10],
    xycoords="data")
annotate("Lattice Load",
    fontsize = 18.0,
    xy=[0,0],
    xytext=[0.1,t_load + 15],
    xycoords="data")
annotate("Lattice Unload",
    fontsize = 18.0,
    xy=[0,0],
    xytext=[0.1,t_load + 15 + t_bloch + t_load],
    xycoords="data")
annotate("Free Propogation",
    fontsize = 18.0,
    xy=[0,0],
    xytext=[0,t_load + 15 + t_bloch + t_load + t_propogate/3],
    xycoords="data")
PyPlot.arrow(0.5,
    16,
    -0.05,
    -10,
    head_width=0.1,
    width=0.03,
    head_length=3,
    overhang=0,
    head_starts_at_zero="true",
    facecolor="black")
PyPlot.arrow(0.5,
    16 + t_bloch + t_load,
    -0.05,
    -10,
    head_width=0.1,
    width=0.03,
    head_length=3,
    overhang=0,
    head_starts_at_zero="true",
    facecolor="black")
axis("off")
fig[:canvas][:draw]() # Update the figure

#function to export the figure. Can put in whatever file extension you want for the image. Note that it's
# very tine consuming to export large image files. 
#@time PyPlot.savefig("BBS_very_high_lat_dep.pdf")

using LinearAlgebra

function integrate_SE_lite(ψ0,V,xs,ts,params)
    #dumb variables
    nTs = length(ts);
    nXs = length(xs);
    reduced_nXs = Int(floor(nXs/params[13]));
    reduced_nTs = Int(floor(nTs/params[14]));
    
    #current and previous state vectors.
    ψ_current = (1+1im).*zeros(nXs);
    ψ_previous = (1+1im).*zeros(nXs);
    ψf = (1+0.0im).*zeros((reduced_nXs, reduced_nTs));
    
    #Initial condition
    ψ_previous= (1+0im).*[ψ0(xs[j],params) for j=1:nXs];
    
    #more dumb varriables
    dx = xs[2] - xs[1];
    dt = ts[2] - ts[1];
    α = 0.5im*dt/dx^2;
    
    #This is where the 💰 is made
    for i=2:nTs
        V_i = [V(xs[j],ts[i],params) for j=1:nXs];
        ξ_i = 1 .+ (0.5im*dt).*(2/dx^2 .+ V_i);
        γ_i = 1 .- (0.5im*dt).*(2/dx^2 .+ V_i);
        U_1 = SymTridiagonal(ξ_i,-α.*ones(nXs-1));
        U_2 = SymTridiagonal(γ_i,α.*ones(nXs-1));
        
        #rewrite state vectors
        ψ_current = U_1\(U_2*ψ_previous);
        ψ_previous = ψ_current;
        
        #only store in output vector if timestep is a multiple of n_reduce_t, and only store spatial values
        #that are a multiple of n_reduce_x
        if i%params[14] == 0 && Int(i/params[14])<=reduced_nTs
            ψf[:,Int(i/params[14])] = [ψ_current[Int(i)] for i=params[13]:params[13]:params[13]*reduced_nXs];
        end
    end
    
    return ψf
end

function ψ0_hot(x,params)
    σ = params[12]    #Heisenberg limited wavefunction, velocity spread σ
    μ = 0;
    a = 1/sqrt(2pi*σ^2)*exp(-(x-μ)^2/(2σ^2));
    return a
end

#Let's do an adiabatic load, then hold for time t_bloch, then adiabatic unload, then free propogation
#Note that there is no force/acceleration in the potential: atoms are just sitting in a still lattice.
function V_lattice_hold(x,t,params)
    bdry = params[1];
    t_propogate = params[2];
    t_load = params[3];
    n_bloch = params[2];
    n_ramp = params[5];
    n_lat_dep = params[6];
    grav_ramp = params[7];
    ramp_rate = params[8];
    bloch_period = params[9];
    t_bloch = params[10];
    dv = params[11];
    dx = params[12];
    
    #Adding a force term into the potential. In our dimensionless equations, force is related to ramp rate
    # via F = r/2 (I've worked thorugh this separately, it's in a Mathematica notebook). I turn on the force 
    #during the lattice hold, turn it off otherwise. Should see Bloch oscillations
    F = ramp_rate/2;
    
    🐒 = F*x*unit_step(t-t_load)*unit_step(t_load + t_bloch - t) + 
        n_lat_dep * cos(2*x)*(t/t_load*unit_step(t)*unit_step(t_load-t) +
        1*unit_step(t-t_load)*unit_step(t_load + t_bloch - t) + 
        (t_load - t + t_load + t_bloch)/t_load*unit_step(t - t_load - t_bloch)*unit_step(2*t_load + t_bloch - t));
    return 🐒
end

bdry = 60;
t_propogate = 15;
t_load = pi;
n_bloch = 2;
n_ramp = 1.0;
n_lat_dep = 3.0;
grav_ramp = 0.8576;      #for ramping due to gravity, bloch period is ~7 in these units
ramp_rate = grav_ramp * n_ramp;
bloch_period = 8/ramp_rate;
t_bloch = n_bloch*bloch_period;
dv = 0.05;
dx = 1/(2*dv);
n_reduce_x = 10;
n_reduce_t = 2;

params = [bdry,t_propogate,t_load,n_bloch,n_ramp,n_lat_dep,grav_ramp,ramp_rate,bloch_period,
    t_bloch,dv,dx,n_reduce_x,n_reduce_t];

xs = -bdry*dx:0.05:bdry*dx;
ts = 0:0.04:(2*t_load + t_bloch + t_propogate);

println("Metrics on solving for wavefunction, including time to solve: ")
@time ψ = integrate_SE_lite(ψ0_hot,V_lattice_hold,xs,ts,params);
println("You just used up ", sizeof(ψ)/10^9, " Gb to store ψ.")

#Plot the hot atom Bloch oscillations

using PyPlot
#import Colors.delta

#define weird x & y vectors for making the timing plot on the right. 
x = 1:length(ψ[:,1])
y = 1:length(ψ[1,:])
dy = .1;
dx = .1;

#density plot of dynamics
fig= figure("Plotting state evolution and timing",figsize=(20,20))

subplot(221)
#main plot
main = pcolormesh(abs.(ψ'),cmap = 🍎)
axis("off")
title("Hot Atom Undulating", fontsize = 24.0)
xlabel("Position")
ylabel("Time")

#subplot on the side for showing timing of sequence
subplot(264)
timing1 = PyPlot.plot([0,1],[0,0],color = :black)
timing2 = PyPlot.plot([0,1],[t_load,t_load],color = :black)
timing3 = PyPlot.plot([0,1],[t_load+t_bloch,t_load+t_bloch],color = :black)
timing4 = PyPlot.plot([0,1],[2*t_load+t_bloch,2*t_load+t_bloch],color = :black)
timing5 = PyPlot.plot([0,1],[2*t_load+t_bloch+t_propogate,2*t_load+t_bloch+t_propogate],color = :black)
#annotate("Look, data!",
#    xy=[x[convert(Int64,floor(length(x)/4.1))];y[convert(Int64,floor(length(y)/4.1))]],
#    xytext=[x[convert(Int64,floor(length(x)/4.1))]+0.1dx;y[convert(Int64,floor(length(y)/4.1))]+0.1dy],
#    xycoords="data")
xlim(0,1)
ylim(0 - 0.1,2*t_load+t_bloch+t_propogate + 0.1)
annotate("Frequency ramping",
    fontsize = 18.0,
    xy=[0,0],
    xytext=[-0.05,t_load + t_propogate/2 + 10],
    xycoords="data")
annotate("Lattice Load",
    fontsize = 18.0,
    xy=[0,0],
    xytext=[0.1,t_load + 15],
    xycoords="data")
annotate("Lattice Unload",
    fontsize = 18.0,
    xy=[0,0],
    xytext=[0.1,t_load + 15 + t_bloch + t_load],
    xycoords="data")
annotate("Free Propogation",
    fontsize = 18.0,
    xy=[0,0],
    xytext=[0,t_load + 15 + t_bloch + t_load + t_propogate/3],
    xycoords="data")
PyPlot.arrow(0.5,
    16,
    -0.05,
    -10,
    head_width=0.1,
    width=0.03,
    head_length=3,
    overhang=0,
    head_starts_at_zero="true",
    facecolor="black")
PyPlot.arrow(0.5,
    16 + t_bloch + t_load,
    -0.05,
    -10,
    head_width=0.1,
    width=0.03,
    head_length=3,
    overhang=0,
    head_starts_at_zero="true",
    facecolor="black")
axis("off")
fig[:canvas][:draw]() # Update the figure

#function to export the figure. Can put in whatever file extension you want for the image. Note that it's
# very tine consuming to export large image files. 
#@time PyPlot.savefig("BBS_very_high_lat_dep.pdf")

#Someday this will plot bandstructure
n = 9;
N_Lat_Dep = 4;
Ω = N_Lat_Dep/4;
kinetic_energy = [2*(k-floor(n/2))^2 for k=1:n];

function ham_BBS(t, kinetic_energy, ramp_rate)
    Ham = Tridigiag(Ω.*exp(-ramp_rate * t^2 * 0.5im).*ones(n), kinetic_energy, Ω.*exp(ramp_rate * t^2 * 0.5im).*ones(n));
    return Ham
end

ham_lattice = Tridiagonal()



M = SymTridiagonal(2.0 .*ones(10),ones(9))
M
