# The tutorial

#### Load the packages
using Quantica
using GLMakie

### 1A. Lattices
#lattice(sublats, ...)

sublat_A = sublat((0.0, -0.5/sqrt(3.0)), name = :A)
sublat_B = sublat((0.0, 0.5/sqrt(3.0)), name = :B)

lat = lattice(sublat_A, sublat_B)

qplot(lat, siteradius = 0.1)

#### Periodic lattices
#lattice(sublats, bravais)

lat0 = lattice(sublat((0.0, -0.5/sqrt(3.0)), name = :A),
            sublat((0.0,  0.5/sqrt(3.0)), name = :B));
lat = lattice(lat0; 
    bravais = SA[cos(pi/3) sin(pi/3); -cos(pi/3) sin(pi/3)]')
qplot(lat, siteradius = 0.1)

##### Supercells out of an unitcell
lat0 = LatticePresets.honeycomb(a0 = 1);
lat = lat0 |> supercell(9);
qplot(lat, siteradius = 0.1)

# 1D 2D and 3D
qplot(LP.hcp(; a0 = 1)|> supercell(5))

### Periodic system along one axis

qplot(lat0, siteradius = 0.1)
rot_lat = lat0 |> supercell((1, -1), (1, 1))
qplot(rot_lat, siteradius = 0.1);

# Define a rotation axis
axis = (0,1);
perpaxis = 
    Quantica.normalize(
        Quantica.bravais(lat0).matrix * 
        SA[1 1; -1 1] * SA[-axis[2], axis[1]])
 arm_lat = rot_lat |> supercell((0,1), 
        region = r -> abs(dot((-1,0), r)) < 10)
qplot(arm_lat, siteradius = 0.03)

zz_lat = rot_lat |> supercell((1,0),
        region = r -> abs(dot((0,1), r)) < 10)

qplot(zz_lat, siteradius = 0.03)

### Bounded Lattices
#missing bravais makes the system bounded

lat = lat0 |> supercell(; region = r -> norm(r) < 30)
qplot(lat, siteradius = 0.1)

circular_with_a_hole(r) = (2<norm(r)< 4);
lat = lat0 |> supercell(;region = r -> circular_with_a_hole(r))
qplot(lat, siteradius = 0.1)

#### A summary: ZZ nanoribbon of bilayer graphene in Bernal stacking
function latBLG_unbounded(a0 = 1)
    dinter =  1.36 * a0
    sAbot = sublat((0.0,-1.0a0/sqrt(3.0), - dinter/2); name = :Ab)
    sBbot = sublat((0.0, 0.0a0/sqrt(3.0), - dinter/2); name = :Bb)
    sAtop = sublat((0.0, 0.0a0/sqrt(3.0), + dinter/2); name = :At)
    sBtop = sublat((0.0, 1.0a0/sqrt(3.0), + dinter/2); name = :Bt)
    br = a0 * SA[cos(pi/3) sin(pi/3) 0; -cos(pi/3) sin(pi/3) 0]'
    lat = lattice(sAtop, sBtop, sAbot, sBbot; bravais = br)
    return lat
end
latBLG_unbounded() 
zz_blg = latBLG_unbounded()  |> 
    supercell((1, -1), (1, 1)) |> 
    supercell((1,0), region = r -> abs(dot([0,1,1], r)) < 10)
zz_blg = latBLG_unbounded() |> supercell((1, -1), (1, 1)) 
    |> supercell((1,0), region = r -> abs(dot([0,1,1], r)) < 10 && -2<abs(dot([1,0,1], r))<2)
zz_blg_bounded = zz_blg |> 
    supercell(region = r -> abs(dot([0,1,1], r)) < 10 && -2<abs(dot([1,0,1], r))<2)


## 1B. The model 

##### Spinless single-orbital system
a0 = 0.246

# lattice
lat_graph = LP.honeycomb(; a0 = a0) |> supercell(3)
# some params
t0 = 2.7
V0 = 5
range = a0/sqrt(3) 
β = 3;

# model
ons = onsite(V0; sublats = :A) + # staggered potential between A and B lattices
    onsite(-V0; sublats = :B) + 
    onsite((r) -> rand(1), sublats = (:A,:B));#random non-magnetic disorder in both lattices
hops = hopping((r, dr) -> t0 * exp(-β*(norm(dr)/a0 - 1)) * I, range = range);
model = ons + hops;
# hamiltonian
ham = hamiltonian(lat_graph, model);
qplot(ham, siteradius = 0.03, inspector = true)

## 1C. Worked out examples
### Magnetic systems: Peierls phases. Parametric Hamiltonians
model = hopping(1); peierls = @hopping!((t, r, dr; A = r -> SA[0,0]) -> t * cis(-dr' * A(r)));
h = LP.honeycomb() |> hamiltonian(model) |> supercell(10) |> hamiltonian(peierls);
h(t = 2, A = Returns(SA[1,1]));

### The Kane-Mele model
SOC(dr) = ifelse(
            iseven(
                round(Int, atan(dr[2], dr[1])/(pi/3))), im, -im); 
# Kane-Mele spin-orbit coupling'
model = hopping(1, range = 1/√3) +
     @hopping((r, dr; α = 0) -> α * SOC(dr); 
    sublats = :A => :A, range = 1) - 
    @hopping((r, dr; α = 0) -> α * SOC(dr); sublats = :B => :B, range = 1);
h = LP.honeycomb(a0 = 1) |> hamiltonian(model);
qplot(h(α = 0.02), inspector = true)

##### Bandstructure
`bands(h::Hamiltonian{...}, ...; kwargs...)`
b = bands(h(α = 0.05), range(0, 2pi, length=60), range(0, 2pi, length = 60))
qplot(b, color = (psi, e, k) -> angle(psi[1] / psi[2]),
     colormap = :cyclic_mrybm_35_75_c68_n256)

qplot(b, color = (psi, e, k) -> angle(psi[1] / psi[2]), 
    colormap = :cyclic_mrybm_35_75_c68_n256, hide = :points)
b = bands(h(α = 0.0), polypath(0,4,101); mapping = (:M, :K, :M, :K´,:M));

qplot(b, color = (psi, e, k) -> angle(psi[1] / psi[2]), hide = :points)

b = bands(h(α = 0.05), polypath(0,4,101); mapping = (:M, :K, :M, :K´,:M));

qplot(b, color = (psi, e, k) -> angle(psi[1] / psi[2]), hide = :points)


### Spectrum
h = HP.graphene() |> supercell(50) |> supercell()
es, phis = spectrum(h,  solver = EigenSolvers.LinearAlgebra());
hist(real.(es), bins = 2000)

### Transport in a Bogoliubov de Gennes system. A SNS junction
const σz = SA[1 0; 0 -1]
const σx = SA[0 1; 1 0]
const L = 10
const W = 10

#Here I build a BdG system consisting of two superconducting leads (S) and a normal 
#(scattering) region (N) SNS junction 


# The scattering or central region
h = LP.square() |> hamiltonian(hopping(σz), orbitals = 2) |> 
    supercell(region = RP.rectangle((2L, 2W)))

# Left and right leads
hl = LP.square() |> hamiltonian(hopping(σz) + onsite(0.1 * σx), orbitals = 2) |> 
    supercell((-1,0), region = r -> -W <= r[2] <= W)
hr = LP.square() |> hamiltonian(hopping(σz) + onsite(0.1 * σx), orbitals = 2) |> 
    supercell((1,0), region = r -> -W <= r[2] <= W)
    
# Green's functions of the lattice
gl = hl |> greenfunction(GS.Schur(boundary = -L))
gr = hl |> greenfunction(GS.Schur(boundary = L))

#Integrating out the leads
g = h |> attach(gl; region = r -> r[1] == -L) |> attach(gr; region = r -> r[1] == L) 

# System's dressed GF
G = g |> greenfunction(GS.SparseLU())

# Plot
qplot(g, inspector = true)

### Observables
# Josephsons
j = josephson(g, 4.5; phases = 100)
@time jlist = j()
scatter(jlist, xlabel = "\phi", ylabel = "J_c(\phi)")



#σ = conductance(g[i,j])
#ldos(g[i,j], kwargs...)
#current(g[i,j], kwargs...)

### Moiré system    
function my_twisted_bilayer_graphene(;
    twistindex = 1, twistindices = (twistindex, 1), a0 = 0.246, interlayerdistance = 1.36a0,
    rangeintralayer = a0/sqrt(3), rangeinterlayer = 4a0/sqrt(3), hopintra = 2.70,
    hopinter = 0.48, modelintra = hopping(hopintra, range = rangeintralayer), kw...)

    # Sublattices
    (m, r) = twistindices
    θ = acos((3m^2 + 3m*r +r^2/2)/(3m^2 + 3m*r + r^2))
    sAbot = sublat((0.0, -0.5a0/sqrt(3.0), - interlayerdistance / 2); name = :Ab)
    sBbot = sublat((0.0,  0.5a0/sqrt(3.0), - interlayerdistance / 2); name = :Bb)
    sAtop = sublat((0.0, -0.5a0/sqrt(3.0),   interlayerdistance / 2); name = :At)
    sBtop = sublat((0.0,  0.5a0/sqrt(3.0),   interlayerdistance / 2); name = :Bt)
    brbot = a0 * SA[ cos(pi/3) sin(pi/3) 0; -cos(pi/3) sin(pi/3) 0]'
    brtop = a0 * SA[ cos(pi/3) sin(pi/3) 0; -cos(pi/3) sin(pi/3) 0]'
    # Supercell matrices sc.
    # The one here is a [1 0; -1 1] rotation of the one in Phys. Rev. B 86, 155449 (2012)
    if gcd(r, 3) == 1
        scbot = SA[m -(m+r); (m+r) 2m+r] * SA[1 0; -1 1]
        sctop = SA[m+r -m; m 2m+r] * SA[1 0; -1 1]
    else
        scbot = SA[m+r÷3 -r÷3; r÷3 m+2r÷3] * SA[1 0; -1 1]
        sctop = SA[m+2r÷3 r÷3; -r÷3 m+r÷3] * SA[1 0; -1 1]
    end
    latbot = lattice(sAbot, sBbot; bravais = brbot)
    lattop = lattice(sAtop, sBtop; bravais = brtop)
    
    # Atomistic Hamiltonian of BLG
    htop = hamiltonian(lattop, modelintra; ) |> supercell(sctop)
    hbot = hamiltonian(latbot, modelintra; ) |> supercell(scbot)
    let R = SA[cos(θ/2) -sin(θ/2) 0; sin(θ/2) cos(θ/2) 0; 0 0 1]
        Quantica.transform!(htop, r -> R * r)
    end
    let R = SA[cos(θ/2) sin(θ/2) 0; -sin(θ/2) cos(θ/2) 0; 0 0 1]
        Quantica.transform!(hbot, r -> R * r)
    end
    modelinter = hopping((r,dr) -> (
        hopintra * exp(-3*(norm(dr)/a0 - 1))  *  dot(dr, SVector(1,1,0))^2/sum(abs2, dr) -
        hopinter * exp(-3*(norm(dr)/a0 - interlayerdistance/a0)) * dr[3]^2/sum(abs2, dr)),
        range = rangeinterlayer)
    return combine(hbot, htop)#; coupling = modelinter)
end

h = my_twisted_bilayer_graphene()

b = bands(h, range(0, 2pi, length=60), range(0, 2pi, length = 60))
qplot(b, color = (psi, e, k) -> angle(psi[1] / psi[2]), 
    colormap = :cyclic_mrybm_35_75_c68_n256)
