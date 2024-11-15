using DrWatson
@quickactivate "Rydberg-Quantum-Simulator"
using MKL
using ITensors
using LinearAlgebra
using ITensorMPS
using Statistics

ITensors.Strided.disable_threads()
BLAS.set_num_threads(2)
ITensors.enable_threaded_blocksparse(false)

@show Threads.nthreads()
@show Sys.CPU_THREADS
@show BLAS.get_num_threads()
@show ITensors.Strided.get_num_threads()
@show ITensors.using_threaded_blocksparse()
println()

function get_ee(ψ::AbstractMPS, si::Int, sites::Vector{Index{Int64}})::Float64
    orthogonalize!(ψ, si)
    _, s = svd(ψ[si], (linkind(ψ, si-1), sites[si]))
    SvN = 0.0
    for n in 1:dim(s, 1)
        p = s[n,n]^2
        SvN -= p*log(p)
    end
    return SvN
end

mutable struct EnergyObserver <: AbstractObserver
   energy_tol::Float64
   last_energy::Float64

   EnergyObserver(energy_tol=0.0) = new(energy_tol,1000.0)
end

function ITensorMPS.checkdone!(o::EnergyObserver;kwargs...)
    sw = kwargs[:sweep]
    energy = kwargs[:energy]
    if abs(energy-o.last_energy)/abs(energy) < o.energy_tol
      println("Stopping DMRG after sweep $sw")
      return true
    else
        # Otherwise, update last_energy and keep going
        println("ratio = $(abs(energy-o.last_energy)/abs(energy))")
        o.last_energy = energy
        return false
    end
end

function run1job(Ω::Float64, Δ::Float64, L::Int64, )
    sites = siteinds("S=1/2", L)
    trunc = 5
    os = OpSum()
    for j = 1:L
        # Single-qubit terms
        os +=  Ω,     "Sx", j
        os += -Δ, "ProjUp", j

        # Interaction terms
        for i = j+1:min(j+trunc, L)
            Vij = -1/(i - j)^6
            os += Vij, "ProjUp", j, "ProjUp", i
        end
    end
    H = MPO(os, sites)
    psi0 = randomMPS(sites, 5)

    nsweeps = 128
    cutoff = 1E-8
    maxdim = [10, 20, 50, 80, 100, 200, 300, 400, 500]

    obs = EnergyObserver(1E-13)
    _, psi = dmrg(H,psi0; nsweeps, observer=obs, cutoff, maxdim)

    ee = get_ee(psi, div(L,2), sites)
    ni = mean(expect(psi, "ProjUp"))

    return ee, ni
end

let
    L = parse(Int64, ARGS[3])
    # lsΩ = LinRange(.4,.8,51)
    lsΩ = 0.45:0.008:1.5
    nlsΩ = length(lsΩ)
    # lsΔ = LinRange(-1.09,-.95,51)
    di = parse(Float64, ARGS[1])
    df = parse(Float64, ARGS[2])
    dd = 0.014
    lsΔ = di:dd:df
    nlsΔ = length(lsΔ)
    matee = zeros(Float64, nlsΩ, nlsΔ)
    matni = zeros(Float64, nlsΩ, nlsΔ)
    for Ωi in eachindex(lsΩ), Δi in eachindex(lsΔ)
        matee[Ωi, Δi], matni[Ωi, Δi] = run1job(lsΩ[Ωi], lsΔ[Δi], L)
    end

    rslt = @strdict matee matni lsΩ lsΔ
    fnstr = @strdict di df L
    fn = savename(fnstr,"jld2")
    wsave(datadir("sims","largeRegion",fn), rslt)
end
