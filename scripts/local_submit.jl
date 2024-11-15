using DrWatson
@quickactivate "Rydberg-Quantum-Simulator"

function submit1job(di,df)
    L = 256
    nthreads = 2
    mem_hint = 8

    runfn = scriptsdir("DMRGmain.jl")
    fnstr = @strdict di df L
    fn = savename(fnstr,"o")
    rsltfdn = datadir("sims","largeRegion")
    mkpath(rsltfdn)

    fullfn = joinpath(rsltfdn, savename(fnstr,"jld2"))
    if !isfile(fullfn)
        fnw = savename(fnstr)
        bashContent = "#!/bin/bash
MKL_NUM_THREADS=$nthreads nohup julia --heap-size-hint=$(mem_hint)G $runfn $di $df $L > $(rsltfdn)/$fn 2>&1 &"
        tmpfdn = projectdir("_research","tmp")
        mkpath(tmpfdn)
        bashfn = projectdir("_research","tmp", fnw*".sh")
        io = open(bashfn, "w")
        println(io, bashContent)
        close(io)

        run(`chmod u+rwx $bashfn`)
        run(`bash $bashfn`)
        rm(bashfn)
        sleep(1)
        println("submitted: $fullfn")
    else
        println("exist: $fullfn")
    end
end

let
    lsΔ = -1.03:0.014:-0.56
    matΔ = reshape(lsΔ, 2, :)
    for idx in axes(matΔ,2)
        submit1job(matΔ[1,idx], matΔ[2,idx])
    end
end
