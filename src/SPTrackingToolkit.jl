module SPTrackingToolkit
    using DefaultArrays
    using LapSolve
    using LinearAlgebra
    export track
    export SPT

    const ENDMARK=-1
    const STARTMARK=1
    const FAKEMARK=-100
    const COMPLETEMARK=0
    const MISSINGFRAME=-1

    Base.@kwdef struct SPT
        maxtimegap::Int
        maxdist::Float64
        dims::Int
        rest::Int=0
        X::UnitRange{Int}=1:dims
        I::UnitRange{Int}=(dims+1):(dims+rest)
        verbose::Bool=false
    end

    function segments(config::SPT,Ma,Mb)
        pa=size(Ma,2)
        pb=size(Mb,2)
        C=DefaultArray(Inf,pa,pb)

        @inbounds for i in 1:pa, j in 1:pb
            cost=norm(Ma[config.X,i] .- Mb[config.X,j])
            cost <= config.maxdist || continue
            C[i,j] = cost
        end

        S=solve_stiff_lap(C,1.05)
        solution=zeros(Int,length(S),3)

        @inbounds for i in eachindex(S)
            ia, ib = S[i]
            solution[i,1]=ia
            solution[i,2]=ib
            solution[i,3]=COMPLETEMARK
            ia<0 && (solution[i,3]=STARTMARK)
            ib<0 && (solution[i,3]=ENDMARK)
            (ia<0 && ib<0) && (solution[i,3]=FAKEMARK)
        end

        solution
    end
    @inline function tC(v::Vector)
        v[1]+length(v)-2
    end
    function buildtracksegments(config::SPT,seg)
        tracks=Vector{Int}[]
        @inbounds for x in eachrow(seg[1])
            x[end]==FAKEMARK && continue
            push!(tracks,[1;x[1:2]])
        end
        @inbounds for t in 2:length(seg)
            @inbounds for x in eachrow(seg[t])
                x[end]==FAKEMARK && continue
                x[end]==STARTMARK && begin push!(tracks,[t;x[1:2]]); continue; end

                for iT in eachindex(tracks)
                    (tC(tracks[iT])!=t) && continue;
                    if tracks[iT][end]==x[1]
                        push!(tracks[iT],x[2])
                        break
                    end
                end
            end
        end
        tracks
    end

    function buildcost_closegaps(config::SPT,frames,tracks)
        startseg=findall(x->x[2]==MISSINGFRAME,tracks)
        endseg=findall(x->x[end]==MISSINGFRAME,tracks)

        n=length(startseg)
        m=length(endseg)
        C=DefaultArray(Inf,m,n)
        @inbounds for i in 1:m, j in 1:n  ## GAP filling block
            ti=endseg[i]
            tj=startseg[j]

            if ti==tj
                continue
            end
            frameendidx = findlast(x->x != MISSINGFRAME,tracks[ti][2:end])+1
            framestartidx = findfirst(x->x != MISSINGFRAME,tracks[tj][2:end])+1
            Tend = frameendidx-2+tracks[ti][1]
            Tstart= framestartidx-2+tracks[tj][1]

            Δt=Tstart-Tend-1
            if 0< Δt <= config.maxtimegap
                pend=frames[Tend][config.X,tracks[ti][frameendidx]]
                pstart=frames[Tstart][config.X,tracks[tj][framestartidx]]
                cost=norm(pend.-pstart)/sqrt(Δt+1)
                cost <= config.maxdist || continue
                C[i,j] = cost
            end
        end
        return C
    end
    function closegaps(config::SPT,frames,origtracks)
        tracks=origtracks
        @label restartclosing
        startseg=findall(x->x[2]==MISSINGFRAME,tracks)
        endseg=findall(x->x[end]==MISSINGFRAME,tracks)
        n=length(startseg)
        m=length(endseg)
        C=buildcost_closegaps(config::SPT,frames,tracks)
        S=solve_stiff_lap(C)
        newlinks=[(endseg[a],startseg[b]) for (a,b) in S if (0<a<=m && 0<b<=n)]
        config.verbose && println("Started GAP closing: ", length(newlinks)," gaps identified")
        linktr=Vector{Int}[]
        forbidden=Set{Int}()
        if length(newlinks)>0
            @inbounds for l in newlinks
                ti,tj=l
                if ti in forbidden || tj in forbidden
                    continue
                end
                frameendidx = findlast(x->x != MISSINGFRAME,tracks[ti][2:end])+1
                framestartidx = findfirst(x->x != MISSINGFRAME,tracks[tj][2:end])+1
                Tend = frameendidx-2+tracks[ti][1]
                Tstart= framestartidx-2+tracks[tj][1]
                Δt=Tstart-Tend-1
                A=tracks[ti][1:frameendidx]
                B=tracks[tj][framestartidx:end]
                newtrack=[A;fill(MISSINGFRAME,Δt);B]
                push!(linktr,newtrack)
                push!(forbidden,ti)
                push!(forbidden,tj)
            end
        else
            config.verbose && println("GAP closing terminated")
            return tracks
        end
        for i in 1:length(tracks)
            if i in forbidden
                continue;
            end
            push!(linktr,tracks[i])
        end
        tracks = linktr
        config.verbose && println("GAP closing is being restarted: closed ",div(length(forbidden),2), " gaps")
        @goto restartclosing
    end
    function get_tracked(specs::SPT,frames,links)
        allP=Array{Float64,2}[]
        @inbounds for l in links
            t=l[1]
            idx=1
            n_particleinfo=length(specs.X)+length(specs.I)
            P=zeros(n_particleinfo+1,sum(x>0 for x in l)-1)
            for i in 2:length(l)
                f=l[i]
                if f>0
                    P[:,idx].=[t;frames[t][:,f]]
                    idx+=1
                end
                t+=1

            end
            push!(allP,P)
        end
        sort!(allP,by=x->-size(x,2))
        allP
    end
    function buildlinks(config::SPT,frames)
        l=length(frames)
        seg=[segments(config,frames[fa],frames[fb]) for (fa,fb) in zip(1:(l-1),2:l)]
        closegaps(config,frames,buildtracksegments(config,seg))
    end
    function track(config::SPT,frames)
        links=buildlinks(config,frames)
        get_tracked(config,frames,links)
    end
end # module
