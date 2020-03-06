module SPTrackingToolkit
    using DefaultArrays
    using SolveLAP
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
        C=DefaultArray(Inf,pa+pb,pa+pb)
        mincost=Inf
        maxcost=0.0

        for i in 1:pa, j in 1:pb
            cost=norm(Ma[config.X,i] .- Mb[config.X,j])
            cost <= config.maxdist || continue
            C[i,j] = cost
            mincost=min(cost,mincost)
            maxcost=max(cost,maxcost)
        end

        for i in (pa+1):(pa+pb), j in (pb+1):(pa+pb)
            iT=j-pb
            jT=i-pa
            if C[iT,jT] != C.default
                C[i,j]=mincost
            end
        end
        for i in (pb+1):(pa+pb)
            C[i-pb,i]=maxcost*1.05
        end
        for i in (pa+1):(pa+pb)
            C[i,i-pa]=maxcost*1.05
        end
        solution=[1:(pa+pb) solve_lap(C)[1] zeros(Int,pa+pb)]
        for i in 1:(pa+pb)
            ia,ib,_=solution[i,:]
            if ia<=pa && ib<=pb
                solution[i,3]=COMPLETEMARK
            end
            if ia<=pa && ib>pb
                solution[i,3]=ENDMARK
                solution[i,2]=MISSINGFRAME
            end
            if ia>pa && ib<=pb
                solution[i,3]=STARTMARK
                solution[i,1]=MISSINGFRAME
            end
            if ia>pa && ib>pb
                solution[i,3]=FAKEMARK
                solution[i,2]=FAKEMARK
                solution[i,1]=FAKEMARK
            end
        end
        solution
    end
    @inline function tC(v::Vector)
        v[1]+length(v)-2
    end
    function buildtracksegments(config::SPT,seg)
        tracks=Vector{Int}[]
        for x in eachrow(seg[1])
            x[end]==FAKEMARK && continue
            push!(tracks,[1;x[1:2]])
        end
        for t in 2:length(seg)
            for x in eachrow(seg[t])
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
        C=DefaultArray(Inf,m+n,m+n)
        maxG=0.0
        minG=Inf

        for i in 1:m, j in 1:n  ## GAP filling block
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
                cost=norm(pend.-pstart)/sqrt(Δt)
                cost <= config.maxdist || continue
                C[i,j] = cost
                maxG=max(maxG,C[i,j])
                minG=min(minG,C[i,j])
            end
        end

        for i in 1:m, j in 1:n  ## no GAP filling block
            if isinf(C[i,j])
               continue
            end
            C[j+m,i+n]=minG
        end

        for i in 1:n  ## no GAP filling block
            C[i+m,i]=maxG*1.05
        end
        for i in 1:m  ## no GAP filling block
            C[i,i+n]=maxG*1.05
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
        C= buildcost_closegaps(config::SPT,frames,tracks)
        S=solve_lap(C)
        newlinks=[(endseg[r[1]],startseg[r[2]]) for r in eachrow([1:(n+m) S[1]]) if all(r.<=(m,n))]
        linktr=Vector{Int}[]
        forbidden=Set{Int}()
        if length(newlinks)>0
            for l in newlinks
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
            println("GAP closing terminated")
            return tracks
        end
        for i in 1:length(tracks)
            if i in forbidden
                continue;
            end
            push!(linktr,tracks[i])
        end
        tracks = linktr
        println("GAP closing is being restarted: closed ",div(length(forbidden),2), " gaps")
        @goto restartclosing
    end
    function get_tracked(specs::SPT,frames,links)
        allP=Array{Float64,2}[]
        for l in links
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
