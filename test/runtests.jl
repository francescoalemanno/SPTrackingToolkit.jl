using Test
using Random
using SPTrackingToolkit
using SPTrackingToolkit: buildlinks
function synthframe(t::Int,N::Int,ndims::Int)
    rng=MersenneTwister(215234)
    P=kron(1:ndims,(1:5:(5*N))').*1.0
    Δ=0.4
    for f in 1:t
        P=P.+randn(rng,ndims,N)*Δ
    end
    P
end
function TEST(ndims)
    M=[[1, -1, 1, -1, 1], [1, 1, 2, 1, 2], [1, 2, 3, 2, 3], [1, 3, 4, 3, 4], [1, 4, 5, 4, -1]]
    specs = SPT(maxtimegap=10,maxdist=10,dims=ndims,verbose=false)
    frames=[
            synthframe(10,5,ndims)[:,2:end],
            synthframe(11,5,ndims),
            synthframe(12,5,ndims)[:,2:end],
            synthframe(13,5,ndims)[:,1:(end-1)]
            ]
    @test length(track(specs,frames)) == 5
    @test buildlinks(specs,frames) == M
end
@testset "Simple Low dimensional Test" begin
    for i in 2:4
        TEST(i)
    end
end
