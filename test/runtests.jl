using Aqua
using GeometricPreconditioner
using Test

#
@testset verbose = true showtiming = true "Base" begin
    for name ∈ (:aqua, :default)
        @testset "$(name)" begin
            test_name = Symbol(:test_, name)
            include("$(test_name).jl")
            @eval $test_name()
        end
    end
end
