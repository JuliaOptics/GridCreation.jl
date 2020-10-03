using Test
using GridCreation



@testset "fft alignment" begin
    xfft, yfft = mkCartVecs(1., 4, center=CFFT);
    # alignment
    @test xfft[3] ≈ 0
    @test yfft[3] ≈ 0

    # spacing
    @test xfft[1] ≈ -2
    @test yfft[1] ≈ -2
end


@testset "intersample alignment" begin
    xinter, yinter = mkCartVecs(1., 4, center=CInterSample);
    @test xinter[2] ≈ -0.5
    @test xinter[3] ≈  0.5
    @test yinter[2] ≈ -0.5
    @test yinter[3] ≈  0.5
end

@testset "first sample centering" begin
# test sample index centering is properly aligned
    x0, y0 = mkCartVecs(1., 4, center=1);
    @test x0[1] ≈ 0
    @test y0[1] ≈ 0
end

@testset "polar grid" begin
    xfft, yfft = mkCartVecs(1., 4, center=CFFT);
    ρ, θ = cartVecsToPolarGrid(xfft, yfft);
    @test ρ[3,3] ≈ 0
    @test ρ[3,2] ≈ 1
    @test ρ[2,3] ≈ 1
    @test ρ[3,4] ≈ 1
    @test ρ[4,3] ≈ 1
    @test ρ[4,4] ≈ √2
    @test ρ[2,2] ≈ √2
    @test ρ[2,4] ≈ √2
    @test ρ[4,2] ≈ √2
    @test θ[3,3] ≈ 0
    @test θ[3,4] ≈ 0
    @test θ[3,2] ≈ π
    @test θ[3,1] ≈ π
    @test θ[2,3] ≈ -π/2
    @test θ[4,3] ≈ π/2
end
