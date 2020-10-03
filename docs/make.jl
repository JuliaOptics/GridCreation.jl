using Documenter, GridCreation

makedocs(modules = [GridCreation], doctest = false, sitename = "GridCreation")
deploydocs(repo = "github.com/JuliaOptics/GridCreation.jl.git")
