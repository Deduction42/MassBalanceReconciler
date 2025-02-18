using MassBalanceReconciler
using Test

#For coverage, run with
#julia --startup-file=no --depwarn=yes --threads=auto -e 'using Coverage; clean_folder(\"src\"); clean_folder(\"test\"); clean_folder(\"ext\")'
#julia --startup-file=no --depwarn=yes --threads=auto --code-coverage=user --project=. -e 'using Pkg; Pkg.test(coverage=true)'
#julia --startup-file=no --depwarn=yes --threads=auto coverage.jl

@testset "MassBalanceReconciler.jl" begin
    # Write your tests here.
end
