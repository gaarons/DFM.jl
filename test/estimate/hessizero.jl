using DSGE
using Base.Test

# Test `hessizero` in context of Rosenbrock function
function rosenbrock(x::Vector)
    a = 100
    b = 1
    return a*(x[2]-x[1]^2.0)^2.0 + b*(1-x[1])^2.0
end

function rosenbrock_hessian(x::Vector)
    H = zeros(eltype(x), 2, 2)
    H[1,1] = 2.0 - 400.0 * x[2] + 1200.0 * x[1]^2.0
    H[1,2] = -400.0 * x[1]
    H[2,1] = -400.0 * x[1]
    H[2,2] = 200.0
    return H
end

# At the min, ensure no negatives in diagonal
x0 = [1.0, 1.0]
hessian_expected = rosenbrock_hessian(x0)
hessian, _ = DSGE.hessizero(rosenbrock, x0; check_neg_diag=true)
@test_matrix_approx_eq hessian_expected hessian

# Not at the min (indeed, we are at max), ensure throws error
x1 = [1.0, 1.0]
rosenbrock_neg(x) = -rosenbrock(x)
@test_throws Exception hessian, _ = DSGE.hessizero(rosenbrock_neg, x1; check_neg_diag=true)

# Not at the min, check matches closed form
x1 = Vector[[0.5, 1.5], [-1.0, -1.0]]
for x in x1
    hessian_expected = rosenbrock_hessian(x)
    hessian, _ = DSGE.hessizero(rosenbrock, x)
    @test_matrix_approx_eq hessian_expected hessian
end
