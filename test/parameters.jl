using DSGE

for T in subtypes(Transform)
    u = parameter(:σ_pist, 2.5230, (1e-8, 5.), (1e-8, 5.), T(), fixed=false)
    @test ( transform_to_real_line(u) |> x -> transform_to_model_space(u,x) ) == u.value
    
    if !isa(T,Type{DSGE.Untransformed})
        # check transform_to_real_line and transform_to_model_space to different things if T is not DSGE.Untransformed
        @test transform_to_real_line(u,u.value) != transform_to_model_space(u,u.value)
    end
end

# probability
N = 10^2
u = parameter(:bloop, 2.5230, (1e-8, 5.), (1e-8, 5.), DSGE.SquareRoot(); fixed = true)
v = parameter(:cat, 2.5230, (1e-8, 5.), (1e-8, 5.), DSGE.Exponential(), Gamma(2.00, 0.1))

pvec =  ParameterVector{Float64}(N)
for i in 1:length(pvec)
	pvec[i] = (i%2 == 0) ? u : v
end
@test_approx_eq logpdf(pvec) 50*logpdf(v)
@test_approx_eq pdf(pvec) exp(50*logpdf(v))

updated = update(pvec, ones(length(pvec)))
update!(pvec, ones(length(pvec)))

@test all(updated .== pvec)
@test logpdf(pvec) == logpdf(updated)

# test we only update unfixed parameters
for p in pvec
	if p.fixed
		@test p.value == 2.5230
	elseif isa(p, Parameter)
		@test p.value == one(Float64)
	end
end

# vector of new values must be the same length
@test_throws AssertionError update!(pvec, ones(length(pvec)-1))

for w in [parameter(:moop, 3.0, fixed=false), parameter(:moop, 3.0; scaling = log, fixed=false)]
	# new values must be of the same type
	@test_throws MethodError parameter(w, one(Int))

	# new value is out of bounds
	@test_throws ParamBoundsError parameter(w, -1.)
end

# subspecs
function sstest(m::Model990)

    # Change all the fields of an unfixed parameter
    m <= parameter(:ι_w, 0.000, (0.0, .9999), (0.0,0.9999), DSGE.Untransformed(), Normal(0.0,1.0), fixed=false,
                   description="ι_w: No description available.",
                   tex_label="\\iota_w")


    # Change an unfixed parameter to be fixed
    m <= parameter(:ι_p, 0.000, fixed=true,
                   description= "ι_p: The persistence of last period's inflation in
                   the equation that describes the intertemporal
                   change in prices for intermediate goods producers
                   who cannot adjust prices. The change in prices is a
                   geometric average of steady-state inflation
                   (π_star, with weight (1-ι_p)) and last period's
                   inflation (π_{t-1})).",
                   tex_label="\\iota_p")


    # Change a fixed parameter
    m <= parameter(:δ, 0.02,  fixed=true,
                   description="δ: The capital depreciation rate.", tex_label="\\delta" )     


    # Overwrite a fixed parameter with an unfixed parameter
    m <= parameter(:ϵ_p, 0.750, (1e-5, 10.),   (1e-5, 10.),     DSGE.Exponential(),    GammaAlt(0.75, 0.4),        fixed=false,  scaling = x -> 1 + x/100,
                   description="ϵ_p: No description available.",
                   tex_label="\\varepsilon_{p}")     

    steadystate!(m)
end

m = Model990()
sstest(m)

@test m[:ι_w].value == 0.0
@test m[:ι_w].valuebounds == (0.0, .9999)
@test m[:ι_w].transform == DSGE.Untransformed()
@test m[:ι_w].transform_parameterization == (0.0,0.9999)
@test isa(m[:ι_w].prior.value, Normal)

@test m[:ι_p].value == 0.0
@test m[:ι_p].valuebounds == (0.0, 0.0)
@test isnull(m[:ι_p].prior)
@test m[:ι_p].fixed == true

@test m[:δ].value == 0.02
@test m[:δ].valuebounds == (0.02, 0.02)
@test isnull(m[:δ].prior)
@test m[:δ].fixed == true

@test m[:ϵ_p].value == 0.750
@test m[:ϵ_p].transform == DSGE.Exponential()
@test isa(m[:ϵ_p].prior.value, Gamma)
@test m[:ϵ_p].fixed==false

nothing
