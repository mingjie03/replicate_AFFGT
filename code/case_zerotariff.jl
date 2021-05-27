# Take parameters from the draft, solve for the zero-tariff eqlm.
using JLD, Random, StatsBase, DataFrames
include("solve_eq.jl")

# Parameters taken from the draft
## Fixed values
θ = 4.0; # input varieties
σ = 4.0; # final-good varieties
f_d = 1.0;
f_u = 1.0;

## Values measured from data
α_d = 1.0-0.4517;
α_u = 1.0;
L_US = 0.4531;
L_RoW = 9.5469;

## Estimated values
A_d_RoW = 0.2752;
A_d_US = 1.0;
A_u_RoW = 0.1121;
A_u_US = 1.0;

τ_d = 3.0066;
τ_u = 2.5971;

# Other assumptions
# Consider the zero-tariff eqlm
t_d_RoW_US = 0.0;
t_u_RoW_US = 0.0;

# Assume no export subsidy ν=0
ν = 0.0;

# Values could be directly computed
## Free entry implies that
y_d_US  = free_entry(σ,f_d);
y_d_RoW = free_entry(σ,f_d);
y_u_US  = free_entry(θ,f_u);
y_u_RoW = free_entry(θ,f_u);

## μ = elasticity/(elasticity-1)
μ_d = σ/(σ-1);
μ_u = θ/(θ-1);

ell_u_US  = ell_u_i(f_u,y_u_US,A_u_US);
ell_u_RoW = ell_u_i(f_u,y_u_RoW,A_u_RoW);

# keep solutions where eqlm really hold (error tol = 1e-6)
function verify_solution(rand_seed::Int64)
	# solve the eqlm
	Random.seed!(rand_seed)
	endogenous_value = (nlsolve(f!,randn(8)).zero) .^2
	w_US = endogenous_value[1];
	w_RoW = endogenous_value[2];
	M_d_US = endogenous_value[3];
	M_u_US = endogenous_value[4];
	M_d_RoW = endogenous_value[5];
	M_u_RoW = endogenous_value[6];
	T_US = endogenous_value[7];
	T_RoW = endogenous_value[8];

	# back out relevant stats
	bar_α_u = 1/(α_u^α_u*(1-α_u)^(1-α_u));
	mc_u_US = (bar_α_u*w_US^α_u)/A_u_US;
	p_u_US_US = price_s_ij(μ_u,1.0,mc_u_US,ν); # assume τ_ii=1.0
	p_u_US_RoW  = price_s_ij(μ_u,τ_u,mc_u_US,ν);

	mc_u_RoW    = (bar_α_u*w_RoW^α_u)/A_u_RoW
	p_u_RoW_US  = price_s_ij(μ_u,τ_u,mc_u_RoW,ν)
	P_u_RoW_US  = P_s_ji(M_u_RoW,t_u_RoW_US,p_u_RoW_US,θ)
	P_u_US_US   = P_s_ji(M_u_US,0.0,p_u_US_US,θ)
	P_u_US      = P_s_i(P_u_RoW_US,P_u_US_US,θ)
	mc_d_US     = mc_d_i(α_d,w_US,P_u_US,A_d_US)
	ell_d_US    = ell_d_i(α_d,mc_d_US,f_d,y_d_US,w_US)

	p_u_RoW_RoW = price_s_ij(μ_u,1.0,mc_u_RoW,ν)
	P_u_RoW_RoW = P_s_ji(M_u_RoW,0.0,p_u_RoW_RoW,θ)
	P_u_US_RoW  = P_s_ji(M_u_US,0.0,p_u_US_RoW,θ)
	P_u_RoW     = P_s_i(P_u_US_RoW,P_u_RoW_RoW,θ)
	mc_d_RoW    = mc_d_i(α_d,w_RoW,P_u_RoW,A_d_RoW)
	ell_d_RoW   = ell_d_i(α_d,mc_d_RoW,f_d,y_d_RoW,w_RoW)

	p_d_US_US  = price_s_ij(μ_d,1.0,mc_d_US,ν)
	P_d_US_US  = P_s_ji(M_d_US,0.0,p_d_US_US,σ)
	p_d_RoW_US = price_s_ij(μ_d,τ_d,mc_d_RoW,ν)
	P_d_RoW_US = P_s_ji(M_d_RoW,t_d_RoW_US,p_d_RoW_US,σ)
	P_d_US     = P_s_i(P_d_RoW_US,P_d_US_US,σ)
	c_US_US    = c_ji(w_US,L_US,T_US,0.0,p_d_US_US,σ,P_d_US)

	c_RoW_US  = c_ji(w_US,L_US,T_US,t_d_RoW_US,p_d_RoW_US,σ,P_d_US)

	p_d_US_RoW  = price_s_ij(μ_d,τ_d,mc_d_US,ν)
	p_d_RoW_RoW = price_s_ij(μ_d,1.0,mc_d_RoW,ν)
	P_d_RoW_RoW = P_s_ji(M_d_RoW,0.0,p_d_RoW_RoW,σ)
	P_d_US_RoW  = P_s_ji(M_d_US,0.0,p_d_US_RoW,σ)
	P_d_RoW     = P_s_i(P_d_US_RoW,P_d_RoW_RoW,σ)
	c_RoW_RoW   = c_ji(w_RoW,L_RoW,T_RoW,0.0,p_d_RoW_RoW,σ,P_d_RoW)

	c_US_RoW  = c_ji(w_RoW,L_RoW,T_RoW,0.0,p_d_US_RoW,σ,P_d_RoW)

	Q_u_US_US = Q_u_ji(α_d,mc_d_US,f_d,y_d_US,P_u_US,P_u_US_US,θ)
	x_US_US   = x_ji(Q_u_US_US,0.0,p_u_US_US,P_u_US_US,θ)

	Q_u_US_RoW = Q_u_ji(α_d,mc_d_RoW,f_d,y_d_RoW,P_u_RoW,P_u_US_RoW,θ)
	x_US_RoW   = x_ji(Q_u_US_RoW,0.0,p_u_US_RoW,P_u_US_RoW,θ)

	Q_u_RoW_RoW = Q_u_ji(α_d,mc_d_RoW,f_d,y_d_RoW,P_u_RoW,P_u_RoW_RoW,θ)
	x_RoW_RoW = x_ji(Q_u_RoW_RoW,0.0,p_u_RoW_RoW,P_u_RoW_RoW,θ)

	Q_u_RoW_US  = Q_u_ji(α_d,mc_d_US,f_d,y_d_US,P_u_US,P_u_RoW_US,θ)
	x_RoW_US  = x_ji(Q_u_RoW_US,t_u_RoW_US,p_u_RoW_US,P_u_RoW_US,θ)

	# verify the solution
	# labor market clear
	eq1 = L_US - (M_d_US*ell_d_US + M_u_US*ell_u_US)
	eq2 = L_RoW - (M_d_RoW*ell_d_RoW + M_u_RoW*ell_u_RoW)
	# Good market clearing
	eq3 = y_d_US - (c_US_US + τ_d*c_US_RoW)
	eq4 = y_d_RoW - (c_RoW_RoW + τ_d*c_RoW_US)
	eq5 = y_u_US - (M_d_US*x_US_US + M_d_RoW*τ_u*x_US_RoW)
	eq6 = y_u_RoW - (M_d_RoW*x_RoW_RoW + M_d_US*τ_u*x_RoW_US)
	# Tax revenue
	eq7 = T_US - (- ν*M_d_US*c_US_US  *p_d_US_US   - ν*M_d_US*M_u_US*x_US_US  *p_u_US_US)   - (t_d_RoW_US*M_d_RoW*c_RoW_US*p_d_RoW_US + t_u_RoW_US*M_d_US*M_u_RoW*x_RoW_US*p_u_RoW_US - ν*M_d_US*c_US_RoW*p_d_US_RoW - ν*M_d_RoW*M_u_US*x_US_RoW*p_u_US_RoW)
	eq8 = T_RoW - (- ν*M_d_RoW*c_RoW_RoW*p_d_RoW_RoW - ν*M_d_RoW*M_u_RoW*x_RoW_RoW*p_u_RoW_RoW) - (- ν*M_d_RoW*c_RoW_US*p_d_RoW_US - ν*M_d_US*M_u_RoW*x_RoW_US*p_u_RoW_US)

	res = vcat(eq1,eq2,eq3,eq4,eq5,eq6,eq7,eq8)
	if (res .< fill(1e-6,8)) == ones(8)
		final_res = vcat(w_US,w_RoW,M_d_US,M_u_US,M_d_RoW,M_u_RoW,T_US,T_RoW)
	else
		final_res = fill(NaN,8)
	end

	return final_res
end

# Collect results
simulation_counts = 1000;
results = zeros(8,simulation_counts);
for i in 1:simulation_counts
	results[:,i] = verify_solution(i)
end
df = convert(DataFrame,results')
rename!(df,:x1=>:w_US,:x2=>:w_RoW,:x3=>:M_d_US,:x4=>:M_u_US,:x5=>:M_d_RoW,:x6=>:M_u_RoW,:x7=>:T_US,:x8=>:T_RoW)
df = filter(row -> !isnan(row.w_US),df)
df[!,:rel_wage] = df[!,:w_RoW] ./df[!,:w_US]

# summary
function stats(vec::Array{Float64,1})
	res = zeros(6)
	res[1] = minimum(vec)
	res[2] = percentile(vec,5)
	res[3] = percentile(vec,50)
	res[4] = mean(vec)
	res[5] = percentile(vec,95)
	res[6] = maximum(vec)
	return res
end

stats(df[!,:rel_wage])
stats(df[!,:M_d_US])
stats(df[!,:M_u_US])
stats(df[!,:M_d_RoW])
stats(df[!,:M_u_RoW])
stats(df[!,:T_US])
stats(df[!,:rel_wage])

# Save one simulation
sol = verify_solution(1)
save("../output/zerotariff_endogenous_value.jld",
		"w_US",sol[1],"w_RoW",sol[2],
		"M_d_US",sol[3],"M_u_US",sol[4],
		"M_d_RoW",sol[5],"M_u_RoW",sol[6],
		"T_US",sol[7],"T_RoW",sol[8])

