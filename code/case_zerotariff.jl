# Take parameters from the draft, solve for the zero-tariff eqlm.
using JLD, Random, StatsBase
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

# Does results vary with initial guess?
simulation_counts = 1000;
results = zeros(simulation_counts,7);
for i in 1:simulation_counts
	Random.seed!(i)
	endogenous_value = (nlsolve(f!,randn(8)).zero) .^2
	results[i,1] = endogenous_value[2] ./endogenous_value[1] # relative wage to US
	results[i,2] = endogenous_value[3] # M^d_US
	results[i,3] = endogenous_value[4] # M^u_US
	results[i,4] = endogenous_value[5] # M^d_RoW
	results[i,5] = endogenous_value[6] # M^u_RoW
	results[i,6] = endogenous_value[7] # T_US
	results[i,7] = endogenous_value[8] # T_RoW
end

results_summary = zeros(7,4);
for col in 1:7
	results_summary[col,1]=percentile(results[:,col],5)
	results_summary[col,2]=percentile(results[:,col],50)
	results_summary[col,3]=percentile(results[:,col],75)
	results_summary[col,4]=percentile(results[:,col],90)
end

# When initial guess close to one, results are wrong. T is not small enough to be treated as zero.
nlsolve(f!,ones(8)).zero .^2
nlsolve(f!,fill(0.99,8)).zero .^2

# Save one simulation
Random.seed!(1)
sol = nlsolve(f!,randn(8));
endogenous_value = sol.zero .^2;
w_US = endogenous_value[1];
w_RoW = endogenous_value[2];
M_d_US = endogenous_value[3];
M_u_US = endogenous_value[4];
M_d_RoW = endogenous_value[5];
M_u_RoW = endogenous_value[6];
T_US = endogenous_value[7];
T_RoW = endogenous_value[8];

println("wage,US:", w_US)
println("wage,RoW:", w_RoW)
println("relative wage to US:", w_RoW/w_US)
println("M_d_US:", M_d_US)
println("M_u_US:", M_u_US)
println("M_d_RoW:", M_d_RoW)
println("M_u_RoW:", M_u_RoW)
println("T_US:",T_US)
println("T_RoW:",T_RoW)

save("../output/zerotariff_endogenous_value.jld",
		"w_US",w_US,"w_RoW",w_RoW,
		"M_d_US",M_d_US,"M_u_US",M_u_US,
		"M_d_RoW",M_d_RoW,"M_u_RoW",M_u_RoW,
		"T_US",T_US,"T_RoW",T_RoW)

