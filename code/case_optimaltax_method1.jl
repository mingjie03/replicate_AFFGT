# Case: optimal tariff
# Method: given tariff -> solve for the eqlm -> calculate the implied U -> find the max U

using Optim,JLD,Random,StatsBase

include("solve_eq.jl")
include("revert_solve_eq.jl")

df = load("../output/zerotariff_endogenous_value.jld")

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
## Assume no export subsidy ν=0
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

function cal_util(t_d::Float64,t_u::Float64,sol_guess::Array{Float64,1})

	global t_d_RoW_US = t_d
	global t_u_RoW_US = t_u

	# solve the eqlm
	sol = nlsolve(f!,sol_guess);
	endogenous_value = sol.zero .^2;
	w_US = endogenous_value[1];
	w_RoW = endogenous_value[2];
	M_d_US = endogenous_value[3];
	M_u_US = endogenous_value[4];
	M_d_RoW = endogenous_value[5];
	M_u_RoW = endogenous_value[6];
	T_US = endogenous_value[7];
	T_RoW = endogenous_value[8];

	# calculate price index
	p_d_US_RoW,p_d_RoW_US,p_d_US_US,p_d_RoW_RoW,p_u_US_RoW,p_u_RoW_US,p_u_US_US,p_u_RoW_RoW,x_US_RoW,x_RoW_US,x_US_US,x_RoW_RoW,c_US_RoW,c_RoW_US,c_US_US,c_RoW_RoW=revert_solve_eq(w_US,w_RoW,M_d_US,M_u_US,M_d_RoW,M_u_RoW,T_US,T_RoW,t_d_RoW_US,t_u_RoW_US)
	P_d_RoW_US = P_s_ji(M_d_RoW,t_d_RoW_US,p_d_RoW_US,σ)
	P_d_US_US = P_s_ji(M_d_US,0.0,p_d_US_US,σ)
	P_d_US = P_s_i(P_d_RoW_US,P_d_US_US,σ)

	# calculate util
	U_H = (w_US*L_US+T_US)/P_d_US
	return -U_H
end

function optim_util(x)
	return cal_util(x[1]^2,x[2]^2,guess_opttariff)
end

simulation_count = 1000;
results = zeros(simulation_count,3);
for i in 1:simulation_count
	Random.seed!(i)
	sol = optimize(optim_util,randn(2));
	tariff = Optim.minimizer(sol) .^2
	results[i,1] = tariff[1]
	results[i,2] = tariff[2]
	results[i,3] = -Optim.minimum(sol)
end

results_summary = zeros(3,10);
for col in 1:3
	results_summary[col,1] = minimum(results[:,col])
	results_summary[col,2] = percentile(results[:,col],5)
	results_summary[col,3] = percentile(results[:,col],10)
	results_summary[col,4] = percentile(results[:,col],25)
	results_summary[col,5] = percentile(results[:,col],50)
	results_summary[col,6] = mean(results[:,col])
	results_summary[col,7] = percentile(results[:,col],75)
	results_summary[col,8] = percentile(results[:,col],90)
	results_summary[col,9] = percentile(results[:,col],95)
	results_summary[col,10] = maximum(results[:,col])
end

