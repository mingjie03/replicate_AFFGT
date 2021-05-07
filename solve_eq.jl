# Take parameters from the draft, solve for the zero-tariff eqlm to replicate Table 3.
# MZ

using NLsolve
#using KNITRO

include("functions.jl")

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
t = 0.0;
T = 0.0;

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


function f!(F,x)
	# x = (sqrt(w_US),sqrt(w_RoW),sqrt(M_d_US),sqrt(M_u_US),sqrt(M_d_RoW),sqrt(M_u_RoW))

	bar_α_u = 1/(α_u^α_u*(1-α_u)^(1-α_u));
	mc_u_US = (bar_α_u*x[1]^(2*α_u))/A_u_US;
	p_u_US_US = price_s_ij(μ_u,1.0,mc_u_US,ν); # assume τ_ii=1.0
	p_u_US_RoW  = price_s_ij(μ_u,τ_u,mc_u_US,ν);

	mc_u_RoW    = (bar_α_u*x[2]^(2*α_u))/A_u_RoW
	p_u_RoW_US  = price_s_ij(μ_u,τ_u,mc_u_RoW,ν)
	P_u_RoW_US  = P_s_ji(x[6]^2,t,p_u_RoW_US,θ)
	P_u_US_US   = P_s_ji(x[4]^2,t,p_u_US_US,θ)
	P_u_US      = P_s_i(P_u_RoW_US,P_u_US_US,θ)
	mc_d_US     = mc_d_i(α_d,x[1]^2,P_u_US,A_d_US)
	ell_d_US    = ell_d_i(α_d,mc_d_US,f_d,y_d_US,x[1]^2)

	p_u_RoW_RoW = price_s_ij(μ_u,1.0,mc_u_RoW,ν)
	P_u_RoW_RoW = P_s_ji(x[6]^2,t,p_u_RoW_RoW,θ)
	P_u_US_RoW  = P_s_ji(x[4]^2,t,p_u_US_RoW,θ)
	P_u_RoW     = P_s_i(P_u_US_RoW,P_u_RoW_RoW,θ)
	mc_d_RoW    = mc_d_i(α_d,x[2]^2,P_u_RoW,A_d_RoW)
	ell_d_RoW   = ell_d_i(α_d,mc_d_RoW,f_d,y_d_RoW,x[2]^2)

	# labor market clear
	F[1] = L_US - (x[3]^2*ell_d_US + x[4]^2*ell_u_US)
	F[2] = L_RoW - (x[5]^2*ell_d_RoW + x[6]^2*ell_u_RoW)

	p_d_US_US  = price_s_ij(μ_d,1.0,mc_d_US,ν) # Assume intra-trade iceberg cost = 1
	P_d_US_US  = P_s_ji(x[3]^2,t,p_d_US_US,σ)
	p_d_RoW_US = price_s_ij(μ_d,τ_d,mc_d_RoW,ν)
	P_d_RoW_US = P_s_ji(x[5]^2,t,p_d_RoW_US,σ)
	P_d_US     = P_s_i(P_d_RoW_US,P_d_US_US,σ)
	c_US_US    = c_ji(x[1]^2,L_US,T,t,p_d_US_US,σ,P_d_US)

	c_RoW_US  = c_ji(x[1]^2,L_US,T,t,p_d_RoW_US,σ,P_d_US)

	p_d_US_RoW  = price_s_ij(μ_d,τ_d,mc_d_US,ν)
	p_d_RoW_RoW = price_s_ij(μ_d,1.0,mc_d_RoW,ν) # Assume τ_RoW_RoW=1, is it true?
	P_d_RoW_RoW = P_s_ji(x[5]^2,t,p_d_RoW_RoW,σ)
	P_d_US_RoW  = P_s_ji(x[3]^2,t,p_d_US_RoW,σ)
	P_d_RoW     = P_s_i(P_d_US_RoW,P_d_RoW_RoW,σ)
	c_RoW_RoW   = c_ji(x[2]^2,L_RoW,T,t,p_d_RoW_RoW,σ,P_d_RoW)

	c_US_RoW  = c_ji(x[2]^2,L_RoW,T,t,p_d_US_RoW,σ,P_d_RoW)

	Q_u_US_US = Q_u_ji(α_d,mc_d_US,f_d,y_d_US,P_u_US,P_u_US_US,θ)
	x_US_US   = x_ji(Q_u_US_US,t,p_u_US_US,P_u_US_US,θ)

	Q_u_US_RoW = Q_u_ji(α_d,mc_d_RoW,f_d,y_d_RoW,P_u_RoW,P_u_US_RoW,θ)
	x_US_RoW   = x_ji(Q_u_US_RoW,t,p_u_US_RoW,P_u_US_RoW,θ)

	Q_u_RoW_RoW = Q_u_ji(α_d,mc_d_RoW,f_d,y_d_RoW,P_u_RoW,P_u_RoW_RoW,θ)
	x_RoW_RoW = x_ji(Q_u_RoW_RoW,t,p_u_RoW_RoW,P_u_RoW_RoW,θ)

	Q_u_RoW_US  = Q_u_ji(α_d,mc_d_US,f_d,y_d_US,P_u_US,P_u_RoW_US,θ)
	x_RoW_US  = x_ji(Q_u_RoW_US,t,p_u_RoW_US,P_u_RoW_US,θ)

	# Good market clearing
	F[3] = y_d_US - (c_US_US + τ_d*c_US_RoW)
	F[4] = y_d_RoW - (c_RoW_RoW + τ_d*c_RoW_US)
	F[5] = y_u_US - (x[3]^2*x_US_US + x[5]^2*τ_u*x_US_RoW)
	F[6] = y_u_RoW - (x[5]^2*x_RoW_RoW + x[3]^2*τ_u*x_RoW_US)
end

sol = nlsolve(f!,ones(6));
endogenous_value = sol.zero .^2;
w_US = endogenous_value[1];
w_RoW = endogenous_value[2];
M_d_US = endogenous_value[3];
M_u_US = endogenous_value[4];
M_d_RoW = endogenous_value[5];
M_u_RoW = endogenous_value[6];

println("wage,US:", w_US)
println("wage,RoW:", w_RoW)
println("relative wage to US:", w_RoW/w_US)
println("M_d_US:", M_d_US)
println("M_u_US:", M_u_US)
println("M_d_RoW:", M_d_RoW)
println("M_u_RoW:", M_u_RoW)
