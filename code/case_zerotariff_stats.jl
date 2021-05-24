# Replicate Table 3

using JLD
include("functions.jl")

# Load zero-tariff eqlm values
df = load("../output/zerotariff_endogenous_value.jld")

# Copy from case_zerotariff.jl to calculate relevant stats
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
t_d_US_RoW = 0.0;
t_d_RoW_US = 0.0;
t_d_US_US = 0.0;
t_d_RoW_RoW = 0.0;
t_u_US_RoW = 0.0;
t_u_RoW_US = 0.0;
t_u_US_US = 0.0;
t_u_RoW_RoW = 0.0;

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

bar_α_u = 1/(α_u^α_u*(1-α_u)^(1-α_u));
mc_u_US = (bar_α_u*sqrt(df["w_US"])^(2*α_u))/A_u_US;
p_u_US_US = price_s_ij(μ_u,1.0,mc_u_US,ν); # assume τ_ii=1.0
p_u_US_RoW  = price_s_ij(μ_u,τ_u,mc_u_US,ν);

mc_u_RoW    = (bar_α_u*sqrt(df["w_RoW"])^(2*α_u))/A_u_RoW;
p_u_RoW_US  = price_s_ij(μ_u,τ_u,mc_u_RoW,ν);
P_u_RoW_US  = P_s_ji(sqrt(df["M_u_RoW"])^2,t_u_RoW_US,p_u_RoW_US,θ);
P_u_US_US   = P_s_ji(sqrt(df["M_u_US"])^2,t_u_US_US,p_u_US_US,θ);
P_u_US      = P_s_i(P_u_RoW_US,P_u_US_US,θ);
mc_d_US     = mc_d_i(α_d,sqrt(df["w_US"])^2,P_u_US,A_d_US);
ell_d_US    = ell_d_i(α_d,mc_d_US,f_d,y_d_US,sqrt(df["w_US"])^2);

p_u_RoW_RoW = price_s_ij(μ_u,1.0,mc_u_RoW,ν);
P_u_RoW_RoW = P_s_ji(sqrt(df["M_u_RoW"])^2,t_u_RoW_RoW,p_u_RoW_RoW,θ);
P_u_US_RoW  = P_s_ji(sqrt(df["M_u_US"])^2,t_u_US_RoW,p_u_US_RoW,θ);
P_u_RoW     = P_s_i(P_u_US_RoW,P_u_RoW_RoW,θ);
mc_d_RoW    = mc_d_i(α_d,sqrt(df["w_RoW"])^2,P_u_RoW,A_d_RoW);
ell_d_RoW   = ell_d_i(α_d,mc_d_RoW,f_d,y_d_RoW,sqrt(df["w_RoW"])^2);

p_d_US_US  = price_s_ij(μ_d,1.0,mc_d_US,ν);
P_d_US_US  = P_s_ji(sqrt(df["M_d_US"])^2,t_d_US_US,p_d_US_US,σ);
p_d_RoW_US = price_s_ij(μ_d,τ_d,mc_d_RoW,ν);
P_d_RoW_US = P_s_ji(sqrt(df["M_d_RoW"])^2,t_d_RoW_US,p_d_RoW_US,σ);
P_d_US     = P_s_i(P_d_RoW_US,P_d_US_US,σ);
c_US_US    = c_ji(sqrt(df["w_US"])^2,L_US,sqrt(df["T_US"])^2,t_d_US_US,p_d_US_US,σ,P_d_US);

c_RoW_US  = c_ji(sqrt(df["w_US"])^2,L_US,sqrt(df["T_US"])^2,t_d_RoW_US,p_d_RoW_US,σ,P_d_US);

p_d_US_RoW  = price_s_ij(μ_d,τ_d,mc_d_US,ν);
p_d_RoW_RoW = price_s_ij(μ_d,1.0,mc_d_RoW,ν);
P_d_RoW_RoW = P_s_ji(sqrt(df["M_d_RoW"])^2,t_d_RoW_RoW,p_d_RoW_RoW,σ);
P_d_US_RoW  = P_s_ji(sqrt(df["M_d_US"])^2,t_d_US_RoW,p_d_US_RoW,σ);
P_d_RoW     = P_s_i(P_d_US_RoW,P_d_RoW_RoW,σ);
c_RoW_RoW   = c_ji(sqrt(df["w_RoW"])^2,L_RoW,sqrt(df["T_RoW"])^2,t_d_RoW_RoW,p_d_RoW_RoW,σ,P_d_RoW);

c_US_RoW  = c_ji(sqrt(df["w_RoW"])^2,L_RoW,sqrt(df["T_RoW"])^2,t_d_US_RoW,p_d_US_RoW,σ,P_d_RoW);

Q_u_US_US = Q_u_ji(α_d,mc_d_US,f_d,y_d_US,P_u_US,P_u_US_US,θ);
x_US_US   = x_ji(Q_u_US_US,t_u_US_US,p_u_US_US,P_u_US_US,θ);

Q_u_US_RoW = Q_u_ji(α_d,mc_d_RoW,f_d,y_d_RoW,P_u_RoW,P_u_US_RoW,θ);
x_US_RoW   = x_ji(Q_u_US_RoW,t_u_US_RoW,p_u_US_RoW,P_u_US_RoW,θ);

Q_u_RoW_RoW = Q_u_ji(α_d,mc_d_RoW,f_d,y_d_RoW,P_u_RoW,P_u_RoW_RoW,θ);
x_RoW_RoW = x_ji(Q_u_RoW_RoW,t_u_RoW_RoW,p_u_RoW_RoW,P_u_RoW_RoW,θ);

Q_u_RoW_US  = Q_u_ji(α_d,mc_d_US,f_d,y_d_US,P_u_US,P_u_RoW_US,θ);
x_RoW_US  = x_ji(Q_u_RoW_US,t_u_RoW_US,p_u_RoW_US,P_u_RoW_US,θ);

# Replicate Table 3
Ω_HH = (df["M_u_US"]*df["M_d_US"]*p_u_US_US*x_US_US)/(df["M_d_US"]*(p_d_US_RoW*c_US_RoW+p_d_US_US*c_US_US))
Ω_FH = (df["M_u_RoW"]*df["M_d_US"]*p_u_RoW_US*x_RoW_US)/(df["M_d_US"]*(p_d_US_RoW*c_US_RoW+p_d_US_US*c_US_US))
Ω_FF = (df["M_u_RoW"]*df["M_d_RoW"]*p_u_RoW_RoW*x_RoW_RoW)/(df["M_d_RoW"]*(p_d_RoW_US*c_RoW_US+p_d_RoW_RoW*c_RoW_RoW))
Ω_HF = (df["M_u_US"]*df["M_d_RoW"]*p_u_US_RoW*x_US_RoW)/(df["M_d_RoW"]*(p_d_RoW_US*c_RoW_US+p_d_RoW_RoW*c_RoW_RoW))
b_H_H = (df["M_d_US"]*p_d_US_US*c_US_US)/(df["w_US"]*L_US)
b_H_F = (df["M_d_RoW"]*p_d_RoW_US*c_RoW_US)/(df["w_US"]*L_US)
λ_d_H = (df["M_d_US"]*(p_d_US_RoW*c_US_RoW+p_d_US_US*c_US_US))/(df["w_US"]*L_US)




