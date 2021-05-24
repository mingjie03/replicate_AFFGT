# Replicate Table 3

using JLD
include("functions.jl")
include("revert_solve_eq.jl")

# Load zero-tariff eqlm values and values needed to compute stats
df = load("../output/zerotariff_endogenous_value.jld")
L_US = 0.4531;

# revert relevant stats
p_d_US_RoW,p_d_RoW_US,p_d_US_US,p_d_RoW_RoW,p_u_US_RoW,p_u_RoW_US,p_u_US_US,p_u_RoW_RoW,x_US_RoW,x_RoW_US,x_US_US,x_RoW_RoW,c_US_RoW,c_RoW_US,c_US_US,c_RoW_RoW=revert_solve_eq(df["w_US"],df["w_RoW"],df["M_d_US"],df["M_u_US"],df["M_d_RoW"],df["M_u_RoW"],df["T_US"],df["T_RoW"],0.0,0.0)

# Replicate Table 3
Ω_HH = (df["M_u_US"]*df["M_d_US"]*p_u_US_US*x_US_US)/(df["M_d_US"]*(p_d_US_RoW*c_US_RoW+p_d_US_US*c_US_US))
Ω_FH = (df["M_u_RoW"]*df["M_d_US"]*p_u_RoW_US*x_RoW_US)/(df["M_d_US"]*(p_d_US_RoW*c_US_RoW+p_d_US_US*c_US_US))
Ω_FF = (df["M_u_RoW"]*df["M_d_RoW"]*p_u_RoW_RoW*x_RoW_RoW)/(df["M_d_RoW"]*(p_d_RoW_US*c_RoW_US+p_d_RoW_RoW*c_RoW_RoW))
Ω_HF = (df["M_u_US"]*df["M_d_RoW"]*p_u_US_RoW*x_US_RoW)/(df["M_d_RoW"]*(p_d_RoW_US*c_RoW_US+p_d_RoW_RoW*c_RoW_RoW))
b_H_H = (df["M_d_US"]*p_d_US_US*c_US_US)/(df["w_US"]*L_US)
b_H_F = (df["M_d_RoW"]*p_d_RoW_US*c_RoW_US)/(df["w_US"]*L_US)
λ_d_H = (df["M_d_US"]*(p_d_US_RoW*c_US_RoW+p_d_US_US*c_US_US))/(df["w_US"]*L_US)
