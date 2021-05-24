using NLsolve
#using KNITRO

include("functions.jl")

function f!(F,x)
	# x = (sqrt(w_US),sqrt(w_RoW),sqrt(M_d_US),sqrt(M_u_US),sqrt(M_d_RoW),sqrt(M_u_RoW),sqrt(T_US),sqrt(T_RoW))

	bar_α_u = 1/(α_u^α_u*(1-α_u)^(1-α_u));
	mc_u_US = (bar_α_u*x[1]^(2*α_u))/A_u_US;
	p_u_US_US = price_s_ij(μ_u,1.0,mc_u_US,ν); # assume τ_ii=1.0
	p_u_US_RoW  = price_s_ij(μ_u,τ_u,mc_u_US,ν);

	mc_u_RoW    = (bar_α_u*x[2]^(2*α_u))/A_u_RoW
	p_u_RoW_US  = price_s_ij(μ_u,τ_u,mc_u_RoW,ν)
	P_u_RoW_US  = P_s_ji(x[6]^2,t_u_RoW_US,p_u_RoW_US,θ)
	P_u_US_US   = P_s_ji(x[4]^2,t_u_US_US,p_u_US_US,θ)
	P_u_US      = P_s_i(P_u_RoW_US,P_u_US_US,θ)
	mc_d_US     = mc_d_i(α_d,x[1]^2,P_u_US,A_d_US)
	ell_d_US    = ell_d_i(α_d,mc_d_US,f_d,y_d_US,x[1]^2)

	p_u_RoW_RoW = price_s_ij(μ_u,1.0,mc_u_RoW,ν)
	P_u_RoW_RoW = P_s_ji(x[6]^2,t_u_RoW_RoW,p_u_RoW_RoW,θ)
	P_u_US_RoW  = P_s_ji(x[4]^2,t_u_US_RoW,p_u_US_RoW,θ)
	P_u_RoW     = P_s_i(P_u_US_RoW,P_u_RoW_RoW,θ)
	mc_d_RoW    = mc_d_i(α_d,x[2]^2,P_u_RoW,A_d_RoW)
	ell_d_RoW   = ell_d_i(α_d,mc_d_RoW,f_d,y_d_RoW,x[2]^2)

	p_d_US_US  = price_s_ij(μ_d,1.0,mc_d_US,ν)
	P_d_US_US  = P_s_ji(x[3]^2,t_d_US_US,p_d_US_US,σ)
	p_d_RoW_US = price_s_ij(μ_d,τ_d,mc_d_RoW,ν)
	P_d_RoW_US = P_s_ji(x[5]^2,t_d_RoW_US,p_d_RoW_US,σ)
	P_d_US     = P_s_i(P_d_RoW_US,P_d_US_US,σ)
	c_US_US    = c_ji(x[1]^2,L_US,x[7]^2,t_d_US_US,p_d_US_US,σ,P_d_US)

	c_RoW_US  = c_ji(x[1]^2,L_US,x[7]^2,t_d_RoW_US,p_d_RoW_US,σ,P_d_US)

	p_d_US_RoW  = price_s_ij(μ_d,τ_d,mc_d_US,ν)
	p_d_RoW_RoW = price_s_ij(μ_d,1.0,mc_d_RoW,ν)
	P_d_RoW_RoW = P_s_ji(x[5]^2,t_d_RoW_RoW,p_d_RoW_RoW,σ)
	P_d_US_RoW  = P_s_ji(x[3]^2,t_d_US_RoW,p_d_US_RoW,σ)
	P_d_RoW     = P_s_i(P_d_US_RoW,P_d_RoW_RoW,σ)
	c_RoW_RoW   = c_ji(x[2]^2,L_RoW,x[8]^2,t_d_RoW_RoW,p_d_RoW_RoW,σ,P_d_RoW)

	c_US_RoW  = c_ji(x[2]^2,L_RoW,x[8]^2,t_d_US_RoW,p_d_US_RoW,σ,P_d_RoW)

	Q_u_US_US = Q_u_ji(α_d,mc_d_US,f_d,y_d_US,P_u_US,P_u_US_US,θ)
	x_US_US   = x_ji(Q_u_US_US,t_u_US_US,p_u_US_US,P_u_US_US,θ)

	Q_u_US_RoW = Q_u_ji(α_d,mc_d_RoW,f_d,y_d_RoW,P_u_RoW,P_u_US_RoW,θ)
	x_US_RoW   = x_ji(Q_u_US_RoW,t_u_US_RoW,p_u_US_RoW,P_u_US_RoW,θ)

	Q_u_RoW_RoW = Q_u_ji(α_d,mc_d_RoW,f_d,y_d_RoW,P_u_RoW,P_u_RoW_RoW,θ)
	x_RoW_RoW = x_ji(Q_u_RoW_RoW,t_u_RoW_RoW,p_u_RoW_RoW,P_u_RoW_RoW,θ)

	Q_u_RoW_US  = Q_u_ji(α_d,mc_d_US,f_d,y_d_US,P_u_US,P_u_RoW_US,θ)
	x_RoW_US  = x_ji(Q_u_RoW_US,t_u_RoW_US,p_u_RoW_US,P_u_RoW_US,θ)

	# labor market clear
	F[1] = L_US - (x[3]^2*ell_d_US + x[4]^2*ell_u_US)
	F[2] = L_RoW - (x[5]^2*ell_d_RoW + x[6]^2*ell_u_RoW)

	# Good market clearing
	F[3] = y_d_US - (c_US_US + τ_d*c_US_RoW)
	F[4] = y_d_RoW - (c_RoW_RoW + τ_d*c_RoW_US)
	F[5] = y_u_US - (x[3]^2*x_US_US + x[5]^2*τ_u*x_US_RoW)
	F[6] = y_u_RoW - (x[5]^2*x_RoW_RoW + x[3]^2*τ_u*x_RoW_US)

	# Tax revenue
	F[7] = x[7]^2 - (t_d_US_US  *x[3]^2*c_US_US  *p_d_US_US   + t_u_US_US  *x[3]^2*x[4]^2*x_US_US  *p_u_US_US   - ν*x[3]^2*c_US_US  *p_d_US_US   - ν*x[3]^2*x[4]^2*x_US_US  *p_u_US_US)   - (t_d_RoW_US*x[5]^2*c_RoW_US*p_d_RoW_US + t_u_RoW_US*x[3]^2*x[6]^2*x_RoW_US*p_u_RoW_US - ν*x[3]^2*c_US_RoW*p_d_US_RoW - ν*x[5]^2*x[4]^2*x_US_RoW*p_u_US_RoW)
	F[8] = x[8]^2 - (t_d_RoW_RoW*x[5]^2*c_RoW_RoW*p_d_RoW_RoW + t_u_RoW_RoW*x[5]^2*x[6]^2*x_RoW_RoW*p_u_RoW_RoW - ν*x[5]^2*c_RoW_RoW*p_d_RoW_RoW - ν*x[5]^2*x[6]^2*x_RoW_RoW*p_u_RoW_RoW) - (t_d_US_RoW*x[3]^2*c_US_RoW*p_d_US_RoW + t_u_US_RoW*x[5]^2*x[4]^2*x_US_RoW*p_u_US_RoW - ν*x[5]^2*c_RoW_US*p_d_RoW_US - ν*x[3]^2*x[6]^2*x_RoW_US*p_u_RoW_US)
end
