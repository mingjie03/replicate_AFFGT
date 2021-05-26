# Using eqlm solved by solve_eq.jl, revert the process and return relevant stats

include("functions.jl")

function revert_solve_eq(w_US::Float64,w_RoW::Float64,M_d_US::Float64,M_u_US::Float64,M_d_RoW::Float64,M_u_RoW::Float64,T_US::Float64,T_RoW::Float64,
		t_d_RoW_US::Float64,t_u_RoW_US::Float64,
		t_d_US_RoW::Float64=0.0,t_u_US_RoW::Float64=0.0,
		t_d_US_US::Float64=0.0,t_u_US_US::Float64=0.0,
		t_d_RoW_RoW::Float64=0.0,t_u_RoW_RoW::Float64=0.0,
		σ::Float64=4.0,θ::Float64=4.0,f_d::Float64=1.0,f_u::Float64=1.0,
		α_d::Float64=0.5483,α_u::Float64=1.0,
		L_US::Float64=0.4531,L_RoW::Float64=9.5469,
		A_d_RoW::Float64=0.2752,A_d_US::Float64=1.0,A_u_RoW::Float64=0.1121,A_u_US::Float64=1.0,
		τ_d::Float64=3.0066,τ_u::Float64=2.5971,
		ν::Float64=0.0)

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
	mc_u_US = (bar_α_u*w_US^α_u)/A_u_US;
	p_u_US_US = price_s_ij(μ_u,1.0,mc_u_US,ν);
	p_u_US_RoW  = price_s_ij(μ_u,τ_u,mc_u_US,ν);

	mc_u_RoW    = (bar_α_u*w_RoW^α_u)/A_u_RoW
	p_u_RoW_US  = price_s_ij(μ_u,τ_u,mc_u_RoW,ν)
	P_u_RoW_US  = P_s_ji(M_u_RoW,t_u_RoW_US,p_u_RoW_US,θ)
	P_u_US_US   = P_s_ji(M_u_US,t_u_US_US,p_u_US_US,θ)
	P_u_US      = P_s_i(P_u_RoW_US,P_u_US_US,θ)
	mc_d_US     = mc_d_i(α_d,w_US,P_u_US,A_d_US)
	ell_d_US    = ell_d_i(α_d,mc_d_US,f_d,y_d_US,w_US)

	p_u_RoW_RoW = price_s_ij(μ_u,1.0,mc_u_RoW,ν)
	P_u_RoW_RoW = P_s_ji(M_u_RoW,t_u_RoW_RoW,p_u_RoW_RoW,θ)
	P_u_US_RoW  = P_s_ji(M_u_US,t_u_US_RoW,p_u_US_RoW,θ)
	P_u_RoW     = P_s_i(P_u_US_RoW,P_u_RoW_RoW,θ)
	mc_d_RoW    = mc_d_i(α_d,w_RoW,P_u_RoW,A_d_RoW)
	ell_d_RoW   = ell_d_i(α_d,mc_d_RoW,f_d,y_d_RoW,w_RoW)

	p_d_US_US  = price_s_ij(μ_d,1.0,mc_d_US,ν)
	P_d_US_US  = P_s_ji(M_d_US,t_d_US_US,p_d_US_US,σ)
	p_d_RoW_US = price_s_ij(μ_d,τ_d,mc_d_RoW,ν)
	P_d_RoW_US = P_s_ji(M_d_RoW,t_d_RoW_US,p_d_RoW_US,σ)
	P_d_US     = P_s_i(P_d_RoW_US,P_d_US_US,σ)
	c_US_US    = c_ji(w_US,L_US,T_US,t_d_US_US,p_d_US_US,σ,P_d_US)

	c_RoW_US  = c_ji(w_US,L_US,T_US,t_d_RoW_US,p_d_RoW_US,σ,P_d_US)

	p_d_US_RoW  = price_s_ij(μ_d,τ_d,mc_d_US,ν)
	p_d_RoW_RoW = price_s_ij(μ_d,1.0,mc_d_RoW,ν)
	P_d_RoW_RoW = P_s_ji(M_d_RoW,t_d_RoW_RoW,p_d_RoW_RoW,σ)
	P_d_US_RoW  = P_s_ji(M_d_US,t_d_US_RoW,p_d_US_RoW,σ)
	P_d_RoW     = P_s_i(P_d_US_RoW,P_d_RoW_RoW,σ)
	c_RoW_RoW   = c_ji(w_RoW,L_RoW,T_RoW,t_d_RoW_RoW,p_d_RoW_RoW,σ,P_d_RoW)

	c_US_RoW  = c_ji(w_RoW,L_RoW,T_RoW,t_d_US_RoW,p_d_US_RoW,σ,P_d_RoW)

	Q_u_US_US = Q_u_ji(α_d,mc_d_US,f_d,y_d_US,P_u_US,P_u_US_US,θ)
	x_US_US   = x_ji(Q_u_US_US,t_u_US_US,p_u_US_US,P_u_US_US,θ)

	Q_u_US_RoW = Q_u_ji(α_d,mc_d_RoW,f_d,y_d_RoW,P_u_RoW,P_u_US_RoW,θ)
	x_US_RoW   = x_ji(Q_u_US_RoW,t_u_US_RoW,p_u_US_RoW,P_u_US_RoW,θ)

	Q_u_RoW_RoW = Q_u_ji(α_d,mc_d_RoW,f_d,y_d_RoW,P_u_RoW,P_u_RoW_RoW,θ)
	x_RoW_RoW = x_ji(Q_u_RoW_RoW,t_u_RoW_RoW,p_u_RoW_RoW,P_u_RoW_RoW,θ)

	Q_u_RoW_US  = Q_u_ji(α_d,mc_d_US,f_d,y_d_US,P_u_US,P_u_RoW_US,θ)
	x_RoW_US  = x_ji(Q_u_RoW_US,t_u_RoW_US,p_u_RoW_US,P_u_RoW_US,θ)

	return p_d_US_RoW,p_d_RoW_US,p_d_US_US,p_d_RoW_RoW,p_u_US_RoW,p_u_RoW_US,p_u_US_US,p_u_RoW_RoW,x_US_RoW,x_RoW_US,x_US_US,x_RoW_RoW,c_US_RoW,c_RoW_US,c_US_US,c_RoW_RoW
end
