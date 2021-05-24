# Functions
## eq26.1
function price_s_ij(μ_s::Float64,τ_s_ij::Float64,mc_s_i::Float64,ν_s_ij::Float64)
	price_s_ij = (μ_s*τ_s_ij*mc_s_i)/(1+ν_s_ij)
	return price_s_ij
end
## eq26.2
function mc_d_i(α_d::Float64,w_i::Float64,P_u_i::Float64,A_d_i::Float64)
	bar_α_d = 1/(α_d^α_d*(1-α_d)^(1-α_d))
	mc_d_i = (bar_α_d*w_i^α_d*P_u_i^(1-α_d))/A_d_i
	return mc_d_i
end
## eq27.1
function c_ji(w_i::Float64,L_i::Float64,T_i::Float64,t_d_ji::Float64,p_d_ji::Float64,σ::Float64,P_d_i::Float64)
	c_ji = (w_i*L_i+T_i) * P_d_i^(σ-1) * (1/((1+t_d_ji)*p_d_ji)^σ)
	return c_ji
end
## eq27.2
function x_ji(Q_u_ji::Float64,t_u_ji::Float64,p_u_ji::Float64,P_u_ji::Float64,θ::Float64)
	x_ji = Q_u_ji * P_u_ji^θ * (1/((1+t_u_ji)*p_u_ji)^θ)
	return x_ji
end
## eq27.3
function Q_u_ji(α::Float64,mc_d_i::Float64,f_d_i::Float64,y_d_i::Float64,P_u_i::Float64,P_u_ji::Float64,θ::Float64)
	Q_u_ji = (1-α)*(mc_d_i*(f_d_i+y_d_i)/P_u_i)*(P_u_i/P_u_ji)^θ
	return Q_u_ji
end
## eq28.1
function ell_d_i(α::Float64,mc_d_i::Float64,f_d_i::Float64,y_d_i::Float64,w_i::Float64)
	ell_d_i = (α*mc_d_i*(f_d_i+y_d_i))/w_i
	return ell_d_i
end
## eq28.2
function ell_u_i(f_u_i::Float64,y_u_i::Float64,A_u_i::Float64)
	ell_u_i = (f_u_i+y_u_i)/A_u_i
	return ell_u_i
end
## eq29.1
function P_s_ji(M_s_j::Float64,t_s_ji::Float64,p_s_ji::Float64,elasticity::Float64)
	P_s_ji = M_s_j^(1/(1-elasticity)) * (1+t_s_ji) * p_s_ji
	return P_s_ji 
end
## eq29.2
function P_s_i(P_s_ji::Float64,P_s_ii::Float64,elasticity::Float64)
	P_s_i = (P_s_ji^(1-elasticity)+P_s_ii^(1-elasticity))^(1/(1-elasticity))
	return P_s_i
end
## eq.30
function free_entry(elasticity::Float64,entrycost::Float64)
	y = (elasticity-1)*entrycost
	return y
end

