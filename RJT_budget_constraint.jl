using JuMP, Gurobi
n_stages = 3
model = Model()
    optimizer = optimizer_with_attributes(
        () -> Gurobi.Optimizer(Gurobi.Env()),
        "IntFeasTol"      => 1e-9,
        "TimeLimit"       => 3600,
        # "DualReductions"  => 0,
    )
    set_optimizer(model, optimizer)

    # Index sets
    N = 1:n_stages
    H = 1:2 # Ill/healthy
    T = 1:2 # +/-
    D = 1:2 # Treat/pass
    # The reason we introduce states for value nodes is that they could be stochastic
    U_treat = 1:2 # Treat/pass 
    U_sell = 1:2 # Ill/healthy

    # Treatment costs 
    cost_treat = zeros(length(D),length(U_treat))
    cost_treat[1,1] = -100

    # Profit from selling the pig in the end
    profit_sell = zeros(length(H),length(U_sell))
    profit_sell[1,1] = 300
    profit_sell[2,2] = 1000

    # Initial probability of pig being ill/healthy
    p_init = [0.1 0.9]

    # Transition probabilities of the health states
    p_trans = zeros(length(H),length(D),length(H))
    p_trans[1,1,:] = [0.5 0.5] # Ill, treated
    p_trans[1,2,:] = [0.9 0.1] # Ill, not treated
    p_trans[2,1,:] = [0.1 0.9] # Healthy, treated
    p_trans[2,2,:] = [0.2 0.8] # Healthy, not treated

    # Probabilities of getting positive/negative result from the test
    p_res = zeros(length(H),length(T))
    p_res[1,:] = [0.8 0.2] # Ill
    p_res[2,:] = [0.1 0.9] # Healthy

    # Variables corresponding to the nodes in the RJT
    @variable(model, μ_h0[H] >= 0)
    @variable(model, μ_ht[H,T] >= 0)
    @variable(model, μ_htd[H,T,D] >= 0)
    @variable(model, μ_hdh[1:2,H,D,H] >= 0)
    @variable(model, μ_dht[D,H,T] >= 0)
    @variable(model, μ_dhtd[D,H,T,D] >= 0)
    @variable(model, μ_dhdh[D,H,D,H] >= 0)
    @variable(model, μ_ddht[D,D,H,T] >= 0)
    @variable(model, μ_ddhtd[D,D,H,T,D] >= 0)
    #@variable(model, μ_hdhl[H,D,H] >= 0)

    @variable(model, μ_du[N,D,U_treat] >= 0)
    @variable(model, μ_hu[H,U_sell] >= 0)
    # Decision strategy variable
    @variable(model, δ[N,T,D], Bin)

    # Objective function
    @objective(model, Max, sum(μ_du[n,d,u]*cost_treat[d,u] for n in N, d in D, u in U_treat) + sum(μ_hu[h,u]*profit_sell[h,u] for h in H, u in U_sell))

    # Probability distributions μ sum to 1
    @constraint(model, sum(μ_h0) == 1)
    @constraint(model, sum(μ_ht[h,t] for h in H, t in T) == 1)
    @constraint(model, sum(μ_dht[d,h,t] for d in D, h in H, t in T) == 1)
    @constraint(model, sum(μ_ddht[d,d2,h,t] for d in D, d2 in D, h in H, t in T) == 1)
    @constraint(model,  sum(μ_htd[h,t,d] for h in H, t in T, d in D) == 1)
    @constraint(model,  sum(μ_dhtd[d,h,t,dd] for d in D, h in H, t in T, dd in D) == 1)
    @constraint(model,  sum(μ_ddhtd[d,d2,h,t,dd] for d in D, d2 in D, h in H, t in T, dd in D) == 1)
    @constraint(model, [n in 1:2], sum(μ_hdh[n,h,d,hh] for h in H, d in D, hh in H) == 1)
    @constraint(model,  sum(μ_dhdh[d1,h,d,hh] for d1 in D, h in H, d in D, hh in H) == 1)
    @constraint(model, [n in N], sum(μ_du[n,d,u] for d in D, u in U_treat) == 1)
    @constraint(model, sum(μ_hu[h,u] for h in H, u in U_sell) == 1)


    # Local consistency constraints
    @constraint(model, [h in H], μ_h0[h] == sum(μ_ht[h,t] for t in T))
    @constraint(model, [h in H, t in T], μ_ht[h,t] == sum(μ_htd[h,t,d] for d in D))
    @constraint(model, [d in D], sum(μ_htd[h,t,d] for h in H, t in T) == sum(μ_du[1,d, u] for u in U_treat))
    @constraint(model, [d in D, h in H], sum(μ_htd[h,t,d] for t in T) == sum(μ_hdh[1,h,d,hh] for hh in H))
    @constraint(model, [h in H, d in D], sum(μ_hdh[1,hh,d,h] for hh in H) == sum(μ_dht[d,h,t] for t in T))
    @constraint(model, [h in H, d in D, t in T], μ_dht[d,h,t] == sum(μ_dhtd[d,h,t,dd] for dd in D))
    @constraint(model, [dd in D], sum(μ_dhtd[d,h,t,dd] for d in D, h in H, t in T) == sum(μ_du[2,dd, u] for u in U_treat))
    @constraint(model, [d in D, h in H, dd in D], sum(μ_dhtd[d,h,t,dd] for t in T) == sum(μ_dhdh[d,h,dd,hh] for hh in H))
    @constraint(model, [d in D, hh in H, dd in D], sum(μ_dhdh[d,h,dd,hh] for h in H) == sum(μ_ddht[d,dd,hh,t] for t in T))
    @constraint(model, [d in D, hh in H, dd in D,t in T], μ_ddht[d,dd,hh,t] == sum(μ_ddhtd[d,dd,hh,t,ddd] for ddd in D))
    @constraint(model, [ddd in D], sum(μ_ddhtd[d,dd,hh,t,ddd] for d in D, hh in H, dd in D,t in T) == sum(μ_du[3,ddd, u] for u in U_treat))
    @constraint(model, [ddd in D, hh in H], sum(μ_ddhtd[d,dd,hh,t,ddd] for d in D, dd in D,t in T) == sum(μ_hdh[2,hh,ddd, h] for h in H))
    @constraint(model, [h in H], sum(μ_hdh[2,hh,ddd, h] for hh in H,ddd in D) == sum(μ_hu[h,u] for u in U_sell))




    @constraint(model, [n in N, t in T], sum(δ[n,t,d] for d in D) == 1)

    # Moments μ_{\breve{C}_v} (the moments from above, but with the last variable dropped out)
    # Some such moments were already created, for example μ_htd becomes μ_ht
    
    @variable(model, cμ_h[1:2,H] >= 0)
    @variable(model, cμ_ht[H,T] >= 0)
    @variable(model, cμ_d[N,D] >= 0)
    @variable(model, cμ_hd[1:2,H,D] >= 0)
    @variable(model, cμ_dh[D,H] >= 0)
    @variable(model, cμ_dht[D,H,T] >= 0)
    @variable(model, cμ_dhd[D,H,D] >= 0)
    @variable(model, cμ_ddh[D,D,H] >= 0)
    @variable(model, cμ_ddht[D,D,H,T] >= 0)

    # μ_{\breve{C}_v} = ∑_{x_v} μ_{C_v}
    @constraint(model,[h in H], sum(μ_ht[h,t] for t in T) == cμ_h[1,h])
    @constraint(model,[h in H, t in T], sum(μ_htd[h,t,d] for d in D) == cμ_ht[h,t])
    @constraint(model,[n in N,d in D], sum(μ_du[n,d,u] for u in U_treat) == cμ_d[n,d])
    @constraint(model,[h in H, d in D], sum(μ_hdh[1,h,d,hh] for hh in H) == cμ_hd[1,h,d])
    @constraint(model,[h in H, d in D], sum(μ_dht[d,h,t] for t in T) == cμ_dh[d,h])
    @constraint(model,[h in H, d in D, t in T], sum(μ_dhtd[d,h,t,dd] for dd in D) == cμ_dht[d,h,t])
    @constraint(model,[h in H, d in D, dd in D], sum(μ_dhdh[d,h,dd,hh] for hh in H) == cμ_dhd[d,h,dd])
    @constraint(model,[h in H, d in D, dd in D], sum(μ_ddht[d,dd,h,t] for t in T) == cμ_ddh[d,dd,h])
    @constraint(model,[h in H, d in D, dd in D,t in T], sum(μ_ddhtd[d,dd,h,t,ddd] for ddd in D) == cμ_ddht[d,dd,h,t])
    @constraint(model,[h in H, d in D], sum(μ_hdh[2,h,d,hh] for hh in H) == cμ_hd[2,h,d])
    @constraint(model,[h in H], sum(μ_hu[h,u] for u in U_sell) == cμ_h[2,h])

    # Factorization constraints (Corollary 3 in Parmentier et al.)
    @constraint(model, [h in H], μ_h0[h] == p_init[h])
    @constraint(model,[h in H, t in T], μ_ht[h,t]  == μ_h0[h]*p_res[h,t])
    @constraint(model,[h in H, t in T, d in D], μ_htd[h,t,d] == cμ_ht[h,t]*δ[1,t,d])
    @constraint(model, [n in N, d in D, u in U_treat], μ_du[n,d,u] == cμ_d[n,d]*((d==u) ? 1 : 0))
    @constraint(model,[h in H, d in D, hh in H], μ_hdh[1,h,d,hh] == cμ_hd[1,h,d]*p_trans[h,d,hh])
    @constraint(model,[h in H, d in D, t in T], μ_dht[d,h,t] == cμ_dh[d,h]*p_res[h,t])
    @constraint(model,[h in H, d in D, t in T, dd in D], μ_dhtd[d,h,t,dd] == cμ_dht[d,h,t]*δ[2,t,dd])
    @constraint(model,[h in H, d in D, dd in D, hh in H], μ_dhdh[d,h,dd,hh] == cμ_dhd[d,h,dd]*p_trans[h,dd,hh])
    @constraint(model,[h in H, d in D, dd in D, t in T], μ_ddht[d,dd,h,t] == cμ_ddh[d,dd,h]*p_res[h,t])
    @constraint(model,[h in H, d in D, dd in D, t in T, ddd in D], μ_ddhtd[d,dd,h,t,ddd] == cμ_ddht[d,dd,h,t]*δ[3,t,ddd])
    @constraint(model,[h in H, d in D, hh in H], μ_hdh[2,h,d,hh] == cμ_hd[2,h,d]*p_trans[h,d,hh])
    @constraint(model, [h in H, u in U_sell], μ_hu[h,u] == cμ_h[2,h]*((h==u) ? 1 : 0))

    #budget/logical constraint
    @constraint(model,sum(μ_ddhtd[1,1,h,t,1] + μ_ddhtd[2,1,h,t,1] + μ_ddhtd[1,2,h,t,1] +μ_ddhtd[1,1,h,t,2] for h in H, t in T) <= 0)


    optimize!(model)

    solution_summary(model)





    
