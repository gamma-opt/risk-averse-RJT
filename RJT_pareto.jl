using JuMP, Gurobi, Plots, Random,  LaTeXStrings
rng = MersenneTwister(1234)
n_stages = 5
n = 1
EU_opt = zeros(4)
Cvar_opt = zeros(4)
EU_non_opt = zeros(30)
Cvar_non_opt = zeros(30)
Cvar_cur = zeros(4)
Cvar_cur[1] = 100
strategies = Dict()
for n in 1:4
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
    @variable(model, μ_hdh[H,D,H] >= 0)
    @variable(model, μ_dht[D,H,T] >= 0)
    @variable(model, μ_dhtd[D,H,T,D] >= 0)
    @variable(model, μ_dhdh[D,H,D,H] >= 0)
    @variable(model, μ_ddht[D,D,H,T] >= 0)
    @variable(model, μ_ddhtd[D,D,H,T,D] >= 0)
    @variable(model, μ_ddhdh[D,D,H,D,H] >= 0)
    @variable(model, μ_dddht[D,D,D,H,T] >= 0)
    @variable(model, μ_dddhtd[D,D,D,H,T,D] >= 0)
    @variable(model, μ_dddhdh[D,D,D,H,D,H] >= 0)
    @variable(model, μ_ddddht[D,D,D,D,H,T] >= 0)
    @variable(model, μ_ddddhtd[D,D,D,D,H,T,D] >= 0)
    @variable(model, μ_ddddhdh[D,D,D,D,H,D,H] >= 0)
    @variable(model, μ_dddddh[D,D,D,D,D,H] >= 0)
    #@variable(model, μ_ddddh[D,D,D,D,H] >= 0)

    # Decision strategy variable
    @variable(model, δ[N,T,D], Bin)

    # Objective function
    @objective(model, Max, sum(μ_dddddh[d,dd,ddd,d4,d5,h]*(cost_treat[d,d]+cost_treat[dd,dd]+cost_treat[ddd,ddd]+cost_treat[d4,d4]+cost_treat[d5,d5]+profit_sell[h,h]) for d in D, dd in D, ddd in D,d4 in D,d5 in D, h in H))
    #@objective(model, Max, sum(μ_ddddh[d,dd,ddd,d4,h]*(cost_treat[d,d]+cost_treat[dd,dd]+cost_treat[ddd,ddd]+cost_treat[d4,d4]+profit_sell[h,h]) for d in D, dd in D, ddd in D,d4 in D, h in H))

    profit = zeros(length(D),length(D),length(D),length(D),length(D),length(H))
    for d in D, dd in D, ddd in D,dddd in D,ddddd in D, h in H
        profit[d,dd,ddd,dddd,ddddd,h] = cost_treat[d,d]+cost_treat[dd,dd]+cost_treat[ddd,ddd]+cost_treat[dddd,dddd]+cost_treat[ddddd,ddddd]+profit_sell[h,h]
    end


    # Probability distributions μ sum to 1
    @constraint(model, sum(μ_h0) == 1)
    @constraint(model, sum(μ_ht[h,t] for h in H, t in T) == 1)
    @constraint(model, sum(μ_dht[d,h,t] for d in D, h in H, t in T) == 1)
    @constraint(model, sum(μ_ddht[d,d2,h,t] for d in D, d2 in D, h in H, t in T) == 1)
    @constraint(model,  sum(μ_htd[h,t,d] for h in H, t in T, d in D) == 1)
    @constraint(model,  sum(μ_dhtd[d,h,t,dd] for d in D, h in H, t in T, dd in D) == 1)
    @constraint(model,  sum(μ_ddhtd[d,d2,h,t,dd] for d in D, d2 in D, h in H, t in T, dd in D) == 1)
    @constraint(model,  sum(μ_hdh[h,d,hh] for h in H, d in D, hh in H) == 1)
    @constraint(model,  sum(μ_dhdh[d1,h,d,hh] for d1 in D, h in H, d in D, hh in H) == 1)
    @constraint(model, sum(μ_ddhdh[d,dd,h,ddd,hh] for d in D, dd in D, ddd in D, h in H, hh in H) == 1)
    @constraint(model, sum(μ_dddht[d,dd,ddd,h,t] for d in D, dd in D, ddd in D, h in H, t in T) == 1)
    @constraint(model, sum(μ_dddhtd[d,dd,ddd,h,t,d4] for d in D, dd in D, ddd in D, h in H, t in T, d4 in D) == 1)
    @constraint(model, sum(μ_dddhdh[d,dd,ddd,h,d4,hh] for d in D, dd in D, ddd in D, h in H, hh in H, d4 in D) == 1)
    @constraint(model, sum(μ_ddddht[d,dd,ddd,d4,h,t] for d in D, dd in D, ddd in D, h in H, t in T, d4 in D) == 1)
    @constraint(model, sum(μ_ddddhtd[d,dd,ddd,d4,h,t,d5] for d in D, dd in D, ddd in D, h in H, t in T, d4 in D, d5 in D) == 1)
    @constraint(model, sum(μ_ddddhdh[d,dd,ddd,d4,h,d5,hh]  for d in D, dd in D, ddd in D, h in H, hh in H, d4 in D, d5 in D) == 1)
    @constraint(model, sum(μ_dddddh[d,dd,ddd,d4,d5,h] for d in D, dd in D, ddd in D, d4 in D, d5 in D, h in H) == 1)
    #@constraint(model, sum(μ_ddddh[d,dd,ddd,d4,h] for d in D, dd in D, ddd in D,d4 in D, h in H) == 1)

    


    # Local consistency constraints
    @constraint(model, [h in H], μ_h0[h] == sum(μ_ht[h,t] for t in T))
    @constraint(model, [h in H, t in T], μ_ht[h,t] == sum(μ_htd[h,t,d] for d in D))
    @constraint(model, [d in D, h in H], sum(μ_htd[h,t,d] for t in T) == sum(μ_hdh[h,d,hh] for hh in H))
    @constraint(model, [h in H, d in D], sum(μ_hdh[hh,d,h] for hh in H) == sum(μ_dht[d,h,t] for t in T))
    @constraint(model, [h in H, d in D, t in T], μ_dht[d,h,t] == sum(μ_dhtd[d,h,t,dd] for dd in D))
    @constraint(model, [d in D, h in H, dd in D], sum(μ_dhtd[d,h,t,dd] for t in T) == sum(μ_dhdh[d,h,dd,hh] for hh in H))
    @constraint(model, [d in D, hh in H, dd in D], sum(μ_dhdh[d,h,dd,hh] for h in H) == sum(μ_ddht[d,dd,hh,t] for t in T))
    @constraint(model, [d in D, hh in H, dd in D,t in T], μ_ddht[d,dd,hh,t] == sum(μ_ddhtd[d,dd,hh,t,ddd] for ddd in D))
    @constraint(model, [d in D, dd in D, ddd in D, hh in H], sum(μ_ddhtd[d,dd,hh,t,ddd] for t in T) == sum(μ_ddhdh[d,dd,hh,ddd, h] for h in H))
    @constraint(model, [d in D, dd in D, ddd in D, h in H], sum(μ_ddhdh[d, dd, hh,ddd, h] for hh in H) == sum(μ_dddht[d,dd,ddd,h, t] for t in T))
    @constraint(model, [d in D, dd in D, ddd in D, h in H, t in T], μ_dddht[d, dd,ddd, h,t] == sum(μ_dddhtd[d,dd,ddd,h, t,d4] for d4 in D))
    @constraint(model, [d in D, dd in D, ddd in D, h in H, d4 in D], sum(μ_dddhtd[d, dd,ddd, h,t,d4] for t in T) == sum(μ_dddhdh[d,dd,ddd,h,d4,hh] for hh in H))
    @constraint(model, [d in D, dd in D, ddd in D, h in H, d4 in D], sum(μ_ddddht[d, dd,ddd,d4, h,t] for t in T) == sum(μ_dddhdh[d,dd,ddd,hh,d4,h] for hh in H))
    @constraint(model, [d in D, dd in D, ddd in D, h in H, d4 in D, t in T], μ_ddddht[d, dd,ddd,d4, h,t] == sum(μ_ddddhtd[d,dd,ddd,d4,h,t,d5] for d5 in D))
    @constraint(model, [d in D, dd in D, ddd in D, h in H, d4 in D, d5 in D], sum(μ_ddddhdh[d, dd,ddd,d4, h,d5,hh] for hh in H) == sum(μ_ddddhtd[d,dd,ddd,d4,h,t,d5] for t in T))
    @constraint(model, [d in D, dd in D, ddd in D, h in H, d4 in D, d5 in D], sum(μ_ddddhdh[d, dd,ddd,d4, hh,d5,h] for hh in H) == μ_dddddh[d,dd,ddd,d4,d5,h])


    @constraint(model, [n in N, t in T], sum(δ[n,t,d] for d in D) == 1)

    # Moments μ_{\breve{C}_v} (the moments from above, but with the last variable dropped out)
    # Some such moments were already created, for example μ_htd becomes μ_ht
    
    @variable(model, cμ_h[H] >= 0)
    @variable(model, cμ_ht[H,T] >= 0)
    @variable(model, cμ_hd[H,D] >= 0)
    @variable(model, cμ_dh[D,H] >= 0)
    @variable(model, cμ_dht[D,H,T] >= 0)
    @variable(model, cμ_dhd[D,H,D] >= 0)
    @variable(model, cμ_ddh[D,D,H] >= 0)
    @variable(model, cμ_ddht[D,D,H,T] >= 0)
    @variable(model, cμ_ddhd[D,D,H,D] >= 0)
    @variable(model, cμ_dddh[D,D,D,H] >= 0)
    @variable(model, cμ_dddht[D,D,D,H,T] >= 0)
    @variable(model, cμ_dddhd[D,D,D,H,D] >= 0)
    @variable(model, cμ_ddddh[D,D,D,D,H] >= 0)
    @variable(model, cμ_ddddht[D,D,D,D,H,T] >= 0)
    @variable(model, cμ_ddddhd[D,D,D,D,H,D] >= 0)

    # μ_{\breve{C}_v} = ∑_{x_v} μ_{C_v}
    @constraint(model,[h in H], sum(μ_ht[h,t] for t in T) == cμ_h[h])
    @constraint(model,[h in H, t in T], sum(μ_htd[h,t,d] for d in D) == cμ_ht[h,t])
    @constraint(model,[h in H, d in D], sum(μ_hdh[h,d,hh] for hh in H) == cμ_hd[h,d])
    @constraint(model,[h in H, d in D], sum(μ_dht[d,h,t] for t in T) == cμ_dh[d,h])
    @constraint(model,[h in H, d in D, t in T], sum(μ_dhtd[d,h,t,dd] for dd in D) == cμ_dht[d,h,t])
    @constraint(model,[h in H, d in D, dd in D], sum(μ_dhdh[d,h,dd,hh] for hh in H) == cμ_dhd[d,h,dd])
    @constraint(model,[h in H, d in D, dd in D], sum(μ_ddht[d,dd,h,t] for t in T) == cμ_ddh[d,dd,h])
    @constraint(model,[h in H, d in D, dd in D,t in T], sum(μ_ddhtd[d,dd,h,t,ddd] for ddd in D) == cμ_ddht[d,dd,h,t])
    @constraint(model,[h in H, d in D,dd in D, ddd in D], sum(μ_ddhdh[d,dd,h,ddd,hh] for hh in H) == cμ_ddhd[d,dd,h,ddd])
    @constraint(model,[h in H, d in D,dd in D, ddd in D], sum(μ_dddht[d,dd,ddd,h,t] for t in T) == cμ_dddh[d,dd,ddd,h])
    @constraint(model,[h in H, d in D,dd in D, ddd in D,t in T], sum(μ_dddhtd[d,dd,ddd,h,t,d4] for d4 in D) == cμ_dddht[d,dd,ddd,h,t] )
    @constraint(model,[h in H, d in D,dd in D, ddd in D,d4 in D], sum(μ_dddhdh[d,dd,ddd,h,d4,hh] for hh in H) == cμ_dddhd[d,dd,ddd,h,d4])
    @constraint(model,[h in H, d in D,dd in D, ddd in D,d4 in D], sum(μ_ddddht[d,dd,ddd,d4,h,t] for t in T) == cμ_ddddh[d,dd,ddd,d4,h])
    @constraint(model,[h in H, d in D,dd in D, ddd in D,d4 in D, t in T], sum(μ_ddddhtd[d,dd,ddd,d4,h,t,d5] for d5 in D) == cμ_ddddht[d,dd,ddd,d4,h,t])
    @constraint(model,[h in H, d in D,dd in D, ddd in D,d4 in D, d5 in D], sum(μ_ddddhdh[d,dd,ddd,d4,h,d5,hh] for hh in H) == cμ_ddddhd[d,dd,ddd,d4,h,d5])

    # Factorization constraints (Corollary 3 in Parmentier et al.)
    @constraint(model, [h in H], μ_h0[h] == p_init[h])
    @constraint(model,[h in H, t in T], μ_ht[h,t]  == μ_h0[h]*p_res[h,t])
    @constraint(model,[h in H, t in T, d in D], μ_htd[h,t,d] == cμ_ht[h,t]*δ[1,t,d])
    @constraint(model,[h in H, d in D, hh in H], μ_hdh[h,d,hh] == cμ_hd[h,d]*p_trans[h,d,hh])
    @constraint(model,[h in H, d in D, t in T], μ_dht[d,h,t] == cμ_dh[d,h]*p_res[h,t])
    @constraint(model,[h in H, d in D, t in T, dd in D], μ_dhtd[d,h,t,dd] == cμ_dht[d,h,t]*δ[2,t,dd])
    @constraint(model,[h in H, d in D, dd in D, hh in H], μ_dhdh[d,h,dd,hh] == cμ_dhd[d,h,dd]*p_trans[h,dd,hh])
    @constraint(model,[h in H, d in D, dd in D, t in T], μ_ddht[d,dd,h,t] == cμ_ddh[d,dd,h]*p_res[h,t])
    @constraint(model,[h in H, d in D, dd in D, t in T, ddd in D], μ_ddhtd[d,dd,h,t,ddd] == cμ_ddht[d,dd,h,t]*δ[3,t,ddd])
    @constraint(model,[h in H, d in D, hh in H, dd in D, ddd in D], μ_ddhdh[d,dd,h,ddd,hh] == cμ_ddhd[d,dd,h,ddd]*p_trans[h,ddd,hh])
    @constraint(model,[h in H, d in D,dd in D, ddd in D, t in T], μ_dddht[d,dd,ddd,h,t] == cμ_dddh[d,dd,ddd,h]*p_res[h,t])
    @constraint(model,[h in H, d in D,dd in D, ddd in D,t in T,d4 in D], μ_dddhtd[d,dd,ddd,h,t,d4] == cμ_dddht[d,dd,ddd,h,t]*δ[4,t,d4])
    @constraint(model,[h in H, d in D,dd in D, ddd in D,d4 in D, hh in H], μ_dddhdh[d,dd,ddd,h,d4,hh] == cμ_dddhd[d,dd,ddd,h,d4]*p_trans[h,d4,hh])
    @constraint(model,[h in H, d in D,dd in D, ddd in D,d4 in D, t in T], μ_ddddht[d,dd,ddd,d4,h,t]  == cμ_ddddh[d,dd,ddd,d4,h]*p_res[h,t])
    @constraint(model,[h in H, d in D,dd in D, ddd in D,d4 in D, t in T,d5 in D], μ_ddddhtd[d,dd,ddd,d4,h,t,d5] == cμ_ddddht[d,dd,ddd,d4,h,t]*δ[5,t,d5])
    @constraint(model,[h in H, d in D,dd in D, ddd in D,d4 in D, d5 in D, hh in H], μ_ddddhdh[d,dd,ddd,d4,h,d5,hh] == cμ_ddddhd[d,dd,ddd,d4,h,d5]*p_trans[h,d5,hh])

    M = maximum(profit) - minimum(profit)
    ϵ = minimum(diff(unique(profit))) / 2
    η = @variable(model)
    ρ′_s = Dict{Int64, VariableRef}()
    α = 0.2
    for u in unique(profit)
        λ = @variable(model, binary=true)
        λ′ = @variable(model, binary=true)
        ρ = @variable(model)
        ρ′ = @variable(model)
        @constraint(model, η - u ≤ M * λ)
        @constraint(model, η - u ≥ (M + ϵ) * λ - M)
        @constraint(model, η - u ≤ (M + ϵ) * λ′ - ϵ)
        @constraint(model, η - u ≥ M * (λ′ - 1))
        @constraint(model, 0 ≤ ρ)
        @constraint(model, 0 ≤ ρ′)
        @constraint(model, ρ ≤ λ)
        @constraint(model, ρ′ ≤ λ′)
        @constraint(model, ρ ≤ ρ′)
        @constraint(model, ρ′ ≤ sum(μ_dddddh[d,dd,ddd,dddd,ddddd,h] for d in D, dd in D, ddd in D, dddd in D, ddddd in D, h in H if profit[d,dd,ddd,dddd,ddddd,h] == u))
        @constraint(model, (sum(μ_dddddh[d,dd,ddd,dddd,ddddd,h] for d in D, dd in D, ddd in D,dddd in D, ddddd in D, h in H if profit[d,dd,ddd,dddd,ddddd,h] == u) - (1 - λ)) ≤ ρ)
        ρ′_s[u] = ρ′
    end
    @constraint(model, sum(values(ρ′_s)) == α)

    @constraint(model, (sum(ρ_bar * u for (u, ρ_bar) in ρ′_s)/α) >= Cvar_cur[n]) 



    optimize!(model)

    solution_summary(model)
    println(termination_status(model))

    EU_opt[n] = objective_value(model)
    Cvar_opt[n] = sum(value.(ρ_bar) * u for (u, ρ_bar) in ρ′_s)/α
    if n < 4
        Cvar_cur[n+1] = sum(value.(ρ_bar) * u for (u, ρ_bar) in ρ′_s)/α + 0.5
    end
    strategies[n] = δ
end

for n in 1:30
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
    @variable(model, μ_hdh[H,D,H] >= 0)
    @variable(model, μ_dht[D,H,T] >= 0)
    @variable(model, μ_dhtd[D,H,T,D] >= 0)
    @variable(model, μ_dhdh[D,H,D,H] >= 0)
    @variable(model, μ_ddht[D,D,H,T] >= 0)
    @variable(model, μ_ddhtd[D,D,H,T,D] >= 0)
    @variable(model, μ_ddhdh[D,D,H,D,H] >= 0)
    @variable(model, μ_dddht[D,D,D,H,T] >= 0)
    @variable(model, μ_dddhtd[D,D,D,H,T,D] >= 0)
    @variable(model, μ_dddhdh[D,D,D,H,D,H] >= 0)
    @variable(model, μ_ddddht[D,D,D,D,H,T] >= 0)
    @variable(model, μ_ddddhtd[D,D,D,D,H,T,D] >= 0)
    @variable(model, μ_ddddhdh[D,D,D,D,H,D,H] >= 0)
    @variable(model, μ_dddddh[D,D,D,D,D,H] >= 0)
    #@variable(model, μ_ddddh[D,D,D,D,H] >= 0)

    # Decision strategy variable
    mo = zeros(n_stages,length(T),2)
    for n in N, t in T
        bit = bitrand(rng, 1)
        if bit[1] == 1
            mo[n,t,1] = 1
            mo[n,t,2] = 0
        else
            mo[n,t,1] = 0
            mo[n,t,2] = 1
        end
    end

    # Objective function
    @objective(model, Max, sum(μ_dddddh[d,dd,ddd,d4,d5,h]*(cost_treat[d,d]+cost_treat[dd,dd]+cost_treat[ddd,ddd]+cost_treat[d4,d4]+cost_treat[d5,d5]+profit_sell[h,h]) for d in D, dd in D, ddd in D,d4 in D,d5 in D, h in H))
    #@objective(model, Max, sum(μ_ddddh[d,dd,ddd,d4,h]*(cost_treat[d,d]+cost_treat[dd,dd]+cost_treat[ddd,ddd]+cost_treat[d4,d4]+profit_sell[h,h]) for d in D, dd in D, ddd in D,d4 in D, h in H))

    profit = zeros(length(D),length(D),length(D),length(D),length(D),length(H))
    for d in D, dd in D, ddd in D,dddd in D,ddddd in D, h in H
        profit[d,dd,ddd,dddd,ddddd,h] = cost_treat[d,d]+cost_treat[dd,dd]+cost_treat[ddd,ddd]+cost_treat[dddd,dddd]+cost_treat[ddddd,ddddd]+profit_sell[h,h]
    end


    # Probability distributions μ sum to 1
    @constraint(model, sum(μ_h0) == 1)
    @constraint(model, sum(μ_ht[h,t] for h in H, t in T) == 1)
    @constraint(model, sum(μ_dht[d,h,t] for d in D, h in H, t in T) == 1)
    @constraint(model, sum(μ_ddht[d,d2,h,t] for d in D, d2 in D, h in H, t in T) == 1)
    @constraint(model,  sum(μ_htd[h,t,d] for h in H, t in T, d in D) == 1)
    @constraint(model,  sum(μ_dhtd[d,h,t,dd] for d in D, h in H, t in T, dd in D) == 1)
    @constraint(model,  sum(μ_ddhtd[d,d2,h,t,dd] for d in D, d2 in D, h in H, t in T, dd in D) == 1)
    @constraint(model,  sum(μ_hdh[h,d,hh] for h in H, d in D, hh in H) == 1)
    @constraint(model,  sum(μ_dhdh[d1,h,d,hh] for d1 in D, h in H, d in D, hh in H) == 1)
    @constraint(model, sum(μ_ddhdh[d,dd,h,ddd,hh] for d in D, dd in D, ddd in D, h in H, hh in H) == 1)
    @constraint(model, sum(μ_dddht[d,dd,ddd,h,t] for d in D, dd in D, ddd in D, h in H, t in T) == 1)
    @constraint(model, sum(μ_dddhtd[d,dd,ddd,h,t,d4] for d in D, dd in D, ddd in D, h in H, t in T, d4 in D) == 1)
    @constraint(model, sum(μ_dddhdh[d,dd,ddd,h,d4,hh] for d in D, dd in D, ddd in D, h in H, hh in H, d4 in D) == 1)
    @constraint(model, sum(μ_ddddht[d,dd,ddd,d4,h,t] for d in D, dd in D, ddd in D, h in H, t in T, d4 in D) == 1)
    @constraint(model, sum(μ_ddddhtd[d,dd,ddd,d4,h,t,d5] for d in D, dd in D, ddd in D, h in H, t in T, d4 in D, d5 in D) == 1)
    @constraint(model, sum(μ_ddddhdh[d,dd,ddd,d4,h,d5,hh]  for d in D, dd in D, ddd in D, h in H, hh in H, d4 in D, d5 in D) == 1)
    @constraint(model, sum(μ_dddddh[d,dd,ddd,d4,d5,h] for d in D, dd in D, ddd in D, d4 in D, d5 in D, h in H) == 1)

    


    # Local consistency constraints
    @constraint(model, [h in H], μ_h0[h] == sum(μ_ht[h,t] for t in T))
    @constraint(model, [h in H, t in T], μ_ht[h,t] == sum(μ_htd[h,t,d] for d in D))
    @constraint(model, [d in D, h in H], sum(μ_htd[h,t,d] for t in T) == sum(μ_hdh[h,d,hh] for hh in H))
    @constraint(model, [h in H, d in D], sum(μ_hdh[hh,d,h] for hh in H) == sum(μ_dht[d,h,t] for t in T))
    @constraint(model, [h in H, d in D, t in T], μ_dht[d,h,t] == sum(μ_dhtd[d,h,t,dd] for dd in D))
    @constraint(model, [d in D, h in H, dd in D], sum(μ_dhtd[d,h,t,dd] for t in T) == sum(μ_dhdh[d,h,dd,hh] for hh in H))
    @constraint(model, [d in D, hh in H, dd in D], sum(μ_dhdh[d,h,dd,hh] for h in H) == sum(μ_ddht[d,dd,hh,t] for t in T))
    @constraint(model, [d in D, hh in H, dd in D,t in T], μ_ddht[d,dd,hh,t] == sum(μ_ddhtd[d,dd,hh,t,ddd] for ddd in D))
    @constraint(model, [d in D, dd in D, ddd in D, hh in H], sum(μ_ddhtd[d,dd,hh,t,ddd] for t in T) == sum(μ_ddhdh[d,dd,hh,ddd, h] for h in H))
    @constraint(model, [d in D, dd in D, ddd in D, h in H], sum(μ_ddhdh[d, dd, hh,ddd, h] for hh in H) == sum(μ_dddht[d,dd,ddd,h, t] for t in T))
    @constraint(model, [d in D, dd in D, ddd in D, h in H, t in T], μ_dddht[d, dd,ddd, h,t] == sum(μ_dddhtd[d,dd,ddd,h, t,d4] for d4 in D))
    @constraint(model, [d in D, dd in D, ddd in D, h in H, d4 in D], sum(μ_dddhtd[d, dd,ddd, h,t,d4] for t in T) == sum(μ_dddhdh[d,dd,ddd,h,d4,hh] for hh in H))
    @constraint(model, [d in D, dd in D, ddd in D, h in H, d4 in D], sum(μ_ddddht[d, dd,ddd,d4, h,t] for t in T) == sum(μ_dddhdh[d,dd,ddd,hh,d4,h] for hh in H))
    @constraint(model, [d in D, dd in D, ddd in D, h in H, d4 in D, t in T], μ_ddddht[d, dd,ddd,d4, h,t] == sum(μ_ddddhtd[d,dd,ddd,d4,h,t,d5] for d5 in D))
    @constraint(model, [d in D, dd in D, ddd in D, h in H, d4 in D, d5 in D], sum(μ_ddddhdh[d, dd,ddd,d4, h,d5,hh] for hh in H) == sum(μ_ddddhtd[d,dd,ddd,d4,h,t,d5] for t in T))
    @constraint(model, [d in D, dd in D, ddd in D, h in H, d4 in D, d5 in D], sum(μ_ddddhdh[d, dd,ddd,d4, hh,d5,h] for hh in H) == μ_dddddh[d,dd,ddd,d4,d5,h])


    # Moments μ_{\breve{C}_v} (the moments from above, but with the last variable dropped out)
    # Some such moments were already created, for example μ_htd becomes μ_ht
    
    @variable(model, cμ_h[H] >= 0)
    @variable(model, cμ_ht[H,T] >= 0)
    @variable(model, cμ_hd[H,D] >= 0)
    @variable(model, cμ_dh[D,H] >= 0)
    @variable(model, cμ_dht[D,H,T] >= 0)
    @variable(model, cμ_dhd[D,H,D] >= 0)
    @variable(model, cμ_ddh[D,D,H] >= 0)
    @variable(model, cμ_ddht[D,D,H,T] >= 0)
    @variable(model, cμ_ddhd[D,D,H,D] >= 0)
    @variable(model, cμ_dddh[D,D,D,H] >= 0)
    @variable(model, cμ_dddht[D,D,D,H,T] >= 0)
    @variable(model, cμ_dddhd[D,D,D,H,D] >= 0)
    @variable(model, cμ_ddddh[D,D,D,D,H] >= 0)
    @variable(model, cμ_ddddht[D,D,D,D,H,T] >= 0)
    @variable(model, cμ_ddddhd[D,D,D,D,H,D] >= 0)

    # μ_{\breve{C}_v} = ∑_{x_v} μ_{C_v}
    @constraint(model,[h in H], sum(μ_ht[h,t] for t in T) == cμ_h[h])
    @constraint(model,[h in H, t in T], sum(μ_htd[h,t,d] for d in D) == cμ_ht[h,t])
    @constraint(model,[h in H, d in D], sum(μ_hdh[h,d,hh] for hh in H) == cμ_hd[h,d])
    @constraint(model,[h in H, d in D], sum(μ_dht[d,h,t] for t in T) == cμ_dh[d,h])
    @constraint(model,[h in H, d in D, t in T], sum(μ_dhtd[d,h,t,dd] for dd in D) == cμ_dht[d,h,t])
    @constraint(model,[h in H, d in D, dd in D], sum(μ_dhdh[d,h,dd,hh] for hh in H) == cμ_dhd[d,h,dd])
    @constraint(model,[h in H, d in D, dd in D], sum(μ_ddht[d,dd,h,t] for t in T) == cμ_ddh[d,dd,h])
    @constraint(model,[h in H, d in D, dd in D,t in T], sum(μ_ddhtd[d,dd,h,t,ddd] for ddd in D) == cμ_ddht[d,dd,h,t])
    @constraint(model,[h in H, d in D,dd in D, ddd in D], sum(μ_ddhdh[d,dd,h,ddd,hh] for hh in H) == cμ_ddhd[d,dd,h,ddd])
    @constraint(model,[h in H, d in D,dd in D, ddd in D], sum(μ_dddht[d,dd,ddd,h,t] for t in T) == cμ_dddh[d,dd,ddd,h])
    @constraint(model,[h in H, d in D,dd in D, ddd in D,t in T], sum(μ_dddhtd[d,dd,ddd,h,t,d4] for d4 in D) == cμ_dddht[d,dd,ddd,h,t] )
    @constraint(model,[h in H, d in D,dd in D, ddd in D,d4 in D], sum(μ_dddhdh[d,dd,ddd,h,d4,hh] for hh in H) == cμ_dddhd[d,dd,ddd,h,d4])
    @constraint(model,[h in H, d in D,dd in D, ddd in D,d4 in D], sum(μ_ddddht[d,dd,ddd,d4,h,t] for t in T) == cμ_ddddh[d,dd,ddd,d4,h])
    @constraint(model,[h in H, d in D,dd in D, ddd in D,d4 in D, t in T], sum(μ_ddddhtd[d,dd,ddd,d4,h,t,d5] for d5 in D) == cμ_ddddht[d,dd,ddd,d4,h,t])
    @constraint(model,[h in H, d in D,dd in D, ddd in D,d4 in D, d5 in D], sum(μ_ddddhdh[d,dd,ddd,d4,h,d5,hh] for hh in H) == cμ_ddddhd[d,dd,ddd,d4,h,d5])

    # Factorization constraints (Corollary 3 in Parmentier et al.)
    @constraint(model, [h in H], μ_h0[h] == p_init[h])
    @constraint(model,[h in H, t in T], μ_ht[h,t]  == μ_h0[h]*p_res[h,t])
    @constraint(model,[h in H, t in T, d in D], μ_htd[h,t,d] == cμ_ht[h,t]*mo[1,t,d])
    @constraint(model,[h in H, d in D, hh in H], μ_hdh[h,d,hh] == cμ_hd[h,d]*p_trans[h,d,hh])
    @constraint(model,[h in H, d in D, t in T], μ_dht[d,h,t] == cμ_dh[d,h]*p_res[h,t])
    @constraint(model,[h in H, d in D, t in T, dd in D], μ_dhtd[d,h,t,dd] == cμ_dht[d,h,t]*mo[2,t,dd])
    @constraint(model,[h in H, d in D, dd in D, hh in H], μ_dhdh[d,h,dd,hh] == cμ_dhd[d,h,dd]*p_trans[h,dd,hh])
    @constraint(model,[h in H, d in D, dd in D, t in T], μ_ddht[d,dd,h,t] == cμ_ddh[d,dd,h]*p_res[h,t])
    @constraint(model,[h in H, d in D, dd in D, t in T, ddd in D], μ_ddhtd[d,dd,h,t,ddd] == cμ_ddht[d,dd,h,t]*mo[3,t,ddd])
    @constraint(model,[h in H, d in D, hh in H, dd in D, ddd in D], μ_ddhdh[d,dd,h,ddd,hh] == cμ_ddhd[d,dd,h,ddd]*p_trans[h,ddd,hh])
    @constraint(model,[h in H, d in D,dd in D, ddd in D, t in T], μ_dddht[d,dd,ddd,h,t] == cμ_dddh[d,dd,ddd,h]*p_res[h,t])
    @constraint(model,[h in H, d in D,dd in D, ddd in D,t in T,d4 in D], μ_dddhtd[d,dd,ddd,h,t,d4] == cμ_dddht[d,dd,ddd,h,t]*mo[4,t,d4])
    @constraint(model,[h in H, d in D,dd in D, ddd in D,d4 in D, hh in H], μ_dddhdh[d,dd,ddd,h,d4,hh] == cμ_dddhd[d,dd,ddd,h,d4]*p_trans[h,d4,hh])
    @constraint(model,[h in H, d in D,dd in D, ddd in D,d4 in D, t in T], μ_ddddht[d,dd,ddd,d4,h,t]  == cμ_ddddh[d,dd,ddd,d4,h]*p_res[h,t])
    @constraint(model,[h in H, d in D,dd in D, ddd in D,d4 in D, t in T,d5 in D], μ_ddddhtd[d,dd,ddd,d4,h,t,d5] == cμ_ddddht[d,dd,ddd,d4,h,t]*mo[5,t,d5])
    @constraint(model,[h in H, d in D,dd in D, ddd in D,d4 in D, d5 in D, hh in H], μ_ddddhdh[d,dd,ddd,d4,h,d5,hh] == cμ_ddddhd[d,dd,ddd,d4,h,d5]*p_trans[h,d5,hh])

    M = maximum(profit) - minimum(profit)
    ϵ = minimum(diff(unique(profit))) / 2
    η = @variable(model)
    ρ′_s = Dict{Int64, VariableRef}()
    α = 0.2
    for u in unique(profit)
        λ = @variable(model, binary=true)
        λ′ = @variable(model, binary=true)
        ρ = @variable(model)
        ρ′ = @variable(model)
        @constraint(model, η - u ≤ M * λ)
        @constraint(model, η - u ≥ (M + ϵ) * λ - M)
        @constraint(model, η - u ≤ (M + ϵ) * λ′ - ϵ)
        @constraint(model, η - u ≥ M * (λ′ - 1))
        @constraint(model, 0 ≤ ρ)
        @constraint(model, 0 ≤ ρ′)
        @constraint(model, ρ ≤ λ)
        @constraint(model, ρ′ ≤ λ′)
        @constraint(model, ρ ≤ ρ′)
        @constraint(model, ρ′ ≤ sum(μ_dddddh[d,dd,ddd,dddd,ddddd,h] for d in D, dd in D, ddd in D, dddd in D, ddddd in D, h in H if profit[d,dd,ddd,dddd,ddddd,h] == u))
        @constraint(model, (sum(μ_dddddh[d,dd,ddd,dddd,ddddd,h] for d in D, dd in D, ddd in D,dddd in D, ddddd in D, h in H if profit[d,dd,ddd,dddd,ddddd,h] == u) - (1 - λ)) ≤ ρ)
        ρ′_s[u] = ρ′
    end
    @constraint(model, sum(values(ρ′_s)) == α)


    optimize!(model)
    EU_non_opt[n] = objective_value(model)
    println(objective_value(model))
    Cvar_non_opt[n] = sum(value.(ρ_bar) * u for (u, ρ_bar) in ρ′_s)/α

end

p = plot(Cvar_non_opt, EU_non_opt, seriestype=:scatter, mc=:blue, label="dominated strategies")
plot!(Cvar_opt, EU_opt, seriestype=:scatter, mc=:orange, label="nondominated strategies")
plot!(legend=:outerbottom)
xlabel!(L"Cvar at $\alpha$=20%")
ylabel!("EU")

savefig("pareto.svg")







    
