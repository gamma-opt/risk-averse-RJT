using JuMP, Gurobi, DelimitedFiles, DecisionProgramming

j = parse(Int64, ARGS[1])


# Calculate the utility of the path
function profit(c,a,f,t,N)
    if f != t
        return 0
    else
        utility = sum(-c[k]*(a[k]-1) for k in 1:N)
        if f == 1 # No failure
            utility += 100 
        end
        return utility
    end
end

# Calculate probability of failure 
function p_fail(b,c,l,a,f,p_low,p_high)
    p = l==1 ? p_low : p_high
    prob_fail = p/exp(b*sum(c.*(a.-1)))
    if f == 1
        return 1-prob_fail
    else
        return prob_fail
    end
end

# A helper function for returning an array of possible paths (cartesian indices) 
# for nodes whose number of states is defined by the vector arr
function cart(arr)
    cart_arr = CartesianIndices(zeros(Tuple(arr)))
    index_arr = []
    for index in cart_arr
        push!(index_arr, [elem for elem in Tuple(index)])
    end
    return index_arr
end

function solve_nmonitoring(N; verbose=false)

    c_k = rand(N)
    b = 0.03
    fortification(k, a) = [c_k[k], 0][a]

    @info("Creating the influence diagram.")
    diagram = InfluenceDiagram()

    add_node!(diagram, ChanceNode("L", [], ["high", "low"]))

    for i in 1:N
        add_node!(diagram, ChanceNode("R$i", ["L"], ["high", "low"]))
        add_node!(diagram, DecisionNode("A$i", ["R$i"], ["yes", "no"]))
    end

    add_node!(diagram, ChanceNode("F", ["L", ["A$i" for i in 1:N]...], ["failure", "success"]))

    add_node!(diagram, ValueNode("T", ["F", ["A$i" for i in 1:N]...]))

    generate_arcs!(diagram)

    X_L = [rand(), 0]
    X_L[2] = 1.0 - X_L[1]
    add_probabilities!(diagram, "L", X_L)

    for i in 1:N
        x_R, y_R = rand(2)
        X_R = ProbabilityMatrix(diagram, "R$i")
        X_R["high", "high"] = max(x_R, 1-x_R)
        X_R["high", "low"] = 1 - max(x_R, 1-x_R)
        X_R["low", "low"] = max(y_R, 1-y_R)
        X_R["low", "high"] = 1-max(y_R, 1-y_R)
        add_probabilities!(diagram, "R$i", X_R)
    end

    X_F = ProbabilityMatrix(diagram, "F")
    x_F, y_F = rand(2)
    for s in paths([State(2) for i in 1:N])
        denominator = exp(b * sum(fortification(k, a) for (k, a) in enumerate(s)))
        X_F[1, s..., 1] = max(x_F, 1-x_F) / denominator
        X_F[1, s..., 2] = 1.0 - X_F[1, s..., 1]
        X_F[2, s..., 1] = min(y_F, 1-y_F) / denominator
        X_F[2, s..., 2] = 1.0 - X_F[2, s..., 1]
    end
    add_probabilities!(diagram, "F", X_F)

    Y_T = UtilityMatrix(diagram, "T")
    for s in paths([State(2) for i in 1:N])
        cost = sum(-fortification(k, a) for (k, a) in enumerate(s))
        Y_T[1, s...] = 0 + cost
        Y_T[2, s...] = 100 + cost
    end
    add_utilities!(diagram, "T", Y_T)


    @info("Creating the decision model.")
    generate_diagram!(diagram)
    model = Model()
    z = DecisionVariables(model, diagram, names=true)
    variables = RJTVariables(model, diagram, z, names=true)
    α = 0.15
    CVaR = conditional_value_at_risk(model, diagram, variables, α)
    @objective(model, Max, CVaR)

    model2 = Model()
    z2 = DecisionVariables(model2, diagram, names=true)
    variables2 = PathCompatibilityVariables(model2, diagram, z2)
    CVaR = conditional_value_at_risk(model2, diagram, variables2, α)
    @objective(model2, Max, CVaR)

    @info("Starting the optimization process.")
    optimizer = optimizer_with_attributes(
        () -> Gurobi.Optimizer()
    )
    set_optimizer(model, optimizer)
    set_optimizer(model2, optimizer)

    optimize!(model)
    optimize!(model2)
    
    return solve_time(model), solve_time(model2)
end

N_arr = collect(1:6)
sol_times = zeros(2,length(N_arr))
for (i, n_stages) in enumerate(N_arr)
    t1, t2 = solve_nmonitoring(n_stages)
    sol_times[1,n_stages] = t1 
    sol_times[2,n_stages] = t2 
end

writedlm("results/nmonitoring_indicator_"*string(j)*".csv", sol_times, ',')
