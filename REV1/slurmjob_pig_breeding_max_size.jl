using Pkg
using JuMP, Gurobi, Random, DecisionProgramming, DelimitedFiles


function main()
    pp = 0
    Random.seed!(8)
    sol_times = zeros(2,NN-1)
    solve_times_DP = zeros(2,NN-1)
    solve_times_RJT = zeros(NN-1)
    N = 2
    while pp <= 10

        p_ill = rand()*0.2                  # Initial probability of pig being ill

        spec = 0.5 + rand()*0.5             # Specificity of the test
        sens = 0.5 + rand()*0.5             # Sensitivity of the test

        p_vals = sort(rand(4))              # Health state update probabilities
        p_hpi = p_vals[1]                   # Probability of ill untreated pig to become healthy, assumed to be the smallest of the random probabilities
        p_hti = p_vals[2]                   # Probability of ill treated pig to become healthy
        p_hph = p_vals[3]                   # Probability of healthy untreated pig to remain healthy
        p_hth = p_vals[4]                   # Probability of healthy treated pig to remain healthy, assumed to be the largest of the random probabilities

        c_treat = -rand()*100               # Cost of treatment
        u_sell = sort(rand(2).*1000).+100   # Utility of selling a pig, ill or healthy (healthy is more valuable)


        @info("Creating the influence diagram.")
        diagram = InfluenceDiagram()

        add_node!(diagram, ChanceNode("H1", [], ["ill", "healthy"]))
        for i in 1:N-1
            # Testing result
            add_node!(diagram, ChanceNode("T$i", ["H$i"], ["positive", "negative"]))
            # Decision to treat
            add_node!(diagram, DecisionNode("D$i", ["T$i"], ["treat", "pass"]))
            # Cost of treatment
            # Health of next period
            add_node!(diagram, ChanceNode("H$(i+1)", ["H$(i)", "D$(i)"], ["ill", "healthy"]))
        end
        add_node!(diagram, ValueNode("V4", ["H$N", ["D$i" for i in 1:N-1]...]))

        generate_arcs!(diagram)

        # Add probabilities for node H1
        add_probabilities!(diagram, "H1", [p_ill, 1-p_ill])

        # Declare probability matrix for health nodes H_2, ... H_N-1, which have identical information sets and states
        X_H = ProbabilityMatrix(diagram, "H2")
        X_H["healthy", "pass", :] = [1-p_hph, p_hph]
        X_H["healthy", "treat", :] = [1-p_hth, p_hth]
        X_H["ill", "pass", :] = [1-p_hpi, p_hpi]
        X_H["ill", "treat", :] = [1-p_hti, p_hti]

        # Declare probability matrix for test result nodes T_1...T_N
        X_T = ProbabilityMatrix(diagram, "T1")
        X_T["ill", "positive"] = 1-spec
        X_T["ill", "negative"] = spec
        X_T["healthy", "negative"] = 1-sens
        X_T["healthy", "positive"] = sens

        for i in 1:N-1
            add_probabilities!(diagram, "T$i", X_T)
            add_probabilities!(diagram, "H$(i+1)", X_H)
        end

        Y_T = UtilityMatrix(diagram, "V4")
        cost_of_treatment = [c_treat 0.0]
        for s in paths([State(2) for i in 1:N-1])
            cost = sum(cost_of_treatment[a] for a in s)
            Y_T[1, s...] = u_sell[1] + cost
            Y_T[2, s...] = u_sell[2] + cost
        end


        add_utilities!(diagram, "V4", Y_T)

        @info("Creating the decision model.")
        generate_diagram!(diagram)
        model = Model()
        z = DecisionVariables(model, diagram, names=true)
        variables = RJTVariables(model, diagram, z, names=true)
        α = 0.15
        CVaR = conditional_value_at_risk(model, diagram, variables, α)
        @objective(model, Max, CVaR)

        @info("Starting the optimization process.")
        optimizer = optimizer_with_attributes(
            () -> Gurobi.Optimizer()
        )
        set_optimizer(model, optimizer)

        optimize!(model)
        pp = solve_time(model)
        N = N +1
    end
   println(pp)
   println(N)
end

main()