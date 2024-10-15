using Logging
using JuMP, HiGHS
using DecisionProgramming

const N = 4

# Replace the RJT_generation algorithm from DP to create the RJT
function ID_to_RJT(diagram::InfluenceDiagram)
    C_rjt = Dict{Name, Vector{Name}}()
    A_rjt = []
    names = get_keys(diagram.Nodes)
    for j in length(diagram.Nodes):-1:1
        C_j = copy(get_values(diagram.I_j)[j])
        push!(C_j, names[j])
        for a in A_rjt 
            if a[1] == names[j]
                push!(C_j, setdiff(C_rjt[a[2]], [a[2]])...)
            end
        end
        C_j = unique(C_j)
        C_j_aux = sort([(elem, findfirst(isequal(elem), names)) for elem in C_j], by = last)
        C_j = [C_j_tuple[1] for C_j_tuple in C_j_aux]
        C_rjt[names[j]] = C_j
        if length(C_rjt[names[j]]) > 1
            u = maximum([findfirst(isequal(name), names) for name in setdiff(C_j, [names[j]])])
            push!(A_rjt, (names[u], names[j]))
        end
    end
    prepend!(C_rjt["T2"],["H1"])
    prepend!(C_rjt["D2"],["H1"])
    prepend!(C_rjt["H3"],["H1"])
    prepend!(C_rjt["T3"],["H1"])
    prepend!(C_rjt["T3"],["H2"])
    prepend!(C_rjt["D3"],["H1"])
    prepend!(C_rjt["D3"],["H2"])
    prepend!(C_rjt["H4"],["H1"])
    prepend!(C_rjt["H4"],["H2"])
    return C_rjt, A_rjt
end

@info("Creating the influence diagram.")
diagram = InfluenceDiagram()

add_node!(diagram, ChanceNode("H1", [], ["ill", "healthy"]))
for i in 1:N-1
    # Testing result
    add_node!(diagram, ChanceNode("T$i", ["H$i"], ["positive", "negative"]))
    # Decision to treat
    add_node!(diagram, DecisionNode("D$i", ["T$i"], ["treat", "pass"]))
    # Cost of treatment
    add_node!(diagram, ValueNode("V$i", ["D$i"]))
    # Health of next period
    add_node!(diagram, ChanceNode("H$(i+1)", ["H$(i)", "D$(i)"], ["ill", "healthy"]))
end
add_node!(diagram, ValueNode("V4", ["H$N"]))

generate_arcs!(diagram)

# Add probabilities for node H1
add_probabilities!(diagram, "H1", [0.1, 0.9])

# Declare probability matrix for health nodes H_2, ... H_N-1, which have identical information sets and states
X_H = ProbabilityMatrix(diagram, "H2")
X_H["healthy", "pass", :] = [0.2, 0.8]
X_H["healthy", "treat", :] = [0.1, 0.9]
X_H["ill", "pass", :] = [0.9, 0.1]
X_H["ill", "treat", :] = [0.5, 0.5]

# Declare probability matrix for test result nodes T_1...T_N
X_T = ProbabilityMatrix(diagram, "T1")
X_T["ill", "positive"] = 0.8
X_T["ill", "negative"] = 0.2
X_T["healthy", "negative"] = 0.9
X_T["healthy", "positive"] = 0.1

for i in 1:N-1
    add_probabilities!(diagram, "T$i", X_T)
    add_probabilities!(diagram, "H$(i+1)", X_H)
end

for i in 1:N-1
    add_utilities!(diagram, "V$i", [-100.0, 0.0])
end

add_utilities!(diagram, "V4", [300.0, 1000.0])

@info("Creating the decision model.")
model, z, x_s = generate_model(diagram, model_type="RJT")
elements = filter( element -> count(x -> occursin("ill",x),split(string(element))) >= 1,x_s.data["H4"].statevars)
# println(size(x_s.data["H4"].statevars))
# println(size(elements))
@constraint(model, sum(elements) <= 0.35)


@info("Starting the optimization process.")
optimizer = optimizer_with_attributes(
    () -> HiGHS.Optimizer()
)
set_optimizer(model, optimizer)

optimize!(model)

@info("Extracting results.")

Z = DecisionStrategy(diagram, z)
S_probabilities = StateProbabilities(diagram, Z)
U_distribution = UtilityDistribution(diagram, Z)

@info("Printing decision strategy:")
print_decision_strategy(diagram, Z, S_probabilities)

@info("Printing utility distribution.")
print_utility_distribution(U_distribution)

@info("Printing statistics")
print_statistics(U_distribution)

@info("State probabilities:")
print_state_probabilities(diagram, S_probabilities, [["H$i" for i in 1:N]...])
print_state_probabilities(diagram, S_probabilities, [["T$i" for i in 1:N-1]...])
print_state_probabilities(diagram, S_probabilities, [["D$i" for i in 1:N-1]...])

@info("Conditional state probabilities")
for state in ["ill", "healthy"]
    S_probabilities2 = StateProbabilities(diagram, Z, "H1", state, S_probabilities)
    print_state_probabilities(diagram, S_probabilities2, [["H$i" for i in 1:N]...])
    print_state_probabilities(diagram, S_probabilities2, [["T$i" for i in 1:N-1]...])
    print_state_probabilities(diagram, S_probabilities2, [["D$i" for i in 1:N-1]...])
end