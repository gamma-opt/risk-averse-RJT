using Logging
using JuMP, Gurobi
using DecisionProgramming, Plots, DelimitedFiles, Statistics

NN = 5
n_sample = 50
 
times_pigfarm_RJT = zeros(n_sample,4)
times_pigfarm_DP = zeros(n_sample,4)
 
for i in 1:n_sample
    result_pigfarm = readdlm("results/pigfarm_"*string(i)*".csv", ',')
    times_pigfarm_RJT[i,:] .= result_pigfarm[1,:]
    times_pigfarm_DP[i,:] .= result_pigfarm[2,:]
end

solve_times_RJT = Array(mean(times_pigfarm_RJT, dims=1))
solve_times_DP = Array(mean(times_pigfarm_DP, dims=1))
std_RJT = Array(std(times_pigfarm_RJT, dims=1))
std_DP = Array(std(times_pigfarm_DP, dims=1))


scatter(collect(2:NN), mean(times_pigfarm_RJT, dims=1)', xlabel="Number of breeding periods", ylabel="Solution time (s)", label="RJT", color=palette(:tab10)[1], legend=:topleft, yaxis=:log, yticks=[10.0^k for k in -3:3], xticks=[2,3,4,5])
scatter!(collect(2:NN), mean(times_pigfarm_DP, dims=1)', label="DecisionProgramming", color=palette(:tab10)[2])
plot!(collect(2:NN), mean(times_pigfarm_RJT, dims=1)', yerr = std_RJT, msc = palette(:tab10)[1], label=false, ls=:dash, color=palette(:tab10)[1])
plot!(collect(2:NN), mean(times_pigfarm_DP, dims=1)', yerr = std_DP, msc = palette(:tab10)[2], label=false, ls=:dash, color=palette(:tab10)[2])
Plots.pdf("sol_times_DP_vs_RJT_pig_farm.pdf")
