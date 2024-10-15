using Logging
using JuMP, Gurobi
using DecisionProgramming, Plots, DelimitedFiles, Statistics

NN = 6
n_sample = 50
 
times_nmonitoring_RJT = zeros(47,6)
times_nmonitoring_DP = zeros(47,6)
 
for i in 1:n_sample
    if !(i in [11, 12, 27])
        result_nmonitoring = readdlm("results_nmonitoring/nmonitoring_indicator_"*string(i)*".csv", ',')
        if i > 27
            times_nmonitoring_RJT[i-3,:] .= result_nmonitoring[1,:]
            times_nmonitoring_DP[i-3,:] .= result_nmonitoring[2,:]
        elseif i > 12
            times_nmonitoring_RJT[i-2,:] .= result_nmonitoring[1,:]
            times_nmonitoring_DP[i-2,:] .= result_nmonitoring[2,:]
        else
            times_nmonitoring_RJT[i,:] .= result_nmonitoring[1,:]
            times_nmonitoring_DP[i,:] .= result_nmonitoring[2,:]
        end
    end
end

solve_times_RJT = Array(mean(times_nmonitoring_RJT, dims=1))
solve_times_DP = Array(mean(times_nmonitoring_DP, dims=1))
std_RJT = Array(std(times_nmonitoring_RJT, dims=1))
std_DP = Array(std(times_nmonitoring_DP, dims=1))


scatter(collect(1:NN), mean(times_nmonitoring_RJT, dims=1)', xlabel="N", ylabel="Solution time (s)", label="RJT", color=palette(:tab10)[1], legend=:topleft, yaxis=:log, yticks=[10.0^k for k in -3:3], xticks=[1,2,3,4,5,6])
scatter!(collect(1:NN), mean(times_nmonitoring_DP, dims=1)', label="DecisionProgramming", color=palette(:tab10)[2])
plot!(collect(1:NN), mean(times_nmonitoring_RJT, dims=1)', yerr = std_RJT, msc = palette(:tab10)[1], label=false, ls=:dash, color=palette(:tab10)[1])
plot!(collect(1:NN), mean(times_nmonitoring_DP, dims=1)', yerr = std_DP, msc = palette(:tab10)[2], label=false, ls=:dash, color=palette(:tab10)[2])
Plots.pdf("sol_times_DP_vs_RJT_nmonitoring.pdf")
