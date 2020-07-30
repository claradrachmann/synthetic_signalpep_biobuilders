using JuMP
using GLPK

println("big imports done")

function to_seq(mat)
    String(mapslices(row -> to_letter(row), mat, dims=2)[:,1])
end

function to_letter(row)
    if row[1] == 1
        'A'
    elseif row[2] == 1
        'C'
    elseif row[3] == 1
        'G'
    elseif row[4] == 1
        'T'
    else
        throw(ArgumentError("Must have one defined state"))
    end
end

function from_seq(seq)
    n = length(seq)
    mat = zeros(n, 4)
    for (i, c) in enumerate(seq)
        mat[i, letter_to_index[c]] = 1
    end
    return mat
end

letter_to_index = Dict('A' => 1, 'C' => 2, 'G' => 3, 'T' => 4)

function domesticate_for_var_seq(seq, var, model)
    n = length(seq)
    mat = from_seq(uppercase(seq))
    @constraint(model, [i = 1:(N-n)], sum(mat .* x[i:(i+n-1), :]) <= n - 1 + var)
end

function domesticate_for_seq(seq, model)
    domesticate_for_var_seq(seq, 0, model)
end

function domesticate_for_var_site(seq, var, model)
    domesticate_for_var_seq(seq, var, model)
    domesticate_for_var_seq(reverse_complement(seq), var, model)
end

function domesticate_for_site(seq, model)
    domesticate_for_var_site(seq, 0, model)
end

function reverse_complement(seq)
    reverse(map(x -> complements[x], uppercase(seq)))
end

complements = Dict('A' => 'T', 'C' => 'G', 'G' => 'C', 'T' => 'A')

println("reading started")

using CSV
pwm_raw = CSV.read(ARGS[1]*"-hom-upseqs.final.csv")
N = size(pwm_raw,1) - 3
#N = 716
#fkT = 1
pwm = convert(Matrix, pwm_raw)[1:N,:]
avg_entropy = sum(exp.(-pwm) .* pwm)/N/log(2)

println("reading done")

m = Model(with_optimizer(GLPK.Optimizer))

@variable(m, x[1:N,1:4], Bin)
#@variable(m, 0 <= x[1:N,1:4] <= 1)
@objective(m, Min, sum(pwm .* x)/log(2));

@constraint(m, [i = 1:N], sum(x[i, :]) == 1);

#println(to_seq(value.(x)))
#println("Objective value: ",objective_value(m))
#println()
#original = objective_value(m)

# MoClo promoter constraint
@constraint(m, x[N, 1] == 1)

# Global GC content
@constraint(m, 0.25 <= sum(x[:, [2,3]])/N <= 0.65);

println("basic constraints done")

# A/T homopolymers
max_homopoly_w = 10
domesticate_for_seq(repeat('A', max_homopoly_w), m)
domesticate_for_seq(repeat('T', max_homopoly_w), m)

# G/C homopolymers
max_homopoly_s = 10
domesticate_for_seq(repeat('C', max_homopoly_s), m)
domesticate_for_seq(repeat('G', max_homopoly_s), m)

println("homopolymer constraints done")

# Windowed GC content
window_n = 50
gc_min_window = 0.24
gc_max_window = 0.76
@constraint(m, [i = 1:(N-window_n)], gc_min_window <= sum(x[i:(i+window_n), [2,3]])/window_n <= gc_max_window)

println("synthability constraints done")

#optimize!(m)
#println(to_seq(value.(x)))
#println("Objective value: ",objective_value(m))
#println()
#synthable = objective_value(m)

### Domesticate for Golden Gate
## Domesticate for Loop
# BsaI
domesticate_for_site("GGTCTC", m)
# SapI
domesticate_for_site("GCTCTTC", m)
## Domesticate for MoClo
@variable(m, MoClo, Bin)
# BsaI
# BpiI
domesticate_for_var_site("GAAGAC", MoClo, m)
# Esp3I
domesticate_for_var_site("CGTCTC", MoClo, m)
## Domesticate for Mobius Assembly (like Loop but seemingly better)
@variable(m, Mobius, Bin)
# BsaI
# AarI
domesticate_for_var_site("CACCTGC", Mobius, m)
## Domesticate for GoldenBraid
@variable(m, GoldenBraid, Bin)
# BsaI
# Esp3I
domesticate_for_var_site("CGTCTC", GoldenBraid, m)
# BtgZI
domesticate_for_var_site("GCGATG", GoldenBraid, m)

### Domesticate for Chi recombination hotspots (E. coli)?
@variable(m, Chi, Bin)
domesticate_for_var_site("GCTGGTGG", Chi, m)
domesticate_for_var_seq(repeat('A', 8), Chi, m)

### Domesticate for linearization
@variable(m, Linearize, Bin)
# SwaI
domesticate_for_var_site("ATTTAAAT", Linearize, m)


### Domesticate for cloning into pSB1C3
# EcoRI
domesticate_for_site("GAATTC", m)
# PstI
domesticate_for_site("CTGCAG", m) # Was "CGTCAG"

#optimize!(m)
#println(to_seq(value.(x)))
#println("Objective value: ",objective_value(m))
#println()

println("domestication constraints done")

@objective(m, Min, sum(pwm .* x)/log(2) + 2*Linearize + 1.5*MoClo + 1*Mobius + 0.3*GoldenBraid + 0.2*Chi);
println("Expected numbers of changes at this level of conservation")
println("If full domestication: ", 3/avg_entropy)

output_file = open(ARGS[1]*"-promoter-samples.nt.fasta","w")

#optimize!(m)
optimize!(m)
S_0 = objective_value(m)

println("With actual domestication: ", value(2*Linearize + 1.5*(1-MoClo) + 1*(1-Mobius) + 0.3*(1-GoldenBraid) + 0.2*(1-Chi))/avg_entropy)
println()

println("Closest synthesizable consensus generated")
println("Objective value: ",value(sum(pwm .* x)/log(2)))
write(output_file,">"*ARGS[1]*"-synthesizable-consensus\n")
write(output_file,to_seq(value.(x)))
write(output_file,"\n")
println()

for v in all_variables(m)
    unset_binary(v)
end

function gumbel(x)
    x - log(-log(rand(Float64)))
end

function truncGumbel(x, ub)
    x - log(exp(-ub + x) - log(rand(Float64)))
end

function uniform_bits_like(xs)
    uniform_bits_n(length(xs))
end

function uniform_bits_n(n)
    map(round,rand(Float64, n))
end

dom_vars = [Linearize, MoClo, Mobius, GoldenBraid, Chi]
D = length(dom_vars)

using StatsBase

dom_logweights = [2.0, 1.5, 1.0, 0.3, 0.2]
dom_weights = 2 .^ (-dom_logweights .- 1.0)

println("Generating noisy samples")
for fkT in [1.0, 1.5, 2.0]
    for sample_n in 1:parse(Int,ARGS[2])
        println("fkT: ", fkT)

        attempts = 0

        pwm = convert(Matrix, pwm_raw)[1:N,:]/fkT
        pos_weights = exp.(-pwm)

        s_star = nothing
        s_last = (nothing, nothing)
        while s_star == nothing
            s_dom = map(x -> x ? 1 : 0, rand(Float64, D) .< dom_weights)
            s_pos = map(x -> sample(Weights(pos_weights[x,:])),1:N)
            s = (s_dom,s_pos)

            if attempts == 0
                dom_changes = 1:D
                pos_changes = 1:N
            else
                dom_changes = findall(!, s_last[1] .== s[1])
                pos_changes = findall(!, s_last[2] .== s[2])
            end
            for i in vcat(dom_changes, map(x -> x+5, pos_changes))
                if i <= D
                    fix(dom_vars[i], s[1][i], force=true)
                else
                    i = i-D
                    for j in 1:4
                        fix(x[i,j], s[2][i] == j ? 1 : 0, force=true)
                    end
                end
            end
            s_last = s

            optimize!(m)
            if termination_status(m) == JuMP.MOI.OPTIMAL
                s_star = s
            end
            attempts += 1
            if attempts % 100 == 0
                println("Attempts: ", attempts)
            end
        end

        println("Attempts: ", attempts)
        #println(s_star[1])
        println("")
        write(output_file,">"*ARGS[1]*"-noisy-fkT-"*string(fkT)*"-"*string(sample_n)*"\n")
        write(output_file,String(map(x -> "ACGT"[x], s_star[2])))
        write(output_file,"\n")

    end
end

close(output_file)
