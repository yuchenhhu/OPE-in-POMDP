using Random, Distributions, LinearAlgebra, RollingFunctions, Statistics, JLD

function y_func(x,h)
    d = Normal()
    eps = rand(d, 1)[1]
    x*h+0.1*eps
end

function hmm_sim(pi_on,ui,t,burn)
    
    x0 = 0; h0 = 0
    t = t + burn
    
    # create matrices for storage
    x = fill(99, t); h = fill(99, t)
    w = fill(99, t)
    y = fill(99.99, t)
    
    # get initial state & treatment
    p_trans = (ui+0+h0)/4
    x[1] = x0 + rand(Binomial(1,p_trans), 1)[1]
    h[1] = h0
    w[1] = rand(Binomial(1,pi_on[x[1]+1]), 1)[1]
    y[1] = y_func(x[1],h[1])
    
    # transition
    for i = 2:t
        p_trans1 = (ui+w[i-1]+h[i-1])/4
        p_trans2 = (ui+1-w[i-1]+h[i-1])/4
        
        if x[i-1] == 0
            h[i] = h[i-1]
            x[i] = x[i-1] + rand(Binomial(1,p_trans1), 1)[1]
        else
            if rand(Binomial(1,p_trans2), 1)[1]==1
                h[i] = 1
                x[i] = 1
            else
                h[i] = 0
                x[i] = 0
            end
        end
        w[i] = rand(Binomial(1,pi_on[x[i]+1]), 1)[1]
        y[i] = y_func(x[i],h[i])
    end
    return [x[(burn+1):end], w[(burn+1):end], y[(burn+1):end]]
end

function hmm_sim_n(pi_on,n,t,burn = 50)
    # generate chains of different types
    types = rand(Binomial(1,0.5), n)
    chain_x = Array{Int64}(undef,n,t); chain_w= Array{Int64}(undef,n,t); chain_y = zeros(n,t)
    for i = 1:n
        chain_on = hmm_sim(pi_on,types[i],t,burn)
        chain_x[i,:] = chain_on[1]
        chain_w[i,:] = chain_on[2]
        chain_y[i,:] = chain_on[3]
    end
    return [chain_x, chain_w, chain_y]
end

function off_est(chain_x,chain_w,chain_y,pi_on,pi_off,t,k)
    
    if k==-1
        est = mean(chain_y)
        return est
    else
        all_one = fill(1,t)
        weights_ind = chain_w.*pi_off[chain_x.+1]./pi_on[chain_x.+1] .+
            (all_one-chain_w).*(all_one-pi_off[chain_x.+1])./(all_one-pi_on[chain_x.+1])
        weights_roll = rolling(prod,weights_ind,(k+1))
        chain_y_weighted = weights_roll.*chain_y[(k+1):end]
        est = mean(chain_y_weighted)
        return est
    end
end

function all_k(n,t,k_all,pi_on,ey_off,C0,C00)
    # generate the chain
    chain_on = hmm_sim_n(pi_on,n,t)
    chain_x = Int.(chain_on[1])
    chain_w = Int.(chain_on[2])
    chain_y = chain_on[3]
    
    # prepare for estimation
    ey_off_est = zeros(length(k_all)) # point estimate
    sd_est = zeros(length(k_all)) # variance estimate 
    # prepare for lepski's selection
    analyt_low = fill(99.99,length(k_all)); analyt_high = fill(99.99,length(k_all))
    analyt_mse = 99.99; marker = false
    k_all = reverse(k_all)
  
    # k_max
    ey_off_est_temp = zeros(n)
    for j = 1:n
        ey_off_est_temp[j] = off_est(chain_x[j,:],chain_w[j,:],chain_y[j,:],
            pi_on,pi_off,t,k_all[1])
    end
    ey_off_est[1] = mean(ey_off_est_temp); sd_est[1] = std(ey_off_est_temp)
    analyt_low[1] = ey_off_est[1] - sd_est[1]/sqrt(n)
    analyt_high[1] = ey_off_est[1] + sd_est[1]/sqrt(n)
    slope_low = analyt_low[1]; slope_high = analyt_high[1]
  
    for i = 2:length(k_all)
        ey_off_est_temp = zeros(n)
        for j = 1:n
            ey_off_est_temp[j] = off_est(chain_x[j,:],chain_w[j,:],chain_y[j,:],
                pi_on,pi_off,t,k_all[i])
        end
        ey_off_est[i] = mean(ey_off_est_temp); sd_est[i] = std(ey_off_est_temp)
        analyt_low[i] = ey_off_est[i] - sd_est[i]/sqrt(n)
        analyt_high[i] = ey_off_est[i] + sd_est[i]/sqrt(n)
        if (((analyt_high[i]<slope_low)|(analyt_low[i]>slope_high))&(marker==false))
            analyt_mse = (ey_off_est[i-1]-ey_off)^2
            marker = true
			best_k = k_all[i-1]
        end
        slope_low = max(slope_low,analyt_low[i]); slope_high = min(slope_high,analyt_high[i])
    end
    
    if (analyt_mse==99.99) 
        analyt_mse = (ey_off_est[length(k_all)]-ey_off)^2
		best_k = k_all[length(k_all)]
    end
    
    # k chosen by formula (11)
    all_mse = reverse((ey_off_est[:].-ey_off).^2)
    theo_k = round(Int,C00*log(C0*n*t))
    theo_mse = all_mse[theo_k+2]
    
    return([reverse(ey_off_est),reverse(sd_est), # ests and sd ests
        push!(push!(all_mse,analyt_mse),theo_mse), # MSE of lepski
        best_k,theo_k]) # k selected by lepski
end

seed_num = rand(1:1000000, 1)[1]
Random.seed!(seed_num)

# policies of interest
pi_on = [0.5; 0.5]
pi_off = [1; 0]

# calculate the optimal constant C0 based on (11)
t0 = 1/(log(4/3))
M1_2 = 1
M2 = 1
zeta = 2
C0 = 4*M1_2/(M2+2*M1_2*(exp(zeta)/(exp(zeta)-1)))
C00 = t0/(t0*zeta+2)

# start simulation
ey_off = 0.24448574
k_all = -1:8
t = 300
n_all = collect(10:5:30).^2
B = 100

res1 = zeros(0)
res2 = zeros(0)
res3 = zeros(0)
res4 = zeros(0)
for b in 1:B
    for n in n_all
        results = all_k(n,t,k_all,pi_on,ey_off,C0,C00)
        append!(res1, results[1])
        append!(res2, results[2])
        append!(res3, results[3])
        append!(res4, results[4])
    end
end

save("res/out_" * string(seed_num) * ".jld", "data", [res1,res2,res3,res4])
