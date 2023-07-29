using Random, Distributions, LinearAlgebra, RollingFunctions, Statistics, JLD

function reward(s1, a)
    w1 = -2; w2 = -3 
    w3 = -1; w4 = -1  
    R = 0
    if (s1 < 70) 
        R = R + w2 
    end
    if (s1 > 150) 
        R = R + w1
    end
    if (s1 > 70) & (s1 < 80) 
        R = R + w4 
    end
    if (s1 > 120) & (s1 < 150) 
        R = R +  w3
    end
    return(R)
end

function genmod_on(n, capT, u = reward, burn = 50, init_states = 999)
  
    capT = capT + burn; p = 8
    S = zeros(n, p, capT); A = zeros(n, capT); R = zeros(n, capT)
    S_gl = zeros(n, capT); S_ex = zeros(n, capT)
  
    beta = [0.9, 0.1, 0.1, -0.01, -0.01, -2, -4]
  
    mu_food = 78; sigma_food = 10; d_food = Normal(mu_food,sigma_food)
    mu_act = 819; sigma_act = 10; d_act = Normal(mu_act,sigma_act)
    mu_act_mod = 31; sigma_act_mod = 5; d_act_mod = Normal(mu_act_mod,sigma_act_mod)
    sigma_er = 5; d_er = Normal(0,sigma_er)
    sigma_id = 10; d_id = Normal(0,sigma_id)
    du = Uniform()
    
    for i = 1:n 
        eps_id = rand(d_id, 1)[1]
        if init_states==999
            if rand(du,1)[1] < 0.2
                S[i, 2, 1] = S[i, 2, 1] + max(rand(d_food, 1)[1], 0)
            end
            if rand(du,1)[1] < 0.2
                S[i, 3, 1] = S[i, 3, 1] + max(rand(d_act, 1)[1], 0)
            end
            if rand(du,1)[1] < 0.4
                S[i, 3, 1] = S[i, 3, 1] + max(rand(d_act_mod, 1)[1], 0)
            end
            
            S[i, 1, 1] = rand(Normal(100,25), 1)[1]
            if S[i, 1, 1] > 250 # prevent simulated glucose from going too high
                S[i, 1, 1] = 250
            end
            if S[i, 1, 1] < 50 # prevent simulated glucose from going too low
                S[i, 1, 1] = 50
            end
            
            S[i, 4, 1] = 100; S[i, 5, 1] = 0
            S[i, 6, 1] = 0; S[i, 7, 1] = 0; S[i, 8, 1] = 0
            A[i, 1] = rand(Binomial(1,0.3), 1)[1]
            S_gl[i, 1] = S[i, 4, 1]
            S_ex[i, 1] = S[i, 3, 1] + S[i, 7, 1]
        end
        
        if init_states!=999
            S[i, :, 1] = init.states[i, ]
            A[i, 1] = rand(Binomial(1,0.3), 1)[1]
        end
        
        if rand(du,1)[1] < 0.2
            S[i, 2, 2] = S[i, 2, 2] + max(rand(d_food, 1)[1], 0)
        end
        if rand(du,1)[1] < 0.2
            S[i, 3, 2] = S[i, 3, 2] + max(rand(d_act, 1)[1], 0)
        end
        if rand(du,1)[1] < 0.4
            S[i, 3, 2] = S[i, 3, 2] + max(rand(d_act_mod, 1)[1], 0)
        end
        
        S[i, 4, 2] = S[i, 1, 1] # glucose t-1
        S[i, 5, 2] = S[i, 2, 1] # food t-1
        S[i, 6, 2] = S[i, 5, 1] # food t-2
        S[i, 7, 2] = S[i, 3, 1] # activity t-1
        S[i, 8, 2] = S[i, 7, 1] # activity t-2
        
        eps = rand(d_er, 1)[1]
        S[i, 1, 2] = (1 - beta[1]) * 100 + beta[1] * S[i, 1, 1] + beta[2] * S[i, 2, 1] + beta[3] * S[i, 2, 1] + beta[4] * S[i, 3, 1] + beta[5] * S[i, 3, 1] + beta[6] * A[i, 1] + beta[7] * A[i, 1] + eps + eps_id
        if S[i, 1, 2] > 250
            S[i, 1, 2] = 250  # prevent simulated glucose from going too high
        end
        if (S[i, 1, 2] < 50) 
            S[i, 1, 2] = 50  # prevent simulated glucose from going too low
        end
        A[i, 2] = rand(Binomial(1,0.3), 1)[1]
        S_gl[i, 2] = S[i, 4, 2]
        S_ex[i, 2] = S[i, 3, 2] + S[i, 7, 2]
        
        for t=3:capT
            if rand(du,1)[1] < 0.2
                S[i, 2, t] = S[i, 2, t] + max(rand(d_food, 1)[1], 0)
            end
            if rand(du,1)[1] < 0.2
                S[i, 3, t] = S[i, 3, t] + max(rand(d_act, 1)[1], 0)
            end
            if rand(du,1)[1] < 0.4
                S[i, 3, t] = S[i, 3, t] + max(rand(d_act_mod, 1)[1], 0)
            end
            eps = rand(d_er, 1)[1]
            S[i, 1, t] = (1 - beta[1]) * 100 + beta[1] * S[i, 1, t - 1] + beta[2] * S[i, 2, t - 1] + beta[3] * S[i, 2, t - 2] + beta[4] * S[i, 3, t - 1] + beta[5] * S[i, 3, t - 2] + beta[6] * A[i, t - 1] + beta[7] * A[i, t - 2] + eps + eps_id
            if S[i, 1, t] > 250
                S[i, 1, t] = 250  # prevent simulated glucose from going too high
            end
            if S[i, 1, t] < 50
                S[i, 1, t] = 50  # prevent simulated glucose from going too low
            end
            S[i, 4, t] = S[i, 1, t - 1]
            S[i, 5, t] = S[i, 2, t - 1]
            S[i, 6, t] = S[i, 5, t - 1]
            S[i, 7, t] = S[i, 3, t - 1]
            S[i, 8, t] = S[i, 7, t - 1]  
            
            A[i, t] = rand(Binomial(1,0.3), 1)[1]
            S_gl[i, t] = S[i, 4, t]
            S_ex[i, t] = S[i, 3, t] + S[i, 7, t]
        end
        
        for t = 1:(capT - 1)
            R[i, t] = u(S[i, 1, t], A[i, t])
        end
    end
    
    if init_states==999
        S = S[:,:, (burn+1):end]
        S_gl = S_gl[:, (burn+1):end]
        S_ex = S_ex[:, (burn+1):end]
        A = A[:, (burn+1):end]
        R = R[:, (burn+1):end]
    else
        S = S[:,:, 1:(capT - burn)]
        S_gl = S_gl[:, 1:(capT - burn)]
        S_ex = S_ex[:, 1:(capT - burn)]
        A = A[:, 1:(capT - burn)]
        R = R[:, 1:(capT - burn)]
    end
    
    [S_gl, S_ex, A, R]
end

function off_est(chain_on,n,t,k)
    
    chain_gl = chain_on[1]; chain_ex = chain_on[2]
    chain_a = chain_on[3]
    chain_r = chain_on[4]
    weights_roll = zeros(n, t-k)
    
    if k==-1
        est = mean(chain_r,dims=2)
        return est
    else
        p = (chain_gl.>=110).*(chain_ex.<=100)
        all_one = fill(1,n,t)
        weights_ind = chain_a.*p./fill(0.3,n,t) .+
            (all_one-chain_a).*(all_one-p)./fill(0.7,n,t)
        for i in 1:n
            weights_roll[i,:] = rolling(prod,weights_ind[i,:],(k+1))
        end
        chain_r_weighted = weights_roll.*chain_r[:,(k+1):end]
        est = mean(chain_r_weighted,dims=2)
        return est
    end
end

function all_k(n,t,k_all,ey_off)
    # generate observed chain
    chain_on = genmod_on(n,t)
    # prepare for estimation
    ey_off_est = zeros(length(k_all)) # point estimate
    sd_est = zeros(length(k_all)) # variance estimate 
    # prepare for lepski's selection
    analyt_low = zeros(length(k_all)); analyt_high = zeros(length(k_all))
    analyt_mse = 99.99; marker = false
    k_all = reverse(k_all)
    
    # k_max
    results = off_est(chain_on,n,t,k_all[1])
    ey_off_est[1] = mean(results); sd_est[1] = std(results)
    analyt_low[1] = ey_off_est[1] - 1.96*sd_est[1]/sqrt(n)
    analyt_high[1] = ey_off_est[1] + 1.96*sd_est[1]/sqrt(n)
    slope_low = analyt_low[1]; slope_high = analyt_high[1]
    
    # k < k_max
    for i = 2:length(k_all)
        results = off_est(chain_on,n,t,k_all[i])
        ey_off_est[i] = mean(results); sd_est[i] = std(results)
        analyt_low[i] = ey_off_est[i] - 1.96*sd_est[i]/sqrt(n)
        analyt_high[i] = ey_off_est[i] + 1.96*sd_est[i]/sqrt(n)
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

    return([reverse(ey_off_est),reverse(sd_est), # ests and sd ests
        push!(reverse((ey_off_est[:].-ey_off).^2),analyt_mse), # MSE of lepski
        best_k]) # k selected by lepski
end

seed_num = rand(1:1000000, 1)[1]
Random.seed!(seed_num)

# start simulation
ey_off = -1.69938
k_all = -1:6
t_all = collect(50:50:150)
n_all = collect(50:50:150)
B = 30

res1 = zeros(0)
res2 = zeros(0)
res3 = zeros(0)
res4 = zeros(0)
for b in 1:B
    for n in n_all
        for t in t_all
            results = all_k(n,t,k_all,ey_off)
            append!(res1, results[1])
            append!(res2, results[2])
            append!(res3, results[3])
            append!(res4, results[4])
        end
    end
end

save("/res/out_" * string(seed_num) * ".jld", "data", [res1,res2,res3,res4])
