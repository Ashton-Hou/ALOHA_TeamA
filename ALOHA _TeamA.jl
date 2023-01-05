using Random, Plots, Distributions

# G = 每 Tfr 產生的數量
# S = 每 Tfr 成功的數量

# 隨機種子碼
rng = MersenneTwister()

function ALOHA(;Tfr = 1, λ = 0.005, STA_n = 5, frame_n = 10000, mod = "exp", slotted = false)
    # 隨機抽取發送每個frame的間格時間
    STA_lst = []
    for i in (1:STA_n)
        # rng = MersenneTwister(i)
        if mod == "exp"
            push!(STA_lst, cumsum(rand(rng, Exponential(1/λ), frame_n)))
        end
        if mod == "tri"
            push!(STA_lst, cumsum(rand(rng, TriangularDist(1, 200), frame_n)))
        end
        if mod == "nor"
            push!(STA_lst, cumsum(rand(rng, Normal(100, 2), frame_n)))
        end
    end

    # 整合並排序發送frame的時間
    sort_lst = []
    for i in (1:STA_n)
        for j in (1:frame_n)
            push!(sort_lst, (STA_lst[i][j] + (j-1))*Tfr)
        end
    end
    sort!(sort_lst)
    
    # 是否slotted
    slotted_ = false
    if slotted == true
        global slotted_
        slotted_ = true
        for i in (1:length(sort_lst))
            sort_lst[i] = ceil(sort_lst[i])
        end
    end
    
    # 計算總時長
    T = (maximum(sort_lst) + frame_n)*Tfr
    
    # 計算G的大小
    G = (STA_n*frame_n) / (T/Tfr)
    
    # pure偵測碰撞
    collide_frame_lst = []
    for i in (1:(length(sort_lst)-1))
        send_time = getindex(sort_lst, i)
        if i == length(sort_lst)
            break
        end
        if (sort_lst[i+1] - sort_lst[i]) < 1
            push!(collide_frame_lst, i)
            push!(collide_frame_lst, i+1)
        end
    end
    
    # 將資料型態轉成set 以此將重複的data剔除
    collide_frame_set = Set(collide_frame_lst)
    
    # 計算S的大小
    S = (STA_n*frame_n - length(collide_frame_set)) / (T/Tfr)

    P = (STA_n*frame_n - length(collide_frame_set)) / (STA_n*frame_n)
    
    return G, S, P
end

G, S = ALOHA(slotted = true)
println("G = ", G)
println("S = ", S)
println("pass = ", S/G)

res_lst = []
X_lst = []
Y_lst = []
Y_lst_ = []
# 呼叫函式
for i in (1:2:500)
    print(i)
    push!(res_lst, ALOHA(STA_n = i, mod = "exp", slotted = false))
end

# 繪圖用的list
X_lst = [0.0]
Y_lst = [0.0]
Y_lst_ = [0.0]
for i in (1:length(res_lst))
    push!(X_lst, res_lst[i][1])
    push!(Y_lst, res_lst[i][2])
    G_ = res_lst[i][1]
    if slotted_ == true
        push!(Y_lst_, G_*exp(-G_))
    else 
        push!(Y_lst_, G_*exp(-2G_))
    end
end

# 尋找S的最大值和出現的位置
println()
S_max = maximum(Y_lst)
println("S (max) = ", S_max)
S_max_index = findall(x -> x == S_max, Y_lst)
println("S (max) index = ", S_max_index)

# 繪製 gif
chart = @animate for i ∈ 1:length(res_lst)
    plot(X_lst[1:i], Y_lst[1:i], seriestype=:scatter, label="Simulated", xlabel="G", ylabel="S")
    plot!(X_lst[1:i], Y_lst_[1:i], label="Ideal", title = S_max)
end every 1

display(gif(chart, "chart.gif", fps = 30));



