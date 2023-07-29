using JLD, DelimitedFiles

searchdir(path,key) = filter(x->occursin(key,x), readdir(path))
file_list = searchdir("/res/",".jld")

combinded = load("/res/"*file_list[1])["data"]
result1 = combinded[1]
result2 = combinded[2]
result3 = combinded[3]
result4 = combinded[4]
for i = 2:length(file_list)
    int = load("/res/"*file_list[i])["data"]
    result1 = vcat(result1,int[1])
    result2 = vcat(result2,int[2])
    result3 = vcat(result3,int[3])
    result4 = vcat(result4,int[4])
end

writedlm("result1.csv",  result1, ',')
writedlm("result2.csv",  result2, ',')
writedlm("result3.csv",  result3, ',')
writedlm("result4.csv",  result4, ',')
