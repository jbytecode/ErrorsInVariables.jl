using Eive 


n = 500
x = rand(n)
otherx = randn(n)
betas = [5.0, 5]
e = randn(n)
u = randn(n)

dirtyx = x .+ u 

y = 5 .+ 5 .* x .+ 5 .* otherx .+ e 

X = hcat(ones(n), x, otherx)
println("Clean Regression:")
println(X \ y)

Xd = hcat(ones(n), dirtyx, otherx)
println("Dirty Regression:")
println(Xd \ y)

println("EIVE Regression")
result = eive(dirtyx, y, otherx)
println(result)