
function sphereFunc(x::AbstractArray)
    return sum(x.^2)
end

function rosenbrockfunc(x::AbstractVector)
    sum = 0.0
    for i in 1:length(x) - 1
        sum += 100*(x[i + 1] - x[i]^2)^2 + (x[i] - 1)^2
    end
    return sum
end

function rastriginfunc(x::AbstractVector)
    d = length(x)

    sum = 0.0
    for i in 1:d
        sum += (x[i]^2 - 10*cos(2*π*x[i]))
    end
    return 10*d + sum
end

function ackleyfunc(x::AbstractVector)
    a = 20
    b = 0.2
    c = 2*π
    d = length(x)

    sum1 = 0.0
    sum2 = 0.0
    for i in 1:d
        sum1 += x[i]^2
        sum2 += cos(c*x[i])
    end
    return -a*exp(-b*sqrt(sum1/d))  - exp(sum2/d) + a + exp(1)
end