using JuMP, Ipopt


m = Model(solver = IpoptSolver(print_level=0))

@variable(m, x[1:n, 1:n, 1:a], Bin)    # departures from i to j at time t
@variable(m, y[1:n, 1:a], Bin)    # arrivals at i at time t

@constraint(m, c1[for i in 1:n], sum(sum(x[i,j,t] for j in 1:n) for t in 1:slots) = 1)    # can only leave ride once
@constraint(m, c2[for i in 1:n], sum(sum(x[j,i,t] for j in 1:n) for t in 1:slots) = 1)    # can only come from one ride
@constraint(m, c3[for i in 2:n], sum(y[i,t] for t in 2:slots) = 1)    # can only ride each ride once

# difference between departure and arrival time must be at least the wait time at i, if arriving at i at time t
@constraint(m, c5[for i in 1:n, for t in 1:slots], sum(t*x[i,j,t] for j in 1:n) - t*y[i,t] >= y[i,t] * w[i,t])

@Expression(m, walk, sum(sum(sum(d[i,j] * x[i,j,t] for j in 1:j) for i in 1:n) for t in 1:slots))    # walk time
@Expression(m, wait, sum(sum(sum(t*x[i,j,t] j in 1:n) for t in 1:slots) - sum(t*y[i,t] for t in 1:slots) for i in 1:n))    # wait time
@objective(m, Min, walk + wait)    # minimize walk and wait time
solve(m)
