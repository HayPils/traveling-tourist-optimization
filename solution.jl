using JuMP, Ipopt

A = [:hulk, :drdoom, :spiderman, :kong, :jurassic, :hp]
n = length(A)
slots = 44

d = zeros(n, slots)

w = zeros(n, slots)
T = collect(1:1:slots)

raw = readcsv("stigler.csv")
(m,n) = size(raw)

n_waitTimes = 2:n      # columns containing waitTimes
n_rides = 3:m          # rows containing rides names

waitTimes = raw[1,n_waitTimes][:]   # the list of waitTimes (convert to 1-D array)
rides = raw[n_rides,1][:]           # the list of rides (convert to 1-D array)

data = NamedArray( raw[n_rides,n_waitTimes], (rides,waitTimes), ("Rides","Wait Times") );

m = Model(solver = IpoptSolver(print_level=0))

@variable(m, x[1:n, 1:n, 1:slots], Bin)    # departures from i to j at time t
@variable(m, y[1:n, 1:slots], Bin)    # arrivals at i at time t
@constraint(m, y[1,1] = 1)    # start at the entrance of the park

@constraint(m, c1[for i in 1:n], sum(sum(x[i,j,t] for j in 1:n) for t in 1:slots) = 1)    # can only leave ride once
@constraint(m, c2[for i in 1:n], sum(sum(x[j,i,t] for j in 1:n) for t in 1:slots) = 1)    # can only come from one ride
@constraint(m, c3[for i in 2:n], sum(y[i,t] for t in 2:slots) = 1)    # can only ride each ride once
# if not arriving at time t at ride i, you cannot leave any ride at a time t - dij and expect to walk in time to arrive at time t
@constraint(m, c4[for i in 1:n, for t in 1:slots], sum(x[j,i,t-d[j,i]] for j in 1:n) <= y[i, t])
# you cannot leave i before you arrive at i
@constraint(m, c5[for i in 1:n, for t in 1:slots], sum(x[i,j,t] for j in 1:n) <= sum(y[i,k] for k in 1:t))
# difference between departure and arrival time must be at least the wait time at i, if arriving at i at time t
@constraint(m, c6[for i in 1:n, for t in 1:slots], sum(t*x[i,j,t] for j in 1:n) - t*y[i,t] >= y[i,t] * w[i,t])

@Expression(m, walk, sum(sum(sum(d[i,j] * x[i,j,t] for j in 1:j) for i in 1:n) for t in 1:slots))    # walk time
@Expression(m, wait, sum(sum(sum(t*x[i,j,t] j in 1:n) for t in 1:slots) - sum(t*y[i,t] for t in 1:slots) for i in 1:n))    # wait time
@objective(m, Min, walk + wait)    # minimize walk and wait time
solve(m)
