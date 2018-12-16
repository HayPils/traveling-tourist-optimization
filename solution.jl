using JuMP, Gurobi, NamedArrays
s = 40    # number of time s
z = 15    # interval of time s
# import wait times
raw_wait = readcsv("wait_times.csv")
(u,v) = size(raw_wait)

n_waitTimes = 2:v      # columns containing waitTimes
n_rides = 2:u          # rows containing rides names

waitTimes = raw_wait[1,n_waitTimes][:]   
rides = raw_wait[n_rides,1][:]           

w = NamedArray( raw_wait[n_rides,n_waitTimes], (rides,waitTimes), ("Rides","Wait Times") );
r = u - 1    # number of rides

m = Model(solver=GurobiSolver(OutputFlag=0))

@variable(m, x[1:r, 1:r, 1:s], Bin)    # departures from i to j at time t
@variable(m, y[1:r, 1:s], Bin)    # arrivals at i at time t

@constraint(m, c1[i in 1:r], sum(sum(x[i,j,t] for j in 1:r) for t in 1:s) == 1)    # can only leave ride once
@constraint(m, c2[i in 1:r], sum(sum(x[j,i,t] for j in 1:r) for t in 1:s) == 1)    # can only come from one ride
@constraint(m, c3[i in 1:r], sum(y[i,t] for t in 1:s) == 1)    # ride each ride once and only once
@constraint(m, c4[t in 1:s], sum(y[i,t] for i in 1:r) <= 1)   # can only arrive at one ride at one time
# you cannot leave i before you arrive at i
@constraint(m, c5[i in 1:r, t in 1:s], sum(x[i,j,t] for j in 1:r) <= sum(y[i,k] for k in 1:t))
# difference between departure and arrival time must be at least the wait time at i, if arriving at i at time t
@constraint(m, c6[i in 1:r, t in 1:s], z*(sum(sum(k*x[i,j,k] for j in 1:r) for k in 1:s) - t*y[i,t]) >= y[i,t] * w[i,t])

# wait time
@expression(m, wait, sum(sum(sum(z*t*x[i,j,t] for j in 1:r) for t in 1:s) - sum(z*t*y[i,t] for t in 1:s) for i in 1:r))
@objective(m, Min, wait)    # minimize wait time
solve(m)
println("Min wait time: ", getobjectivevalue(m))
start_time = waitTimes[1]
end_time = waitTimes[n-1]
sched_data = Array{String}(r, 2)
sched_headers = ["Arrival Time", "Wait (min)"]

opty = getvalue(y)
optx = getvalue(x)
order = 1

ride_order = copy(rides)
for i in 1:s
    for j in 1:r
        if (opty[j,i] == 1.0)
            ride_order[order] = rides[j]
            sched_data[order,1] = waitTimes[i]
            sched_data[order,2] = string(w[j,i])
            order += 1
        end
    end
end

schedule = NamedArray( sched_data[1:r,1:2], (ride_order, sched_headers), ("Rides","") );
println(schedule)
