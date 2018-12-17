# import wait times
raw_walk = readcsv("walk_times.csv")
r = size(raw_walk)[1] - 1 # number of rides

rides = raw_walk[r,1][:]           
c = raw_walk[2:r+1, 2:r+1]
# import wait times
raw_wait = readcsv("wait_times.csv")
slots = size(raw_wait[1, :])[1] - 1
w = Matrix(0,slots*15) # asymmetric temporal wait time matrix (wait time at ride i at minute t)
interval_len = 15    # interval length in minutes
raw_wait = raw_wait[2:r+1,2:slots+1]
for i in 1:r
    w_row = []
    for j in 1:slots
        for k in 1:interval_len
            push!(w_row, raw_wait[i,j])
        end
    end
    w = [w ; transpose(w_row)]
end

T = 615    # minutes in the day

m = Model(solver=GurobiSolver(OutputFlag=0))
@variable(m, x[1:r, 1:r], Bin)                                      # decision matrix going from ride i to j
@variable(m, a[1:r, 1:T], Bin)
# TSP constraints
@constraint(m, c1[j in 1:r], sum( x[i,j] for i in 1:r ) == 1)       # one out-edge
@constraint(m, c2[i in 1:r], sum( x[i,j] for j in 1:r ) == 1)       # one in-edge
@constraint(m, c3[i in 1:r], x[i,i] == 0 )                          # no self-loops

@constraint(m, c4[i in 1:r, t in 1:T], a[i,t] + w[i,t] <= sum(sum(k*a[j,k]*x[i,j] for k in t:T) for j in 1:r))
@constraint(m, c5[i in 1:r], sum(a[i,t] for t in 1:T) == 1)
@constraint(m, c6[t in 1:T], sum(a[i,t] for i in 1:r) <= 1)

# Miller-Tucker-Zemlin variables and constraints to eliminate subtours
@variable(m, u[1:r])
@constraint(m, c7[i in 1:r, j in 2:r], u[i] - u[j] + r*x[i,j] <= r-1 )

@expression(m, walk, sum(x[i,j]*c[i,j] for i in 1:r,ca j in 1:r))
@expression(m, wait, sum(sum(a[i,t]*w[i,t] for t in T) for i in 1:r))
@objective(m, Min, walk + wait)   # minimize total cost

solve(m)
println(getobjectivevalue(m))
optx = getvalue(x)

opta = getvalue(a)
for i in 1:r
    print(sum(t*opta[i,t] for t in 1:T), " ")
end
