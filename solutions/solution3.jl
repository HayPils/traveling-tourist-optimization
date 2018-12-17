using JuMP, Gurobi
function import_walk(filename)        # import walk times
    raw_walk = readcsv(filename)
    r = size(raw_walk)[1] - 1 # number of rides
    return raw_walk[2:r+1, 2:r+1]
end

function import_wait(filename, interval_len)    # import wait times
    raw_wait = readcsv(filename)
    slots = size(raw_wait[1, :])[1] - 1
    w = Matrix(0,slots*interval_len)
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
    return w
end

function optimal_theme_park_tour(walk_λ, wait_λ)
    c = import_walk("walk_times.csv")    # symmetric walk time matrix
    w = import_wait("wait_times.csv", 15)     # asymmetric temporal wait time matrix (wait time at ride i at minute t)
    β = 15     # max grace period in between rides
    T = 150    # minutes in the day

    m = Model(solver=GurobiSolver(OutputFlag=0))
    @variable(m, x[1:r, 1:r], Bin)                                      # decision matrix going from ride i to j
    @variable(m, a[1:r, 1:T], Bin)                                      # decision matrix riding ride i at minute t
    # TSP constraints
    @constraint(m, c1[j in 1:r], sum( x[i,j] for i in 1:r ) == 1)       # one out-edge
    @constraint(m, c2[i in 1:r], sum( x[i,j] for j in 1:r ) == 1)       # one in-edge
    @constraint(m, c3[i in 1:r], x[i,i] == 0 )                          # no self-loops

    if (wait_λ > 0)
        #@constraint(m, c4[i in 1:r-1], sum(a[i,t]*w[i,t] for t in 1:T) + sum(x[i,j]*c[i,j] for j in 1:r) == sum(t*sum(a[j,t]*x[i,j] for j in 1:r) for t in 1:T) - sum(t*a[i,t] for t in 1:T))
        @constraint(m, c4[i in 1:r], (1 - x[i,1])*sum(a[i,t]*(t + w[i,t]) for t in 1:T) + sum(x[i,j]*c[i,j] for j in 2:r) <= sum(t*sum(a[j,t]*x[i,j] for j in 1:r) for t in 1:T))
        @constraint(m, c5[i in 1:r], sum(a[i,t]*(β + t + w[i,t]) for t in 1:T) + sum(x[i,j]*c[i,j] for j in 2:r) >= sum(t*sum(a[j,t]*x[i,j] for j in 2:r) for t in 1:T))
        @constraint(m, c6[i in 1:r], sum(a[i,t] for t in 1:T) == 1)
        @constraint(m, c7[t in 1:T], sum(a[i,t] for i in 1:r) <= 1)
        @constraint(m, a[1,1] == 1)
    end

    # Miller-Tucker-Zemlin variables and constraints to eliminate multiple subtours
    @variable(m, u[1:r])
    @constraint(m, c8[i in 1:r, j in 2:r], u[i] - u[j] + r*x[i,j] <= r-1 )

    @expression(m, walk, sum(x[i,j]*c[i,j] for i in 1:r, j in 1:r))    # walk time
    @expression(m, wait, sum(sum(a[i,t]*w[i,t] for t in 1:T) for i in 1:r))    # wait time
    @objective(m, Min, walk_λ*walk + wait_λ*wait)   # minimize total walk and wait time

    solve(m)
    println(getobjectivevalue(m))
    optx = getvalue(x)
    opta = getvalue(a)
    for i in 1:r
        print(sum(t*opta[i,t] for t in 1:T), " ")
    end
end

optimal_theme_park_tour(1,1)
