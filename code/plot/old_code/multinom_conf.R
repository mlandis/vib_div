# How many samples, N, are needed to guarantee that
#   Pr( |f_m - pi_m| < d for all m in M ) > 1 - a?

# For a time bin, you might only have K samples. Use K/N to plot
# alpha channel for state-lineage bins


# significance level
alpha = 0.05

# desired precision
d = 0.05

# observed frequencies
f = c(1,2,3,4); f = f/sum(f)

# number of categories
k = length(f)

N = 510
D = 0.05
n = c()
y = c()
for (m in 1:k) {
    z = qnorm(p=(1-alpha/(2*m)),mean=0,sd=1)
    cat("z =", z, "\n")
    n[m] = z^2*(1/m) * (1 - 1/m) / D^2
    cat("n_m = ", n[m], "\n")
    d[m] = sqrt( z^2 * (1/m) * (1 - 1/m) / N )
    cat("y_m = ", y[m], "\n")
}


# what is the value of d for the number of samples?
calculate_ske = function(s, k, alpha=0.05, D=0.05) {
    # s     : number of samples
    # k     : number of bins
    # alpha : significance level
    # D     : estimate difference must satisfy | f_i - pi_i | > D for all i in k
    # N     : min num samples needed for
    n = c()
    d = c()

    for (m in 1:k) {
        z = qnorm(p=(1-alpha/(2*m)),mean=0,sd=1)
        n[m] = ceiling(z^2*(1/m) * (1 - 1/m) / D^2)
        d[m] = sqrt( z^2 * (1/m) * (1 - 1/m) / s )
    }
    return( list(n=max(n), d=max(d), alpha=alpha, s=s, D=D) )
}