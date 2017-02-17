polya_urn_model = function(G0, n, a) {
  S = c()
  
  for (i in 1:n) {
    if (runif(1) < a / (a + length(S))) {
      # Add a new value
      new = G0(n=1,sd=1)
      S = c(S, new)
    } else {
      # Pick out a value from the urn, and add back a
      # value from the same cluster
      S1 <- S[sample(1:length(S), 1)]
      S <- c(S, S1)
    }
  }
  
  S
}
