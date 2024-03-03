def GreedyMaxCoverage(S, F, n, k):
    # Sort the sets by their finish times
    pairs = [(S[i], F[i]) for i in range(n)]
    pairs.sort(key=lambda x: x[1])
    # Initialize the result
    result = []
    # Initialize the current index
    i = 0
    # Iterate through the sets
    while i < n:
        # Add the current set to the result
        result.append(pairs[i])
        # Find the next set that is disjoint from the current set
        while i < n and pairs[i][0] < pairs[i-1][1]:
            i += 1
    return result