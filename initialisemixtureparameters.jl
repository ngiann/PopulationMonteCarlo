function initialisemixtureparameters(D, K)

    # initialise mixture components
    α = ones(K) / K

    # initialise means
    μ = [randn(D)*10 for k in 1:K]

    # initialise covariance matrices
    Σ = [Matrix(I, D, D)*5 for k in 1:K]

    MixtureModel([MvNormal(μ[i], Σ[i]) for i in 1:K], α)

end
