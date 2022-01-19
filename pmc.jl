function pmc(D, logp; K=3, maxiter=1, numsamples=100)

    mixture = initialisemixtureparameters(D, K)

    for iteration in 1:maxiter

        # draw samples from mixture of size D × numsamples

        X = drawsamples(mixture, numsamples)

        # compute important weights

        ω = zeros(numsamples)

        for n in 1:numsamples
            ω[n] = exp(logp(X[:,n]) - logpdf(mixture, X[:,n]))
        end

        # normalise weights

        w = ω / sum(ω)

        # compute log-responsibilities

        logρ = zeros(numsamples, K)

        for k in 1:K
            logρ[:, k] = logpdf(mixture.components[k], X)
        end


        # compute responsibilities
        ρ = zeros(numsamples, K)

        for n in 1:numsamples, k in 1:K
            ρ[n,k] = exp(logρ[n,k] - logsumexp(logρ[n,:]))
        end

        @show sum(ρ, dims=1)

        # calculate new mixture coefficients
        α = zeros(K)

        for k in 1:K
            for n in 1:numsamples
                α[k] += w[n] * ρ[n,k]
            end
        end

        @show α

        # calculate new mixture means
        μ = [zeros(D) for k in 1:K]

        for k in 1:K
            for n in 1:numsamples
                μ[k] += w[n] * ρ[n,k] * X[:,n] / α[k]
            end
        end

        # calculate new mixture covariances
        Σ = [zeros(D, D) for k in 1:K]

        for k in 1:K

            for n in 1:numsamples
                Σ[k] += w[n] * ρ[n,k] * (μ[k]-X[:,n])*(μ[k]-X[:,n])' / α[k]
            end

            Σ[k] = (Σ[k] + Σ[k]')/2 + 1e-1I

        end

        # update mixture
        mixture = MixtureModel([MvNormal(μ[i], Σ[i]) for i in 1:K], α)


    end


    # Plot result
    # Plot target density
    if D == 2
        figure(1)
        cla()
        θval = collect(LinRange(-6.0, 6.0, 800))
        title("target unnormalised posterior")
        contourf( repeat(θval,1,length(θval)),  repeat(θval',length(θval),1), map(θ -> exp(logp(θ)), [[x;y] for x in θval, y in θval]), cmap=plt.cm.binary)
        ax = axis()
    end

    for comp in mixture.components
        plot_ellipse(comp.μ, comp.Σ)
    end

    return mixture

end
