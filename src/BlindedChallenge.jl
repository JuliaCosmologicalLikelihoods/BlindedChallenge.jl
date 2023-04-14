module BlindedChallenge

using Artifacts
using NPZ

function __init__()

    global k_grid = npzread(joinpath(artifact"simulations_data", "k_grid.npy"))
    global covariance = npzread(joinpath(artifact"simulations_data", "covariance.npy"))
    global mono = npzread(joinpath(artifact"simulations_data", "monopole.npy"))
    global quad = npzread(joinpath(artifact"simulations_data", "quadrupole.npy"))

end

function mask_cov(cov, n)
    elements_to_remove_left = 0
    elements_to_remove_right = 100-n
    datavec_len = 100-elements_to_remove_left-elements_to_remove_right
    intermediate_Cov = zeros(300, 2*datavec_len)
    first_Cov = zeros(datavec_len*2,datavec_len*2)

    for i in 1:datavec_len
        intermediate_Cov[:,i] = cov[:, i+elements_to_remove_left]
        intermediate_Cov[:,i+datavec_len] = cov[:, i+elements_to_remove_left+100]
    end

    for i in 1:100-elements_to_remove_left-elements_to_remove_right
        first_Cov[i,:] = intermediate_Cov[i+elements_to_remove_left,:]
        first_Cov[i+100-elements_to_remove_left-elements_to_remove_right,:] = intermediate_Cov[i+elements_to_remove_left+100,:]
    end

    yerror_Mono = ones(20)
    yerror_Quad = ones(20)

    for i in 1:20
        yerror_Mono[i] = sqrt(cov[i,i])
        yerror_Quad[i] = sqrt(cov[100+i,100+i])
    end

    #let us enforce the covariance to be symmetric
    Covariance = (first_Cov .+ first_Cov')./2

    return Covariance, yerror_Mono, yerror_Quad
end

function create_data(n)

    cov, yerror_Mono, yerror_Quad = mask_cov(covariance, n)
    data = vcat(mono[1:n], quad[1:n])
    k = k_grid[1:n]

    return data, k, cov, yerror_Mono, yerror_Quad

end

end # module BlindedChallenge
