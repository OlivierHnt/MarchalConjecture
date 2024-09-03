module MarchalConjecture

    using RadiiPolynomial

include("symmetry.jl")
    export EvenCos, OddCos, EvenSin, OddSin, EvenFourier, project_sym!, norm1_sym, opnorm1_sym

include("marchal.jl")
    export F_marchal!, DF_marchal!

# inverse discret Fourier transform

dft_grid2cheb(x_grid::Vector{<:Sequence}) =
    Sequence(space(x_grid[1]), [idft!(getindex.(x_grid, i)) for i ∈ indices(space(x_grid[1]))])

dft_grid2cheb(A_grid::Vector{<:LinearOperator}) =
    LinearOperator(domain(A_grid[1]), codomain(A_grid[1]),
        [idft!(getindex.(A_grid, i, j)) for i ∈ indices(codomain(A_grid[1])), j ∈ indices(domain(A_grid[1]))])

function idft!(v_grid::Vector)
    N = length(v_grid) - 1
    N_dft = 2*N
    c = zeros(eltype(v_grid), Chebyshev(N))
    for n ∈ 0:N
        c[n] = v_grid[N+1] + (-1)^n * v_grid[1]
        for k ∈ 1:N-1
            c[n] += 2 * v_grid[N+1-k] * cospi(k*n/N)
        end
        c[n] /= N_dft
    end
    c[N] /= 2
    return c
end
function idft!(v_grid::Vector{Complex{Interval{T}}}) where {T}
    N = length(v_grid) - 1
    N_dft = 2*N
    c = zeros(Complex{Interval{T}}, Chebyshev(N))
    for n ∈ 0:N
        c[n] = v_grid[N+1] + ExactReal((-1)^n) * v_grid[1]
        for k ∈ 1:N-1
            c[n] += ExactReal(2) * v_grid[N+1-k] * setprecision(BigFloat, 256) do
                cospi(interval(BigFloat, k*n // N))
            end
        end
        c[n] /= ExactReal(N_dft)
    end
    c[N] /= ExactReal(2)
    return c
end

    export dft_grid2cheb

# fast Fourier transform

function grid2cheb(x_grid::Vector{<:Sequence}, N)
    x_fft = [x_grid[end:-1:1] ; x_grid[2:end-1]]
    return Sequence(space(x_grid[1]),
        [ifft!(getindex.(x_fft, i), Chebyshev(N)) for i ∈ indices(space(x_grid[1]))])
end

function grid2cheb(A_grid::Vector{<:LinearOperator}, N)
    A_fft = [A_grid[end:-1:1] ; A_grid[2:end-1]]
    return LinearOperator(domain(A_grid[1]), codomain(A_grid[1]),
        [ifft!(getindex.(A_fft, i, j), Chebyshev(N)) for i ∈ indices(codomain(A_grid[1])), j ∈ indices(domain(A_grid[1]))])
end

function cheb2grid(x_cheb, N_fft)
    npts = N_fft ÷ 2 + 1
    x_fft = fft.(x_cheb, N_fft)
    return [getindex.(x_fft, i) for i ∈ npts:-1:1]
end

    export grid2cheb, cheb2grid

#

function get_tail(x, K)
    y = copy(x)
    y[(:,-K:K)] .= interval(0)
    return y
end

    export get_tail

end
