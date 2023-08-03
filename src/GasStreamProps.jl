"""
    GasStreamProperties{T}

- `C::SVector{6, Float64}`
- `PO::Float64`
- `mmean::Float64`
- `Ta::Float64`
"""
struct GasStreamProperties{T}
    C::SVector{6,T}           #gas concentratiosn: C = [O H N O2 N2 He]
    PO::T          #[Pa] for 300km (5×10−6;10−6 for 400km and 500km, respectively) computed with NRLMSISE-00
    mmean::T       #[g/mol] mean molecular mass
    Ta::T          #[K] ambient temperature
end

include("EnvironmentalCalcs.jl")
