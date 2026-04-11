const WireIndex = Int

abstract type ReversibleGate end

struct NOTGate <: ReversibleGate
    target::WireIndex
end

struct CNOTGate <: ReversibleGate
    control::WireIndex
    target::WireIndex
end

struct ToffoliGate <: ReversibleGate
    control1::WireIndex
    control2::WireIndex
    target::WireIndex
end

struct ReversibleCircuit
    n_wires::Int
    gates::Vector{ReversibleGate}
    input_wires::Vector{WireIndex}
    output_wires::Vector{WireIndex}
    ancilla_wires::Vector{WireIndex}
    input_widths::Vector{Int}
    output_elem_widths::Vector{Int}  # e.g. [8] for Int8, [8,8] for Tuple{Int8,Int8}
end
