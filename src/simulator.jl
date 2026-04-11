@inline apply!(b::Vector{Bool}, g::NOTGate)     = (b[g.target] ⊻= true; nothing)
@inline apply!(b::Vector{Bool}, g::CNOTGate)    = (b[g.target] ⊻= b[g.control]; nothing)
@inline apply!(b::Vector{Bool}, g::ToffoliGate) = (b[g.target] ⊻= b[g.control1] & b[g.control2]; nothing)

function simulate(circuit::ReversibleCircuit, input::Integer)
    length(circuit.input_widths) == 1 || error("simulate(circuit, input) requires single-input circuit, got $(length(circuit.input_widths)) inputs")
    return _simulate(circuit, (input,))
end

function simulate(circuit::ReversibleCircuit, inputs::Tuple{Vararg{Integer}})
    return _simulate(circuit, inputs)
end

function _simulate(circuit::ReversibleCircuit, inputs::Tuple)
    bits = zeros(Bool, circuit.n_wires)

    offset = 0
    for (k, w) in enumerate(circuit.input_widths)
        v = inputs[k]
        for i in 1:w
            bits[circuit.input_wires[offset + i]] = (v >> (i - 1)) & 1 == 1
        end
        offset += w
    end

    for gate in circuit.gates
        apply!(bits, gate)
    end

    for w in circuit.ancilla_wires
        bits[w] && error("Ancilla wire $w not zero — Bennett construction bug")
    end

    return _read_output(bits, circuit.output_wires, circuit.output_elem_widths)
end

"""
Read the output value from the simulation bit vector. Returns Int8/16/32/64
for single-element outputs, or a Tuple for multi-element (insertvalue) outputs.
Note: return type is inherently unstable (depends on circuit's output_elem_widths).
"""
function _read_output(bits, output_wires, elem_widths)
    if length(elem_widths) == 1
        return _read_int(bits, output_wires, 1, elem_widths[1])
    else
        # Multi-element return → tuple
        vals = Vector{Int64}(undef, length(elem_widths))
        off = 0
        for (k, ew) in enumerate(elem_widths)
            vals[k] = _read_int(bits, output_wires, off + 1, ew)
            off += ew
        end
        return Tuple(vals)
    end
end

function _read_int(bits, wires, start, width)
    raw = UInt64(0)
    for i in 0:width-1
        raw |= UInt64(bits[wires[start + i]]) << i
    end
    if width == 8;      reinterpret(Int8, UInt8(raw & 0xFF))
    elseif width == 16; reinterpret(Int16, UInt16(raw & 0xFFFF))
    elseif width == 32; reinterpret(Int32, UInt32(raw & 0xFFFFFFFF))
    elseif width == 64; reinterpret(Int64, raw)
    else                Int(raw & ((UInt64(1) << width) - 1))
    end
end
