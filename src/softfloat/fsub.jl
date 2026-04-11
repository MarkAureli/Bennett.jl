"""
    soft_fsub(a::UInt64, b::UInt64) -> UInt64

IEEE 754 double-precision subtraction on raw bit patterns.
Implemented as a - b = a + (-b) via soft_fadd and soft_fneg.
"""
function soft_fsub(a::UInt64, b::UInt64)::UInt64
    return soft_fadd(a, soft_fneg(b))
end
