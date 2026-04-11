# --- Operand: SSA variable or integer constant ---

struct IROperand
    kind::Symbol       # :ssa or :const
    name::Symbol       # SSA name (if :ssa)
    value::Int         # constant value (if :const)
end

ssa(name::Symbol)    = IROperand(:ssa, name, 0)
iconst(value::Int)   = IROperand(:const, Symbol(""), value)

# --- Instructions ---

abstract type IRInst end

struct IRBinOp <: IRInst
    dest::Symbol
    op::Symbol         # :add, :sub, :mul, :and, :or, :xor, :shl, :lshr, :ashr
    op1::IROperand
    op2::IROperand
    width::Int
end

struct IRICmp <: IRInst
    dest::Symbol
    predicate::Symbol  # :eq, :ne, :ult, :slt, :ugt, :sgt, :ule, :sle, :uge, :sge
    op1::IROperand
    op2::IROperand
    width::Int         # width of operands (result is always i1)
end

struct IRSelect <: IRInst
    dest::Symbol
    cond::IROperand    # i1
    op1::IROperand     # true value
    op2::IROperand     # false value
    width::Int         # width of result
end

struct IRRet <: IRInst
    op::IROperand
    width::Int
end

struct IRInsertValue <: IRInst
    dest::Symbol
    agg::IROperand       # aggregate operand (or :zero for zeroinitializer)
    val::IROperand       # value to insert
    index::Int           # 0-based element index
    elem_width::Int      # bit width of each element
    n_elems::Int         # number of elements in the aggregate
end

# --- v0.3: branch, phi, basic blocks ---

struct IRCast <: IRInst
    dest::Symbol
    op::Symbol         # :sext, :zext, :trunc
    operand::IROperand
    from_width::Int
    to_width::Int
end

struct IRPtrOffset <: IRInst
    dest::Symbol
    base::IROperand     # pointer SSA name
    offset_bytes::Int   # byte offset from base
end

struct IRVarGEP <: IRInst
    dest::Symbol
    base::IROperand     # pointer SSA name (flat wire array)
    index::IROperand    # 0-based element index (runtime SSA)
    elem_width::Int     # bit width per element
end

struct IRLoad <: IRInst
    dest::Symbol
    ptr::IROperand      # pointer (or ptr+offset SSA name)
    width::Int          # load width in bits
end

struct IRExtractValue <: IRInst
    dest::Symbol
    agg::IROperand       # aggregate operand
    index::Int           # 0-based element index
    elem_width::Int      # bit width of each element
    n_elems::Int         # number of elements in the aggregate
end

struct IRCall <: IRInst
    dest::Symbol
    callee::Function       # Julia function to compile and inline
    args::Vector{IROperand}
    arg_widths::Vector{Int}
    ret_width::Int
end

struct IRBranch <: IRInst
    cond::Union{IROperand, Nothing}      # nothing for unconditional
    true_label::Symbol
    false_label::Union{Symbol, Nothing}  # nothing for unconditional
end

struct IRSwitch <: IRInst
    cond::IROperand                                   # value being switched on
    cond_width::Int                                    # bit width of condition
    default_label::Symbol                              # default target
    cases::Vector{Tuple{IROperand, Symbol}}            # (case_val, target_label)
end

struct IRPhi <: IRInst
    dest::Symbol
    width::Int
    incoming::Vector{Tuple{IROperand, Symbol}}  # (value, from_block)
end

struct IRBasicBlock
    label::Symbol
    instructions::Vector{IRInst}  # non-terminator instructions
    terminator::IRInst             # IRBranch or IRRet
end

# --- Parsed IR bundle ---

struct ParsedIR
    ret_width::Int
    args::Vector{Tuple{Symbol, Int}}
    blocks::Vector{IRBasicBlock}
    ret_elem_widths::Vector{Int}   # [8] for i8, [8,8] for [2 x i8]
end

# Backward compat: parsed.instructions flattens all blocks
function Base.getproperty(p::ParsedIR, name::Symbol)
    if name === :instructions
        insts = IRInst[]
        for block in getfield(p, :blocks)
            append!(insts, block.instructions)
            push!(insts, block.terminator)
        end
        return insts
    else
        return getfield(p, name)
    end
end

Base.propertynames(::ParsedIR) = (:ret_width, :args, :blocks, :instructions, :ret_elem_widths)
