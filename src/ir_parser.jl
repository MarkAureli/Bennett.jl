# ---------- helpers ----------

function parse_ssa_name(s::AbstractString)
    s = strip(s)
    if startswith(s, "%\"") && endswith(s, '"')
        return Symbol(s[3:end-1])
    elseif startswith(s, '%')
        return Symbol(s[2:end])
    else
        error("Not an SSA name: $s")
    end
end

function parse_operand(s::AbstractString)
    s = strip(s)
    startswith(s, '%') && return ssa(parse_ssa_name(s))
    s == "true"  && return iconst(1)
    s == "false" && return iconst(0)
    return iconst(parse(Int, s))
end

# ---------- function header ----------

const RE_FUNCDEF = r"^define\s+i(\d+)\s+@[\w.]+\((.*?)\)"
const RE_ARG     = r"i(\d+)\s+(?:(?:signext|zeroext)\s+)?(%(?:\"[^\"]+\"|[\w.]+))"

function parse_function_header(line::AbstractString)
    m = match(RE_FUNCDEF, line)
    m === nothing && error("Cannot parse function header: $line")
    ret_width = parse(Int, m[1])
    args = Tuple{Symbol,Int}[]
    for am in eachmatch(RE_ARG, m[2])
        push!(args, (parse_ssa_name(am[2]), parse(Int, am[1])))
    end
    return (ret_width, args)
end

# ---------- instruction line ----------

const RE_DEST    = r"^(%(?:\"[^\"]+\"|[\w.]+))\s*=\s*(.+)$"
const RE_BINOP   = r"^(add|sub|mul|and|or|xor|shl|lshr|ashr)(?:\s+(?:nsw|nuw|exact))*\s+i(\d+)\s+(.+?),\s*(.+)$"
const RE_ICMP    = r"^icmp\s+(\w+)\s+i(\d+)\s+(.+?),\s*(.+)$"
const RE_SELECT  = r"^select\s+i1\s+(.+?),\s*i(\d+)\s+(.+?),\s*i(\d+)\s+(.+)$"
const RE_RET     = r"^ret\s+i(\d+)\s+(.+)$"
const RE_BR_COND = r"^br\s+i1\s+(.+?),\s*label\s+(%(?:\"[^\"]+\"|[\w.]+)),\s*label\s+(%(?:\"[^\"]+\"|[\w.]+))$"
const RE_BR_UNC  = r"^br\s+label\s+(%(?:\"[^\"]+\"|[\w.]+))$"
const RE_CAST    = r"^(sext|zext|trunc)\s+i(\d+)\s+(.+?)\s+to\s+i(\d+)$"
const RE_PHI     = r"^phi\s+i(\d+)\s+(.+)$"
const RE_PHI_INC = r"\[\s*(.+?),\s*(%(?:\"[^\"]+\"|[\w.]+))\s*\]"
const RE_LABEL   = r"^([\w.]+):\s*(?:;.*)?$"

function parse_instruction(line::AbstractString)
    line = strip(line)

    # skip
    (isempty(line) || startswith(line, ';') || startswith(line, '!') ||
     line == "}" || startswith(line, "define") || startswith(line, "attributes") ||
     startswith(line, "source_filename") || startswith(line, "target ")) &&
        return nothing

    # br (conditional)
    m = match(RE_BR_COND, line)
    if m !== nothing
        return IRBranch(parse_operand(m[1]), parse_ssa_name(m[2]), parse_ssa_name(m[3]))
    end

    # br (unconditional)
    m = match(RE_BR_UNC, line)
    if m !== nothing
        return IRBranch(nothing, parse_ssa_name(m[1]), nothing)
    end

    # ret
    m = match(RE_RET, line)
    if m !== nothing
        return IRRet(parse_operand(m[2]), parse(Int, m[1]))
    end

    # assignment: %dest = rhs
    m = match(RE_DEST, line)
    m === nothing && return nothing
    dest = parse_ssa_name(m[1])
    rhs  = strip(m[2])

    # cast (sext, zext, trunc)
    m2 = match(RE_CAST, rhs)
    if m2 !== nothing
        return IRCast(dest, Symbol(m2[1]), parse_operand(m2[3]),
                      parse(Int, m2[2]), parse(Int, m2[4]))
    end

    # phi
    m2 = match(RE_PHI, rhs)
    if m2 !== nothing
        width = parse(Int, m2[1])
        incoming = Tuple{IROperand, Symbol}[]
        for pm in eachmatch(RE_PHI_INC, m2[2])
            push!(incoming, (parse_operand(pm[1]), parse_ssa_name(pm[2])))
        end
        return IRPhi(dest, width, incoming)
    end

    # icmp
    m2 = match(RE_ICMP, rhs)
    if m2 !== nothing
        return IRICmp(dest, Symbol(m2[1]), parse_operand(m2[3]),
                      parse_operand(m2[4]), parse(Int, m2[2]))
    end

    # select
    m2 = match(RE_SELECT, rhs)
    if m2 !== nothing
        return IRSelect(dest, parse_operand(m2[1]),
                        parse_operand(m2[3]), parse_operand(m2[5]),
                        parse(Int, m2[2]))
    end

    # binary op
    m2 = match(RE_BINOP, rhs)
    if m2 !== nothing
        return IRBinOp(dest, Symbol(m2[1]), parse_operand(m2[3]),
                       parse_operand(m2[4]), parse(Int, m2[2]))
    end

    error("Unrecognised LLVM IR instruction: $line")
end

# ---------- full IR string → block-structured ParsedIR ----------

function parse_ir(ir::AbstractString)
    lines = split(ir, '\n')
    ret_width = 0
    args = Tuple{Symbol,Int}[]
    blocks = IRBasicBlock[]

    current_label = :entry
    current_insts = IRInst[]

    for line in lines
        sline = strip(line)

        if startswith(sline, "define")
            (ret_width, args) = parse_function_header(sline)
            continue
        end

        # block label
        m = match(RE_LABEL, sline)
        if m !== nothing
            current_label = Symbol(m[1])
            current_insts = IRInst[]
            continue
        end

        inst = parse_instruction(sline)
        inst === nothing && continue

        if inst isa IRRet || inst isa IRBranch
            push!(blocks, IRBasicBlock(current_label, current_insts, inst))
            current_label = Symbol("")
            current_insts = IRInst[]
        else
            push!(current_insts, inst)
        end
    end

    return ParsedIR(ret_width, args, blocks, [ret_width])
end
