mutable struct WireAllocator
    next_wire::Int
    free_list::Vector{Int}   # sorted DESCENDING — pop! gets minimum, O(1)
end

WireAllocator() = WireAllocator(1, Int[])

function allocate!(wa::WireAllocator, n::Int)
    wires = Int[]
    for _ in 1:n
        if !isempty(wa.free_list)
            # Reuse freed wire (pop min from descending-sorted list = pop last element)
            push!(wires, pop!(wa.free_list))
        else
            push!(wires, wa.next_wire)
            wa.next_wire += 1
        end
    end
    return wires
end

"""Return wires to the allocator for reuse. Wires MUST be in zero state."""
function free!(wa::WireAllocator, wires::Vector{Int})
    for w in wires
        # Insert into descending-sorted list
        idx = searchsortedlast(wa.free_list, w; rev=true) + 1
        insert!(wa.free_list, idx, w)
    end
end

wire_count(wa::WireAllocator) = wa.next_wire - 1
