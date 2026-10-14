@eval VortexPasta.Reconnections begin
function _reconnect_from_cache_serial!(cache::ReconnectFastCache)
    (; crit, cl, nodes, velocities, node_prev, node_next, Ls,) = cache
    @assert crit.nthreads == 1
    reconnection_count = Ref(0)
    reconnection_length_loss = Ref(zero(number_type(nodes)))
    CellLists.foreach_pair(cl, nodes; folded = Val(true)) do x⃗, i, j
        info = should_reconnect(crit, nodes, velocities, i, j; Ls, node_prev, node_next)
        info === nothing && return
        (; is, js, length_before, length_after) = info
        i⁻, i⁺ = is
        j⁻, j⁺ = js
        reconnection_count[] += 1
        reconnection_length_loss[] += length_before - length_after
        Main.record_event(nodes, is, js, Ls, length_before, length_after)
        # Reconnect i⁻ -> j⁺ and j⁻ -> i⁺
        @inbounds begin
            node_next[i⁻] = j⁺
            node_prev[j⁺] = i⁻
            node_next[j⁻] = i⁺
            node_prev[i⁺] = j⁻
        end
    end
    (; reconnection_count = reconnection_count[], reconnection_length_loss = reconnection_length_loss[],)
end


end
