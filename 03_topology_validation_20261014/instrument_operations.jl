@eval VortexPasta.Timestepping begin
function after_advection!(iter::VortexFilamentSolver)
    (; fs, prob, nstep, refinement, fold_periodic, stats, to) = iter
    Main.audit_measure(iter,0)
    if nstep % iter.step_reconnect == 0
        # Perform reconnections, possibly changing the number of filaments.
        @timeit to "Reconnect" rec = reconnect!(iter)
        @debug lazy"Number of reconnections: $(rec.reconnection_count)"
        stats.reconnection_count += rec.reconnection_count
        stats.reconnection_length_loss += rec.reconnection_length_loss
        stats.filaments_removed_count += rec.filaments_removed_count
        stats.filaments_removed_length += rec.filaments_removed_length
        stats.reconnection_passes += rec.npasses
    end
    Main.audit_measure(iter,1)
    if nstep % iter.step_refinement == 0
        L_fold = periods(prob.p)  # box size (periodicity)
        fields = fields_to_resize(iter)
        @timeit to "Fold and refine" let
            ℓ_lost = fold_and_refine!(fs, fields; L_fold, refinement, fold_periodic)
            stats.refinement_length_loss += ℓ_lost
        end
    end
    Main.audit_measure(iter,2)
    nothing
end


end
