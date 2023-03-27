# generate interaction
#  there are full form or some reduced form based on symmetry

"""
generate interaction, without spin flip term
interaction=gene_interaction(1.0,0.1,2)
interaction=gene_interaction(1.0,0.1,4)
# this does not assume any symmetry
"""
function gene_interaction(U,J,N_spin_orbital)
    N_orbital=trunc(Int,N_spin_orbital/2)
    interaction=[]
    # introband
    for orb in 1:N_orbital
        push!(interaction,(get_idx(orb,1),get_idx(orb,2),U))
    end
    # interband, opposite spin
    for orb1 in 1:(N_orbital-1)
        for orb2 in (orb1+1):(N_orbital)
            push!(interaction,(get_idx(orb1,1),get_idx(orb2,2),U-2*J))
            push!(interaction,(get_idx(orb1,2),get_idx(orb2,1),U-2*J))
        end
    end
    # interband, same spin
    for orb1 in 1:(N_orbital-1)
        for orb2 in (orb1+1):(N_orbital)
            push!(interaction,(get_idx(orb1,1),get_idx(orb2,1),U-3*J))
            push!(interaction,(get_idx(orb1,2),get_idx(orb2,2),U-3*J))
        end
    end
    interaction
end

"""
we can simplify the number of calculation when all orbital are same.
(We manually encode the information for each case)
Assuming N_oribtal>1
For example
U=1.0
J=0.1
gene_interaction_degenerate(U,J,N_spin_orbital)
"""
function gene_interaction_degenerate(U,J,N_spin_orbital)
    N_orbital=trunc(Int,N_spin_orbital/2)
    interaction=Vector{Any}(undef,3)
    # introband
    interaction[1]=(get_idx(1,1),get_idx(1,2),U*N_orbital)
    # interband, opposite spin
    interaction[2]=(get_idx(1,1),get_idx(2,2),(U-2*J)*(N_orbital-1)*N_orbital)
    # interband, same spin
    interaction[3]=(get_idx(1,1),get_idx(2,1),(U-3*J)*(N_orbital-1)*N_orbital)
    interaction
end

# also, there is a special case for SU(N), without J
"""
for degenerate case,
"""
function gene_interaction_degenerate_J_0(U,N_spin_orbital)
    [(1,2,U*N_spin_orbital*(N_spin_orbital-1)/2)]
end

