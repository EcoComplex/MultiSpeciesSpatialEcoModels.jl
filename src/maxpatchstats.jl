

"""
    function max_patch_highlight(m)

Takes a landscape `m` and highlights the maximum patch, 
assigns it with value 2 and all the other patches with value 1
and return a new landscape

"""
function max_patch_highlight(m)
    labels = label_components(m)            
    patch_size = component_lengths(labels)
    lmax = findmax(patch_size[2:end])[2]
    @inbounds for t in eachindex(labels)
        if labels[t] == lmax
            labels[t] = 2
        elseif labels[t]>0
            labels[t] = 1
        end
    end
    return labels
end

"""
    function max_patch_highlight(m)

Takes a landscape `m` and returns the size of the maximum patch

"""
function max_patch_size(m)
    labels = label_components(m)            
    patch_size = component_lengths(labels)
    return maximum(patch_size[2:end])
end
