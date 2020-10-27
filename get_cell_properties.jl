"""
Calculate or update cell-related properties
"""
##
function get_center_pos(cellid, model)::Tuple{Float64,Float64}
    pos = [0.0, 0.0]
    count = 0
    
    @inbounds for (_,px) in model.agents
        if px.cellid==cellid
            @. pos .+= px.pos
            count += 1
        end
    end
    
    return pos[1]/count, pos[2]/count
end


function get_area(cellid, model)::Int64
    A = 0
    @inbounds for (_,ag) in model.agents
        if ag.cellid==cellid
            A += 1
        end
    end
    return A
end


function get_perimeter(cellid, model)::Int64
    P = 0
    @inbounds for (_,px) in model.agents
        if px.cellid==cellid && (checkNeighbors(px,model)!=0)
            P += 1
        end
    end
    return P
end


function checkNeighbors(pixel, model)::Int
    # Checks if a pixel is an edge
    # Returns the number of neighboring pixels of different cellid (non-zero if edge)
    N = 0
    @inbounds for nc in node_neighbors(pixel, model)
        npx = model.agents[Agents.coord2vertex((nc[1],nc[2]), model)]
    # for npx in pixel_neighbors(pixel, model)
        if npx.cellid != pixel.cellid
            N += 1
        end
    end
    return N
end


function get_perimeter_change(target_pixel, neighbor_pixel, model)::Tuple{Int,Int}
    ΔPₜ = -1 # Currently, center should be an edge. By copying neighbor_pixel to the center, edge count decreses by one.
    ΔPₙ = 0
    N = 0 # a counter
    @inbounds for nc in node_neighbors(target_pixel, model)
        npx = model.agents[Agents.coord2vertex((nc[1],nc[2]), model)]
    # for px in pixel_neighbors(target_pixel, model)
        if npx.cellid==target_pixel.cellid
            if npx.neighbors==0
                ΔPₜ += 1
            end
        elseif npx.cellid==neighbor_pixel.cellid
            if npx.neighbors==1
                ΔPₙ -= 1 # This npx will no longer be an edge, thus decrease by 1
            end
            N += 1
        end
    end
    if N < 8 # At least one of the neighbor pixels is not neighbor_pixel.cellid, meaning after copying, the center is an edge.
        ΔPₙ += 1
    end
    return ΔPₜ, ΔPₙ
end


function update_neighbors!(target_pixel, neighbor_pixel, model)
    # Update the Pixels.neighbors of all of the neighboring pixels
    @inbounds for nc in node_neighbors(target_pixel, model)
        npx = model.agents[Agents.coord2vertex((nc[1],nc[2]), model)]
    # for npx in pixel_neighbors(target_pixel, model)
        if npx.cellid==neighbor_pixel.cellid
            npx.neighbors -= 1
            target_pixel.neighbors -= 1
        elseif npx.cellid==target_pixel.cellid
            npx.neighbors += 1
            target_pixel.neighbors += 1
        end
    end
end


function pbc(pos, box_size)::Tuple{Float64,Float64}
    # Resets the position opposite side of the box (due to periodic bc)
    return @. pos - floor(pos / box_size) * box_size
end


function checkFrag(pixel, model)::Tuple{Float64, Float64}
    # Checks if a cell is fragmented due to periodic bc
    # Returns adjusted pixels position close to cell center

    # c_pos = model.Cells[pixel.cellid].centerPos
    new_p_pos = @. pixel.pos -
    ((pixel.pos - model.Cells[pixel.cellid].centerPos) > model.griddims/2) * model.griddims +
    ((model.Cells[pixel.cellid].centerPos - pixel.pos) > model.griddims/2) * model.griddims

    return new_p_pos
end


function get_new_pos(center_pos, pixel_pos, area, a)::Tuple{Float64,Float64}
    return @. (center_pos * area + a*pixel_pos) / (area + a)
end


function update_center_pos!(target_pixel, neighbor_pixel, model)
    # For some reason... this does not update correctly
    t_cell = model.Cells[target_pixel.cellid]
    n_cell = model.Cells[neighbor_pixel.cellid]

    # t_center_pos = get_new_pos(t_cell.centerPos, checkFrag(target_pixel,model), t_cell.area, -1)
    # n_center_pos = get_new_pos(n_cell.centerPos, checkFrag(neighbor_pixel,model), n_cell.area, 1)
    t_px_pos = checkFrag(target_pixel,model)
    n_px_pos = checkFrag(neighbor_pixel,model)

    t_center_pos = @. (t_cell.centerPos * t_cell.area - t_px_pos) / (t_cell.area - 1)
    n_center_pos = @. (n_cell.centerPos * n_cell.area + n_px_pos) / (n_cell.area + 1)

    model.Cells[target_pixel.cellid].centerPos = pbc(t_center_pos, model.griddims)
    model.Cells[neighbor_pixel.cellid].centerPos = pbc(n_center_pos, model.griddims)

    # model.Cells[target_cellid].centerPos = pbc(get_center_pos(target_cellid, model), model.griddims)
    # model.Cells[neighbor_cellid].centerPos = pbc(get_center_pos(neighbor_cellid, model), model.griddims)

end


function update_polarity!(model)
    # polarity is a property which depends upon the history of the cell's motion
    # polarity updates only once every 1 MCS step
    @inbounds for cell in model.Cells
        Δpos = cell.centerPos .- cell.prev_pos
        Δpos = @. Δpos - 
        (Δpos > model.griddims/2)*model.griddims +
        (-Δpos > model.griddims/2)*model.griddims 

        cell.polarity = @. cell.polarity + Δpos - cell.polarity / model.τ
        cell.prev_pos = cell.centerPos
    end
end


