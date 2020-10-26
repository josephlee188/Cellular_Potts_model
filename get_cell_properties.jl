"""
Calculate or update cell-related properties
"""
##
function get_center_pos(cellid, model)::Tuple{Float64,Float64}
    x = 0
    y = 0
    count = 0
    @inbounds for (_,px) in model.agents
        if px.cellid==cellid
            x += px.pos[1]
            y += px.pos[2]
            count += 1
        end
    end

    return x/count, y/count
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


function get_perimeter_change(target_pixel, neighbor_pixel, model)::Tuple{Int,Int}
    ΔPₜ = -1 # Currently, center should be an edge. By copying neighbor_pixel to the center, edge count decreses by one.
    ΔPₙ = 0
    N = 0 # A counter
    for nc in node_neighbors(target_pixel, model)
        npx = model.agents[Agents.coord2vertex((nc[1],nc[2]), model)]
    # for px in pixel_neighbors(target_pixel, model)
        if npx.cellid==target_pixel.cellid
            if npx.neighbors==0
                ΔPₜ += 1
            end
        elseif npx.cellid==neighbor_pixel.cellid
            if npx.neighbors==1
                ΔPₙ -= 1
                N += 1
            end
        end
    end
    if N < 8 # At least one of the neighbor pixels is not neighbor_pixel.cellid, meaning after copying, the center is an edge.
        ΔPₙ += 1
    end
    return ΔPₜ, ΔPₙ
end


function pbc(x, box_size)::Float64
    # Resets the position opposite side of the box (due to periodic bc)
    return x - floor(x / box_size) * box_size
end


function update_center_pos!(target_pixel, neighbor_pixel, model)

    t_cell = model.Cells[target_pixel.cellid]
    n_cell = model.Cells[neighbor_pixel.cellid]

    t_x, t_y = checkFrag(target_pixel, model)
    n_x, n_y = checkFrag(neighbor_pixel, model)

    new_tc_x = (t_cell.centerPos[1] * t_cell.area - t_x) / (t_cell.area - 1)
    new_tc_y = (t_cell.centerPos[2] * t_cell.area - t_y) / (t_cell.area - 1)
    new_nc_x = (n_cell.centerPos[1] * n_cell.area + n_x) / (n_cell.area + 1)
    new_nc_y = (n_cell.centerPos[2] * n_cell.area + n_y) / (n_cell.area + 1)

    new_tc_x = pbc(new_tc_x, model.griddims[1])
    new_tc_y = pbc(new_tc_y, model.griddims[2])
    new_nc_x = pbc(new_nc_x, model.griddims[1])
    new_nc_y = pbc(new_nc_y, model.griddims[2])

    model.Cells[target_pixel.cellid].centerPos = (new_tc_x, new_tc_y)
    model.Cells[neighbor_pixel.cellid].centerPos = (new_nc_x, new_nc_y)

end


function update_polarity!(model)
    # polarity is a property which depends upon the history of the cell's motion
    # polarity updates only once every 1 MCS step
    for cell in model.Cells
        @. cell.polarity += (cell.centerPos - cell.prev_pos) - cell.polarity / model.τ
        @. cell.prev_pos = cell.centerPos
    end
end


