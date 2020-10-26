
function get_conserv_ΔE(target_pixel, neighbor_pixel, model)::Float64
    
    target_cell_A = model.Cells[target_pixel.cellid].area
    neighbor_cell_A = model.Cells[neighbor_pixel.cellid].area
    target_cell_P = model.Cells[target_pixel.cellid].perimeter # Current target perimeter
    neighbor_cell_P = model.Cells[neighbor_pixel.cellid].perimeter # Current neighbor perimeter

    ΔP_t, ΔP_n = get_perimeter_change(target_pixel, neighbor_pixel, model)

    ΔEₜ_area = model.K_E * (2*(model.A₀ - target_cell_A) + 1)
    ΔEₙ_area = model.K_E * (2*(neighbor_cell_A - model.A₀) + 1)
    ΔEₜ_peri = model.K_P * (2*ΔP_t*(target_cell_P - model.P₀) + ΔP_t^2)
    ΔEₙ_peri = model.K_P * (2*ΔP_n*(neighbor_cell_P - model.P₀) + ΔP_n^2)
    
    return ΔEₜ_area + ΔEₙ_area + ΔEₜ_peri + ΔEₙ_peri
end


function get_adh_ΔE(target_pixel, neighbor_pixel, model)
    N_old = 0
    N_new = 0
    for nc in node_neighbors(target_pixel, model)
        npx = model.agents[Agents.coord2vertex((nc[1],nc[2]), model)]
    # for npx in pixel_neighbors(target_pixel, model)
        if npx.cellid==target_pixel.cellid
            N_new += 1
        elseif npx.cellid==neighbor_pixel.cellid
            N_old += 1
        end
    end
    return model.J_CC*(N_new - N_old)
end
    

function get_persis_ΔE(target_pixel, neighbor_pixel, model)
    target_cell = model.Cells[target_pixel.cellid]
    neighbor_cell = model.Cells[neighbor_pixel.cellid]

    Δr_target = (target_cell.centerPos .- target_pixel.pos) ./ (target_cell.area - 1)
    Δr_target = Δr_target ./ norm(Δr_target)
    Δr_neighbor = (neighbor_pixel.pos .- neighbor_cell.centerPos) ./ (neighbor_cell.area + 1)
    Δr_neighbor = Δr_neighbor ./ norm(Δr_neighbor)

    ΔE = dot(target_cell.polarity, Δr_target) / norm(target_cell.polarity) +
    dot(neighbor_cell.polarity, Δr_neighbor) / norm(neighbor_cell.polarity)

    return - model.S * ΔE
end


function get_total_ΔE(target_pixel, neighbor_pixel, model)
    ΔE = 
    get_conserv_ΔE(target_pixel, neighbor_pixel, model) + 
    get_adh_ΔE(target_pixel, neighbor_pixel, model) + 
    get_persis_ΔE(target_pixel, neighbor_pixel, model)
    return ΔE
end