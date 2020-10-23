function get_conserv_ΔE(target_pixel, neighbor_pixel, model, newPeri_target, newPeri_neighbor)
    
    target_cell_A = model.Cells[target_pixel.cellid].area
    neighbor_cell_A = model.Cells[neighbor_pixel.cellid].area
    target_cell_P = model.Cells[target_pixel.cellid].perimeter # Current target perimeter
    neighbor_cell_P = model.Cells[neighbor_pixel.cellid].perimeter # Current neighbor perimeter

    ΔP_t = newPeri_target - target_cell_P
    ΔP_n = newPeri_neighbor - neighbor_cell_P

    ΔEₜ_area = model.K_E * (2(model.A₀ - target_cell_A) + 1)
    ΔEₙ_area = model.K_E * (2(neighbor_cell_A - model.A₀) + 1)
    ΔEₜ_peri = model.K_P * (2ΔP_t*(target_cell_P - model.P₀) + ΔP_t^2)
    ΔEₙ_peri = model.K_P * (2ΔP_n*(neighbor_cell_P - model.P₀) + ΔP_n^2)
    
    return ΔEₜ_area + ΔEₙ_area + ΔEₜ_peri + ΔEₙ_peri
end

function get_adh_ΔE(target_pixel, neighbor_pixel, model)
    return 0
end
    
function get_persis_ΔE(target_pixel, neighbor_pixel, model)
    return 0
end

function get_total_ΔE(target_pixel, neighbor_pixel, model)
    ΔE = 
    get_conserv_ΔE(target_pixel, neighbor_pixel, model) + 
    get_adh_ΔE(target_pixel, neighbor_pixel, model) + 
    get_persis_ΔE(target_pixel, neighbor_pixel, model)
    return ΔE
end