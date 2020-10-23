"""
Julia implementation of Cellular Potts Model
2020-10-16~
"""
##
using Agents, AgentsPlots

mutable struct Pixel <: AbstractAgent
    id::Int
    pos::Tuple{Int,Int}
    cellid::Int
    isedge::Bool
end

mutable struct Cell
    # cellid::Int
    centerPos::Tuple{Float64,Float64}
    area::Int
    perimeter::Int
end

function checkFrag(pixel, model)::Tuple{Float64, Float64}
    # Returns the shifted x, y position of a pixel if it is far away from the cell's center pos
    c_x, c_y = model.Cells[pixel.cellid].centerPos
    p_x, p_y = pixel.pos
    
    new_px = p_x + 
    - ((p_x - c_x) > (model.griddims[1]/2)) * model.griddims[1] + 
    ((p_x - c_x) < (model.griddims[1]/2)) * model.griddims[1]

    new_py = p_y +
    - ((p_y - c_y) > (model.griddims[1]/2)) * model.griddims[2] + 
    ((p_y - c_y) < (model.griddims[1]/2)) * model.griddims[2]

    return new_px, new_py
end

function checkEdge(pixel, model)::Bool
    # Check if a pixel is edges of a cell
    neighbor_coords = node_neighbors(pixel, model)
    for (i,nc) in enumerate(neighbor_coords)
        if !(i in (1,3,6,8)) # Diagonals are not even looked at
            if model.agents[Agents.coord2vertex((nc[1],nc[2]), model)].cellid != pixel.cellid
                return true
            end
        end
    end
    return false
end


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
        if px.cellid==cellid && checkEdge(px, model)
            P += 1
        end
    end
    return P
end


function pbc(x, box_size)::Float64
    # Resets the position opposite side of the box
    return x - floor(x / box_size) * box_size
end

function update_center_pos!(target_pixel, neighbor_pixel, model)
    # Update the center position of the cells when target pixel copies neighbor pixel
    # Should consider periodic boundary condistion
    tc_x, tc_y = model.Cells[target_pixel.cellid].centerPos
    nc_x, nc_y = model.Cells[neighbor_pixel.cellid].centerPos

    t_x, t_y = checkFrag(target_pixel, model)
    n_x, n_y = checkFrag(neighbor_pixel, model)

    tc_area = model.Cells[target_pixel.cellid].area
    nc_area = model.Cells[neighbor_pixel.cellid].area

    new_tc_x = (tc_x * tc_area - t_x) / (tc_area - 1)
    new_tc_y = (tc_y * tc_area - t_y) / (tc_area - 1)

    new_nc_x = (nc_x * nc_area + n_x) / (nc_area + 1)
    new_nc_y = (nc_y * nc_area + n_y) / (nc_area + 1)

    model.Cells[target_pixel.cellid].centerPos = (new_tc_x, new_tc_y)
    model.Cells[neighbor_pixel.cellid].centerPos = (new_nc_x, new_nc_y)

end

function update_perimeter!(target_pixel, neighbor_pixel, model)
    # Target pixel must be an edge. But by copying neighbor pixel to the position of target pixel,
    # new edges can be made...
    model.Cells[target_pixel.cellid].perimeter -= 1
    model.Cells[neighbor_pixel.cellid].perimeter += 1

    target_cell_edges = 0
    neighbor_cell_edges = 0
    neighbor_coords = node_neighbors(target_pixel, model)
    for nc in neighbor_coords
        n_pixel = model.agents[Agents.coord2vertex((nc[1],nc[2]), model)]
        if n_pixel.isedge 
            if n_pixel.cellid==target_pixel.cellid
                target_cell_edges += 1
                
            elseif n_pixel.cellid==neighbor_pixel.cellid
                neighbor_cell_edges += 1
            end
        end
    end
       
    target_pixel.cellid = neighbor_pixel.cellid

    target_pixel.isedge = checkEdge(target_pixel, model)

    new_target_cell_edges = 0
    new_neighbor_cell_edges = 0
    neigbor_coords = node_neighbors(target_pixel, model)
    for nc in neighbor_coords
        n_pixel = model.agents[Agents.coord2vertex((nc[1],nc[2]), model)]
        if checkEdge(n_pixel, model)
            if n_pixel.cellid==target_pixel.cellid
                new_target_cell_edges += 1
            elseif n_pixel.cellid==neighbor_pixel.cellid
                new_neighbor_cell_edges += 1
            end
            n_pixel.isedge = true
        else
            n_pixel.isedge = false
        end
    end

    model.Cells[target_pixel.cellid].perimeter += (new_target_cell_edges - target_cell_edges)
    model.Cells[neighbor_pixel.cellid].perimeter += (new_neighbor_cell_edges - neighbor_cell_edges)

end


function get_conserv_ΔE()
    return 0
end

function get_adh_ΔE()
    return 0
end
    
function get_persis_ΔE()
    return 0
end

function get_total_ΔE(target_pixel, neighbor_pixel, model)
    ΔE = 
    get_conserv_ΔE(target_pixel, neighbor_pixel, model) + 
    get_adh_ΔE(target_pixel, neighbor_pixel, model) + 
    get_persis_ΔE(target_pixel, neighbor_pixel, model)
    return ΔE
end
    

######### Initilization ########
function initialize(; numCells=25, griddims=(50,50), Temperature=30)::AgentBasedModel
    space = GridSpace(griddims, periodic=true, moore=true)
    matCellIdx = reshape(1:numCells, Int.((√numCells, √numCells)))
    matCellIdx = repeat(matCellIdx, inner=(10,10))
    Cells = [Cell((0,0), 0, 0) for _ in 1:numCells]
    properties = Dict(:T=>Temperature, :Cells=>Cells, :griddims=>griddims, :numCells=>numCells)
    model = ABM(Pixel, space; properties = properties)

    # Initialize the Pixels
    node_idx = 1
    for i in 1:griddims[1]
        for j in 1:griddims[2]
            pixel = Pixel(node_idx, (j, i), matCellIdx[j,i], false)
            add_agent_pos!(pixel, model)
            node_idx += 1
        end
    end

    for (_, ag) in model.agents
        ag.isedge = checkEdge(ag, model)
    end

    # Initialize the cell properties
    for cellid in 1:numCells
        model.Cells[cellid].area = get_area(cellid, model)
        model.Cells[cellid].perimeter = get_perimeter(cellid, model)
        model.Cells[cellid].centerPos = get_center_pos(cellid, model)
    end

    return model
end

########## Model step! ##########
function cpm_step!(model)
    # This ONE cpm_step consists of many Monte-Carlo steps accross the whole grid
    n_pixels = length(model.agents)
    # I could use @distribute these MC-steps
    for _ in 1:n_pixels
        target_pixel = model.agents[rand(1:n_pixels)]
        if target_pixel.isedge
            ncoord = node_neighbors(target_pixel, model)[rand(1:8)]
            neighbor_pixel = model.agents[Agents.coord2vertex((ncoord[1],ncoord[2]), model)]

            if target_pixel.cellid == neighbor_pixel.cellid
                continue
            else
                # Calc Boltzmann probability
                # ΔE = get_total_ΔE(target_pixel, neighbor_pixel, model)
                ΔE = 2*rand()-1
                if ΔE < 0
                    p = 1
                else
                    p = exp(-ΔE / model.T)
                end

                # Update the pixel's cell id
                if p > rand()
                    # Update cell center pos (both target and neighbor)
                    update_center_pos!(target_pixel, neighbor_pixel, model)
                    # Increment the target cell-area & decrement neighbor cell-area
                    model.Cells[target_pixel.cellid].area -= 1
                    model.Cells[neighbor_pixel.cellid].area += 1
                    # Update the perimeters
                    update_perimeter!(target_pixel, neighbor_pixel, model)
                    # target_pixel.cellid = neighbor_pixel.cellid # This is done in `update_perimeter!`
                end
            end
        end
    end
end


##################################################################################
##
ac_edges_black(x) = x.isedge ? x.cellid : :black
ac_edges(x) = Int64(x.isedge)

model = initialize(numCells=100, griddims=(100,100))

@time step!(model, dummystep, cpm_step!, 1)

plotabm(
    model; 
    as=3.6, 
    ac=ac_edges_black, 
    am=:square, 
    size=(800,800),
    showaxis=false,
    grid=false,
    markerstrokewidth=0,
    markerstrokealpha=nothing,
    )


plotabm(
    model; 
    as=4.4, 
    ac=ac, 
    am=:square, 
    size=(500,500),
    showaxis=false,
    grid=false,
    markerstrokewidth=0,
    markerstrokealpha=nothing,
    )


