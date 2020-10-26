"""
Julia implementation of Cellular Potts Model
2020-10-16~
"""
##
using Agents, AgentsPlots
using LinearAlgebra
include("./get_energy_diff.jl")
include("./get_cell_properties.jl")
using BenchmarkTools

mutable struct Pixel <: AbstractAgent
    id::Int
    pos::Tuple{Int,Int}
    cellid::Int
    # isedge::Bool
    neighbors::Int # Number of distinct neighbors (only at sides and up down). 
    # It is an edge if this number != 0
end

mutable struct Cell
    centerPos::Tuple{Float64,Float64}
    area::Int
    perimeter::Int
    polarity::Tuple{Float64,Float64}
    prev_pos::Tuple{Float64,Float64}
end

function checkFrag(pixel, model)::Tuple{Float64, Float64}
    # Checks if a cell is fragmented due to periodic bc
    # Returns adjusted pixels position close to cell center

    c_x, c_y = model.Cells[pixel.cellid].centerPos
    p_x, p_y = pixel.pos
    
    new_px = p_x + 
    - ((p_x - c_x) > (model.griddims[1]/2)) * model.griddims[1] + 
    ((c_x - p_x) > (model.griddims[1]/2)) * model.griddims[1]

    new_py = p_y +
    - ((p_y - c_y) > (model.griddims[1]/2)) * model.griddims[2] + 
    ((c_y - p_y) > (model.griddims[1]/2)) * model.griddims[2]

    return new_px, new_py
end


"""
# Using this function (returns neighbor Pixels) slows down the code and doubles the allocations.

function pixel_neighbors(pixel, model)::Array{Pixel}
    # Returns neighboring Pixel (total 8)
    pixels = Array{Pixel}(undef, 8)
    neighbor_coords = node_neighbors(pixel, model)
    for (i,nc) in enumerate(neighbor_coords)
        pixels[i] = model.agents[Agents.coord2vertex((nc[1],nc[2]), model)]
    end
    return pixels
end
"""


function checkNeighbors(pixel, model)::Int
    # Checks if a pixel is an edge
    # Returns the number of neighboring pixels of different cellid (non-zero if edge)
    N = 0
    for nc in node_neighbors(pixel, model)
        npx = model.agents[Agents.coord2vertex((nc[1],nc[2]), model)]
    # for npx in pixel_neighbors(pixel, model)
        if npx.cellid != pixel.cellid
            N += 1
        end
    end
    return N
end


function update_neighbors!(target_pixel, neighbor_pixel, model)
    # Update the Pixels.neighbors of all of the neighboring pixels
    for nc in node_neighbors(target_pixel, model)
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


######### Initilization ########
function initialize(; 
    numCells=25, 
    griddims=(50,50), 
    Temperature=30, 
    A₀=100, 
    P₀=32,
    K_E=1,
    K_P=1,
    J_CC=-10,
    S=0.1,
    τ=5)::AgentBasedModel
    
    
    space = GridSpace(griddims, periodic=true, moore=true)
    matCellIdx = reshape(1:numCells, Int.((√numCells, √numCells)))
    matCellIdx = repeat(matCellIdx, inner=(10,10))
    
    Cells = [Cell((0,0), 
    0, 
    0, 
    (0,0), 
    (0,0)) for _ in 1:numCells]

    properties = Dict(
        :T=>Temperature, 
        :Cells=>Cells, 
        :griddims=>griddims, 
        :numCells=>numCells, 
        :A₀=>A₀, 
        :P₀=>P₀,
        :K_E=>K_E,
        :K_P=>K_P,
        :J_CC=>J_CC,
        :S=>S,
        :τ=>τ,
        )

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
        ag.neighbors = checkNeighbors(ag, model)
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
        if target_pixel.neighbors>0 # This is an edge
            ncoord = node_neighbors(target_pixel, model)[rand(1:8)]
            neighbor_pixel = model.agents[Agents.coord2vertex((ncoord[1],ncoord[2]), model)]
            # neighbor_pixel = pixel_neighbors(target_pixel, model)[rand(1:8)] # This is a copy (not a reference)
            if target_pixel.cellid == neighbor_pixel.cellid
                continue
            else
                # Calc Boltzmann probability
                ΔE = get_total_ΔE(target_pixel, neighbor_pixel, model)
                # ΔE = 2*rand()-1
                if ΔE < 0
                    p = 1
                else
                    p = exp(-ΔE / model.T)
                end

                # Update the pixel's cell id
                if p > rand()
                    # Increment the target cell-area & decrement neighbor cell-area
                    model.Cells[target_pixel.cellid].area -= 1
                    model.Cells[neighbor_pixel.cellid].area += 1
                    # Update the cell perimeters
                    ΔP_t, ΔP_n = get_perimeter_change(target_pixel, neighbor_pixel, model)
                    model.Cells[target_pixel.cellid].perimeter += ΔP_t
                    model.Cells[neighbor_pixel.cellid].perimeter += ΔP_n
                    # Update the pixel properties
                    update_center_pos!(target_pixel, neighbor_pixel, model)
                    update_neighbors!(target_pixel, neighbor_pixel, model)
                    target_pixel.cellid = neighbor_pixel.cellid
                end
            end
        end
    end
    update_polarity!(model)
end


##################################################################################
##
model = initialize(numCells=100, griddims=(100,100), A₀=100, P₀=31, J_CC=-10);
##
# @benchmark step!(model, dummystep, cpm_step!, 1)
step!(model, dummystep, cpm_step!, 1)
##
plotabm(
    model; 
    # as=1.5, 
    as=1.65,
    # ac=ac_edges_black, 
    ac = ac_cellidx,
    am=:square, 
    size=(400,400),
    showaxis=false,
    grid=false,
    markerstrokewidth=0,
    markerstrokealpha=nothing,
    )

Cells = model.Cells
centerX = [cell.centerPos[1] for cell in Cells]
centerY = [cell.centerPos[2] for cell in Cells]
polX = [cell.polarity[1] for cell in Cells]
polY = [cell.polarity[2] for cell in Cells]
col = [1:225]
scatter!(centerX,centerY, markersize=3, color=:black, leg=false)
quiver!(centerX,centerY, quiver=(polX, polY), color=:black, leg=false)


##
# plotabm(
#     model; 
#     as=4.4, 
#     ac=ac, 
#     am=:square, 
#     size=(500,500),
#     showaxis=false,
#     grid=false,
#     markerstrokewidth=0,
#     markerstrokealpha=nothing,
#     )


ac_edges_black(x) = (x.neighbors>0) ? x.cellid : :black
ac_cellidx(x) = Int64(x.cellid)
    


## Test
target_pixel = model.agents[1]
neighbor_pixel = model.agents[10000]
get_total_ΔE(target_pixel, neighbor_pixel, model)
