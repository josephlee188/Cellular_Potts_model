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
    @inbounds for _ in 1:n_pixels
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
                if ΔE < 0
                    p = 1
                else
                    p = exp(-ΔE / model.T)
                end

                # Update the system
                if p > rand()
                    # Update the cell perimeters
                    ΔP_t, ΔP_n = get_perimeter_change(target_pixel, neighbor_pixel, model)
                    model.Cells[target_pixel.cellid].perimeter += ΔP_t
                    model.Cells[neighbor_pixel.cellid].perimeter += ΔP_n
                    
                    # Update the pixel properties
                    update_neighbors!(target_pixel, neighbor_pixel, model)
                    update_center_pos!(target_pixel, neighbor_pixel, model)
                    
                    # Increment the target cell-area & decrement neighbor cell-area
                    model.Cells[target_pixel.cellid].area -= 1
                    model.Cells[neighbor_pixel.cellid].area += 1
                    
                    # target_cellid = target_pixel.cellid
                    # neighbor_cellid = neighbor_pixel.cellid

                    target_pixel.cellid = neighbor_pixel.cellid
                end
            end
        end
    end
    update_polarity!(model)
end


#####################################################################################################################
##
model = initialize(
    numCells=100, 
    griddims=(100,100), 
    Temperature=30, 
    A₀=100, 
    P₀=32,
    K_E=1,
    K_P=1,
    J_CC=-10,
    S=20,
    τ=50,
    );
##
# @benchmark step!(model, dummystep, cpm_step!, 1)
@time step!(model, dummystep, cpm_step!, 100)

##
ac_edges_black(x) = (x.neighbors>0) ? x.cellid : :black
ac_cellidx(x) = Int64(x.cellid)
plotabm(
    model; 
    # as=1.5, 
    as=2.7,
    # ac=ac_edges_black, 
    ac = ac_cellidx,
    am=:square, 
    size=(600,600),
    showaxis=false,
    grid=false,
    markerstrokewidth=0,
    markerstrokealpha=nothing,
    )

Cells = model.Cells
centerX = [cell.centerPos[1] for cell in Cells]
centerY = [cell.centerPos[2] for cell in Cells]
# vX = [cell.centerPos[1]-cell.prev_pos[1] for cell in Cells]
# vY = [cell.centerPos[2]-cell.prev_pos[2] for cell in Cells]
# v = [norm((vX[i],vY[i])) for i in 1:length(Cells)]
# vX ./= v
# vY ./= v
polX = [cell.polarity[1] for cell in Cells]
polY = [cell.polarity[2] for cell in Cells]

cenX = [get_center_pos(i,model)[1] for i in 1:length(Cells)]
cenY = [get_center_pos(i,model)[2] for i in 1:length(Cells)]

scatter!(centerX,centerY, markersize=3, color=:black, leg=false)
scatter!(cenX,cenY, markersize=3, color=:red, leg=false)
# quiver!(centerX,centerY, quiver=(polX, polY).*0.1, color=:black, leg=false)

## Test
target_pixel = model.agents[1]
neighbor_pixel = model.agents[10000]
get_total_ΔE(target_pixel, neighbor_pixel, model)

## ################ Generate gif ################
model = initialize(numCells=100, griddims=(100,100), A₀=100, P₀=31, J_CC=-10, S=10, τ=20);
Time = 1000
anim = @time @animate for i ∈ 1:Time
    plotabm(
    model; 
    # as=1.5, 
    as=2.7,
    # ac=ac_edges_black, 
    ac = ac_cellidx,
    am=:square, 
    size=(600,600),
    showaxis=false,
    grid=false,
    markerstrokewidth=0,
    markerstrokealpha=nothing,
    )
    centerX = [cell.centerPos[1] for cell in model.Cells]
    centerY = [cell.centerPos[2] for cell in model.Cells]
    polX = [cell.polarity[1] for cell in model.Cells]
    polY = [cell.polarity[2] for cell in model.Cells]
    cenX = [get_center_pos(i,model)[1] for i in 1:length(model.Cells)]
    cenY = [get_center_pos(i,model)[2] for i in 1:length(model.Cells)]
    scatter!(centerX,centerY, markersize=3, color=:black, leg=false)
    scatter!(cenX,cenY, markersize=3, color=:red, leg=false)
    # quiver!(centerX,centerY, quiver=(polX, polY).*0.1, color=:black, leg=false)
    step!(model, dummystep, cpm_step!, 1)
    print(i, " ")
end
gif(anim, "test5.gif", fps = 15)








## Test
Cells = [cell for cell in model.Cells];
A = [get_area(i,model)==cell.area for (i,cell) in enumerate(Cells)];
sum(A)
P = [(get_perimeter(i,model)==cell.perimeter) for (i,cell) in enumerate(Cells)];
sum(P)

N = Array{Int}(undef, 10000)
for i in 1:length(model.agents)
    N[i] = model.agents[i].neighbors
end

NN = Array{Int}(undef, 10000);
for i in 1:length(model.agents)
    ag = model.agents[i]
    NN[i] = checkNeighbors(ag, model)
end

sum(NN .== N)

Cx = [get_center_pos(i, model)[1]==cell.centerPos[1] for (i,cell) in enumerate(Cells)]
sum(Cx)
Cy = [get_center_pos(i, model)[2]==cell.centerPos[2] for (i,cell) in enumerate(Cells)]
sum(Cy)
C = [[get_center_pos(i, model),cell.centerPos] for (i,cell) in enumerate(Cells)]