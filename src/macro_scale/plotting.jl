
"""
    animate_macro_result(res::NamedTuple, setup::SolveSetup{dim}; file_name ="Myresult.mp4", kwargs...) where {dim}

Return an mp4 file with name `file_name`, which by default is "Myresult.mp4". 

This shows the evolution of the solution `res` for a macro scale simulation of a multi scale problem based on input `SolveSetup`. 

Further arguments like `scale` for visibel deformation can be used.
    
"""
function animate_macro_result(res::NamedTuple, setup::SolveSetup{dim}; file_name ="Myresult.mp4", kwargs...) where {dim}
    (; t, a) = res
    fig = Makie.Figure(size=(1200,800))
    tᵒᵇˢ, aᵒᵇˢ = _prepare_plots!(fig, res, setup; kwargs...)

	file = joinpath(file_name)
	anim = Makie.record(fig, file, eachindex(t); framerate=1) do i
        tᵒᵇˢ[] = t[i]
        aᵒᵇˢ[] = a[i]
	end
	return file, fig, anim
end

"""
    _prepare_plots!(pos, res::NamedTuple, setup::SolveSetup{dim}; scale::Real=1.0, titlestart::String="macroscale solution") where {dim}

Sets up the necessary plots with a specified layout position `pos` and observables `tᵒᵇˢ` and `aᵒᵇˢ` to visualize the simulation results `res` of a multi scale problem `SolveSetup`:
deformation `u`, chemical potential `μ`.

Optional arguments include a scaling factor `scale` for visualizing the deformation and a title `titlestart`. By default the `scale` is set to `1.0` and `titlestart` to `macroscale solution`.

"""
function _prepare_plots!(pos, res::NamedTuple, setup::SolveSetup{dim};
        scale::Real=1.0, 
        titlestart::String="macroscale solution") where {dim}
    (; grid, dh) = setup
    (; t, a)= res

    μ_all = [evaluate_at_grid_nodes(dh, res.a[i], :μ) for i in eachindex(res.t)]
    μ_min, μ_max = minimum(μ -> minimum(μ), μ_all), maximum(μ -> maximum(μ), μ_all)
         
    mesh = _prepare_plotable_mesh(grid)

    
    tᵒᵇˢ = Makie.Observable(t[1])
    aᵒᵇˢ = Makie.Observable(a[1])

    title = Makie.@lift titlestart*" at t=$( round($(tᵒᵇˢ); sigdigits=4) )"
    Makie.Label(pos[1,1:2], title)

    ax  = Makie.Axis3(pos[3,1], aspect=:equal, title="undeformed grid with RVE location at cell[1] qp[1]")
    Makie.wireframe!(ax, mesh; color=:black)
    for (cellid, qp) in ((1,1),)
        cv = setup.assemblysetup.cv.u
        coords = getcoordinates(grid, cellid)
        reinit!(cv, coords)
        x = spatial_coordinate(cv, qp, coords)

        Makie.meshscatter!(x[1], x[2], x[3], markersize = 10)
        break
    end
    
    u = Makie.@lift evaluate_at_grid_nodes(dh, $(aᵒᵇˢ), :u)
    defmesh = Makie.@lift _prepare_plotable_mesh(grid, ( $(u) .* scale ))
    ax = Makie.Axis3(pos[2,1], aspect=:equal, title="deformed grid")
    Makie.mesh!(ax, defmesh; color=Makie.RGB(1.0,0.5,0.5), shading=Makie.NoShading)
    Makie.wireframe!(ax, defmesh; color=:black)

    subpos = pos[2,2]
    μ = Makie.@lift evaluate_at_grid_nodes(dh, $(aᵒᵇˢ), :μ)
    ax = Makie.Axis3(subpos[1,1], aspect=:equal, title="chemical potential")
    colorsettings = (colorrange=(0.0, 1.0), colormap=:viridis)
    Makie.mesh!(ax, mesh; color=μ, colorsettings..., shading=Makie.NoShading)
    Makie.Colorbar(subpos[1,2]; colorsettings...)

    return tᵒᵇˢ, aᵒᵇˢ
end


"""
    animate_combined_result(res::NamedTuple, resᵣᵥₑ::NamedTuple, setup::SolveSetup{dim}; file_name ="Myresult.mp4", scale=1.0) where {dim}

Return an mp4 file with name `file_name`, which by default is "Myresult.mp4". 

This shows the evolution of the macro scale solution `res` and the RVE solution at a corresponding quadrature point `cell[1]qp[1]` for a multi scale simulation based on input `SolveSetup`.

Further arguments like `scale` for visibel deformation can be used, which by default is set to `1.0`.

"""
function animate_combined_result(res::NamedTuple, resᵣᵥₑ::NamedTuple, setup::SolveSetup{dim}; file_name ="Myresult.mp4", scale=1.0) where {dim}
    (; rvesetup) = setup.assemblysetup

    fig = Makie.Figure(size=(1200,800))

    tᵒᵇˢ, aᵒᵇˢ = _prepare_plots!(fig[1:2,1:2], res, setup; scale=scale)
    rveᵒᵇˢ = _prepare_plots!(fig[2,2], resᵣᵥₑ, rvesetup; titlestart="RVE results", scale=scale) 

	file = joinpath(file_name)
	anim = Makie.record(fig, file, eachindex(res.t); framerate=1) do i
        tᵒᵇˢ[] = res.t[i]
        aᵒᵇˢ[] = res.a[i]
        rveᵒᵇˢ[1][] = resᵣᵥₑ.t[i]
        rveᵒᵇˢ[2][] = resᵣᵥₑ.a[i]
	end
	return file, fig, anim
end