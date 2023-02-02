
function plot_result()
    
    # -------------------------------------------------------------------------------------------------------
    # Get the values of the optimal variables.                                                              |
    # -------------------------------------------------------------------------------------------------------

    generation = value.(m[:G])
    storage_gen = value.(m[:CH])
    curtailment = value.(m[:CU])
    storage_level = value.(m[:L])
    lost_load = value.(m[:LL])

    power = value.(m[:PW])[vcat(DISP, NDISP)]
    #storage = value.(m[:PW])[S]

    generation_combined = vcat(generation, value.(m[:feedin]), lost_load.data')
    cum_generation = vcat(zeros(1, length(T)), cumsum(generation_combined, dims=1))

    cum_charge = vcat(
        zeros(1, length(T)),
        -cumsum(vcat(storage_gen.data, curtailment.data'), dims=1)
    )

    water_storage_level = value.(m[:WL])

    # -------------------------------------------------------------------------------------------------------
    # Plot hourly electricity generation and demand by source.                                              |
    # -------------------------------------------------------------------------------------------------------

    fig = Figure(resolution = (1600, 1000), fontsize = 30)
    ax = Axis(fig[2,1], ylabel = "kW", xlabel="Hours")
    colors = ["#754937", "#8e8eb5", "#ffeb3b", "#518696"]
    labels = vcat(string.(P), "Lost load")

    legend_items = []
    legend_items_labels = String[]
    for row in 1:size(cum_generation, 1)-1
        b = band!(
            ax,
            1:length(T),
            cum_generation[row,:],
            cum_generation[row+1,:],
            color=(colors[row], 0.8),
            linewidth=0
        )
        push!(legend_items, b)
        push!(legend_items_labels, labels[row])
    end

    colors2 = ["#8e8eb5", "#ae123a"]
    for row in 1:size(cum_charge, 1)-1
        b = band!(
            ax,
            1:length(T),
            cum_charge[row,:],
            cum_charge[row+1,:],
            color=(colors2[row], 0.8),
            linewidth=0
        )
        row == size(cum_charge, 1)-1 && push!(legend_items, b)
    end

    push!(legend_items_labels, "Curtailment")

    #hlines!(ax, 0, color=:black, linewidth=2)

    demand = value.(m[:elec_demand_irrigation]) .+ OD

    p_demand = lines!(
        ax,
        -demand.data,
        color=:black,
        label="Demand",
        linewidth=2
    )
    push!(legend_items, p_demand)
    push!(legend_items_labels, "Demand")

    Legend(
        fig[1,1],
        legend_items,
        legend_items_labels,
        orientation = :horizontal,
        nbanks = 1
    )

    save("Results/Hourly_generation_and_demand_by_source.png", fig)

    # -------------------------------------------------------------------------------------------------------
    # Plot mean day hourly electricity generation and demand by source.                                     |
    # -------------------------------------------------------------------------------------------------------
	
    meanDay_generation_combined = Matrix{Float64}(undef, size(generation_combined, 1), 24)
	for row in 1:size(meanDay_generation_combined, 1)
        for hour in 1:24
            hourIndices = [hour+(i-1)*24 for i in 1:Int(length(T)/24)]
            meanDay_generation_combined[row, hour] = 1/(length(T)/24) * sum(generation_combined[row, hourIndices])
        end
    end
    meanDay_cum_generation = vcat(zeros(1, 24), cumsum(meanDay_generation_combined, dims=1))

    meanDay_storage_gen = Matrix{Float64}(undef, size(storage_gen.data, 1), 24)
    for row in 1:size(meanDay_storage_gen, 1)
        for hour in 1:24
            hourIndices = [hour+(i-1)*24 for i in 1:Int(length(T)/24)]
            meanDay_storage_gen[row, hour] = 1/(length(T)/24) * sum(storage_gen.data[row, hourIndices])
        end
    end
    meandDay_curtailment = Vector{Float64}(undef, 24)
    for hour in 1:24
        hourIndices = [hour+(i-1)*24 for i in 1:Int(length(T)/24)]
        meandDay_curtailment[hour] = 1/(length(T)/24) * sum(curtailment.data[hourIndices])
    end
    meanDay_cum_charge = vcat(
        zeros(1, 24),
        -cumsum(vcat(meanDay_storage_gen, meandDay_curtailment'), dims=1)
    )
	
    fig = Figure(resolution = (1600, 1000), fontsize = 30)
    ax = Axis(fig[2,1], ylabel = "kW", xlabel="Hours")
    colors = ["#754937", "#8e8eb5", "#ffeb3b", "#518696"]
    labels = vcat(string.(P), "Lost load")

    legend_items = []
    legend_items_labels = String[]
    for row in 1:size(meanDay_cum_generation, 1)-1
        b = band!(
            ax,
            1:24,
            meanDay_cum_generation[row,:],
            meanDay_cum_generation[row+1,:],
            color=(colors[row], 0.8),
            linewidth=0
        )
        push!(legend_items, b)
        push!(legend_items_labels, labels[row])
    end

    colors2 = ["#8e8eb5", "#ae123a"]
    for row in 1:size(meanDay_cum_charge, 1)-1
        b = band!(
            ax,
            1:24,
            meanDay_cum_charge[row,:],
            meanDay_cum_charge[row+1,:],
            color=(colors2[row], 0.8),
            linewidth=0
        )
        row == size(meanDay_cum_charge, 1)-1 && push!(legend_items, b)
    end

    push!(legend_items_labels, "Curtailment")

    hlines!(ax, 0, color=:black, linewidth=2)

    demand = value.(m[:elec_demand_irrigation]) .+ OD
	meanDay_demand = Vector{Float64}(undef, 24)
    for hour in 1:24
        hourIndices = [hour+(i-1)*24 for i in 1:Int(length(T)/24)]
        meanDay_demand[hour] = 1/(length(T)/24) * sum(demand.data[hourIndices])
    end

    p_demand = lines!(
        ax,
        -meanDay_demand,
        color=:black,
        label="Demand",
        linewidth=2
    )
    push!(legend_items, p_demand)
    push!(legend_items_labels, "Demand")

    Legend(
        fig[1,1],
        legend_items,
        legend_items_labels,
        orientation = :horizontal,
        nbanks = 1
    )

    save("Results/Mean_day_hourly_generation_and_demand_by_source.png", fig)

    # -------------------------------------------------------------------------------------------------------
    # Plot hourly electricity generation and demand by source for two example days.                         |
    # -------------------------------------------------------------------------------------------------------
    daysoftheyear = [1, 213]
    for day in daysoftheyear
        dayhours = 1+(day-1)*24:day*24

        fig = Figure(resolution = (1600, 1000), fontsize = 30)
        ax = Axis(fig[2,1], ylabel = "kW", xlabel="Hours")
        colors = ["#754937", "#8e8eb5", "#ffeb3b", "#518696"]
        labels = vcat(string.(P), "Lost load")

        legend_items = []
        legend_items_labels = String[]
        for row in 1:size(cum_generation, 1)-1
            b = band!(
                ax,
                1:24,
                cum_generation[row, dayhours],
                cum_generation[row+1, dayhours],
                color=(colors[row], 0.8),
                linewidth=0
            )
            push!(legend_items, b)
            push!(legend_items_labels, labels[row])
        end

        colors2 = ["#8e8eb5", "#ae123a"]
        for row in 1:size(cum_charge, 1)-1
            b = band!(
                ax,
                1:24,
                cum_charge[row,dayhours],
                cum_charge[row+1,dayhours],
                color=(colors2[row], 0.8),
                linewidth=0
            )
            row == size(cum_charge, 1)-1 && push!(legend_items, b)
        end

        push!(legend_items_labels, "Curtailment")

        #hlines!(ax, 0, color=:black, linewidth=2)

        demand = value.(m[:elec_demand_irrigation]) .+ OD

        p_demand = lines!(
            ax,
            -demand.data[dayhours],
            color=:black,
            label="Demand",
            linewidth=2
        )
        push!(legend_items, p_demand)
        push!(legend_items_labels, "Demand")

        Legend(
            fig[1,1],
            legend_items,
            legend_items_labels,
            orientation = :horizontal,
            nbanks = 1
        )

        save("Results/Day_of_Year_$day"*"_hourly_generation_and_demand_by_source.png", fig)
    end

    # -------------------------------------------------------------------------------------------------------
    # Plot hourly electricity storage levels.                                                               |
    # -------------------------------------------------------------------------------------------------------

    legend_items = []
    legend_items_labels = String[]

    fig = Figure(resolution = (1600, 1000), fontsize = 30)
    ax = Axis(fig[2,1], ylabel = "kWh", title="Storage Levels")#, height=120)
    for (i,s) in enumerate(S)
        storage_lines = lines!(ax, storage_level[s,:].data, color=colors2[i])
        push!(legend_items, storage_lines)
        push!(legend_items_labels, S[i])
    end

    Legend(
        fig[1,1],
        legend_items,
        legend_items_labels,
        orientation = :horizontal,
        nbanks = 1
    )

    save("Results/Hourly_electricity_storage_levels.png", fig)

    # -------------------------------------------------------------------------------------------------------
    # Plot daily electricity storage levels.                                                               |
    # -------------------------------------------------------------------------------------------------------

    dailyStorage_level = DenseAxisArray{Float64}(undef, S, 1:Int(length(T)/24))
    for s in S
        for i in 1:Int(length(T)/24)
            dailyStorage_level[s, i] = 1/24 * sum(storage_level[s, 24*(i-1)+1 : 24*i])
        end
    end

    legend_items = []
    legend_items_labels = String[]

    fig = Figure(resolution = (1600, 1000), fontsize = 30)
    ax = Axis(fig[2,1], ylabel = "kWh", title="Daily Average Storage Levels")#, height=120)
    for (i,s) in enumerate(S)
        storage_lines = lines!(ax, dailyStorage_level[s,:].data, color=colors2[i])
        push!(legend_items, storage_lines)
        push!(legend_items_labels, S[i])
    end

    Legend(
        fig[1,1],
        legend_items,
        legend_items_labels,
        orientation = :horizontal,
        nbanks = 1
    )

    save("Results/Daily_average_electricity_storage_levels.png", fig)

    # -------------------------------------------------------------------------------------------------------
    # Plot target capacity by facility.                                                                     |
    # -------------------------------------------------------------------------------------------------------

    fig = Figure(resolution = (1600, 1000), fontsize = 30)
    power_ticks = axes(power, 1)
    ax = Axis(
        fig[1,1],
        title="Capacity",
        xticks=(1:length(power_ticks), power_ticks),
        ylabel="kW",
        xticklabelrotation=pi/4
    )

    barplot!(
        ax,
        power.data,
        color=["#754937", "#8e8eb5", "#ffeb3b"]#, "#518696"]
    )

    save("Results/Target_capacity_by_facility.png", fig)

    # -------------------------------------------------------------------------------------------------------
    # Plot hourly water storage levels.                                                                     |
    # -------------------------------------------------------------------------------------------------------

    legend_items = []
    legend_items_labels = String[]

    fig = Figure(resolution = (1600, 1000), fontsize = 30)
    ax = Axis(fig[2,1], ylabel = "m^3", title="Water Storage Levels")#, height=120)
    for (i,wta) in enumerate(WTa)
        storage_lines = lines!(ax, water_storage_level[wta,:].data, color=:blue)
        push!(legend_items, storage_lines)
        push!(legend_items_labels, WTa[i])
    end

    Legend(
        fig[1,1],
        legend_items,
        legend_items_labels,
        orientation = :horizontal,
        nbanks = 1
    )

    save("Results/Hourly_water_storage_levels.png", fig)

    # -------------------------------------------------------------------------------------------------------
    # Plot daily water storage levels.                                                                      |
    # -------------------------------------------------------------------------------------------------------
    
    dailyWater_storage_level = DenseAxisArray{Float64}(undef, WTa, 1:Int(length(T)/24))
    for wta in WTa
        for i in 1:Int(length(T)/24)
            dailyWater_storage_level[wta, i] = 1/24 * sum(water_storage_level[wta, 24*(i-1)+1 : 24*i])
        end
    end

    legend_items = []
    legend_items_labels = String[]

    fig = Figure(resolution = (1600, 1000), fontsize = 30)
    ax = Axis(fig[2,1], ylabel = "m^3", title="Daily Average Water Storage Levels")#, height=120)
    for (i,wta) in enumerate(WTa)
        storage_lines = lines!(ax, dailyWater_storage_level[wta,:].data, color=:blue)
        push!(legend_items, storage_lines)
        push!(legend_items_labels, WTa[i])
    end

    Legend(
        fig[1,1],
        legend_items,
        legend_items_labels,
        orientation = :horizontal,
        nbanks = 1
    )

    save("Results/Daily_average_water_storage_levels.png", fig)

    # -------------------------------------------------------------------------------------------------------
    # Plot hourly water efficiency.                                                                         |
    # -------------------------------------------------------------------------------------------------------

    fig = Figure(resolution = (1600, 1000), fontsize = 30)
    ax = Axis(fig[1,1], ylabel = "%", title="Irrigation Efficiency")
    barplot!(ax, WE .* 100, color=:blue)

    save("Results/Hourly_irrigation_efficiency.png", fig)

    # -------------------------------------------------------------------------------------------------------
    # Plot daily water efficiency.                                                                          |
    # -------------------------------------------------------------------------------------------------------

    dailyWE = zeros(Int(length(T)/24))
    for i in 1:Int(length(T)/24)
        dailyWE[i] = 1/24 * sum(WE[24*(i-1)+1 : 24*i])
    end

    fig = Figure(resolution = (1600, 1000), fontsize=30)
    ax = Axis(fig[1,1], ylabel = "%", title="Daily Average Irrigation Efficiency")
    barplot!(ax, dailyWE .* 100, color=:blue)

    save("Results/Daily_average_irrigation_efficiency.png", fig)

    return


end