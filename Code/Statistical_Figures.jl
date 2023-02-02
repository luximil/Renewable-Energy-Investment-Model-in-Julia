function Compute_Model_Statistics()
    pv_facilites = map(x -> startswith(x, "pv") ? x : "", NDISP)
    if sum(pv_facilites .!= "") > 0
        avg_pv_generation = 1/length(T) * sum(value.(m[:feedin])[pv_facilites, :])
        println("Average generation from PV: $avg_pv_generation kWh.")
    end

    wind_facilites = map(x -> startswith(x, "wind") ? x : "", NDISP)
    if sum(wind_facilites .!= "") > 0
        avg_wind_generation = 1/length(T) * sum(value.(m[:feedin])[wind_facilites, :])
        println("Average hourly generation from wind: $avg_wind_generation kWh.")
    end

    avg_disp_non_bat_generation = 1/length(T) * (sum(value.(m[:G])[DISP,:]) - sum(value.(m[:G])[S,:]))
    println("Avg. hourly generation from dispatchable and non-battery facilities: $avg_disp_non_bat_generation kWh.")

    avg_lost_load = 1/length(T) * sum(value.(m[:LL]))
    max_lost_load = maximum(value.(m[:LL]))
    println("Average hourly lost load: $avg_lost_load kWh. Maximum lost load: $max_lost_load kWh.")

    avg_curtailment = 1/length(T) * sum(value.(m[:CU]))
    max_curtailment = maximum(value.(m[:CU]))
    println("Average hourly curtailment: $avg_curtailment kWh. Maximum curtailment: $max_curtailment kWh.")

    avg_irrigation_eff = 1/length(T) * sum(WE) * 100
    println("Average irrigation efficiency: $avg_irrigation_eff %.")

    yearly_elec_cost = value.(m[:yearly_elec_cost])
    yearly_investment_cost = value.(m[:yearly_investment_cost])

    println("Total yearly electricity cost: $yearly_elec_cost DZD.")
    println("Total yearly investment cost (computed as annuity): $yearly_investment_cost DZD.")

    total_irrigation = value.(m[:total_water_irrigated])
    total_irrigation_text = "Total water irrigated: $total_irrigation m^3."
    irrigated_to_target_ratio = total_irrigation / (length(T)/24 * DEW)
    irr_ratio_text = "$irrigated_to_target_ratio times the effective water needed."
    println(total_irrigation_text*irr_ratio_text)
end