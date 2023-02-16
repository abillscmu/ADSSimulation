# Airbus LEO Cycling Dataset

### Raw Data

#### RoomTemperature
#All Files saved as txt (tab-separated). The files are separated by the "first cycle". Within each "First Cycle", the time should be (roughly) correct. Ns represents the cycling phases, which are:
RoomTemperature = Dict(
    "SegmentToExperiment" => Dict(
        1 => (kwargs...) -> "Rest for 15 min",
        2 => (kwargs...) -> "Discharge at 680 mA until 2.5 V",
        3 => (kwargs...) -> "Rest for 15 min",
        4 => (kwargs...) -> "Charge at 1.7 A until 4.2 V",
        5 => (kwargs...) -> "Charge at 4.2 V until 68 mA",
        6 => (kwargs...) -> "Rest for 15 min",
        7 => (kwargs...) -> "GEIS",
        8 => (kwargs...) -> "Rest for 15 min",
        9 => (kwargs...) -> "Discharge at 680 mA until DOD is 20%",
        10 => (kwargs...) -> "Discharge at 68 mA for 30 s",
        11 => (kwargs...) -> "Discharge at 680 mA until DOD is 60%",
        12 => (kwargs...) -> "Discharge at 68 mA for 30 s",
        13 => (kwargs...) -> "Discharge at 680 mA until 2.5 V",
        14 => (kwargs...) -> "Charge at $(kwargs[1]) mA until $(kwargs[2])",
        15 => (kwargs...) -> "Charge at $(kwargs[2]) until 170 mA",
        16 => (kwargs...) -> "Discharge at $(kwargs[3])mA $(kwargs[4])"),
    "CELLS" => Dict(
        "SN-01" => "Baseline",
        "SN-02" => "Baseline",
        "SN-03" => "IncreasedDOD",
        "SN-04" => "IncreasedDOD",
        "SN-05" => "IncreasedEOCV",
        "SN-06" => "IncreasedEOCV",
        "SN-07" => "IncreasedDRate",
        "SN-08" => "IncreasedDRateAndEOCV",
        "SN-17" => "Baseline",
        "SN-18" => "IncreasedDOD",
        "SN-19" => "IncreasedEOCV",
        "SN-20" => "IncreasedDRate",
        "SN-21" => "IncreasedDRate"),
    "CYCLE_PARAMETERS" => Dict(
        "Baseline" => Dict(
            "CHARGE_RATE" => 1133, #(mA)
            "VEOC" => 4.1, #(V)
            "DISCHARGE_RATE" => 680, #(mA)
            "EOD_CRITERIA" => "for 30 min"),
        "IncreasedDOD" => Dict(
            "CHARGE_RATE" => 1133, #(mA)
            "VEOC" => 4.1, #(V)
            "DISCHARGE_RATE" => 680, #(mA)
            "EOD_CRITERIA" => "for 60 min"),
        "IncreasedEOCV" => Dict(
            "CHARGE_RATE" => 1133, #(mA)
            "VEOC" => 4.2, #(V)
            "DISCHARGE_RATE" => 680, #(mA)
            "EOD_CRITERIA" => "for 30 min"),
        "IncreasedDRate" =>Dict(
            "CHARGE_RATE" => 1133, #(mA)
            "VEOC" => 4.1, #(V)
            "DISCHARGE_RATE" => 1133, #(mA)
            "EOD_CRITERIA" => "for 18 min"),
        "IncreasedDRateAndEOCV" => Dict(
            "CHARGE_RATE" => 1133, #(mA)
            "VEOC" => 4.2, #(V)
            "DISCHARGE_RATE" => 1133, #(mA)
            "EOD_CRITERIA" => "for 18 min")))


#### 25deg
#Very similar to room temperature. All Files saved as txt (tab-separated). The files are separated by the "first cycle". Within each "First Cycle", the time should be (roughly) correct. Ns represents the cycling phases, which are:
#=
    1 => "Rest for 15 min"
    2 => "Discharge at 680 mA until 2.5 V"
    3 => "Rest for 15 min"
    4 => "Charge at 1.7 A until 4.2 V"
    5 => "Charge at 4.2 V until 68 mA"
    6 => "Rest for 15 min"
    7 => "GEIS"
    8 => "Rest for 15 min"
    9 => "Discharge at 680 mA until DOD is 20%"
    10 => "Discharge at 68 mA for 30 s"
    11 => "Discharge at 680 mA until DOD is 60%"
    12 => "Discharge at 68 mA for 30 s"
    13 => "Discharge at 680 mA until 2.5 V"
    14 => "Charge at $CHARGE_RATE mA until $VEOC"
    15 => "Charge at $VEOC until 170 mA"
    16 => "Discharge at $DISCHARGE_RATE mA $EOD_CRITERIA"


CELLS: 
    SN-22 => DecreasedEOCV
    SN-23 => DecreasedEOCV
    SN-24 => DecreasedEOCV
    SN-25 => MoreDecreasedEOCV
    SN-26 => MoreDecreasedEOCV
    SN-27 => MoreDecreasedEOCV
    SN-28 => IncreasedCRate
    SN-29 => IncreasedCRate
    SN-30 => IncreasedCRate


PARAMETERS:
    DecreasedEOCV =>
        CHARGE_RATE => 1133 #(mA)
        VEOC => 4 #(V)
        DISCHARGE_RATE => 680 #(mA)
        EOD_CRITERIA => "for 30 min"
    MoreDecreasedEOCV =>
        CHARGE_RATE => 1133 #(mA)
        VEOC => 3.9 #(V)
        DISCHARGE_RATE => 680 #(mA)
        EOD_CRITERIA => "for 30 min"
    IncreasedCRate =>
        CHARGE_RATE => 1700 #(mA)
        VEOC => 4.1 #(V)
        DISCHARGE_RATE => 1133 #(mA)
        EOD_CRITERIA => "for 18 min"
=#