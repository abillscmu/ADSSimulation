function cell_geometry()
    #Calculate Cell Volumes Based on Capacity
    T⁻ = 86.7e-6
    T⁺ = 66.2e-6
    Tˢ = 12e-6

    W⁻ = 2*61.5e-2
    W⁺ = 2*61.5e-2

    L⁺ = 5.8e-2
    L⁻ = 5.8e-2

    Vₛ⁻ = W⁻*L⁻*T⁻
    Vₛ⁺ = W⁺*L⁺*T⁺
    Vₑ⁻ = W⁻*L⁻*T⁻
    Vₑˢ = W⁺*L⁺*Tˢ
    Vₑ⁺ = W⁺*L⁺*T⁺

    

    cellgeometry = ComponentArray(
        Vₛ⁻ = W⁻*L⁻*T⁻/2,
        Vₛ⁺ = W⁺*L⁺*T⁺/2,
        Vₑ⁻ = W⁻*L⁻*T⁻,
        Vₑˢ = W⁺*L⁺*Tˢ,
        Vₑ⁺ = W⁺*L⁺*T⁺,
        T⁺ = T⁺,
        T⁻ = T⁻,
)

end