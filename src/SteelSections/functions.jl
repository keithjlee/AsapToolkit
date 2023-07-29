#get all tabulated sections
allW() = [W(name) for name in names[Wrange]]
allC() = [C(name) for name in names[Crange]]
allL() = [L(name) for name in names[Lrange]]
allLL() = [LL(name) for name in names[LLrange]]
allWT() = [WT(name) for name in names[WTrange]]
allHSSRect() = [HSSRect(name) for name in names[HSSRectrange]]
allHSSRound() = [HSSRound(name) for name in names[HSSRoundrange]]

#get all names
Wnames = names[Wrange]
Cnames = names[Crange]
Lnames = names[Lrange]
LLnames = names[LLrange]
WTnames = names[WTrange]
HSSRectnames = names[HSSRectrange]
HSSRoundnames = names[HSSRoundrange]

unitfactors = Dict(:mm => [1, 1, 1, 1, 1, 1],
    :m => [1e-6, 1e6, 1e6, 1e-12, 1e-12, 1e-12],
    :in => [1/25.4^2, 25.4^2, 25.4^2, 1/25.4^4, 1/25.4^4, 1/25.4^4])

#convert units
function toASAPframe(section::TorsionAllowed, E::Real, G::Real; unit = :mm)
    println("E and G should have distance units in [mm]")
    @assert in(unit, keys(unitfactors))

    vals = [section.A,
        E,
        G,
        section.Ix,
        section.Iy,
        section.J]

    factor = unitfactors[unit]

    return Section((vals .* factor)...)

end

trussunitfactors = Dict(:mm => [1., 1.],
    :m => [1e-6, 1e6],
    :in => [1/25.4^2, 25.4^2])

function toASAPtruss(section::AbstractSection, E::Real; unit = :mm)
    println("E should have distance unit in [mm]")
    @assert in(unit, keys(trussunitfactors))

    vals = [section.A, E] .* trussunitfactors[unit]

    return TrussSection(vals...)

end