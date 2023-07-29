abstract type AbstractSection end
abstract type TorsionAllowed <: AbstractSection end

"""
Wide-flange sections
"""
struct W <: TorsionAllowed
    name::String
    A #area [mm²]
    d #depth [mm]
    bf #flange width [mm]
    tw #web thickness [mm]
    tf #flange thickness [mm]
    Ix #strong moment of Inertia [mm⁴]
    Zx #strong plastic modulus [mm³]
    Sx #strong section modulus [mm³]
    rx #strong radius of gyration [mm]
    Iy #weak moment of Inertia [mm⁴]
    Zy #weak plastic modulus [mm³]
    Sy #weak section modulus [mm³]
    ry #weak radius of gyration [mm]
    J #torsional constant [mm⁴]
    Cw #warping constant [mm⁶]

    function W(name::String)

        @assert in(name, data[Wrange, 1]) "Name not found"

        irow = findfirst(data[:,1] .== name)
        values = [data[irow, colDict[i]] for i in Wfields]
        values[2:end] .*= Wfactors

        return new(values...)
    end
end

"""
Channels
"""
struct C <: TorsionAllowed
    name::String
    A #area [mm²]
    d #depth [mm]
    bf #flange width [mm]
    tw #web thickness [mm]
    tf #flange thickness [mm]
    Ix #strong moment of Inertia [mm⁴]
    Zx #strong plastic modulus [mm³]
    Sx #strong section modulus [mm³]
    rx #strong radius of gyration [mm]
    Iy #weak moment of Inertia [mm⁴]
    Zy #weak plastic modulus [mm³]
    Sy #weak section modulus [mm³]
    ry #weak radius of gyration [mm]
    J #torsional constant [mm⁴]
    Cw #warping constant [mm⁶]

    function C(name::String)

        @assert in(name, data[Crange, 1]) "Name not found"

        irow = findfirst(data[:,1] .== name)
        values = [data[irow, colDict[i]] for i in Wfields]
        values[2:end] .*= Cfactors

        return new(values...)
    end
end

"""
Angles
"""
struct L <: TorsionAllowed
    name::String
    A #area [mm²]
    d #leg length 1 [mm]
    b #leg length 2 [mm]
    t #thickness [mm]
    Ix #strong moment of Inertia [mm⁴]
    Zx #strong plastic modulus [mm³]
    Sx #strong section modulus [mm³]
    rx #strong radius of gyration [mm]
    Iy #weak moment of Inertia [mm⁴]
    Zy #weak plastic modulus [mm³]
    Sy #weak section modulus [mm³]
    ry #weak radius of gyration [mm]
    J #torsional constant [mm⁴]
    Cw #warping constant [mm⁶]

    function L(name::String)
        @assert in(name, data[Lrange, 1]) "Name not found"
        irow = findfirst(data[:,1] .== name)
        values = [data[irow, colDict[i]] for i in Lfields]
        values[2:end] .*= Lfactors

        return new(values...)
    end


end

"""
Double angles
"""
struct LL <: AbstractSection
    name::String
    A #area [mm²]
    d #leg length 1 [mm]
    b #leg length 2 [mm]
    t #thickness [mm]
    Ix #strong moment of Inertia [mm⁴]
    Zx #strong plastic modulus [mm³]
    Sx #strong section modulus [mm³]
    rx #strong radius of gyration [mm]
    Iy #weak moment of Inertia [mm⁴]
    Zy #weak plastic modulus [mm³]
    Sy #weak section modulus [mm³]
    ry #weak radius of gyration [mm]

    function LL(name::String)
        @assert in(name, data[LLrange, 1]) "Name not found"
        irow = findfirst(data[:,1] .== name)
        values = [data[irow, colDict[i]] for i in LLfields]
        values[2:end] .*= LLfactors

        return new(values...)
    end
end

"""
WT sections
"""
struct WT <: TorsionAllowed
    name::String
    A #area [mm²]
    d #depth [mm]
    bf #flange width [mm]
    tw #web thickness [mm]
    tf #flange thickness [mm]
    Ix #strong moment of Inertia [mm⁴]
    Zx #strong plastic modulus [mm³]
    Sx #strong section modulus [mm³]
    rx #strong radius of gyration [mm]
    Iy #weak moment of Inertia [mm⁴]
    Zy #weak plastic modulus [mm³]
    Sy #weak section modulus [mm³]
    ry #weak radius of gyration [mm]
    J #torsional constant [mm⁴]
    Cw #warping constant [mm⁶]

    function WT(name::String)
        @assert in(name, data[WTrange, 1]) "Name not found"
        irow = findfirst(data[:,1] .== name)
        values = [data[irow, colDict[i]] for i in Wfields]
        values[2:end] .*= Wfactors

        return new(values...)
    end
end

"""
HSS rectangular
"""
struct HSSRect <: TorsionAllowed
    name::String
    A #area [mm²]
    d #depth [mm]
    b #width [mm]
    t #thickness [mm]
    Ix #strong moment of Inertia [mm⁴]
    Zx #strong plastic modulus [mm³]
    Sx #strong section modulus [mm³]
    rx #strong radius of gyration [mm]
    Iy #weak moment of Inertia [mm⁴]
    Zy #weak plastic modulus [mm³]
    Sy #weak section modulus [mm³]
    ry #weak radius of gyration [mm]
    J #torsional constant [mm⁴]

    function HSSRect(name::String)
        @assert in(name, data[HSSRectrange, 1]) "Name not found"
        irow = findfirst(data[:,1] .== name)
        values = [data[irow, colDict[i]] for i in HSSRectfields]
        values[2:end] .*= HSSRectfactors

        return new(values...)
    end
end

"""
HSS round
"""
struct HSSRound <: TorsionAllowed
    name::String
    A #area [mm²]
    OD #depth [mm]
    t #flange thickness [mm]
    Ix #strong moment of Inertia [mm⁴]
    Zx #strong plastic modulus [mm³]
    Sx #strong section modulus [mm³]
    rx #strong radius of gyration [mm]
    Iy #weak moment of Inertia [mm⁴]
    Zy #weak plastic modulus [mm³]
    Sy #weak section modulus [mm³]
    ry #weak radius of gyration [mm]
    J #torsional constant [mm⁴]

    function HSSRound(name::String)
        @assert in(name, data[HSSRoundrange, 1]) "Name not found"
        irow = findfirst(data[:,1] .== name)
        values = [data[irow, colDict[i]] for i in HSSRoundfields]
        values[2:end] .*= HSSRoundfactors

        return new(values...)
    end
end

