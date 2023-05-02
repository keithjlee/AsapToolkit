using Asap, AsapToolkit

#sections
begin
    W = allW()
    iStiffness = sortperm(getproperty.(W, :Ix))
    nsecs = length(W)

    mat = Steel_Nmm
    girder = toASAPframe((W[iStiffness])[Int(round(nsecs * 0.5))], mat.E, mat.G)
    joist = toASAPframe((W[iStiffness])[Int(round(nsecs * .2))], mat.E, mat.G)

    H = allHSSRound()
    iArea = sortperm(getproperty.(H, :A))
    nh = length(H)
    tube = toASAPframe((H[iArea])[Int(round(nh * 0.75))], mat.E, mat.G)

    girder.ρ = joist.ρ = tube.ρ = mat.ρ
end

nx = 4
ny = 3
nz = 10
dx = 4000
dy = 3200
dz = 3500

model = AsapToolkit.generateFrame(nx,
    dx,
    ny,
    dy,
    nz,
    dz,
    1000,
    girder,
    girder,
    joist,
    tube)

solve!(model)