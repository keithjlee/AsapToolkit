section = toASAPframe(rand(allW()), Steel_Nmm.E, Steel_Nmm.G)

n1 = Node([0., 0., 0.], :fixed)
n2 = Node([6000., 0., 0.], :fixed)
nodes = [n1, n2]

element = Element(n1, n2, section)
elements = [element]

load1 = LineLoad(element, [0., 0., -30])
pointloads = [PointLoad(element, rand(), 25e3 .* [0., 0., -rand()]) for _ = 1:5]
loads = [load1; pointloads]

model = Model(nodes, elements, loads)
solve!(model)

fanalysis = InternalForces(element, model)
danalysis = ElementDisplacements(element, model)