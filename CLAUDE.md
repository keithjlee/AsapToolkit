# AsapToolkit.jl

Utility layer over Asap.jl (`../Asap`): structure generators, steel section databases, internal-force diagrams, geometry/displacement sampling, Grasshopper IO. Vendored `AsapSections` submodule for polygon section properties.

**Asap is undergoing a v1.0 modernization (`../Asap/docs/MODERNIZATION.md`); this package migrates in lockstep at Phase 5b.** Key change: `src/ForceAnalysis/` (the `InternalForces` diagram machinery) is being **absorbed into Asap core** (Phase 3, equilibrium-based recovery) and will be deleted here, replaced by a mechanical rename port of consumers. Generators, SteelSections, AsapSections, and IO stay.

## Commands

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Structure

- `Generation/` — parametric model generators (`Frame`, `SpaceFrame`, `Warren2D`, `Pratt2D`, `TrussFrame`, `GridFrame`, `BakerTruss`, ground structures)
- `ForceAnalysis/` — `InternalForces` (samples P, My, Vy, Mz, Vz along elements by superposing per-(load × release) analytic formulas), `load_envelopes`. **Being absorbed into Asap core.**
- `Geometry/` — `ElementDisplacements` sampling
- `SteelSections/` — AISC-style sections from XLSX (`W`, `C`, `L`, `HSSRect`, …) + `toASAPframe`/`toASAPtruss` bridges
- `AsapSections/` — vendored polygon section-property solver
- `IO/` — Grasshopper/JSON exchange

## Traps and coupling

- **Axis-naming trap**: `InternalForces.My` is built from `Flocal[6]` and pairs with `Vy` — it is the moment about local **z**; `.Mz` is the moment about local y. The modern core API fixes this (port map: `old.My → new.Mz`, `old.Mz → new.My`, `P → N`).
- `Flocal` here omits `element.Q`, so diagrams for `GravityLoad`ed elements are silently wrong.
- Heavy direct field access into Asap structs (`.id`, `.elements`, `.section`, `.R`, `.LCS`, `.forces`, …) and non-exported reach-ins (`Asap.nodeids`, release types) — this is the de-facto API contract the modernization must re-provide via accessors.
