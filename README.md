# ATR-MIR Ray Path Simulator

A Python GUI to simulate MIR ray paths in ATR geometries (rectangle/trapezoid). It computes intersections, refraction/reflection with TIR checks, and shows per-hit incident/critical angles and coordinates.


## Use
- Choose Rectangle or Trapezoid; for trapezoid set Base, Top, and constrain by Height or Angle.
- Set crystal 
 and the per-boundary media (0=Bottom, 1=Left-slanted, 2=Top, 3=Right-slanted).
- Define ray origin and angle. Run Simulation.
- Interactions panel shows: boundary, hit X/Y, incident angle, critical angle (if applicable), n_incn_trn, and direction (or TIR).
