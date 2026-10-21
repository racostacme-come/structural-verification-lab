"""Run with python examples/cantilever.py after installing the package."""

from shearscope import Beam, exact_tip, modes, solve_static

beam = Beam(length=1.2, width=0.025, height=0.012)
for integration in ("full", "reduced"):
    result = solve_static(beam, 32, integration, force=1.0)
    frequencies, _ = modes(beam, 32, integration)
    print(f"{integration}: tip = {result.displacement[-1]:.8g} m; f = {frequencies} Hz")
print(f"Exact static tip: {exact_tip(beam):.8g} m")
