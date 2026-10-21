"""Run from the repository root after installation: python examples/two_layer.py."""

from yieldledger import Bar, Material, solve_step

bar = Bar([40, 60], (Material(200000, 200, 2000, 4000), Material(100000, 400, 5000, 5000)), 10)
step = solve_step(bar, 0.4)
print(f"Interface displacement: {step.displacement[1]:.9f} mm")
print(f"Reaction: {bar.area * step.stress[-1]:.9f} N")
print(f"Newton corrections: {step.iterations}")
