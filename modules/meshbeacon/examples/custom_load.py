"""A uniformly loaded bar with prescribed displacement at both endpoints."""

from meshbeacon.core import adapt, errors

result = adapt(lambda x: 2.0, boundary=(1.0, 3.0), tolerance=0.05)
l2, energy = errors(result.solution, lambda x: 1 + 3 * x - x**2, lambda x: 3 - 2 * x)
print(f"stop={result.stop_reason}; elements={result.solution.elements}")
print(f"L2 error={l2:.6e}; energy error={energy:.6e}")
