"""Run after installation: python examples/through_the_fold.py."""

from foldtrace import Truss, trace

model = Truss(ratio=0.2)
path = trace(model, q_stop=2.2)
negative = [point for point in path.points if point.tangent < 0]
print(f"{len(path.points)} equilibria; {len(negative)} have negative symmetric stiffness")
print(f"Stop: {path.stop_reason}; q={path.points[-1].q:.6f}")
print(f"Maximum scaled residual: {max(abs(p.residual) for p in path.points):.3e}")
