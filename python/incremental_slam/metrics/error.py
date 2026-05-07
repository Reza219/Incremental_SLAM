from __future__ import annotations

from dataclasses import dataclass

from incremental_slam.linearization.edge_jacobian import jacobian_edge_jr
from incremental_slam.types import Graph


@dataclass(frozen=True)
class ChiSquaredSummary:
    chi2: float
    num_scalar_measurements: int

    @property
    def normalized_chi2(self) -> float:
        return self.chi2 / max(1, self.num_scalar_measurements)


def compute_chi2_summary(g: Graph) -> ChiSquaredSummary:
    total = 0.0
    m_total = 0
    for edge in g.edges:
        _, re = jacobian_edge_jr(edge, g, g.state_size)
        total += float(re @ re)
        m_total += int(re.size)
    return ChiSquaredSummary(chi2=total, num_scalar_measurements=m_total)


def compute_global_error(g: Graph) -> float:
    """Paper's normalized chi-squared error: Nchi^2_t = 2 c(x_t) / M_t.

    The linearization code returns whitened residuals, so sum(re^T re) equals
    2 c(x_t). M_t is the number of scalar residual equations, not the number of
    graph edges.
    """
    return compute_chi2_summary(g).normalized_chi2
