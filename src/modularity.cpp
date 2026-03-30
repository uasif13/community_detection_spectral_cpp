#include "modularity.hpp"
#include <algorithm>
#include <stdexcept>

ModResult compute_modularity(
    const Graph& g,
    const std::vector<bool>* keep,
    const std::vector<int>& member,
    int n1,
    int n2)
{
    // ── Pass 1: count m_sparse, intra, inter ─────────────────────────────────
    int m_sparse = 0;
    int intra_kept = 0;
    int inter_kept = 0;

    for (int i = 0; i < g.m; i++) {
        if (keep && !(*keep)[i]) continue;
        m_sparse++;
        auto [u, v] = g.edges[i];
        if (member[u] == member[v])
            intra_kept++;
        else
            inter_kept++;
    }

    ModResult res;
    res.m_sparse        = m_sparse;
    res.preserved_intra = intra_kept;
    res.preserved_inter = inter_kept;
    res.intra_rate      = (n1 > 0) ? static_cast<double>(intra_kept) / n1 : 0.0;
    res.inter_rate      = (n2 > 0) ? static_cast<double>(inter_kept) / n2 : 0.0;
    res.ratio_obs       = (res.intra_rate > 0.0) ? res.inter_rate / res.intra_rate : 0.0;

    if (m_sparse == 0) {
        res.F = 0.0; res.G = 0.0; res.Q = 0.0;
        return res;
    }

    res.F = static_cast<double>(intra_kept) / m_sparse;

    // ── Pass 2: sparse degrees → community volumes → G ──────────────────────
    // Determine number of communities
    int k = 0;
    for (int v = 0; v < g.n; v++)
        if (member[v] + 1 > k) k = member[v] + 1;

    // Accumulate vol_c using sparse degrees (double to avoid overflow)
    std::vector<double> vol(k, 0.0);
    for (int i = 0; i < g.m; i++) {
        if (keep && !(*keep)[i]) continue;
        auto [u, v] = g.edges[i];
        vol[member[u]] += 1.0;
        vol[member[v]] += 1.0;
    }

    // G = sum(vol_c^2) / (4 * m_sparse^2)  — Kahan summation
    double g_sum = 0.0, g_comp = 0.0;
    for (int c = 0; c < k; c++) {
        double term = vol[c] * vol[c];
        double y    = term - g_comp;
        double t    = g_sum + y;
        g_comp      = (t - g_sum) - y;
        g_sum       = t;
    }
    double denom = 4.0 * static_cast<double>(m_sparse) * static_cast<double>(m_sparse);
    res.G = g_sum / denom;
    res.Q = res.F - res.G;

    return res;
}
