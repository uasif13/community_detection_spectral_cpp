# CLAUDE.md — DSpar Experiment 1: C++ Implementation Guide

This document is the complete specification for implementing Experiment 1 from the paper
**"Less is More: Sparsification Clarifies Community Structure"** in C++.
Read this file in full before writing any code.

---

## 1. What This Project Does

The experiment validates a theoretical claim about graph sparsification and community
detection modularity. Concretely:

- Load a real-world graph (SNAP format)
- Run the Leiden community detection algorithm on the original graph → get a **fixed partition**
- For a range of **retention levels α** (e.g. 0.2, 0.4, …, 0.9), run **N independent trials**:
  - Sample a sparse subgraph by keeping each edge with probability `p(e) = α · s(e) / s̄`
    where `s(e) = 1/d_u + 1/d_v` is the DSpar score
  - Evaluate **modularity Q = F − G** on the sparse graph using the **frozen partition**
  - Re-run Leiden on the sparse graph → evaluate Q again (pipeline metric)
  - Record all decomposition terms: F, G, ΔQ_fixed, ΔQ_Leiden, ΔF, ΔG
  - Verify `|ΔQ_fixed − (ΔF − ΔG)| < 1e-12` (must hold to machine precision)
- Write results to CSV files that match the Python reference output column-for-column

The experiment runs on 11 SNAP datasets, each independently. Results must be reproducible
given the same seed.

---

## 2. Mathematical Definitions

These are the exact formulas. Implement them precisely.

### 2.1 DSpar Score

For an undirected edge e = (u, v):

```
s(e) = 1/deg(u) + 1/deg(v)
```

All degrees are from the **original, unsparsified graph**. Scores are computed once and
reused across all trials and all α values.

```
s̄ = (1/m) · Σ_{e ∈ E} s(e)      (global mean score)
```

### 2.2 Edge Sampling (Bernoulli / probabilistic_no_replace method)

For each edge e, independently retain it with probability:

```
p(e) = α · s(e) / s̄
```

This is the **probabilistic_no_replace** method from the Python codebase — each edge is
an independent Bernoulli trial. The expected number of retained edges is α·m, but the
actual count varies. Clamp p(e) to ≤ 1 (enforce as assertion in debug builds).

Note: The Python codebase also has a "paper" method (sampling with replacement + reweighting).
**Do not implement the paper method** — the experiment uses `probabilistic_no_replace`.
Check `exp1_2_theoretical_predictions.py` line: `dspar_sparsify(G, retention=retention,
method="paper", seed=seed)` — but in the Leiden evaluation the graph is converted to
unweighted topology anyway. The Bernoulli approach is equivalent for the modularity
decomposition since we only use graph topology (not weights) in F and G.

### 2.3 Modularity Decomposition

Given a sparse graph G' and a **fixed** membership vector (from the original graph):

```
m' = number of edges in G'  (undirected count)

F = (number of edges in G' whose endpoints share the same community) / m'

For each community c:
  vol_c = Σ_{v ∈ c} deg_{G'}(v)      (degree in the SPARSE graph, not original)

G = Σ_c vol_c² / (4 · m'²)

Q = F − G
```

**Critical detail:** `vol_c` uses the **sparse graph degrees**, not original degrees.
The degree of a node in the sparse graph is just how many retained edges it has.

### 2.4 The Identity (must hold to 1e-12)

```
ΔQ_fixed = Q_sparse_fixed − Q_orig_fixed
ΔF       = F_sparse − F_orig
ΔG       = G_sparse − G_orig

Identity: ΔQ_fixed = ΔF − ΔG    (algebraically exact, not approximate)
```

Verify `|ΔQ_fixed − (ΔF − ΔG)| < 1e-12` for every trial. If it ever fails, there is a
bug in the F or G computation.

### 2.5 DSpar Structural Metrics (also report these)

```
mu_intra = mean DSpar score of intra-community edges (both endpoints in same community)
mu_inter = mean DSpar score of inter-community edges (endpoints in different communities)
delta = mu_intra − mu_inter    (DSpar separation)
```

These use the **original graph** scores and the **fixed partition**.

---

## 3. Reference Python Implementation

The Python experiment lives in `PAPER_EXPERIMENTS/exp1_2_theoretical_predictions.py`.
Key functions and their C++ equivalents:

### Python → C++ mapping

| Python function | C++ equivalent | Notes |
|---|---|---|
| `compute_dspar_scores(G)` | `compute_scores(g)` → `DsparScores` | Returns per-edge scores and mean |
| `classify_edges(G, membership)` | Inline in modularity pass | Intra vs inter classification |
| `compute_G_term_unweighted(G, membership)` | `compute_modularity(g, keep, member).G` | Uses sparse graph degrees |
| `compute_observed_fixed_membership(G, G_sparse, membership)` | `compute_modularity(g, keep, member)` | Fixed membership on sparse graph |
| `run_single_trial(G, membership, Q_orig, retention, seed)` | `run_trial(...)` | Core per-trial function |
| `dspar_sparsify(G, retention, method="paper", seed=seed)` | `sample_edges(scores, alpha, rng)` | Returns keep[] mask |
| `leiden_partition(G)` | `run_leiden(g, keep)` | igraph C API |

### Python experiment loop (translate this exactly)

```python
# One-time setup per dataset
G = load_graph(dataset_name)
membership_fixed, Q_orig_leiden = leiden_partition(G)
Q_fixed_original = modularity_fixed_membership(G, membership_fixed)
theory = compute_theory(G, membership_fixed)  # scores, mu_intra, mu_inter, delta

# Per-alpha, per-trial loop
for retention in RETENTIONS:          # e.g. [0.5, 0.6, 0.7, 0.8, 0.9]
    for trial in range(N_TRIALS):     # e.g. 10
        seed = BASE_SEED + trial      # Python uses: hash-based or sequential
        G_sparse = dspar_sparsify(G, retention=retention, method="paper", seed=seed)
        # Convert to unweighted (paper method returns weighted graph)
        G_sparse_unweighted = nx.Graph(G_sparse)  # topology only

        obs = compute_observed_fixed_membership(G, G_sparse_unweighted, membership_fixed)
        # obs contains: F_observed, modularity_fixed_sparse, preserved_intra/inter counts

        G_sparse_ig = nx_to_igraph(G_sparse_unweighted)
        membership_sparse, Q_leiden_sparse = leiden_partition(G_sparse_unweighted)

        G_sparse_observed = compute_G_term_unweighted(G_sparse_unweighted, membership_fixed)
        dG_observed = G_sparse_observed - theory["G_original"]
        dF_observed = obs["F_observed"] - theory["F_original"]

        dQ_fixed = obs["modularity_fixed_sparse"] - Q_fixed_original
        dQ_reconstructed = dF_observed - dG_observed
        # assert abs(dQ_fixed - dQ_reconstructed) < 1e-12
```

### Output CSV columns (must match exactly)

Raw trial CSV (`<dataset>_theoretical_validation_FIXED.csv`):

| Column | Type | Description |
|---|---|---|
| `retention` | float | α value for this trial |
| `seed` | int | RNG seed used |
| `m` | int | Edge count of original graph |
| `n1` | int | Intra-community edge count (original) |
| `n2` | int | Inter-community edge count (original) |
| `mu_intra` | float | Mean DSpar score of intra edges |
| `mu_inter` | float | Mean DSpar score of inter edges |
| `delta` | float | mu_intra − mu_inter |
| `F_original` | float | F on original graph |
| `G_original` | float | G on original graph |
| `modularity_fixed_original` | float | Q_orig = F_orig − G_orig |
| `m_sparse` | int | Edge count of sparse graph |
| `F_observed` | float | F on sparse graph (fixed partition) |
| `G_sparse_observed` | float | G on sparse graph (fixed partition) |
| `dG_observed` | float | G_sparse − G_orig |
| `F_improvement_observed` | float | F_sparse − F_orig |
| `modularity_fixed_sparse` | float | Q on sparse graph (fixed partition) |
| `modularity_fixed_change` | float | Q_sparse − Q_orig (ΔQ_fixed) |
| `modularity_leiden_sparse` | float | Q on sparse graph (Leiden re-run) |
| `modularity_leiden_change` | float | Q_leiden_sparse − Q_orig |
| `dQ_reconstructed` | float | ΔF − ΔG (must equal modularity_fixed_change) |
| `preserved_intra` | int | Intra edges kept in sparse graph |
| `preserved_inter` | int | Inter edges kept in sparse graph |
| `intra_rate` | float | preserved_intra / n1 |
| `inter_rate` | float | preserved_inter / n2 |
| `ratio_observed` | float | inter_rate / intra_rate |

Summary CSV (`<dataset>_summary_FIXED.csv`) — aggregated over trials per retention:

| Column | Description |
|---|---|
| `retention` | α |
| `modularity_fixed_change_mean` | mean ΔQ_fixed |
| `modularity_fixed_change_std` | std ΔQ_fixed |
| `modularity_leiden_change_mean` | mean ΔQ_Leiden |
| `modularity_leiden_change_std` | std ΔQ_Leiden |
| `F_improvement_observed_mean` | mean ΔF |
| `F_improvement_observed_std` | std ΔF |
| `dG_observed_mean` | mean ΔG |
| `dG_observed_std` | std ΔG |
| `delta_mean` | mean δ (same across trials for one dataset) |

---

## 4. Project Structure

Create this exact directory layout in the empty repo:

```
dspar_exp1/
├── CLAUDE.md                    ← this file
├── CMakeLists.txt
├── configs/
│   └── exp1.toml               ← default run configuration
├── extern/
│   └── toml11/                 ← header-only TOML parser (clone from github)
├── src/
│   ├── main.cpp
│   ├── graph.hpp / graph.cpp
│   ├── rng.hpp                 ← xoshiro256++ header-only
│   ├── dspar.hpp / dspar.cpp
│   ├── sampler.hpp / sampler.cpp
│   ├── modularity.hpp / modularity.cpp
│   ├── leiden.hpp / leiden.cpp
│   ├── experiment.hpp / experiment.cpp
│   ├── config.hpp / config.cpp
│   ├── downloader.hpp / downloader.cpp
│   └── csv_writer.hpp
├── tests/
│   ├── CMakeLists.txt
│   ├── test_main.cpp           ← Catch2 main
│   ├── test_graph.cpp
│   ├── test_dspar.cpp
│   ├── test_sampler.cpp
│   ├── test_modularity.cpp
│   ├── test_leiden.cpp
│   └── test_experiment.cpp
├── data/                       ← auto-created, gitignored
└── results/                    ← auto-created, gitignored
```

---

## 5. Implementation Parts (implement in order)

### Part 1: Graph loading and CSR representation

**Files:** `src/graph.hpp`, `src/graph.cpp`, `tests/test_graph.cpp`

**Data structure:**

```cpp
struct Graph {
    int n;                           // number of nodes
    int m;                           // number of undirected edges
    std::vector<int> deg;            // deg[v] = degree of node v
    std::vector<int> adj_off;        // CSR offsets, size n+1
    std::vector<int> adj;            // CSR neighbor list, size 2*m
    std::vector<std::pair<int,int>> edges; // canonical edge list, u < v, size m
};

Graph load_snap(const std::string& path);
```

**Parsing rules (match Python `load_snap_dataset`):**
- Skip lines starting with `#`
- Parse two integers per line: src dst
- Build set of undirected edges: for each (src,dst), store (min,max)
- Remap all node IDs to contiguous 0-indexed range (sort unique node set, map old→new)
- Remove self-loops (src == dst after remapping)
- Build CSR: for each node v, adj[adj_off[v]..adj_off[v+1]] = all neighbors
- `edges` vector: sorted canonical pairs (u < v), length exactly m
- `deg[v]` = adj_off[v+1] - adj_off[v]

**Tests in `test_graph.cpp`:**

```
TEST_CASE("triangle graph") {
  // 3 nodes, 3 edges: 0-1, 1-2, 0-2
  // Write to tmp file, load, verify:
  // n=3, m=3, deg={2,2,2}
  // adj is symmetric (each neighbor appears in both directions)
  // edges.size() == 3, all u < v
}

TEST_CASE("star K1_10") {
  // Hub=0, leaves=1..10
  // n=11, m=10
  // deg[0]=10, deg[1..10]=1
  // edges: {(0,1),(0,2),...,(0,10)}, all canonical
}

TEST_CASE("self-loop rejection") {
  // Input contains line "3 3" — should not appear in output edges
}

TEST_CASE("degree sum invariant") {
  // For any loaded graph: sum(deg) == 2*m
}

TEST_CASE("SNAP file ca-GrQc") {
  // Load actual file (skip if not present)
  // n == 5242, m == 14496  (known values for ca-GrQc)
}
```

---

### Part 2: DSpar score engine

**Files:** `src/dspar.hpp`, `src/dspar.cpp`, `tests/test_dspar.cpp`

```cpp
struct DsparScores {
    std::vector<double> score;   // score[i] = s(edges[i])
    double mean;                 // s̄ = mean of all scores
    // Structural metrics (computed w.r.t. a partition)
    double mu_intra = 0.0;
    double mu_inter = 0.0;
    double delta = 0.0;          // mu_intra - mu_inter
    int n1 = 0;                  // intra-community edge count
    int n2 = 0;                  // inter-community edge count
};

// Compute scores from graph degrees only (no partition needed)
DsparScores compute_scores(const Graph& g);

// Fill structural metrics given a partition
void compute_structural_metrics(
    DsparScores& scores,
    const Graph& g,
    const std::vector<int>& membership);
```

**Implementation:**
- `score[i] = 1.0/g.deg[g.edges[i].first] + 1.0/g.deg[g.edges[i].second]`
- `mean = accumulate(score) / m` — use Kahan summation for accuracy
- `mu_intra/mu_inter` — iterate over all edges, classify by membership, accumulate

**Tests:**

```
TEST_CASE("star graph edge score") {
  // Hub deg=10, leaf deg=1
  // Score of hub-leaf edge = 1/10 + 1/1 = 1.1  (exactly)
}

TEST_CASE("complete graph K5 uniform scores") {
  // All nodes have deg=4, all edges score = 1/4+1/4 = 0.5
  // mean = 0.5
}

TEST_CASE("mean invariant") {
  // sum(scores) / m == mean (to 1e-14)
}

TEST_CASE("reference match ca-GrQc") {
  // Compare mean score to Python output: 0.07845... (load reference value)
  // Tolerance: 1e-10
}

TEST_CASE("mu_intra > mu_inter for ca-GrQc") {
  // This is the DSpar separation property that drives the experiment
  // Load graph, run Leiden, compute structural metrics
  // Assert delta > 0
}
```

---

### Part 3: RNG and edge sampler

**Files:** `src/rng.hpp` (header-only), `src/sampler.hpp`, `src/sampler.cpp`,
`tests/test_sampler.cpp`

**xoshiro256++ RNG (header-only):**

```cpp
// src/rng.hpp
struct Xoshiro256pp {
    uint64_t s[4];

    explicit Xoshiro256pp(uint64_t seed) {
        // SplitMix64 to initialize state from seed
        // (standard initialization for xoshiro)
        uint64_t x = seed;
        for (int i = 0; i < 4; i++) {
            x += 0x9e3779b97f4a7c15ULL;
            uint64_t z = x;
            z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
            z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
            s[i] = z ^ (z >> 31);
        }
    }

    uint64_t next() {
        // xoshiro256++ step
        const uint64_t result = rotl(s[0] + s[3], 23) + s[0];
        const uint64_t t = s[1] << 17;
        s[2] ^= s[0]; s[3] ^= s[1];
        s[1] ^= s[2]; s[0] ^= s[3];
        s[2] ^= t;
        s[3] = rotl(s[3], 45);
        return result;
    }

    // uniform double in [0, 1)
    double next_double() {
        return (next() >> 11) * (1.0 / (1ULL << 53));
    }

private:
    static uint64_t rotl(uint64_t x, int k) {
        return (x << k) | (x >> (64 - k));
    }
};
```

**Sampler:**

```cpp
// src/sampler.hpp
// Returns keep[i] = true means edge i (canonical) is retained
// Probability: p(i) = alpha * scores.score[i] / scores.mean
// Clamped to [0, 1]: if alpha * score[i]/mean > 1, keep with probability 1
std::vector<bool> sample_edges(
    const DsparScores& scores,
    double alpha,
    Xoshiro256pp& rng);
```

**Seed scheme for reproducibility:**
Trial t at retention index j: `seed = base_seed + (uint64_t)t * 1000003ULL + (uint64_t)j`
This avoids correlation between trials and retention levels.

**Tests:**

```
TEST_CASE("alpha=1.0 keeps all edges") {
  // For uniform score graph (K5): all p(e) = 1.0
  // All keep[i] == true
}

TEST_CASE("expected retention rate") {
  // Run 10000 samples on a 100-edge graph
  // Mean |kept|/100 should be within 3σ of alpha
}

TEST_CASE("seed reproducibility") {
  // Same Xoshiro256pp state → same keep[] array
}

TEST_CASE("hub edges retained less often") {
  // Star graph K1_10: run 1000 trials at alpha=0.5
  // Hub edges (score = 1/10 + 1/1 = 1.1) vs peripheral...
  // Actually there are no peripheral-peripheral edges in star
  // Verify: mean retention of all edges ≈ alpha
}

TEST_CASE("p(e) <= 1 assertion") {
  // Construct graph where alpha * s(e)/mean > 1 for some edge
  // Verify no assertion fires (clamping works)
}
```

---

### Part 4: Modularity F and G computation

**Files:** `src/modularity.hpp`, `src/modularity.cpp`, `tests/test_modularity.cpp`

```cpp
struct ModResult {
    double F;             // intra-edge fraction
    double G;             // null-model penalty
    double Q;             // F - G
    int m_sparse;         // edge count in sparse graph
    int preserved_intra;  // intra edges kept
    int preserved_inter;  // inter edges kept
    double intra_rate;    // preserved_intra / n1
    double inter_rate;    // preserved_inter / n2
    double ratio_obs;     // inter_rate / intra_rate
};

// Compute modularity with fixed membership on (possibly sparse) graph.
// keep[i] == true means edge i is present (keep=nullptr means all edges kept)
ModResult compute_modularity(
    const Graph& g,
    const std::vector<bool>* keep,   // nullptr = use all edges
    const std::vector<int>& member,  // community ID for each node (0-indexed)
    int n1,                          // from DsparScores (intra count in original)
    int n2);                         // from DsparScores (inter count in original)
```

**Implementation (two-pass over edge list):**

Pass 1 — count F numerator and classify preserved edges:
```
intra_kept = 0, inter_kept = 0, m_sparse = 0
for i in 0..m-1:
    if keep == nullptr or keep[i]:
        m_sparse++
        u, v = edges[i]
        if member[u] == member[v]:
            intra_kept++
        else:
            inter_kept++
F = intra_kept / m_sparse   (or 0 if m_sparse == 0)
```

Pass 2 — accumulate community volumes using sparse degrees:
```
// Compute sparse degree of each node
vector<int> deg_sparse(n, 0)
for i in 0..m-1:
    if keep == nullptr or keep[i]:
        u, v = edges[i]
        deg_sparse[u]++
        deg_sparse[v]++

// Accumulate vol_c
int k = max(member) + 1   (number of communities)
vector<double> vol(k, 0.0)
for v in 0..n-1:
    vol[member[v]] += deg_sparse[v]

// Compute G
G = sum(vol[c]^2 for c in 0..k-1) / (4.0 * m_sparse * m_sparse)
```

**Critical implementation notes:**
- Use `double` accumulation for vol[] to avoid integer overflow on large graphs
- `m_sparse` must be computed as edge count (not directed edge count)
- When `keep == nullptr`, `deg_sparse == g.deg` — but compute it fresh for correctness

**Tests:**

```
TEST_CASE("two perfect cliques") {
  // Graph: K3 ∪ K3 (6 nodes, 2 communities of 3)
  // No inter-community edges
  // member = {0,0,0,1,1,1}
  // F = 1.0 (all edges intra)
  // vol[0] = vol[1] = 2*3 = 6  (each node has deg 2 in own clique)
  // G = (36+36)/(4*9) = 72/36 = 2.0  -- Wait, m=6 (3+3 edges)
  // Actually: vol[0] = sum deg in sparse = 2+2+2=6, vol[1]=6
  // G = (6^2 + 6^2)/(4*6^2) = 72/144 = 0.5
  // Q = 1.0 - 0.5 = 0.5
}

TEST_CASE("identity holds for 100 random masks") {
  // Load karate graph, run Leiden, get membership
  // For 100 random keep[] masks:
  //   orig = compute_modularity(g, nullptr, member)
  //   sparse = compute_modularity(g, &keep, member)
  //   dQ = sparse.Q - orig.Q
  //   dF = sparse.F - orig.F
  //   dG = sparse.G - orig.G
  //   assert |dQ - (dF - dG)| < 1e-12
}

TEST_CASE("match Python reference ca-GrQc") {
  // Load graph + saved Python membership + saved Python F/G/Q values
  // Compare to 1e-10 tolerance
}

TEST_CASE("keep=nullptr matches keep=all-true") {
  // compute_modularity(g, nullptr, m) == compute_modularity(g, &all_true, m)
}
```

---

### Part 5: Leiden community detection via igraph C

**Files:** `src/leiden.hpp`, `src/leiden.cpp`, `tests/test_leiden.cpp`

igraph ships an ANSI C library. Link with `find_package(igraph REQUIRED)`.

```cpp
struct LeidenResult {
    std::vector<int> membership;  // membership[v] = community ID, 0-indexed
    double modularity;
    int n_communities;
};

// Build igraph from graph with optional edge mask, run Leiden
LeidenResult run_leiden(
    const Graph& g,
    const std::vector<bool>* keep = nullptr);
    // keep=nullptr → full graph
    // keep[i]=true → edge i included
```

**Implementation:**

```cpp
LeidenResult run_leiden(const Graph& g, const std::vector<bool>* keep) {
    // 1. Count edges to include
    std::vector<igraph_integer_t> edge_list;
    for (int i = 0; i < g.m; i++) {
        if (!keep || (*keep)[i]) {
            edge_list.push_back(g.edges[i].first);
            edge_list.push_back(g.edges[i].second);
        }
    }

    // 2. Build igraph graph
    igraph_t ig;
    igraph_vector_int_t edges_vec;
    igraph_vector_int_init_array(&edges_vec,
        edge_list.data(), edge_list.size());
    igraph_create(&ig, &edges_vec, g.n, IGRAPH_UNDIRECTED);
    igraph_vector_int_destroy(&edges_vec);

    // 3. Run Leiden
    igraph_vector_int_t membership_vec;
    igraph_vector_int_init(&membership_vec, 0);
    igraph_real_t modularity_val;
    igraph_community_leiden(
        &ig,
        /*edge_weights=*/nullptr,
        /*node_weights=*/nullptr,
        /*resolution_parameter=*/1.0,
        /*beta=*/0.01,
        /*start=*/false,
        /*n_iterations=*/2,
        &membership_vec,
        /*nb_clusters=*/nullptr,
        &modularity_val);

    // 4. Extract result
    LeidenResult result;
    result.n_communities = 0;
    result.membership.resize(g.n);
    for (int v = 0; v < g.n; v++) {
        result.membership[v] = VECTOR(membership_vec)[v];
        result.n_communities = std::max(
            result.n_communities, result.membership[v] + 1);
    }
    result.modularity = (double)modularity_val;

    igraph_vector_int_destroy(&membership_vec);
    igraph_destroy(&ig);
    return result;
}
```

**igraph version note:** Target igraph ≥ 0.10. The `igraph_community_leiden` signature
changed between 0.9 and 0.10 — check installed version with `igraph_version(nullptr, nullptr, nullptr)`.

**Tests:**

```
TEST_CASE("karate club community count and modularity") {
  // Q > 0.35, 2 communities (well-known result)
}

TEST_CASE("barbell graph perfect split") {
  // Two cliques of 5 connected by single bridge edge
  // Leiden should find exactly 2 communities
  // Run 10 times (Leiden is randomized), all should give 2 communities
}

TEST_CASE("igraph Q matches our Q") {
  // run_leiden on karate → membership
  // compute_modularity(g, nullptr, membership) → our Q
  // leiden_result.modularity → igraph Q
  // |our_Q - igraph_Q| < 1e-6  (small numerical diff expected)
}

TEST_CASE("all nodes assigned") {
  // For any graph + keep mask, all membership[v] in [0, n_communities)
}
```

---

### Part 6: Single-dataset experiment loop and CSV output

**Files:** `src/experiment.hpp`, `src/experiment.cpp`, `src/csv_writer.hpp`,
`tests/test_experiment.cpp`

```cpp
struct ExperimentConfig {
    std::vector<double> alphas = {0.5, 0.6, 0.7, 0.8, 0.9};
    int n_trials = 10;
    uint64_t base_seed = 42;
    bool run_leiden_resample = true;
};

struct TrialResult {
    // From config
    double retention;
    int trial;
    uint64_t seed;
    // Graph stats
    int m, n1, n2;
    // DSpar structural metrics (constant per dataset)
    double mu_intra, mu_inter, delta;
    // Original graph modularity
    double F_original, G_original, modularity_fixed_original;
    // Sparse graph (fixed membership)
    int m_sparse, preserved_intra, preserved_inter;
    double F_observed, G_sparse_observed;
    double intra_rate, inter_rate, ratio_observed;
    double modularity_fixed_sparse;
    // Derived changes
    double F_improvement_observed;  // ΔF
    double dG_observed;             // ΔG
    double modularity_fixed_change; // ΔQ_fixed
    double dQ_reconstructed;        // ΔF - ΔG (must == modularity_fixed_change)
    // Leiden re-optimization
    double modularity_leiden_sparse;
    double modularity_leiden_change;
    // Verification
    double identity_error;          // |dQ_fixed - dQ_reconstructed|
};

// Run complete experiment for one dataset
std::vector<TrialResult> run_experiment(
    const Graph& g,
    const std::string& dataset_name,
    const ExperimentConfig& cfg);

// Write raw CSV (one row per trial)
void write_raw_csv(
    const std::vector<TrialResult>& results,
    const std::string& path);

// Write summary CSV (mean/std aggregated per retention)
void write_summary_csv(
    const std::vector<TrialResult>& results,
    const std::string& path);
```

**Experiment loop (implement exactly this order):**

```
1. Compute DsparScores for full graph (once)
2. Run Leiden on full graph → membership_fixed (once)
3. Compute ModResult on full graph with membership_fixed → orig (once)
4. Compute structural metrics (mu_intra, mu_inter, delta) (once)
5. For each alpha in cfg.alphas:
     For each trial t in 0..n_trials-1:
       seed = cfg.base_seed + (uint64_t)t * 1000003ULL + alpha_index * 7919ULL
       rng = Xoshiro256pp(seed)
       keep = sample_edges(scores, alpha, rng)
       sparse_mod = compute_modularity(g, &keep, membership_fixed, n1, n2)
       leiden_mod = 0.0
       if cfg.run_leiden_resample:
           leiden_result = run_leiden(g, &keep)
           leiden_mod = leiden_result.modularity
       record TrialResult
```

**CSV writer (`src/csv_writer.hpp`):**

Simple header-only CSV writer. Write header on first call, then one row per trial.
Use `std::fixed << std::setprecision(15)` for all floating-point values.

**Tests:**

```
TEST_CASE("identity holds for all karate trials") {
  // Run experiment on karate graph (34 nodes, 78 edges)
  // For all TrialResult: identity_error < 1e-12
}

TEST_CASE("alpha=1.0 gives zero change") {
  // Add alpha=1.0 to config
  // All dQ_fixed == 0.0 exactly
  // All m_sparse == m
}

TEST_CASE("CSV columns match Python reference") {
  // Load Python reference CSV for ca-GrQc
  // Compare header line exactly (column names as string)
}

TEST_CASE("seed reproducibility") {
  // Run twice with same config → identical TrialResult vectors
}

TEST_CASE("dQ_leiden >= dQ_fixed majority") {
  // Over all karate trials: assert at least 80% satisfy this
  // (Paper claim: always holds empirically)
}
```

---

### Part 7: Multi-dataset runner, CLI, TOML config

**Files:** `src/main.cpp`, `src/config.hpp`, `src/config.cpp`,
`src/downloader.hpp`, `src/downloader.cpp`, `configs/exp1.toml`

**TOML config format:**

```toml
# configs/exp1.toml
[experiment]
n_trials    = 10
base_seed   = 42
alphas      = [0.5, 0.6, 0.7, 0.8, 0.9]
output_dir  = "results/exp1"
run_leiden_resample = true

[[datasets]]
name = "ca-GrQc"
url  = "https://snap.stanford.edu/data/ca-GrQc.txt.gz"

[[datasets]]
name = "ca-HepTh"
url  = "https://snap.stanford.edu/data/ca-HepTh.txt.gz"

[[datasets]]
name = "ca-HepPh"
url  = "https://snap.stanford.edu/data/ca-HepPh.txt.gz"

[[datasets]]
name = "ca-AstroPh"
url  = "https://snap.stanford.edu/data/ca-AstroPh.txt.gz"

[[datasets]]
name = "ca-CondMat"
url  = "https://snap.stanford.edu/data/ca-CondMat.txt.gz"

[[datasets]]
name = "cit-HepTh"
url  = "https://snap.stanford.edu/data/cit-HepTh.txt.gz"

[[datasets]]
name = "cit-HepPh"
url  = "https://snap.stanford.edu/data/cit-HepPh.txt.gz"

[[datasets]]
name = "email-Enron"
url  = "https://snap.stanford.edu/data/email-Enron.txt.gz"

[[datasets]]
name = "facebook-combined"
url  = "https://snap.stanford.edu/data/facebook_combined.txt.gz"

[[datasets]]
name = "wiki-Vote"
url  = "https://snap.stanford.edu/data/wiki-Vote.txt.gz"

[[datasets]]
name = "email-Eu-core"
url  = "https://snap.stanford.edu/data/email-Eu-core.txt.gz"
```

**CLI interface:**

```
./exp1 [options]

Options:
  --config PATH         Path to TOML config [default: configs/exp1.toml]
  --datasets D1,D2,...  Run only these datasets (comma-separated)
  --seed N              Override base_seed from config
  --output DIR          Override output_dir from config
  --trials N            Override n_trials from config
  --resume              Skip datasets that already have output files
  --dry-run             Print config and exit without running
  --no-leiden           Skip Leiden re-optimization (faster, missing leiden columns)
  --help                Show this help
```

**Downloader:**

```cpp
// src/downloader.hpp
// Download URL to local path using libcurl
// Decompress .gz files using zlib
// Return true on success
bool download_and_decompress(
    const std::string& url,
    const std::string& local_path);  // local_path = path to .txt (decompressed)

// Check if dataset exists locally, download if not
// Returns path to decompressed .txt file
std::string ensure_dataset(
    const std::string& name,
    const std::string& url,
    const std::string& data_dir = "data");
```

**Main loop:**

```
parse CLI args
load TOML config
apply CLI overrides
if dry-run: print config, exit 0

create output_dir if not exists

for each dataset in config:
    if --datasets filter and dataset not in filter: skip
    if --resume and output file exists: skip

    print "[dataset] Downloading..."
    path = ensure_dataset(name, url, "data")

    print "[dataset] Loading graph..."
    g = load_snap(path)
    print "  n=%d m=%d\n", g.n, g.m

    print "[dataset] Running experiment..."
    results = run_experiment(g, name, cfg)

    write_raw_csv(results, output_dir + "/" + name + "_theoretical_validation_FIXED.csv")
    write_summary_csv(results, output_dir + "/" + name + "_summary_FIXED.csv")
    print "[dataset] Done. Wrote CSV files.\n"
```

**Progress reporting:** print one line per trial or per alpha block. Include wall-clock
time for each dataset.

---

### Part 8: OpenMP parallelism

Modify `run_experiment` in `src/experiment.cpp`.

The trial loop `for t in 0..n_trials` is embarrassingly parallel. Each trial uses:
- Its own `Xoshiro256pp` RNG (thread-local by construction)
- Its own `std::vector<bool> keep` buffer
- Read-only access to `g`, `scores`, `membership_fixed`

igraph is **not thread-safe for graph construction**. Protect `run_leiden` calls with
a mutex.

```cpp
// In run_experiment():
std::vector<TrialResult> results(cfg.n_trials * cfg.alphas.size());
std::mutex leiden_mtx;

#pragma omp parallel for schedule(dynamic) collapse(2)
for (int ai = 0; ai < (int)cfg.alphas.size(); ai++) {
    for (int t = 0; t < cfg.n_trials; t++) {
        double alpha = cfg.alphas[ai];
        uint64_t seed = cfg.base_seed
            + (uint64_t)t * 1000003ULL
            + (uint64_t)ai * 7919ULL;

        Xoshiro256pp rng(seed);
        auto keep = sample_edges(scores, alpha, rng);
        auto mod = compute_modularity(g, &keep, membership_fixed,
                                      scores.n1, scores.n2);

        double leiden_q = 0.0;
        if (cfg.run_leiden_resample) {
            std::lock_guard<std::mutex> lk(leiden_mtx);
            leiden_q = run_leiden(g, &keep).modularity;
        }

        int idx = ai * cfg.n_trials + t;
        results[idx] = build_trial_result(alpha, t, seed, orig, mod, leiden_q, scores);
    }
}
```

**Compile flags:** Add `-fopenmp` to CMakeLists.txt. Respect `OMP_NUM_THREADS` at runtime.

**Tests:**

```
TEST_CASE("parallel matches serial") {
  // Run experiment with OMP_NUM_THREADS=1 and OMP_NUM_THREADS=4
  // All TrialResult fields identical (same seed scheme)
}

TEST_CASE("no data races") {
  // Instruction: compile with -fsanitize=thread and run under TSAN
  // All tests must pass with zero TSAN reports
}
```

---

### Part 9: Performance optimizations

Implement these after Part 8 passes all tests. Each optimization must leave all tests
passing.

**Opt A — thread-local keep buffer (eliminate per-trial allocation)**

Instead of `auto keep = sample_edges(...)` allocating a new vector each time:

```cpp
// In run_experiment, before parallel region:
// (thread_local storage is per-thread, initialized once)
// Use a per-thread pre-allocated buffer, reset before each use

// Alternative: pass buffer by reference, let caller manage
void sample_edges_inplace(
    const DsparScores& scores,
    double alpha,
    Xoshiro256pp& rng,
    std::vector<bool>& keep_buf);  // caller pre-allocated to size m
```

Preallocate one buffer per thread in the parallel region setup.

**Opt B — precompute per-alpha thresholds (eliminate per-trial division)**

For each alpha, precompute `threshold[i] = alpha * scores.score[i] / scores.mean`
once outside the trial loop. The inner sampling loop becomes:

```cpp
// Precomputed once per alpha:
std::vector<double> threshold(m);
for (int i = 0; i < m; i++)
    threshold[i] = alpha * scores.score[i] / scores.mean;

// Per trial (no division):
for (int i = 0; i < m; i++)
    keep[i] = rng.next_double() < threshold[i];
```

**Opt C — Kahan summation in G computation**

The G term sums vol_c² values which can be large. Use Kahan summation:

```cpp
double g_sum = 0.0, comp = 0.0;
for (int c = 0; c < k; c++) {
    double term = vol[c] * vol[c];
    double y = term - comp;
    double t2 = g_sum + y;
    comp = (t2 - g_sum) - y;
    g_sum = t2;
}
double G = g_sum / (4.0 * (double)m_sparse * (double)m_sparse);
```

**Opt D — prefix-sum sampler for aggressive sparsification (α < 0.3)**

Build CDF once per graph: `cdf[i] = Σ_{j≤i} score[j]`.
To retain exactly `k = round(alpha * m)` edges, sample k values uniformly from
`[0, total_score]` and binary-search in CDF. This is faster than m Bernoulli tests
when most edges are rejected.

Switch at runtime:
```cpp
if (alpha < 0.3 && m > 100000)
    sample_edges_prefix_sum(scores, alpha, rng, keep);
else
    sample_edges_bernoulli(threshold, rng, keep);
```

---

## 6. CMakeLists.txt

```cmake
cmake_minimum_required(VERSION 3.20)
project(dspar_exp1 CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Release build by default
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Compiler flags
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native")
set(CMAKE_CXX_FLAGS_DEBUG   "-g -O0 -fsanitize=address,undefined")

# Dependencies
find_package(igraph REQUIRED)   # sudo apt install libigraph-dev
find_package(OpenMP REQUIRED)   # sudo apt install libomp-dev
find_package(CURL REQUIRED)     # sudo apt install libcurl4-openssl-dev
find_package(ZLIB REQUIRED)     # for .gz decompression

# toml11: header-only, must be in extern/toml11
# Clone: git clone https://github.com/ToruNiina/toml11.git extern/toml11
include_directories(extern/toml11)

# Main binary
add_executable(exp1
    src/main.cpp
    src/graph.cpp
    src/dspar.cpp
    src/sampler.cpp
    src/modularity.cpp
    src/leiden.cpp
    src/experiment.cpp
    src/config.cpp
    src/downloader.cpp)

target_link_libraries(exp1
    PRIVATE
    igraph::igraph
    OpenMP::OpenMP_CXX
    CURL::libcurl
    ZLIB::ZLIB)

target_include_directories(exp1 PRIVATE src)

# Tests (Catch2)
include(FetchContent)
FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.5.0)
FetchContent_MakeAvailable(Catch2)

add_executable(tests
    tests/test_main.cpp
    tests/test_graph.cpp
    tests/test_dspar.cpp
    tests/test_sampler.cpp
    tests/test_modularity.cpp
    tests/test_leiden.cpp
    tests/test_experiment.cpp
    # Source files needed by tests
    src/graph.cpp
    src/dspar.cpp
    src/sampler.cpp
    src/modularity.cpp
    src/leiden.cpp
    src/experiment.cpp)

target_link_libraries(tests
    PRIVATE
    Catch2::Catch2WithMain
    igraph::igraph
    OpenMP::OpenMP_CXX)

target_include_directories(tests PRIVATE src)

enable_testing()
add_test(NAME unit_tests COMMAND tests)
```

---

## 7. Datasets Reference

All 11 datasets from the SNAP repository. The table shows expected graph size after
loading (deduplication, remapping, removing self-loops):

| Dataset | n (nodes) | m (edges) | Notes |
|---|---|---|---|
| ca-GrQc | 5,242 | 14,496 | Smallest, fastest for iteration |
| ca-HepTh | 9,877 | 25,973 | |
| ca-HepPh | 12,008 | 118,489 | |
| ca-AstroPh | 18,772 | 198,050 | |
| ca-CondMat | 23,133 | 93,439 | |
| cit-HepTh | 27,770 | 352,285 | Directed in SNAP → make undirected |
| cit-HepPh | 34,546 | 420,877 | Directed in SNAP → make undirected |
| email-Enron | 36,692 | 183,831 | |
| facebook-combined | 4,039 | 88,234 | |
| wiki-Vote | 7,115 | 100,762 | Directed → make undirected |
| email-Eu-core | 1,005 | 16,064 | Smallest number of nodes |

**Directed graphs:** cit-HepTh, cit-HepPh, wiki-Vote are directed in SNAP format.
The loader already handles this: treating all edges as undirected (for each line
"src dst", add edge (min,max)), deduplicating. This matches the Python `load_snap_dataset`
behavior.

---

## 8. Validation Against Python Reference

After implementing Part 6, run this validation on `ca-GrQc` (the smallest dataset):

1. Run Python experiment (from the original repo):
   ```bash
   cd PAPER_EXPERIMENTS
   python exp1_2_theoretical_predictions.py ca-GrQc
   ```
   This writes `results/exp1_2_theoretical/ca-GrQc_theoretical_validation_FIXED.csv`

2. Run C++ experiment:
   ```bash
   ./exp1 --config configs/exp1.toml --datasets ca-GrQc --seed 42
   ```

3. Compare CSVs with the validation script `scripts/validate.py`:

```python
#!/usr/bin/env python3
"""scripts/validate.py — compare C++ output to Python reference"""
import pandas as pd
import sys

py_csv  = sys.argv[1]  # Python reference
cpp_csv = sys.argv[2]  # C++ output

py  = pd.read_csv(py_csv)
cpp = pd.read_csv(cpp_csv)

# Check columns present
for col in ['retention','F_original','G_original','modularity_fixed_change',
            'dQ_reconstructed','F_improvement_observed','dG_observed',
            'modularity_leiden_change']:
    assert col in cpp.columns, f"Missing column: {col}"

# Aggregate by retention and compare means
for alpha in sorted(py['retention'].unique()):
    py_row  = py[abs(py['retention']-alpha) < 1e-9]
    cpp_row = cpp[abs(cpp['retention']-alpha) < 1e-9]

    for metric in ['modularity_fixed_change', 'F_improvement_observed', 'dG_observed']:
        py_mean  = py_row[metric].mean()
        cpp_mean = cpp_row[metric].mean()
        diff = abs(py_mean - cpp_mean)
        status = "PASS" if diff < 0.005 else "FAIL"
        print(f"alpha={alpha:.1f} {metric}: py={py_mean:.6f} cpp={cpp_mean:.6f} diff={diff:.2e} [{status}]")

# Identity check (must be EXACT)
max_err = cpp['dQ_reconstructed'].sub(cpp['modularity_fixed_change']).abs().max()
print(f"\nMax identity error: {max_err:.2e} (must be < 1e-12)")
assert max_err < 1e-12, "IDENTITY ERROR TOO LARGE"
print("All checks passed.")
```

**Acceptable tolerances:**
- `modularity_fixed_change` mean per alpha: within ±0.003 (Leiden randomness)
- `F_improvement_observed` mean per alpha: within ±0.002
- `dG_observed` mean per alpha: within ±0.002
- `|dQ_reconstructed - modularity_fixed_change|`: must be < 1e-12 for ALL rows

---

## 9. Build and Run Instructions

```bash
# 0. Install system dependencies
sudo apt install libigraph-dev libcurl4-openssl-dev zlib1g-dev libomp-dev

# 1. Clone toml11 (header-only TOML parser)
git clone https://github.com/ToruNiina/toml11.git extern/toml11

# 2. Build
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
cd ..

# 3. Run tests
./build/tests

# 4. Quick smoke test (single small dataset)
./build/exp1 --config configs/exp1.toml --datasets ca-GrQc --seed 42

# 5. Full 11-dataset run (matching paper)
./build/exp1 --config configs/exp1.toml

# 6. Re-run with different seed for independent replication
./build/exp1 --config configs/exp1.toml --seed 123 --output results/run_seed123

# 7. Resume interrupted run
./build/exp1 --config configs/exp1.toml --resume

# 8. Parallel with 4 cores
OMP_NUM_THREADS=4 ./build/exp1 --config configs/exp1.toml

# 9. Debug build (with ASAN/UBSAN)
mkdir build_debug && cd build_debug
cmake .. -DCMAKE_BUILD_TYPE=Debug
make -j$(nproc)
./tests
```

---

## 10. Common Pitfalls to Avoid

**Pitfall 1 — Degree from original graph for scores, degree from sparse graph for G:**
`score[i] = 1/deg_original[u] + 1/deg_original[v]` always uses original degrees.
But `vol_c` for the G term uses `deg_sparse[v]` (how many retained edges touch v).
Mixing these is the most common bug.

**Pitfall 2 — igraph community IDs are not guaranteed to be 0..k-1 contiguous:**
In igraph ≥ 0.10 they are renumbered contiguous. In older versions they may not be.
After extracting membership, always renumber to 0..k-1:
```cpp
std::map<int,int> remap; int next=0;
for (auto& c : result.membership)
    if (!remap.count(c)) remap[c] = next++;
for (auto& c : result.membership) c = remap[c];
result.n_communities = next;
```

**Pitfall 3 — The identity test fails with error ~1e-16 to 1e-14:**
This is fine — that is machine precision. The Python code asserts `< 1e-10`. Use
`1e-12` as the C++ threshold (doubles are exact in the algebraic identity; small
errors come from floating-point order of operations).

**Pitfall 4 — F uses edge count, not directed edge count:**
`m_sparse` is the number of undirected edges kept. F = intra_kept / m_sparse.
If you accidentally double-count edges (both (u,v) and (v,u)), F will be wrong
but Q may coincidentally be close — the identity check will catch this.

**Pitfall 5 — igraph thread safety:**
`igraph_create` and `igraph_community_leiden` share global state for random number
generation. Wrap all igraph calls in a mutex even if they look read-only.

**Pitfall 6 — SNAP directed graph files:**
cit-HepTh, cit-HepPh, wiki-Vote have lines like "FromNodeId ToNodeId". Treat them
as undirected: insert (min, max) into the edge set. The Python code does this via
`edge_set.add((min(s,d), max(s,d)))`.

**Pitfall 7 — CSV floating-point precision:**
Write all floats with `std::setprecision(15)`. The Python reference uses ~15 significant
digits. Lower precision will cause the validation script to show larger diffs than the
true algorithmic differences.

---

## 11. Implementation Order and Testing Checkpoints

Implement in exactly this order. Do not proceed to the next part until all tests for
the current part pass.

```
Part 1 → Part 2 → Part 3 → Part 4 → Part 5 → Part 6 → Part 7 → Part 8 → Part 9
  ↓         ↓         ↓         ↓         ↓         ↓
test_graph  test_dspar test_sampler test_modularity test_leiden test_experiment
```

**After Part 6:** Run the Python validation script on ca-GrQc. Fix any discrepancies
before proceeding.

**After Part 7:** Do a full 11-dataset dry-run to verify config parsing and dataset
discovery work correctly.

**After Part 8:** Run tests under ThreadSanitizer (`cmake -DCMAKE_BUILD_TYPE=Debug
-DCMAKE_CXX_FLAGS="-fsanitize=thread"`). Fix any races.

**After Part 9:** Benchmark email-Enron (183K edges, 10 trials, 5 alphas) with
`OMP_NUM_THREADS=1` vs `OMP_NUM_THREADS=4`. Document the speedup.
