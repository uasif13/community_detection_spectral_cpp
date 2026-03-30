# community_detection_spectral_cpp

### C++ Implementation

We are interested in writing a sequential and parallel version of DSpar. The four experiments are listed below.

All have compute_dspar_scores and dspar_sparsify

1. Modularity Decomposition

exp1_2_theoretical_predictions.py

load graph
run leiden on original G
for [0.1 ... 1.0]: (run_single_trial)
    dspar_sparsify
    compute q_sparse on g_sparse with membership_fixed (= delta Q_fixed)
    compute F and G terms separately -> verify delta Q = delta F - delta G
    rerun leiden on G_sparse -> Q_leiden (=delta Q_Leiden)

2. LFR Benchmark

exp_1_3_lfr_analysis.py

3. Ground-truth recovery

exp2_ground_truth_recovery.py

4. Scalability

exp3_scalability.py
