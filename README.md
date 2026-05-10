# ⚡ Monophasic Power Flow — Newton-Raphson Method

> Implementation of a single-phase power flow solver (balanced power system) in Python using Jupyter Notebook, validated against CEPEL's ANAREDE software.

---

## 📘 Overview

This project implements a **monophasic (single-phase) power flow** algorithm for balanced electrical power systems (SEP), developed as part of an academic assignment. The solver is written in **Python** within a `.ipynb` (Jupyter Notebook) environment and validated against results from the industry-standard **ANAREDE** software by CEPEL.

Two benchmark systems are simulated:

- **4-Bus System** — with full step-by-step analytical solution
- **9-Bus System** — used for broader numerical validation

---

## 🎯 Objectives

- Implement a power flow algorithm capable of solving balanced monophasic power systems
- Simulate and validate the **4-bus** and **9-bus** test systems
- Compare computed results with ANAREDE reference outputs (maximum and minimum errors)
- Discuss and explain sources of numerical discrepancy
- Draw conclusions about solver accuracy and performance

---

## 🗂️ Repository Structure

```
.
├── power_flow.ipynb        # Main Jupyter Notebook with full implementation
├── data/
│   ├── bus4.py / .csv      # 4-bus system input data
│   └── bus9.py / .csv      # 9-bus system input data
├── results/
│   ├── bus4_results.csv    # Computed results for the 4-bus system
│   └── bus9_results.csv    # Computed results for the 9-bus system
├── report/
│   └── report.pdf          # Technical report (analytical solution + comparisons)
└── README.md
```

---

## ⚙️ Methodology

The power flow problem is solved using the **Newton-Raphson** iterative method applied to the nodal power equations:

$$P_i = \sum_{k=1}^{n} |V_i||V_k|\left(G_{ik}\cos\theta_{ik} + B_{ik}\sin\theta_{ik}\right)$$

$$Q_i = \sum_{k=1}^{n} |V_i||V_k|\left(G_{ik}\sin\theta_{ik} - B_{ik}\cos\theta_{ik}\right)$$

Key steps:
1. Build the **admittance matrix** Y-bus from network data
2. Initialize bus voltages (flat start: |V| = 1.0 pu, θ = 0°)
3. Compute power mismatches ΔP and ΔQ
4. Assemble and solve the **Jacobian** matrix
5. Update voltage magnitudes and angles
6. Repeat until convergence (mismatch < tolerance)

---

## 🖥️ Requirements

- Python 3.8+
- Jupyter Notebook or JupyterLab
- NumPy
- Pandas
- Matplotlib (for plots)

Install dependencies with:

```bash
pip install numpy pandas matplotlib jupyter
```

---

## 🚀 Running the Notebook

```bash
git clone https://github.com/your-username/your-repo-name.git
cd your-repo-name
jupyter notebook power_flow.ipynb
```

Run all cells sequentially. Results and comparison tables are generated inline within the notebook.

---

## 📊 Results Summary

### 4-Bus System

| Bus | V (pu) — Computed | V (pu) — ANAREDE | Error (%) |
|-----|-------------------|------------------|-----------|
| 1   | —                 | —                | —         |
| 2   | —                 | —                | —         |
| 3   | —                 | —                | —         |
| 4   | —                 | —                | —         |

### 9-Bus System

| Bus | V (pu) — Computed | V (pu) — ANAREDE | Error (%) |
|-----|-------------------|------------------|-----------|
| ...   | —               | —                | —         |

> Fill in the result tables with your computed values after running the notebook.

---

## 📋 Report

The technical report (located in `/report/`) covers:

1. **Step-by-step analytical solution** of the 4-bus system (hand calculations)
2. **Comparison tables** for both systems — maximum and minimum errors vs. ANAREDE
3. **Discussion** of error sources (convergence tolerance, data rounding, modeling differences)
4. **Conclusions** on solver accuracy and applicability

---

## 🔍 Error Analysis

Differences between computed results and ANAREDE outputs are expected due to:

- Numerical precision and convergence threshold settings
- Rounding in input impedance/admittance data
- Differences in transformer and shunt modeling conventions
- Iteration count limits

---

## 📚 References

- **ANAREDE** — Power Flow Analysis Program, CEPEL (Centro de Pesquisas de Energia Elétrica)
- Stevenson, W. D. — *Elements of Power System Analysis*
- Glover, Sarma & Overbye — *Power Systems Analysis and Design*
- Bergen & Vittal — *Power Systems Analysis*

---

## 👤 Author

**Renan LARRIEU**
Power Systems Course — [Rio de Janeiro State University]
[2022] · [Power Systems Analysis II]

---

## 📜 License

This project is intended for academic use only.
