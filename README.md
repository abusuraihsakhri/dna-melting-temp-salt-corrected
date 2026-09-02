# DNA Melting Temp Salt Corrected

> **Domain:** Clinical Decision Support & Biomedical Computing  
> **Reference Guidelines & Standards:** `Standard Clinical Formulations & ISO/IEC Quality Frameworks`

<div align="center">

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/Python-3.10%20%7C%203.11%20%7C%203.12-3776AB.svg?logo=python&logoColor=white)
![FastAPI](https://img.shields.io/badge/FastAPI-0.111-009688.svg?logo=fastapi&logoColor=white)
![Audit Trail](https://img.shields.io/badge/Audit-HMAC--SHA256_Tamper--Evident-brightgreen.svg)
![Zero-PHI Guard](https://img.shields.io/badge/Guard-Zero--PHI_Outbound-blue.svg)
![Docker](https://img.shields.io/badge/Docker-Ready-2496ED.svg?logo=docker&logoColor=white)

</div>

---

## 📖 What It Does

Mg2+ Correction Extension for Oligo Melting Temperature
Implements Owczarzy 2008 formula for divalent cation correction and
Aho formula for combined monovalent + divalent salt effective concentration.

References:
  - Owczarzy et al. (2008) Biochemistry 47(19):5336-5353
  - Aho et al. (2008) Nucleic Acids Res 36(11):e63

Multi-Primer Pair Comparison Engine
Compare and rank multiple primer pairs by Tm matching, GC balance,
and secondary structure absence.

---

## ⚙️ Key Capabilities & Algorithmic Modules

### 🔬 Analytical Functions

- **`compute_effective_salt_concentration()`**: Compute effective monovalent cation concentration using the Aho formula.

When Mg2+ is present, the effective monovalent salt concentration is:
    ratio = sqrt(oligo_conc_nM * 1e-9) / (3.0 * mg_mM * 1e-3)
    if ratio < 0.22:
        Na_eq = 3.79 * sqrt(oligo_conc_nM * 1e-9) / (3.0 * mg_mM * 1e-3)
        Na_eq = Na_eq * Na_eq * 3.0 * mg_mM * 1e-3
    else:
        Na_eq = Na + 140 * Mg (simplified)

Args:
    na_mM: Monovalent cation (Na+/K+) concentration in mM.
    mg_mM: Divalent cation (Mg2+) concentration in mM.
    oligo_conc_nM: Oligonucleotide concentration in nM.

Returns:
    Effective monovalent salt concentration in M.
- **`owczarzy_2004_salt_correction()`**: Apply Owczarzy 2004 monovalent cation salt correction.

ΔTm = (4.29 * f_GC - 3.95) * ln(Na+) + 0.940 * (ln(Na+))^2

Args:
    tm_nn: Nearest-neighbor Tm without salt correction (°C).
    na_M: Monovalent cation concentration in M.
    f_GC: Fraction of G/C bases (0.0 - 1.0).

Returns:
    Salt-corrected Tm in °C.
- **`owczarzy_2008_mg_correction()`**: Apply Owczarzy 2008 Mg2+ correction to nearest-neighbor Tm.

For Mg2+ containing buffers:
    a = 3.92 * (0.843 - 0.352 * sqrt(Na+) * ln(Na+))
    d = 1.42 * (1.279 - 4.03e-3 * ln(Na+) - 8.03e-3 * (ln(Na+))^2)
    g = 8.31 * (0.486 - 0.258 * ln(Na+) + 5.25e-3 * (ln(Na+))^3)

    1/Tm(Mg) = 1/Tm(Na) + a + b*ln(Mg) + f_GC*(c + d*ln(Mg)) + 1/(2*(Nbp-1)) * (e + f*ln(Mg) + g*(ln(Mg))^2)

Args:
    tm_nn: Nearest-neighbor Tm with Na+ correction (°C).
    na_M: Monovalent cation concentration in M.
    mg_M: Divalent Mg2+ concentration in M.
    f_GC: Fraction of G/C bases (0.0 - 1.0).

Returns:
    Mg2+-corrected Tm in °C.
- **`calculate_tm_with_mg2()`**: Calculate Tm with full salt and Mg2+ correction.

Args:
    sequence: DNA primer sequence (5' to 3').
    na_mM: Na+ concentration in mM.
    mg_mM: Mg2+ concentration in mM.
    oligo_conc_nM: Oligonucleotide concentration in nM.
    nn_table: Optional nearest-neighbor parameter table.

Returns:
    Dictionary with Tm values and correction details.
- **`evaluate_primer_pair()`**: Evaluate a single primer pair with all metrics.

Args:
    forward: Forward primer sequence.
    reverse: Reverse primer sequence.
    na_mM: Na+ concentration in mM.
    mg_mM: Mg2+ concentration in mM.
    oligo_conc_nM: Oligo concentration in nM.

Returns:
    Dictionary with all evaluation metrics.

---

## 📐 Mathematical Formulation & Logic

```text
  Implements Owczarzy 2008 formula for divalent cation correction and
  Aho formula for combined monovalent + divalent salt effective concentration.
  """Compute effective monovalent cation concentration using the Aho formula.
  """Calculate Tm with full salt and Mg2+ correction.
  result = calculate_tm_with_mg2("ATCGATCGATCGATCGATCG", na_mM=50, mg_mM=1.5)
```

---

## 💻 CLI Quickstart & Usage

### 1. Guided Interactive Mode
```bash
python cli.py
```

### 2. Direct Parameterized Evaluation
```bash
python cli.py --task-id <value> --target <value> --primary <value> --secondary <value>
```

### Parameter Reference
- `--task-id`: Specifies input measurement or parameter value.
- `--target`: Specifies input measurement or parameter value.
- `--primary`: Specifies input measurement or parameter value.
- `--secondary`: Specifies input measurement or parameter value.
- `--critical`: Specifies input measurement or parameter value.
- `--status`: Specifies input measurement or parameter value.
- `--input`: Specifies input measurement or parameter value.
- `--output`: Specifies input measurement or parameter value.

### Input Data Schema

| Field | Description | Requirement |
|:------|:------------|:------------|
| `Patient_ID` | Parameter / observation metric | Required |
| `v1` | Parameter / observation metric | Required |
| `v2` | Parameter / observation metric | Required |
| `v3` | Parameter / observation metric | Required |

---

## 🛡️ Security & Enterprise Architecture

* **Zero-PHI Outbound Interceptor:** Active AST and regex inspection blocking SSNs, MRNs, phone numbers, and patient identifiers.
* **Tamper-Evident HMAC-SHA256 Audit Trail:** Chained, cryptographically signed logs for every evaluation and state transition.
* **Air-Gapped LLM Reasoning Adapter:** Agnostic integration for local Ollama instances (`llama3`, `mistral`), Claude 3.5 Sonnet, GPT-4o, and deterministic test mocks.
* **Active Learning Bayesian Calibration:** Dynamic tracker updating worker reliability weights and monitoring Brier calibration drift.
* **FastAPI & Prometheus Telemetry:** Exposes OpenAPI 3.1 REST endpoints and operational Prometheus metrics (`/metrics`).

---

## 🧪 Testing & Verification

Run the automated test suite:

```bash
pytest -v
```

Execute high-throughput batch simulation benchmarks:

```bash
python simulator.py --tasks 1000 --concurrency 8
```

---

## 🐳 Container Deployment

```bash
docker build -t dna-melting-temp-salt-corrected .
docker run -p 8000:8000 dna-melting-temp-salt-corrected
```
