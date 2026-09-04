# DNA Melting Temp Salt Corrected

> **Domain:** Bioinformatics & DNA Primer Design

<div align="center">

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/Python-3.10%20%7C%203.11%20%7C%203.12-3776AB.svg?logo=python&logoColor=white)
![FastAPI](https://img.shields.io/badge/FastAPI-0.111-009688.svg?logo=fastapi&logoColor=white)

</div>

---

## What It Does

DNA melting temperature (Tm) calculator with salt corrections for oligonucleotide primers. Implements:

- **Owczarzy 2004** monovalent cation salt correction
- **Owczarzy 2008** divalent cation (Mg2+) correction
- **Aho formula** for combined monovalent + divalent salt effective concentration
- **Multi-primer pair comparison** and ranking by Tm matching, GC balance, and secondary structure

### References
- Owczarzy et al. (2008) Biochemistry 47(19):5336-5353
- Owczarzy et al. (2004) Biochemistry 43(12):3537-3554
- Aho et al. (2008) Nucleic Acids Res 36(11):e63
- SantaLucia (1998) PNAS 95(4):1460-1465

---

## Installation

```bash
# Clone the repository
git clone https://github.com/abusuraihsakhri/dna-melting-temp-salt-corrected.git
cd dna-melting-temp-salt-corrected

# Install dependencies
pip install -e .

# Or install with development dependencies
pip install -e ".[dev]"
```

---

## Usage

### Command Line Interface

```bash
# Calculate Tm for a single sequence
python -c "from mg2_correction import calculate_tm_with_mg2; import json; print(json.dumps(calculate_tm_with_mg2('ATCGATCGATCGATCGATCG', na_mM=50, mg_mM=1.5), indent=2))"

# Run the oligo Tm calculator
python oligo_tm.py single --v1 12.0 --v2 4.0

# Batch process from CSV
python oligo_tm.py batch -i sample.csv -o results.csv
```

### Python API

```python
from mg2_correction import calculate_tm_with_mg2, owczarzy_2004_salt_correction, owczarzy_2008_mg_correction

# Calculate Tm with full salt and Mg2+ correction
result = calculate_tm_with_mg2(
    sequence="ATCGATCGATCGATCGATCG",
    na_mM=50.0,      # Na+ concentration in mM
    mg_mM=1.5,       # Mg2+ concentration in mM
    oligo_conc_nM=250.0  # Oligonucleotide concentration in nM
)
print(result)
# {
#     "sequence": "ATCGATCGATCGATCGATCG",
#     "length": 20,
#     "gc_fraction": 0.4,
#     "tm_basic_C": 64.0,
#     "tm_na_corrected_C": ...,
#     "tm_mg_corrected_C": ...,
#     ...
# }
```

### Multi-Primer Comparison

```python
from multi_primer_comparison import rank_primer_pairs

pairs = [
    {"id": "pair1", "forward": "ATCGATCGATCGATCG", "reverse": "GCTAGCTAGCTAGCTA"},
    {"id": "pair2", "forward": "GCGCGCGCGCGCGCGC", "reverse": "ATATATATATATATAT"},
]
results = rank_primer_pairs(pairs, na_mM=50, mg_mM=1.5)
```

### REST API Server

```bash
# Start the FastAPI server
python cli.py serve --host 127.0.0.1 --port 8000

# Or using Docker
docker build -t dna-melting-temp-salt-corrected .
docker run -p 8000:8000 dna-melting-temp-salt-corrected
```

API Endpoints:
- `GET /health` - Health check
- `GET /metrics` - Prometheus metrics
- `POST /api/audit` - Submit task for evaluation
- `POST /api/chat` - Query the system
- `GET /api/audit/logs` - View audit trail

---

## Core Functions

### `calculate_tm_with_mg2(sequence, na_mM, mg_mM, oligo_conc_nM)`
Calculate Tm with full salt and Mg2+ correction. Returns a dictionary with Tm values and correction details.

### `owczarzy_2004_salt_correction(tm_nn, na_M, f_GC)`
Apply Owczarzy 2004 monovalent cation salt correction.

ΔTm = (4.29 * f_GC - 3.95) * ln(Na+) + 0.940 * (ln(Na+))^2

### `owczarzy_2008_mg_correction(tm_nn, na_M, mg_M, f_GC, nbp)`
Apply Owczarzy 2008 Mg2+ correction to nearest-neighbor Tm.

### `compute_effective_salt_concentration(na_mM, mg_mM, oligo_conc_nM)`
Compute effective monovalent cation concentration using the Aho formula.

---

## Testing

```bash
# Run all tests
pytest -v

# Run specific test modules
pytest tests/test_mg2_correction.py -v
pytest tests/test_enrichment.py -v
pytest test_oligo_tm.py -v
```

---

## Project Structure

```
dna-melting-temp-salt-corrected/
├── mg2_correction.py          # Core Tm calculation with salt corrections
├── oligo_tm.py                # Oligo Tm calculator CLI
├── multi_primer_comparison.py # Multi-primer pair comparison engine
├── secondary_structure.py     # Hairpin and self-dimer prediction
├── cli.py                     # Main CLI interface
├── simulator.py               # High-throughput simulation
├── enrichment.py              # Enrichment feature modules
├── agents/                    # Enterprise agent framework
│   ├── base.py               # Security, PHI guard, audit trail
│   ├── models.py             # Pydantic data models
│   ├── supervisor.py         # Supervisor orchestrator
│   ├── workers.py            # Specialized worker agents
│   ├── api.py                # FastAPI REST server
│   └── ...
├── tests/                     # Test suite
├── pyproject.toml             # Project configuration
├── Dockerfile                 # Container configuration
└── README.md
```

---

## Security

- **PHI Outbound Guard:** AST and regex inspection blocking SSNs, MRNs, phone numbers, and patient identifiers
- **HMAC-SHA256 Audit Trail:** Chained, cryptographically signed logs for every evaluation
- **Input Validation:** All inputs validated before processing

Set the audit secret key via environment variable:
```bash
export AUDIT_SECRET_KEY="your-secure-key-here"
```

---

## License

MIT License - see [LICENSE](LICENSE) for details.
