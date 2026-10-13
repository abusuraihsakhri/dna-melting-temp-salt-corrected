# DNA Melting Temperature — Salt Corrected

A DNA/DNA duplex melting-temperature calculator for sequence-specific oligonucleotide analysis. The primary workflow estimates the melting temperature of a perfectly matched duplex from nearest-neighbor thermodynamics and salt conditions. The browser interface runs locally, without a backend.

## Use the calculator

Open [the browser calculator](web/index.html) directly, or use the GitHub Pages deployment after it has been configured and verified. Paste an A/C/G/T sequence, set Na⁺, **free** Mg²⁺, and total oligo concentrations, and select **Calculate Tm**. Batch mode accepts one sequence per line, FASTA, or CSV files with a `sequence` header (optional `id`); results can be exported to CSV.

Inputs are processed in the browser. The site does not upload sequences, send telemetry, or need server-side dependencies. The application has no external JavaScript, font, or analytics dependencies.

## Scientific method and limits

- Nearest-neighbor DNA/DNA enthalpy and entropy from the SantaLucia/Allawi parameterization, with terminal A/T or G/C initiation contributions.
- Uncorrected Tm is calculated in Kelvin using `ΔH / (ΔS + R ln(Ct/4))`, then converted to °C. The Ct/4 convention assumes a **non-self-complementary** perfect duplex.
- Monovalent correction uses Owczarzy (2004), applied to **inverse Tm in Kelvin**.
- Mixed Na⁺/Mg²⁺ correction uses Owczarzy (2008), including the ion-regime decision tree and its empirical constants.
- `tm_basic_C` is retained as a diagnostic Wallace estimate, **not** the basis of the reported corrected Tm.
- Free Mg²⁺ must be supplied by the caller. dNTP chelation, mismatches, self-complementary concentration correction, hairpins, competing structures, and additives are not included in the melting-temperature model. Results are predictions, not empirical assay measurements.

The `secondary_structure.py` and `multi_primer_comparison.py` modules provide additional approximate exploratory evaluations and should not be treated as experimentally validated structure free-energy predictions.

**References:** [SantaLucia (1998)](https://doi.org/10.1073/pnas.95.4.1460); Owczarzy et al., *Biochemistry* (2004), 43:3537–3554; [Owczarzy et al. (2008)](https://doi.org/10.1021/bi702363u). Numeric regression cases use independently calculated Biopython `Tm_NN` results (DNA_NN3, salt correction 7).

## Python usage

Requires Python 3.10 or newer.

```bash
python -m pip install -e ".[dev]"
pytest -q
node tests/browser_smoke.cjs  # optional; Node.js
```

```bash
# Estimate one DNA sequence
dna-melting-tm single --sequence ATCGATCGATCGATCGATCG --na 50 --mg 1.5 --oligo 250

# Process input.csv with a sequence column
dna-melting-tm batch --input input.csv --output results.csv
```

```python
from mg2_correction import calculate_tm_with_mg2

result = calculate_tm_with_mg2(
    "ATCGATCGATCGATCGATCG",
    na_mM=50.0, mg_mM=1.5, oligo_conc_nM=250.0,
)
print(result["tm_mg_corrected_C"])  # 59.72 °C
```

For an optional FastAPI service, set a private audit signing key in the environment before starting the legacy audit agent:

```bash
export AUDIT_SECRET_KEY="$(python -c 'import secrets; print(secrets.token_hex(32))')"
python cli.py serve --host 127.0.0.1 --port 8000
```

The **Tm API** is `POST /api/tm` with JSON such as `{"sequence":"ATCGATCGATCGATCGATCG","na_mM":50,"mg_mM":1.5,"oligo_conc_nM":250}`. Additional `/api/audit`, `/api/chat`, `/api/audit/logs`, `/health`, and `/metrics` endpoints are retained from the legacy mock-agent prototype; they are **not** substitutes for the thermodynamic calculations. The Docker Compose service also requires `AUDIT_SECRET_KEY` to be supplied externally; no production key is embedded.

## Development and deployment

- `mg2_correction.py`: standalone thermodynamic calculation (Python standard library).
- `oligo_tm.py`: sequence CLI and streaming CSV batch processing.
- `web/index.html`: static HTML/CSS/JavaScript browser implementation, with no build step.
- `tests/` and `test_oligo_tm.py`: Python and browser arithmetic regression checks.
- `.github/workflows/ci.yml`: Python 3.10–3.12 and browser smoke tests.
- `.github/workflows/pages.yml`: publishes the `web/` directory from `master` using GitHub Pages Actions, if Pages is configured to accept Actions deployments. Root `index.html` redirects to `web/` for branch-based Pages hosting.

The browser tool uses common modern JavaScript/HTML features and is intended for current desktop and mobile browsers. JavaScript must be enabled. The server-side FastAPI audit prototype cannot run on GitHub Pages; the static calculator is implemented in JavaScript instead of Pyodide to avoid WebAssembly loading overhead.

## License

[MIT](LICENSE).
