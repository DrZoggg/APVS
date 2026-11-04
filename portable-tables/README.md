# Read the existing APVS repository tables without R

These three CSV files preserve the values and row labels of small R data-frame
objects already public in [DrZoggg/APVS](https://github.com/DrZoggg/APVS) at commit
`c1bd0e5321ff39b5f7b0fa5242164a4c3bf1b2b2`:

| CSV | Original object | Rows | Original data columns |
| --- | --- | ---: | --- |
| `RMSE_residual.csv` | `APVS Generation/RMSE_residual.rda`: `RMSE_residual` | 414 | `model_type`, `RMSE`, `residual` |
| `performance.csv` | `APVS Generation/performance.rda`: `performance` | 5 | `X`, `nnet`, `rf`, `glmboost`, `glmnet`, `treebag`, `bayes`, `knn`, `pls`, `rpart` |
| `IPSG.csv` | `APVSLevel Interpretability/IPS ref.rda`: `IPSG` | 162 | `GENE`, `NAME`, `CLASS`, `WEIGHT` |

`PROVENANCE.json` records original paths, Git blobs and SHA-256, column types and
factor levels, CSV hashes, and the encoding contract. `__rowname__` preserves the
original row label; it does not establish a sample identifier. Numeric columns
use 17 significant digits. The reserved missing token is `__APVS_R_NA__`; a literal
`NA` string remains a string. R and Python readbacks were checked against the
original values, including numeric bits and factor labels.

From this directory, Python 3 and its standard library can read an existing stored
metric row without R or third-party packages:

```python
import csv
import json
from pathlib import Path

base = Path(".")
meta = json.loads((base / "PROVENANCE.json").read_text(encoding="utf-8"))
with (base / "performance.csv").open(encoding="utf-8", newline="") as handle:
    rows = list(csv.DictReader(handle))
missing = meta["encoding"]["na_token"]
rows = [{key: None if value == missing else value for key, value in row.items()}
        for row in rows]
print(next(row for row in rows if row["X"] == "auc"))
```

This reads retained repository values; it does not recalculate performance or
provide new-sample APVS scoring. The repository README associates the research
materials with **Atherosclerotic plaque vulnerability quantification system for
clinical and biological interpretability**, *iScience* (2023),
[DOI 10.1016/j.isci.2023.107587](https://doi.org/10.1016/j.isci.2023.107587).
Use that formal paper citation when reusing the associated research and retain
the repository commit/object provenance. No specific paper table/figure mapping
or unverified sample meaning is assigned here. The fixed source declares no
license; this conversion adds no new license grant.
