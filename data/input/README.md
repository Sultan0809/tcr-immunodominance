# Input data

## `epitopes_predicted.csv`

48 rows, 44 unique epitopes, from three published immunoinformatics vaccine-design
studies. Four epitopes appear twice because they were reported in both the broad
screen and the final selected set. See [../CITATION.md](../CITATION.md) for full
references.

| Column | Meaning |
|---|---|
| `epitope` | Predicted epitope sequence |
| `organism` | SARS-CoV-2 or HSV-2 |
| `mhc_class` | I (CD8+) or II (CD4+) |
| `source_table` | Table in the source publication (A9, B6-final, C7) |
| `position` | Start residue in the parent protein, where reported |
| `iedb_epitope` | Matching IEDB entry; blank if none |
| `iedb_id` | IEDB epitope ID (`iedb.org/epitope/<id>`) |
| `iedb_identity_pct` | BLAST identity to that entry |
| `iedb_verified` | Experimentally verified for T-cell activation |
| `iedb_match_type` | Self-hit / Independent match / Novel (see below) |
| `study` | Source publication key |
| `study_doi` | DOI of that publication |

### Composition

| Study | Epitopes | Target |
|---|---|---|
| `Khan2024_VirusDisease` | 21 | SARS-CoV-2 S/M/N/E, class I |
| `Khan2024_HumanImmunol` | 17 | SARS-CoV-2 spike (Omicron), 6 class I + 11 class II |
| `Khan2025_MolDivers` | 10 | HSV-2 glycoprotein D, 2 class I + 8 class II |

| `iedb_match_type` | n | Meaning |
|---|---|---|
| Self-hit | 32 | 100% identity to an entry already inside the prediction tools' training data — confirms the epitope is known, but **not** independent corroboration |
| Independent match | 14 | Below 100% identity, so genuinely corroborating evidence |
| Novel | 2 | No IEDB match |

All 10 HSV-2 gD epitopes are independent matches; most SARS-CoV-2 entries are
self-hits.

### Using your own epitopes

Replace this file, keeping at least the `epitope` column, and re-run:

```bash
Rscript scripts/07_predicted_epitopes.R
```

`organism` and `mhc_class` are used for grouping in the output; everything else
is carried through for reference.
