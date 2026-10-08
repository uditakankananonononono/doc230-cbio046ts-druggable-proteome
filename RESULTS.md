# CBIO046T-S results - partial: five-fold expansion claim not testable (needs structure-based pocket prediction). Annotation/sequence classifier only.
Data: 20,190 UniProt reviewed human proteins with unique gene symbols; label = has >=1 Hetionet v1.0 binding compound (DrugBank-derived, 2016): 1,633 positives (8.1%). Pharos TDL was not usable (API 10-record pages). 71.5% of proteins have a family string; the rest are singleton groups in the cold-family split.
- D2 random 5-fold AUROC: HGB 0.925, LR 0.920; family-prior baseline 0.855.
- D3 cold-family AUROC: LR 0.886, HGB 0.878 (gate >= 0.75 MET).
- D4 family-revealing keywords removed, cold-family: HGB 0.856 (gate >= 0.70 MET); the signal remains without the obvious class keywords.
- D5 sequence only (aa composition + length), cold-family: HGB 0.758, LR 0.629. D6: keywords only 0.885, location only 0.673, sequence only 0.758.
- D7 shuffle control 0.472 (chance).
- D8 cold-family OOF AUROC by class: kinase 0.575 (n 625; 277 labelled), receptor 0.814, ion channel 0.670, hydrolase/enzyme 0.727, other 0.906. Within-class discrimination is much weaker than the overall number, which is largely driven by separating target-rich classes from the rest.
- D9 calibration: Brier 0.0592 vs 0.0743 base, ECE 0.020.
- D10 18,557 unlabeled proteins; 212 with p > 0.5, 12 with p > 0.8; top: GAPDH, ACADVL, CD36, CYBB, LBR, UQCR11, PRDX6, KCNK5, KARS1, P2RX7. These are model scores from an old, DrugBank-biased label, not validated druggable targets; several (e.g. GAPDH, HBB) look like binder artifacts. No expansion claim.
Limits: label is "has a known drug binder in a 2016 resource", not curated druggability; UniProt keywords are curated with knowledge that may overlap with drug-target annotation (possible annotation leakage, not testable here); no structure features. 0 novel findings.
