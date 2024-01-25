# CAR-Toner

Fine-tuned model based on ESM2-8M for PCP calculation.

### Web interface

The web interface of CAR-Toner consists of two part: single-mode and batch mode. <br />
In single mode, user can input one sequence and calculate its PCP score. Then system show number of mutated sequences and position of tuning mutations (k-q for tuning down, q-k for tuning up). Once user click tuning button, the PCP prediction of mutated sequences is running. After completed, the tuning results are shown as table and the distribution of PCP scores are also drawn. User can select specific range of PCP score in the density plot and corresponding results are listed in table below. <br />
In batch mode, user can input multiple sequence in text box or in `Fasta` file, the predicted results are shown as table.

![single mode](./web_app/pcp1.png)
![single mode](./web_app/pcp3.png)
![single mode](./web_app/pcp4.png)
![batch mode](./web_app/pcp2.png)

### Citation
1. Tuning Charge Density of Chimeric Antigen Receptor Optimizes CAR-T Cell Fitness. ***Cell Research*** 2023 33, 341-354. PMID: 36882513
2. CAR-Toner: An AI-Driven Approach for CAR Tonic Signaling Prediction and Optimization. ***Cell Research*** 2024 (In press).


