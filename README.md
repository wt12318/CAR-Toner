# PCP-AIOptimizer

Fine-tuned model based on ESM2-8M for PCP calculation.

### Dependency

- pytorch 2.0.0 
- numpy 1.26.0
- pandas 1.5.3 
- tqdm 4.65.0 
- transformers 4.27.1

The model weight is saved in data/checkpoint-22176, user can predict PCP score by input `csv` file with `seq` column:

```shell
python pre_pcp.py -i input.csv -o input_pre.csv -m data/checkpoint-checkpoint-22176/
```

The `pre` column in output file is `log2(PCP)`, thus can be convert to PCP by `round(2^(input_pre$pre))`.


