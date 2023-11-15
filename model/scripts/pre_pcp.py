import torch
import numpy as np
import pandas as pd
from tqdm import tqdm
from transformers import AutoTokenizer
from transformers import AutoModelForSequenceClassification, Trainer
from datasets import Dataset
import argparse

##参数
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="input file path, must have seq column")
parser.add_argument("-o", "--output", help="output file")
parser.add_argument("-m", "--model", help="model dir path")
args = parser.parse_args()
seq_file = args.input
pre_file = args.output
model_path = args.model

def get_predict(seq, trainer):
    pre_tokenized = tokenizer(seq,padding="max_length",truncation=True)
    pre_dataset = Dataset.from_dict(pre_tokenized)
    raw_pred, _, _ = trainer.predict(pre_dataset) 
    return(raw_pred.flatten())

model = AutoModelForSequenceClassification.from_pretrained(model_path, num_labels=1, problem_type = "regression", local_files_only=True)
tokenizer = AutoTokenizer.from_pretrained(model_path)

###prediction
trainer = Trainer(model=model)

pcp_val = pd.read_csv(seq_file)
sequences = list(pcp_val.seq)
res = get_predict(sequences,trainer)
pcp_val["pre"] = res
pcp_val.to_csv(pre_file)


