import json
import pandas as pd

# --- 1. Load all-generated log (all molecules: with duplicates) ---
with open("logs/adaptive_log_all_generated_20250519_014408.json") as f:
    all_generated = json.load(f)

all_gen_records = []
for round_str, entries in all_generated["generated"].items():
    for entry in entries:
        entry = dict(entry)
        entry["round"] = int(entry["round"])
        all_gen_records.append(entry)
df_all = pd.DataFrame(all_gen_records)

# --- 2. Load accepted log (accepted, may have duplicates) ---
with open("logs/adaptive_log_rounds_20250519_014408.json") as f:
    accepted = json.load(f)

accepted_set = set()
for round_str, entries in accepted["iterations"].items():
    for entry in entries:
        # Use canonical SMILES if present
        smi = entry.get("smiles")
        round_num = int(entry.get("round", round_str))
        accepted_set.add((smi, round_num))

# --- 3. Load cleaned/canonical log (unique accepted) ---
with open("logs/cleaned/cleaned_log_round_20250519_014408.json") as f:
    cleaned = json.load(f)

# Unique canonical SMILES
cleaned_unique = set()
for round_str, entries in cleaned["iterations"].items():
    for entry in entries:
        smi = entry.get("smiles")
        # Only the first occurrence per canonical SMILES (across all rounds)
        cleaned_unique.add(smi)

# --- 4. Add flags to df_all ---
def is_accepted(row):
    return (row['smiles'], row['round']) in accepted_set

def is_firstunique(row, seen=set()):
    # will mutate seen!
    smi = row['smiles']
    if smi not in seen and smi in cleaned_unique:
        seen.add(smi)
        return True
    else:
        return False

df_all["accepted"] = df_all.apply(is_accepted, axis=1)

# For keepfirstunique, process by SMILES order in file
seen = set()
keep_first_unique_flags = []
for idx, row in df_all.iterrows():
    smi = row['smiles']
    if smi not in seen and smi in cleaned_unique:
        keep_first_unique_flags.append(True)
        seen.add(smi)
    else:
        keep_first_unique_flags.append(False)
df_all["keepfirstunique"] = keep_first_unique_flags

# --- 5. Save master log ---
df_all.to_parquet("logs/adaptive_master_log.parquet", index=False)
print("Saved logs/adaptive_master_log.parquet")
print(df_all.head())