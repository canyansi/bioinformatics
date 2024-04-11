"""
This script was made to call CNVs from probe data. It prints frequencies of different CNVs per ethnicity per breakpoints
Has a hardcoded path to the sample data. (should be in directory)

Run by:
python CNV_identification.py
"""

#############################
# import packages
import numpy as np
import pandas as pd
from scipy.stats import zscore


#############################
# defining functions
def normalize_data(probe_data):
    """Normalize for sample variability using sample means."""
    sample_means = probe_data.mean(axis=1)
    sample_normalized_data = probe_data.div(sample_means, axis=0)
    # Normalize for probe variability using probe median.
    probe_medians = sample_normalized_data.median(axis=0)
    normalized_probe_data = sample_normalized_data.div(probe_medians, axis=1)
    return normalized_probe_data


def ambiguous_aware_round(estCN, lowthreshold=0.49, upperthreshold=0.51):
    """Simple method to call CN by rounding. Ambiguous areas as
    definied will be labeled as NA"""
    ambiguous_CN = lowthreshold < estCN % 1 < upperthreshold
    if ambiguous_CN:
        return np.nan
    else:
        return round(estCN)


def calculate_CN(normalized_data, norm_columns):
    """Estimate copy number of target probes by dividing the deptths in target CNSL probes
    by the median depth in non CNSL CN=2 probes (normalization probes).
    And then multiplying by 2 (the expected CN in non CNSL probes)"""
    normregions_depth = normalized_data.loc[:, norm_columns].median(axis=1)
    est_CNs = (
        normalized_data.drop(columns=norm_columns)
        .div(normregions_depth, axis=0)
        .multiply(2)
    )
    # definitively calls the CN by custom rounding method
    integer_CNs = est_CNs.map(ambiguous_aware_round)
    return integer_CNs


def call_CNVs(sample_row, minseg_length=4):
    """Stores the breakpoints and keeps track of the current segment length
    returns a list of breakpoint coordinates in format: 'probe_start,probe_stop_cnvtype"""
    breakpoints = []
    current_count = 1
    CN_values = sample_row.values
    # iterate through the CNs (starting from 1 not 0)
    for i in range(1, len(CN_values)):
        # check if current CN value is equal to previous
        if CN_values[i] == CN_values[i - 1]:
            current_count += 1
        else:
            # check if the previous segment is long enough to be classified as CNV. if so
            # saves it in breakpoints list along with the CNV type
            if current_count >= minseg_length:
                bp_id = f"{sample_row.index[i-current_count]},{sample_row.index[i-1]}"
                if CN_values[i - 1] == 1:
                    breakpoints.append(bp_id + "_deletion")
                elif CN_values[i - 1] == 3:
                    breakpoints.append(bp_id + "_duplication")
            current_count = 1
    # Accounts for possible tail segment
    if current_count >= minseg_length:
        bp_id = f"{sample_row.index[len(CN_values)-current_count]},{sample_row.index[len(CN_values)-1]}"
        if CN_values[-1] == 1:
            breakpoints.append(bp_id + "_deletion")
        elif CN_values[-1] == 3:
            breakpoints.append(bp_id + "_duplication")
    return breakpoints


def calculate_bp_freqs(cnv_calls, sample_ethnicities):
    """Function to calculate sample based stats on the CNV calls. takes as input
    cnv calls and the corresponding sample ethnicity meta data the output is a dictionary. It
    holds the CNV call_id (bp_id), the number of total times this call was identified (total).
    the samples containing this CNV (sample_idx), and the number of people per ethnicity with this call
    (ethnicity_counts)"""

    bp_freqs = {}
    for idx, cnv_call in enumerate(cnv_calls):
        if len(cnv_call) == 0:
            continue
        sample_ethnicity = sample_ethnicities[idx]
        for call_id in cnv_call:
            if call_id not in bp_freqs:
                bp_freqs[call_id] = {
                    "ethnicity_counts": {
                        ethnicity: 0 for ethnicity in set(sample_ethnicities)
                    },
                    "sample_idx": [],
                    "total": 0,
                    "type": call_id.split("_")[1],
                }
            bp_freqs[call_id]["sample_idx"].append(idx)
            bp_freqs[call_id]["ethnicity_counts"][sample_ethnicity] += 1
            bp_freqs[call_id]["total"] += 1
    return bp_freqs


def calculate_ethnicity_rates(bp_freqs, sample_ethnicities):
    """Function to calculate overall frequency rates per ethnicity"""
    ethnic_pop_total = pd.Series(sample_ethnicities).value_counts().to_dict()
    for bp_id in bp_freqs.keys():
        ethnicity_counts = bp_freqs[bp_id]["ethnicity_counts"]
        ethnicity_rates = {
            ethnicity: ethnicity_counts[ethnicity] / ethnic_pop_total[ethnicity]
            for ethnicity in ethnicity_counts
        }
        bp_freqs[bp_id]["ethnicity_rates"] = ethnicity_rates
    return bp_freqs


##############
# Actualy run script

# read in data
sample_probe_data = pd.read_csv("cnsl_data.csv", index_col=0)
sample_ethnicities = sample_probe_data["ethnicity"]
probe_data = sample_probe_data.drop(columns=["ethnicity"])

# normalize data
normalized_probe_data = normalize_data(probe_data)

# call CN
norm_columns = normalized_probe_data.filter(like="non").columns
cnsl_probe_CNs = calculate_CN(normalized_probe_data, norm_columns)


# check for malfunctioning probes. find probes that have more than 1000 (10%)
# of the samples more than 3 standard deviations away from the average (zscore)
z_scores = cnsl_probe_CNs.apply(zscore, axis=1)
outlier_columns = (z_scores > 3).sum()[(z_scores > 3).sum() > 1000]
cnsl_probe_CNs = cnsl_probe_CNs.drop(columns=outlier_columns.index)


# call CNV detection function on the copy number DF
cnv_calls = cnsl_probe_CNs.apply(call_CNVs, axis=1)

cnv_freqs = calculate_bp_freqs(cnv_calls, sample_ethnicities)
cnv_freqs = calculate_ethnicity_rates(cnv_freqs, sample_ethnicities)

# display results
print(f"{len(cnv_freqs)} different CNV/breakpoint combinations were detected")

# filter to those that occur at least 10 times.
filtered_breakpoints = {
    bp_id: stats for bp_id, stats in cnv_freqs.items() if stats["total"] > 10
}

# View proportions and counts
counts = pd.DataFrame(
    {bp_id: stats["ethnicity_counts"] for bp_id, stats in filtered_breakpoints.items()}
)
print(counts.T)
proportions = pd.DataFrame(
    {bp_id: stats["ethnicity_rates"] for bp_id, stats in filtered_breakpoints.items()}
)
print(proportions.T)
