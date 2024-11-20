import xarray as xr
import numpy as np
from sklearn.metrics import confusion_matrix
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os

def load_and_flatten_data(filename):
    """
    Load NetCDF file and flatten the PTYPE variable

    Parameters:
    filename (str): Path to NetCDF file

    Returns:
    numpy.ndarray: Flattened array of classifications
    """
    ds = xr.open_dataset(filename)
    ptype = ds['PTYPE'].values

    # Flatten all dimensions (time, lat, lon)
    flat_ptype = ptype.reshape(-1)

    return flat_ptype

def compare_classifications(reference_data, comparison_data, class_labels=[1, 2, 4, 8]):
    """
    Compare two classification arrays and compute metrics, excluding areas where reference data is 0

    Parameters:
    reference_data (numpy.ndarray): Flattened reference classifications
    comparison_data (numpy.ndarray): Flattened comparison classifications
    class_labels (list): List of possible classification values

    Returns:
    dict: Dictionary containing comparison metrics
    """
    # Create mask for non-zero and non-NaN values in reference data
    valid_mask = (reference_data != 0) & ~np.isnan(reference_data) & ~np.isnan(comparison_data)

    # Apply mask to both datasets
    ref_clean = reference_data[valid_mask]
    comp_clean = comparison_data[valid_mask]

    # Calculate confusion matrix
    conf_matrix = confusion_matrix(ref_clean, comp_clean, labels=class_labels)

    # Calculate per-class metrics
    metrics = {}
    for idx, class_label in enumerate(class_labels):
        true_pos = conf_matrix[idx, idx]
        false_pos = conf_matrix[:, idx].sum() - true_pos
        false_neg = conf_matrix[idx, :].sum() - true_pos
        true_neg = conf_matrix.sum() - (true_pos + false_pos + false_neg)

        precision = true_pos / (true_pos + false_pos) if (true_pos + false_pos) > 0 else 0
        recall = true_pos / (true_pos + false_neg) if (true_pos + false_neg) > 0 else 0
        f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0

        metrics[class_label] = {
            'precision': precision,
            'recall': recall,
            'f1_score': f1
        }

    # Calculate overall accuracy
    accuracy = np.sum(ref_clean == comp_clean) / len(ref_clean)

    return {
        'confusion_matrix': conf_matrix,
        'per_class_metrics': metrics,
        'overall_accuracy': accuracy,
        'n_samples': len(ref_clean)
    }

def compare_multiple_datasets_years(reference_dataset, comparison_datasets, start_year, end_year, fzra_directory='./', skip_years=[], class_labels=[1, 2, 4, 8]):
    """
    Compare multiple datasets against a reference dataset across multiple years

    Parameters:
    reference_dataset (str): Name of reference dataset (e.g., 'ECMWF')
    comparison_datasets (list): List of dataset names to compare
    start_year (int): Starting year
    end_year (int): Ending year
    skip_years (list): List of years to skip in the analysis
    class_labels (list): List of possible classification values
    """
    yearly_results = {}

    # Initialize storage for aggregate metrics
    dataset_summaries = {dataset: {
        'accuracy': [],
        'class_metrics': {cls: {'precision': [], 'recall': [], 'f1_score': []} for cls in class_labels}
    } for dataset in comparison_datasets}

    for year in range(start_year, end_year + 1):
        if year in skip_years:
            print(f"\nSkipping year {year} (in skip_years list)")
            continue

        print(f"\nProcessing year {year}")

        # Construct reference filename
        reference_file = os.path.join(fzra_directory, f"{reference_dataset}-{year}.nc")

        # Process each comparison dataset
        year_results = {}

        for dataset in comparison_datasets:
            comp_file = os.path.join(fzra_directory, f"{dataset}-{year}.nc")
            print(f"\nProcessing: {comp_file}")

            try:
                # Load data
                reference_data = load_and_flatten_data(reference_file)
                comp_data = load_and_flatten_data(comp_file)

                # Get comparison metrics
                metrics = compare_classifications(reference_data, comp_data, class_labels)

                # Store results
                year_results[dataset] = metrics

                # Update dataset summaries
                dataset_summaries[dataset]['accuracy'].append(metrics['overall_accuracy'])
                for class_label, class_metrics in metrics['per_class_metrics'].items():
                    dataset_summaries[dataset]['class_metrics'][class_label]['precision'].append(class_metrics['precision'])
                    dataset_summaries[dataset]['class_metrics'][class_label]['recall'].append(class_metrics['recall'])
                    dataset_summaries[dataset]['class_metrics'][class_label]['f1_score'].append(class_metrics['f1_score'])

                # Print metrics
                print(f"\nOverall Accuracy: {metrics['overall_accuracy']:.3f}")
                print(f"Number of valid samples: {metrics['n_samples']}")
                print("\nPer-class metrics:")
                for class_label, class_metrics in metrics['per_class_metrics'].items():
                    print(f"\nClass {class_label}:")
                    print(f"Precision: {class_metrics['precision']:.3f}")
                    print(f"Recall: {class_metrics['recall']:.3f}")
                    print(f"F1-score: {class_metrics['f1_score']:.3f}")

            except FileNotFoundError as e:
                print(f"Warning: File not found: {comp_file}")
                continue

        if year_results:
            yearly_results[year] = year_results

    # Print summary statistics
    print("\n=== Summary Statistics ===")
    print(f"Years analyzed: {sorted(set(range(start_year, end_year + 1)) - set(skip_years))}")
    print(f"Years skipped: {sorted(skip_years)}")

    # Print summary statistics for each dataset
    for dataset in comparison_datasets:
        if dataset_summaries[dataset]['accuracy']:  # Only print if we have data
            print(f"\nDataset: {dataset}")
            print(f"Mean Accuracy: {np.mean(dataset_summaries[dataset]['accuracy']):.3f} ± {np.std(dataset_summaries[dataset]['accuracy']):.3f}")

            print("\nClass Metrics (mean ± std):")
            for class_label in class_labels:
                metrics = dataset_summaries[dataset]['class_metrics'][class_label]
                if metrics['precision']:  # Only print if we have data
                    print(f"\nClass {class_label}:")
                    print(f"Precision: {np.mean(metrics['precision']):.3f} ± {np.std(metrics['precision']):.3f}")
                    print(f"Recall: {np.mean(metrics['recall']):.3f} ± {np.std(metrics['recall']):.3f}")
                    print(f"F1-score: {np.mean(metrics['f1_score']):.3f} ± {np.std(metrics['f1_score']):.3f}")

    return yearly_results

# Example usage:
if __name__ == "__main__":
    reference_dataset = "ECMWF"
    comparison_datasets = ["ERA5", "CFSR", "CR20", "JRA"]
    start_year = 1980
    end_year = 2016
    skip_years = []
    #skip_years = [1985, 1991, 1993, 1997]
    FZRAPATH = '/glade/derecho/scratch/zarzycki/FZRA/ptypes/'

    results = compare_multiple_datasets_years(
        reference_dataset,
        comparison_datasets,
        start_year,
        end_year,
        fzra_directory=FZRAPATH,
        skip_years=skip_years
    )
