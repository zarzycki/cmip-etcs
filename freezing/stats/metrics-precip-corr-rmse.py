import xarray as xr
import numpy as np
from scipy import stats
import pandas as pd

def convert_precipitation(data, dataset, datasets_to_convert):
    """
    Convert precipitation to mm/day

    Parameters:
    data : xarray DataArray
    dataset : str, name of dataset
    datasets_to_convert : list of dataset names that need conversion
    """
    if dataset in datasets_to_convert:
        return data * 86400  # Convert from kg/m^2/s to mm/day
    return data

def calculate_metrics(forecast, observed, dataset_name, offset_hours, threshold=1.0):
    """
    Calculate verification metrics between forecast and observed precipitation

    Parameters:
    forecast : numpy array
    observed : numpy array
    dataset_name : str
    offset_hours : int, time offset in hours
    threshold : float, threshold for contingency metrics
    """
    # Mask out missing values and flatten arrays
    mask = ~(np.isnan(forecast) | np.isnan(observed))
    f = forecast[mask]
    o = observed[mask]

    # Check if we have valid data to process
    if len(f) == 0 or len(o) == 0:
        print(f"\nWarning: No valid data for {dataset_name} with {offset_hours} hour offset")
        return {
            'rmse': np.nan,
            'mae': np.nan,
            'bias': np.nan,
            'pearson_r': np.nan,
            'pod': np.nan,
            'far': np.nan,
            'csi': np.nan,
            'freq_bias': np.nan,
            'mean_forecast': np.nan,
            'mean_observed': np.nan,
            'offset_hours': offset_hours,
            'valid_points': 0
        }

    print(f"\nData Summary for {dataset_name} (offset: {offset_hours} hours):")
    print(f"Number of valid points: {len(f)}")
    print(f"Forecast: mean={np.mean(f):.3f}, median={np.median(f):.3f}, min={np.min(f):.3f}, max={np.max(f):.3f}")
    print(f"Observed: mean={np.mean(o):.3f}, median={np.median(o):.3f}, min={np.min(o):.3f}, max={np.max(o):.3f}")

    # Print precipitation distribution
    percentiles = [50, 75, 90, 95, 99]
    print(f"\nForecast percentiles:")
    for p in percentiles:
        print(f"{p}th percentile: {np.percentile(f, p):.3f}")

    print(f"\nObserved percentiles:")
    for p in percentiles:
        print(f"{p}th percentile: {np.percentile(o, p):.3f}")

    # Basic error metrics
    rmse = np.sqrt(np.mean((f - o)**2))
    mae = np.mean(np.abs(f - o))
    bias = np.mean(f - o)

    # Correlation coefficients
    try:
        pearson_r, pearson_p = stats.pearsonr(f, o)
    except:
        pearson_r = np.nan

    # Contingency table metrics
    hits = np.sum((f >= threshold) & (o >= threshold))
    misses = np.sum((f < threshold) & (o >= threshold))
    false_alarms = np.sum((f >= threshold) & (o < threshold))
    correct_negatives = np.sum((f < threshold) & (o < threshold))

    print(f"\nContingency Table (threshold={threshold} mm/day):")
    print(f"Hits: {hits}")
    print(f"Misses: {misses}")
    print(f"False Alarms: {false_alarms}")
    print(f"Correct Negatives: {correct_negatives}")

    # Safe division for contingency metrics
    pod = hits / (hits + misses) if (hits + misses) > 0 else np.nan
    far = false_alarms / (hits + false_alarms) if (hits + false_alarms) > 0 else np.nan
    csi = hits / (hits + misses + false_alarms) if (hits + misses + false_alarms) > 0 else np.nan
    freq_bias = (hits + false_alarms) / (hits + misses) if (hits + misses) > 0 else np.nan

    return {
        'rmse': rmse,
        'mae': mae,
        'bias': bias,
        'pearson_r': pearson_r,
        'pod': pod,
        'far': far,
        'csi': csi,
        'freq_bias': freq_bias,
        'mean_forecast': np.mean(f),
        'mean_observed': np.mean(o),
        'offset_hours': offset_hours,
        'valid_points': len(f)
    }

def process_verification(start_year, end_year, reference_dataset, test_datasets, datasets_to_convert,
                        data_directory='./', time_offsets=[-6, 0, 6], threshold=0.1):
    """
    Process verification for multiple datasets against a reference with time offsets

    Parameters:
    start_year : int, first year to process
    end_year : int, last year to process (inclusive)
    reference_dataset : str, name of reference dataset
    test_datasets : list of datasets to verify against reference
    datasets_to_convert : list of datasets that need unit conversion
    time_offsets : list of ints, time offsets in hours to test
    threshold : float, threshold for contingency metrics
    """
    all_results = {}
    years = range(start_year, end_year + 1)

    for year in years:
        results = {}

        try:
            # Load reference data
            ref_file = f'{data_directory}/pre-{reference_dataset}-{year}.nc'
            print(f"\nLoading reference data from {ref_file}")
            ref = xr.open_dataset(ref_file)['PRECT']
            ref_data = convert_precipitation(ref, reference_dataset, datasets_to_convert)

            # Process each test dataset
            for dataset in test_datasets:
                print(f"\nProcessing {dataset} for {year}")
                ds_file = f'{data_directory}/pre-{dataset}-{year}.nc'

                try:
                    ds = xr.open_dataset(ds_file)['PRECT']
                    ds_data = convert_precipitation(ds, dataset, datasets_to_convert)

                    # Test different time offsets
                    for offset in time_offsets:
                        # Calculate the number of time steps to shift
                        # Assuming 6-hourly data, one time step = 6 hours
                        time_steps = offset // 6

                        if time_steps != 0:
                            # Shift the data and handle edges
                            shifted_data = ds_data.shift(time=time_steps)

                            # For positive shifts, we lose data at the end
                            # For negative shifts, we lose data at the beginning
                            if time_steps > 0:
                                shifted_data = shifted_data.isel(time=slice(None, -time_steps))
                                ref_slice = ref_data.isel(time=slice(None, -time_steps))
                            else:
                                shifted_data = shifted_data.isel(time=slice(-time_steps, None))
                                ref_slice = ref_data.isel(time=slice(-time_steps, None))
                        else:
                            shifted_data = ds_data
                            ref_slice = ref_data

                        # Ensure we have data to process
                        if len(shifted_data.time) > 0 and len(ref_slice.time) > 0:
                            metrics = calculate_metrics(
                                shifted_data.values,
                                ref_slice.values,
                                dataset,
                                offset,
                                threshold=threshold
                            )
                        else:
                            print(f"Warning: No data available for {dataset} with {offset} hour offset")
                            metrics = {
                                'rmse': np.nan,
                                'mae': np.nan,
                                'bias': np.nan,
                                'pearson_r': np.nan,
                                'pod': np.nan,
                                'far': np.nan,
                                'csi': np.nan,
                                'freq_bias': np.nan,
                                'mean_forecast': np.nan,
                                'mean_observed': np.nan,
                                'offset_hours': offset,
                                'valid_points': 0
                            }

                        results[f"{dataset}_offset_{offset}"] = metrics

                except FileNotFoundError:
                    print(f"Warning: Could not find file {ds_file}")
                    continue
                except Exception as e:
                    print(f"Error processing {dataset}: {str(e)}")
                    continue

            all_results[year] = results

        except FileNotFoundError:
            print(f"Warning: Could not find reference file {ref_file}")
            continue
        except Exception as e:
            print(f"Error processing year {year}: {str(e)}")
            continue

    # Create summary DataFrame
    summary = []
    for year in years:
        if year in all_results:
            for dataset in test_datasets:
                for offset in time_offsets:
                    key = f"{dataset}_offset_{offset}"
                    if key in all_results[year]:
                        metrics = all_results[year][key]
                        metrics['year'] = year
                        metrics['dataset'] = dataset
                        summary.append(metrics)

    df = pd.DataFrame(summary)

    return df, all_results

# Example usage
if __name__ == "__main__":

    start_year = 2001
    end_year = 2016
    reference_dataset = 'IMERG'  # Reference dataset
    test_datasets = ['ERA5','JRA', 'CFSR', 'CR20', 'TEST']  # Datasets to verify
    datasets_to_convert = []  # Datasets needing conversion
    time_offsets = [-6, 0, 6]  # Time offsets in hours
    threshold = 1.0  # mm/day
    FZRAPATH = '/glade/derecho/scratch/zarzycki/FZRA/precip/'

    # Run verification
    df, all_results = process_verification(
        start_year=start_year,
        end_year=end_year,
        reference_dataset=reference_dataset,
        test_datasets=test_datasets,
        datasets_to_convert=datasets_to_convert,
        data_directory=FZRAPATH,
        time_offsets=time_offsets,
        threshold=threshold
    )

    # Print summary statistics only if we have results
    if not df.empty:
        print("\nFinal Summary Statistics:")
        print("========================")
        metrics_to_show = ['rmse', 'pearson_r', 'bias', 'pod', 'far', 'csi',
                          'freq_bias', 'mean_forecast', 'mean_observed', 'valid_points']

        # Group by both dataset and offset
        summary_stats = df.groupby(['dataset', 'offset_hours'])[metrics_to_show].mean().round(3)
        print("\nMetrics by dataset and time offset:")
        print(summary_stats)

        # Find best offset for each dataset based on RMSE
        best_offsets = df.groupby('dataset').apply(lambda x: x.loc[x['rmse'].idxmin()])
        print("\nBest time offset for each dataset (based on minimum RMSE):")
        print(best_offsets[['offset_hours', 'rmse', 'pearson_r', 'valid_points']])
    else:
        print("No valid results to display")
