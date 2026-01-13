from dataclasses import dataclass, field
from typing import List

import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import fdrcorrection
from scipy import stats
from sklearn.preprocessing import StandardScaler
import numpy as np
from .vendor import EmpiricalBrownsMethod
from functools import reduce

@dataclass
class ModelParams:
    """
    Define run and model parameters for NCCM enrichment analysis

    Attributes
    ----------
    y_col: str
        Column name for the dependent variable (NCCM count).
    offset_col: str
        Column name for the offset variable (NCC positions).
    covariates_list: List[str]
        List of column names for covariates to include in the regression.
    maxiter: int
        Maximum number of iterations for the regression fitting.
    method: str
        Optimization method for the regression fitting.
    randomized_p: bool
        Whether to compute randomized p-values (strongly recommended).
    drop_raw_p: bool
        Whether to drop raw p-values from the output (for debugging purposes).
    n_bins: int
        Number of bins to split the data into for binned regression.
    binning_col: str
        Column name to use for binning the data.
    fdr_thresh: float
        FDR threshold for multiple testing correction.
    """
    # Regression parameters
    y_col: str = 'nccm_count'
    offset_col: str = 'nccp'
    background_mutation_count_col = 'ncncm_count'
    covariates_list: List[str] = field(default_factory=lambda:
                    ['scaled_log_ncncm_rate','scaled_log_perc_gc','scaled_log_perc_cpg',
                    'scaled_log_perc_tpc', 'scaled_log_perc_tpa','scaled_log_perc_open',
                    'scaled_yeo_rep','scaled_yeo_phylop','scaled_log_expr'])
    maxiter: int = 1000
    method: str = 'bfgs'

    # Post-hoc test parameters
    randomized_p: bool = True
    drop_raw_p: bool = True

    # Binning stragegy
    n_bins: int = 10
    binning_col: str = 'ncncm_rate'

    # Gene-centric approach
    flank_suffixes: List[str] = field(default_factory=lambda: ['_f100', '_f50', '_f10', '_f5'])

    # Multiple testing
    fdr_thresh: float = 0.1

def get_regions_to_exclude(df, min_pass_fraction: float = 0.9, col_overlap: str = 'overlap_fraction'):
    """
    Identify regions to exclude based on fraction of region passing quality filters.

    Parameters
    ----------
    df: pd.DataFrame
        Input dataframe with 'gene' and 'overlap_fraction' columns.
    min_pass_fraction: float
        Minimum fraction of region that must pass quality filters to be included. Default is 0.9.
    col_overlap: str
        Column name for overlap fraction. Default is 'overlap_fraction'.

    Returns
    -------
    exclude_regions: List[str]
        List of regions to exclude.
    """
    exclude_regions = df[df[col_overlap] < min_pass_fraction]['gene'].unique().tolist()

    return exclude_regions

def load_data(file_path: str):
    """
    Load covariate data from a TSV file.

    Parameters
    ----------
    file_path: str
        Path to the TSV file containing covariate data.

    Returns
    -------
    df_covariates: pd.DataFrame
        Dataframe containing the covariate data.
    """
    df_covariates = pd.read_csv(file_path, sep='\t')
    return df_covariates

def merge_dfs_on_gene(*dataframes):
    """
    Merge multiple dataframes on the 'gene' column.

    Parameters
    ----------
    *dataframes : pd.DataFrame
        Variable number of input dataframes to merge.
    
    Returns
    -------
    df_merged: pd.DataFrame
        Merged dataframe.
    """
    df_merged = reduce(
        lambda left, right: pd.merge(left, right, on='gene', how='inner'),
        dataframes
    )
    return df_merged

def run_gamma_poisson_regression(df, params: ModelParams = None, **kwargs):
    """
    Run a gamma-poisson regression for NCCM count.
    
    Parameters
    ----------
    df: pd.DataFrame
        Input dataframe with gene names, NCCM counts, NCNCM counts, NCC positions (offset), NCNCM positions and all covariates.
    params: ModelParams
        Model parameters for the regression.

    Returns
    -------
    results: statsmodels regression results object
        Fitted model results.
    alpha_dispersion: float
        Alpha dispersion parameter from the fitted model.

    Raises
    ------
    Exception
        If the regression fails to converge.
    """
    df_out = df.copy().reset_index(drop=True)

    # Default parameters if none provided
    if params is None:
        params = ModelParams()

    # Overwrite parameters with any extra keywords provided
    for key, value in kwargs.items():
        if hasattr(params, key):
            setattr(params, key, value)
        else:
            raise ValueError(f"Unknown parameter: {key}")

    # Define dependent variable as NCCM count
    y = df_out[params.y_col].astype(np.float64)

    # Generage background mutation rate covariate
    if ('scaled_log_ncncm_rate' in params.covariates_list):
        if ('ncncm_count' in df_out.columns) & ('ncncp' in df_out.columns):
            # Calculate background mutation rate
            df_out['ncncm_rate'] = df_out[params.background_mutation_count_col]/df_out['ncncp']

            # Scaled log transformed background mutation rate
            df_out['log_ncncm_rate'] = np.log1p(df_out['ncncm_rate'])
            scaler = StandardScaler()
            df_out['scaled_log_ncncm_rate']= scaler.fit_transform(df_out[['log_ncncm_rate']])

    # Define the covariates and the intercept to be estimated
    x = pd.DataFrame(index=df_out.index)
    x['intercept'] = 1.0
    x[params.covariates_list] = df_out[params.covariates_list]
    x = x.astype(np.float64)

    # Define offset as the number of NCC positions
    offset = np.log(df_out[params.offset_col]).astype(np.float64)

    # Run the Regression
    model = sm.NegativeBinomial(y, x, offset = offset)
    results = model.fit(maxiter=params.maxiter, method = params.method)

    # Calculate alpha dispersion (for posthoc testing)
    alpha_dispersion = results.params['alpha'] #len(params.covariates_list)+1

    # Check if the model converged
    if results.converged:
        return results, alpha_dispersion
    else:
        raise Exception("Gamma poisson regression failed to converge. Exiting.")

def test_nccm_enrichment(df, results, alpha_dispersion, params: ModelParams = None, **kwargs):
    """
    Test for NCCM enrichment based on the gamma-poisson regression results.

    Parameters
    ----------
    df: pd.DataFrame
        Input dataframe with gene names, NCCM counts, NCC positions (offset) and all covariates.
    results: statsmodels regression results object from run_gamma_poisson_regression.
        Fitted model results.
    alpha_dispersion: float
        Alpha dispersion parameter from the fitted model from run_gamma_poisson_regression.
    params: ModelParams 
        Model parameters for the regression.
    """
    df_out = df.copy()

    # Default parameters if none provided
    if params is None:
        params = ModelParams()
    
    # Overwrite parameters with any extra keywords provided
    for key, value in kwargs.items():
        if hasattr(params, key):
            setattr(params, key, value)
        else:
            raise ValueError(f"Unknown parameter: {key}")

    # Extract predicted number of NCCMs based on the model and alpha dispersion
    predicted_mean = results.predict()
    observed_counts = df_out[params.y_col]

    # Calculate p-values using the Negative Binomial survival function
    n_param = 1 / alpha_dispersion
    p_param = 1 / (1 + alpha_dispersion * predicted_mean)

    # Raw p-values
    raw_p_values = stats.nbinom.sf(observed_counts - 1, n=n_param, p=p_param)

    # Upper and lower bounds for randomized p-values
    p_upper = stats.nbinom.sf(observed_counts - 1, n=n_param, p=p_param)
    p_lower = stats.nbinom.sf(observed_counts, n=n_param, p=p_param)

    # Compile results
    df_results = df_out.copy()
    df_results['predicted_'+params.y_col] = predicted_mean
    df_results['raw_p_value'] = raw_p_values
    df_results['raw_p_upper'] = p_upper
    df_results['raw_p_lower'] = p_lower

    return df_results

def multiple_testing_correction(df, params: ModelParams = None, p_col = 'p_value', **kwargs):
    """
    Perform multiple testing correction using FDR.

    Parameters
    ----------
    df: pd.DataFrame
        Input dataframe with p-values.
    params: ModelParams
        Model parameters including FDR threshold.
    p_col: str
        Column name for the p-values to correct.
    """
    # Default parameters if none provided
    if params is None:
        params = ModelParams()

    # Overwrite parameters with any extra keywords provided
    for key, value in kwargs.items():
        if hasattr(params, key):
            setattr(params, key, value)
        else:
            raise ValueError(f"Unknown parameter: {key}")
    
    # Run FDR correction
    df_out = df.copy()

    q_col = p_col.replace('p','q')
    df_out[q_col] = fdrcorrection(df_out[p_col], alpha=params.fdr_thresh)[1]
    df_out = df_out.sort_values(by=[q_col], ascending=True)

    return df_out

def nccm_enrichment_analysis(df, params: ModelParams = None, **kwargs):
    """
    Perform NCCM enrichment analysis using a binned gamma-poisson regression approach.

    Parameters
    ----------
    df: pd.DataFrame
        Input dataframe with gene names, NCCM counts, NCC positions (offset) and all covariates.
    params: ModelParams
        Model parameters for the regression.
    """
    # Default parameters if none provided
    if params is None:
        params = ModelParams()
    
    # Overwrite parameters with any extra keywords provided
    for key, value in kwargs.items():
        if hasattr(params, key):
            setattr(params, key, value)
        else:
            raise ValueError(f"Unknown parameter: {key}")

    # Split the data into bins based on the binning column
    dfs_bin = np.array_split(df.sort_values([params.binning_col, 'gene']), params.n_bins)

    # Run the analysis for each bin
    df_results_total = pd.DataFrame()
    for i, df_bin in enumerate(dfs_bin):

        # Run regression and post-hoc testing
        results, alpha_dispersion = run_gamma_poisson_regression(df_bin, params)

        df_results = test_nccm_enrichment(df_bin, results, alpha_dispersion, params)

        # Append results
        df_results_total = pd.concat([df_results_total, df_results])

    # To randomize or not to randomize
    if params.randomized_p:
        np.random.seed(410)
        df_results_total['p_value'] = np.random.uniform(df_results_total['raw_p_lower'], 
                                                        df_results_total['raw_p_upper'])
    else:
        df_results_total['p_value'] = df_results_total['raw_p_value']

    # Drop raw p-values if specified
    if params.drop_raw_p:
        df_results_total = df_results_total.drop(columns = ['raw_p_value','raw_p_upper',
                                                            'raw_p_lower']).reset_index(drop=True)

    return df_results_total

def combine_pvalues_browns_method(df, params: ModelParams = None, **kwargs):
    """
    Combine p-values from multiple window sizes using Empirical Brown's Method.

    Parameters
    ----------
    df: pd.DataFrame
        Input dataframe with p-values and covariate data for multiple window sizes.
    params: ModelParams
        Model parameters including covariate names and p-value prefixes.
    
    Returns
    -------
    df_out: pd.DataFrame
        Dataframe with an additional column 'p_brown' containing combined p-values.
    """
    df_out = df.copy()

    # Default parameters if none provided
    if params is None:
        params = ModelParams()

    # Overwrite parameters with any extra keywords provided
    for key, value in kwargs.items():
        if hasattr(params, key):
            setattr(params, key, value)
        else:
            raise ValueError(f"Unknown parameter: {key}")

    # Pre-calculate column names to save time inside the loop
    data_cols_map = {
        suffix: [f"{cov}{suffix}" for cov in params.covariates_list]
        for suffix in params.flank_suffixes
    }
    p_val_cols = [f"p_value{suffix}" for suffix in params.flank_suffixes]

    # Placeholder for results
    brown_p_values = []
    
    for i, row in df_out.iterrows():        
        # Collect the data vectors for each window size
        data_vectors = []
        for suffix in params.flank_suffixes:
            # Extract the specific columns for this suffix
            cols = data_cols_map[suffix]
            vector = row[cols].values.astype(float)
            data_vectors.append(vector)
            
        # Create matrix expected by Brown's method
        data_matrix = np.array(data_vectors)

        # Collect the p-values
        current_p_values = row[p_val_cols].values.astype(float)

        # Handle potential NaNs in p-values
        if np.isnan(current_p_values).any():
            brown_p_values.append(np.nan)
            continue

        # Run Brown's Method
        try:
            p_brown, p_fisher, scale, degf = EmpiricalBrownsMethod(
                data_matrix, 
                current_p_values, 
                extra_info=True
            )
            brown_p_values.append(p_brown)
        except Exception:
            # Fallback if calculation fails (e.g., singular matrix)
            brown_p_values.append(np.nan)

    # Assign results back to the dataframe
    df_out['p_brown'] = brown_p_values

    return df_out

def calculate_lambda_onesided(p_values):
    """
    Calculates the genomic inflation factor (lambda) from one-sided p-values.

    Parameters:
    ----------- 
    p_values: array-like
         A list or numpy array of one-sided p-values.
    
    Returns:
    --------
    float
        The genomic inflation factor (lambda).
    """
    # P-values as array
    p_values = np.asarray(p_values)

    # Convert one-sided p-values to z-scores
    z_scores = stats.norm.ppf(1 - p_values)

    # Square the z-scores to get chi-squared statistics with one degree of freedom
    chisq_stats = z_scores**2

    # Calculate lambda 
    lambda_gc = np.median(chisq_stats) / stats.chi2.ppf(0.5,1)

    return lambda_gc

def apply_genomic_control_zscore(p_values):
    """
    Applies genomic control to combined p-values. 

    Parameters:
    ----------- 
    p_values: array-like
         A list or numpy array of one-sided p-values to be corrected.

    Returns:
    --------
    numpy.ndarray
        A numpy array containing the corrected p-values.
    """
    p_values = np.asarray(p_values)
    # Convert p-values to z-scores and handle edge cases
    eps = np.finfo(float).eps 
    p_values = np.clip(p_values, eps, 1 - eps)    
    z_scores = stats.norm.ppf(1 - p_values)

    # Calculate lambda for one-sided p-values
    lambda_gc = calculate_lambda_onesided(p_values)
    
    # Correct z-scores using lambda
    corrected_z_scores = z_scores / np.sqrt(lambda_gc)

    # Convert corrected z-scores back to p-values
    adjusted_log_p_values = stats.norm.logsf(corrected_z_scores)
    adjusted_p_values = np.exp(adjusted_log_p_values)

    return adjusted_p_values

def window_based_approach(df, params: ModelParams = None, **kwargs):
    """
    Perform NCCM enrichment analysis using a window-based approach.
    
    Parameters
    ----------
    df: pd.DataFrame
        Input dataframe with gene names, NCCM counts, NCC positions (offset) and all covariates.
    params: ModelParams
        Model parameters for the regression.
    """
    # Default parameters if none provided
    if params is None:
        params = ModelParams()

    # Overwrite parameters with any extra keywords provided
    for key, value in kwargs.items():
        if hasattr(params, key):
            setattr(params, key, value)
        else:
            raise ValueError(f"Unknown parameter: {key}")
    
    results = nccm_enrichment_analysis(df, params)
    results = multiple_testing_correction(results, params)
    return results

def gene_centric_approach(*dataframes, params: ModelParams = None, **kwargs):
    """
    Perform NCCM enrichment analysis for any number of flank sizes.
    
    The function maps the input dataframes to the suffixes defined in 
    params.flank_suffixes based on their order.
    
    Parameters
    ----------
    *dataframes : pd.DataFrame
        Variable number of input dataframes (e.g., df100, df50, df10...).
        MUST be passed in the same order as params.flank_suffixes.
    params : ModelParams
        Model parameters containing the list of suffixes.

    Returns
    -------
    merged_results : pd.DataFrame
        Combined results table.
    """
    # Default parameters if none provided
    if params is None:
        params = ModelParams()
    
    # Overwrite parameters with any extra keywords provided
    for key, value in kwargs.items():
        if hasattr(params, key):
            setattr(params, key, value)
        else:
            raise ValueError(f"Unknown parameter: {key}")
    
    # Check if number of dfs matches number of suffixes
    if len(dataframes) != len(params.flank_suffixes):
        raise ValueError(
            f"Mismatch in arguments! "
            f"params.flank_suffixes defines {len(params.flank_suffixes)} items "
            f"({params.flank_suffixes}), but you provided {len(dataframes)} dataframes."
        )

    # Map dataframes to suffixes 
    data_map = dict(zip(params.flank_suffixes, dataframes))
    
    results_list = []

    # Run analysis for each dataframe
    for suffix, df_input in data_map.items():
        # Run core analysis
        df_res = nccm_enrichment_analysis(df_input, params)

        # Rename columns (Everything except 'gene')
        rename_mapping = {
            col: f"{col}{suffix}" 
            for col in df_res.columns 
            if col != 'gene'
        }
        df_res_renamed = df_res.rename(columns=rename_mapping)
        
        results_list.append(df_res_renamed)

    # Merge sequentially
    merged_results = reduce(
        lambda left, right: pd.merge(left, right, on='gene', how='inner'),
        results_list
    )

    # Combine p-values (Brown's method)
    merged_results = combine_pvalues_browns_method(merged_results, params)

    # Genomic control and correction for multiple testing
    merged_results['p_brown_gc'] = apply_genomic_control_zscore(merged_results['p_brown'])
    merged_results = multiple_testing_correction(merged_results, params, p_col='p_brown_gc')

    return merged_results