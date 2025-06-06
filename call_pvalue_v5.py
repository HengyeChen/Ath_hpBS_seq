import sys
import numpy as np
import scipy.stats as stats
from multiprocessing import Pool
from typing import Tuple, List
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def compute_ttest(i: int) -> Tuple[int, float, float, float, float, int]:
    """
    Compute t-test and Mann-Whitney U test for paired data.
    
    Args:
        i: Index value to select data
        
    Returns:
        Tuple containing:
        - index
        - t-test p-value
        - Mann-Whitney U test p-value 
        - mean of first group
        - mean of second group
        - sample size
    """
    try:
        data1 = data[data[:,idx]==i,r1]
        data2 = data[data[:,idx]==i,r2]
        m1=np.mean(data1)
        m2=np.mean(data2)
        sig=0.00000000000000000000000001
        if len(data1) != len(data2):
            return None
            
        if len(data1) == 0:
            return None
            
        if len(data1) < 5:
            return None
        
        if m1==0 and m2==0:
            return None
        
        if m1==0 and m2!=0 and len(data1)>10 and m2>0.3:
            return i, sig, sig, m1, m2, len(data1)

        if m1!=0 and m2==0 and len(data1)>10 and m1>0.3:
            return i, sig, sig, m1, m2, len(data1)

        _, pvalue1 = stats.ttest_rel(data1, data2)
        _, pvalue2 = stats.mannwhitneyu(data1, data2, alternative='two-sided')
        
        if np.isnan(pvalue1) or np.isnan(pvalue2):
            return None
        
        return i, pvalue1, pvalue2, m1, m2, len(data1)
        
    except Exception as e:
        logging.error(f"Error processing index {i}: {str(e)}")
        return None

def process_results(results: List[tuple], output_file: str) -> None:
    """
    Process and save test results.
    
    Args:
        results: List of result tuples from compute_ttest
        output_file: Path to save results
    """
    # Filter out None results
    results = [r for r in results if r is not None]
    
    result_array = np.array(results)
    
    # Calculate adjusted p-values using Benjamini-Hochberg FDR correction
    pvals_ttest = result_array[:,1]
    pvals_mw = result_array[:,2]
    
    adj_pvals_ttest = stats.false_discovery_control(pvals_ttest)
    adj_pvals_mw = stats.false_discovery_control(pvals_mw)
    
    # Create output array with adjusted p-values
    output_array = np.column_stack((
        result_array[:,0],  # index
        adj_pvals_ttest,    # adjusted t-test p-value 
        adj_pvals_mw,       # adjusted Mann-Whitney p-value
        result_array[:,3],  # mean1
        result_array[:,4],  # mean2
        result_array[:,5]   # count
    ))
    
    # Write results with headers
    header = "Index\tAdj_Ttest_pval\tAdj_MW_pval\tMean1\tMean2\tCount"
    np.savetxt(output_file, output_array, fmt='%d\t%.20f\t%.8f\t%.8f\t%.8f\t%d',
               header=header, comments='')

if __name__ == '__main__':
    try:
        # Validate command line arguments
        if len(sys.argv) != 6:
            raise ValueError("Expected 5 arguments: input_file output_file r1 r2 idx")
            
        name, op = sys.argv[1:3]
        r1, r2, idx = map(lambda x: int(x)-1, sys.argv[3:6])
        
        # Load data efficiently using numpy
        logging.info(f"Loading data from {name}")
        data = np.loadtxt(name, dtype=np.float32)
        
        # Input validation
        if data.size == 0:
            raise ValueError("Empty input file")
            
        # Get unique indices
        unique_idx = np.unique(data[:,idx])
        
        # Configure multiprocessing
        n_cores = min(80, len(unique_idx))  # Don't use more processes than necessary
        
        logging.info(f"Processing {len(unique_idx)} unique indices using {n_cores} cores")
        
        # Parallel processing with context manager
        with Pool(processes=n_cores) as pool:
            results = pool.map(compute_ttest, unique_idx)
        
        logging.info("Processing results")
        process_results(results, op)
        logging.info(f"Results saved to {op}")
        
    except Exception as e:
        logging.error(f"Error: {str(e)}")
        sys.exit(1)
