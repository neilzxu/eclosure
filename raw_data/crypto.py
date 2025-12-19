from pathlib import Path
import sys

import pandas as pd


PY_SRC = Path("python") / "src"
if str(PY_SRC) not in sys.path:
    sys.path.insert(0, str(PY_SRC))

def data_dir(kind: str) -> Path:
    return Path("data") / kind

def process_dataset(filepath, return_column, dead_indicator=None, dataset_name=""):
    """
    Process a single dataset file and extract e-values.
    
    Args:
        filepath: Path to the CSV file
        return_column: Name of the column containing returns
        dead_indicator: Function to determine if a coin is dead
        dataset_name: Name for logging purposes
    
    Returns:
        tuple: (evalues_list, live_count, dead_count)
    """
    df = pd.read_csv(filepath)
    evalues = []
    dead_count = 0
    live_count = 0
    
    for _, row in df.iterrows():
        # Skip empty rows (rows where the return column is NaN and no coin name)
        if pd.isna(row.get(return_column)) and pd.isna(row.get('Name')):
            continue
            
        if dead_indicator and dead_indicator(row):
            # Coin is dead, set e-value to 1
            evalues.append(1.0)
            dead_count += 1
        elif pd.notna(row.get(return_column)):
            # Use return value
            evalues.append(row[return_column])
            live_count += 1
    
    print(f"  {dataset_name}: {len(evalues)} e-values (Live: {live_count}, Dead: {dead_count})")
    return evalues, live_count, dead_count

def save_evalues(evalues, output_path):
    """Save e-values to CSV file."""
    output_df = pd.DataFrame({'evalue': evalues})
    output_df.to_csv(output_path, index=False)

def extract_evalues_from_crypto_data():
    """
    Extract e-values from crypto data files and save them as separate CSV files.
    
    Extracts 9 datasets:
    - 3 from buy and hold (bah) files
    - 3 from 50% rebalancing (default rebalance files)  
    - 3 from 30% rebalancing (using 30% columns from rebalance files)
    """
    
    # Create output directory
    output_dir = data_dir("evalues")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # File configurations
    datasets = [
        # Buy and hold datasets
        ('bah_100.csv', 'Return', lambda row: pd.notna(row.get('Dead?')) and row['Dead?'] == 1, 'bah_100_evalues.csv'),
        ('bah_200.csv', 'Return', lambda row: pd.notna(row.get('Dead?')) and row['Dead?'] == 1, 'bah_200_evalues.csv'),
        ('bah_426.csv', 'Return', lambda row: pd.notna(row.get('Dead?')) and row['Dead?'] == 1, 'bah_426_evalues.csv'),
        
        # 50% rebalancing datasets
        ('rebalance_100.csv', 'Return', lambda row: pd.notna(row.get('Return')) and row['Return'] == 0, 'rebalance_50pct_100_evalues.csv'),
        ('rebalance_200.csv', 'Return', lambda row: pd.notna(row.get('Return')) and row['Return'] == 0, 'rebalance_50pct_200_evalues.csv'),
        ('rebalance_426.csv', 'Return', lambda row: pd.notna(row.get('Return')) and row['Return'] == 0, 'rebalance_50pct_426_evalues.csv'),
        
        # 30% rebalancing datasets
        ('rebalance_100.csv', 'Return (30%)', lambda row: pd.notna(row.get('Return (30%)')) and row['Return (30%)'] == 0, 'rebalance_30pct_100_evalues.csv'),
        ('rebalance_200.csv', 'Return (30%)', lambda row: pd.notna(row.get('Return (30%)')) and row['Return (30%)'] == 0, 'rebalance_30pct_200_evalues.csv'),
        ('rebalance_426.csv', 'Return (30%)', lambda row: pd.notna(row.get('Return (30%)')) and row['Return (30%)'] == 0, 'rebalance_30pct_426_evalues.csv'),
    ]
    
    print("Extracting e-values from crypto data...")
    print("\nProcessing datasets:")
    
    # Process each dataset
    crypto_root = Path("raw_data")
    for filename, return_col, dead_func, output_filename in datasets:
        filepath = crypto_root / filename
        
        # Process the dataset
        evalues, live_count, dead_count = process_dataset(
            filepath, return_col, dead_func, output_filename.replace('.csv', '')
        )
        
        # Save to CSV
        output_path = output_dir / output_filename
        save_evalues(evalues, output_path)
    
    print(f"\nExtraction complete! All files saved to '{output_dir}' directory.")
    
    # Print summary
    print("\nExtracted datasets:")
    print("Buy and Hold:")
    print("  - bah_100_evalues.csv")
    print("  - bah_200_evalues.csv") 
    print("  - bah_426_evalues.csv")
    print("50% Rebalancing:")
    print("  - rebalance_50pct_100_evalues.csv")
    print("  - rebalance_50pct_200_evalues.csv")
    print("  - rebalance_50pct_426_evalues.csv")
    print("30% Rebalancing:")
    print("  - rebalance_30pct_100_evalues.csv")
    print("  - rebalance_30pct_200_evalues.csv")
    print("  - rebalance_30pct_426_evalues.csv")

if __name__ == "__main__":
    extract_evalues_from_crypto_data()
