import numpy as np
import pickle

def read_parsec_isochrones(filename):
    """
    Read PARSEC isochrone file with multiple isochrones separated by # lines.
    
    Parameters:
    -----------
    filename : str
        Path to the PARSEC output file
    
    Returns:
    --------
    isochrones : list of dict
        Each dict contains numpy arrays for one isochrone with keys:
        'Gmag', 'G_BPmag', 'G_RPmag', 'BP_RP', 'MH', 'logAge', 'Mass', 
        'logL', 'logTe', 'logR', 'R_Rsun', 'age_Myr'
    """
    
    # Columns we want to extract
    keep_cols = ['Gmag', 'G_BPmag', 'G_RPmag', 'MH', 'logAge', 'Mass', 'logL', 'logTe']
    
    isochrones = []
    
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # Find header line with column names
    header_line = None
    header_idx = None
    for i, line in enumerate(lines):
        if line.startswith('# Zini'):
            header_line = line[2:].strip()  # Remove '# ' prefix
            header_idx = i
            break
    
    if header_line is None:
        raise ValueError("Could not find header line starting with '# Zini'")
    
    # Parse column names
    column_names = header_line.split()
    
    # Get indices of columns we want
    col_indices = {}
    for col in keep_cols:
        try:
            col_indices[col] = column_names.index(col)
        except ValueError:
            print(f"Warning: Column '{col}' not found in file")
            col_indices[col] = None
    
    # Read data after header
    current_iso_data = []
    current_mh = None
    current_logage = None
    
    for i in range(header_idx + 1, len(lines)):
        line = lines[i].strip()
        
        # Skip empty lines
        if not line:
            continue
        
        # New isochrone marker
        if line.startswith('#'):
            # Save previous isochrone if it exists
            if current_iso_data:
                iso_dict = save_isochrone(current_iso_data, col_indices, 
                                         current_mh, current_logage)
                isochrones.append(iso_dict)
                
                # Reset for new isochrone
                current_iso_data = []
            
            continue
        
        # Parse data line
        values = line.split()
        
        # Convert to floats
        row = []
        for val in values:
            try:
                row.append(float(val))
            except ValueError:
                row.append(np.nan)
        
        # Extract MH and logAge from first row
        if not current_iso_data:
            if col_indices['MH'] is not None:
                current_mh = row[col_indices['MH']]
            if col_indices['logAge'] is not None:
                current_logage = row[col_indices['logAge']]
        
        current_iso_data.append(row)
    
    # Don't forget last isochrone
    if current_iso_data:
        iso_dict = save_isochrone(current_iso_data, col_indices, 
                                 current_mh, current_logage)
        isochrones.append(iso_dict)
    
    return isochrones


def save_isochrone(data_rows, col_indices, mh, logage):
    """
    Convert list of rows to dictionary of numpy arrays.
    
    Parameters:
    -----------
    data_rows : list of lists
        Raw data rows
    col_indices : dict
        Column name to index mapping
    mh : float
        Metallicity [M/H]
    logage : float
        log(age/yr)
    
    Returns:
    --------
    iso_dict : dict
        Dictionary with numpy arrays for each column
    """
    # Convert to numpy array
    data = np.array(data_rows)
    
    iso_dict = {}
    
    # Extract each column
    for col_name, idx in col_indices.items():
        if idx is not None:
            iso_dict[col_name] = data[:, idx]
        else:
            # Column not found, create NaN array
            iso_dict[col_name] = np.full(len(data), np.nan)
    
    # Calculate derived quantities
    # BP-RP color
    if 'G_BPmag' in iso_dict and 'G_RPmag' in iso_dict:
        iso_dict['BP_RP'] = iso_dict['G_BPmag'] - iso_dict['G_RPmag']
    
    # Stellar radius: R/R_sun = sqrt(L/L_sun) / (T/T_sun)^2
    # log(R/R_sun) = 0.5 * log(L/L_sun) - 2 * log(T/T_sun)
    # log(T/T_sun) = log(T) - log(5772)
    if 'logL' in iso_dict and 'logTe' in iso_dict:
        T_sun = 5772  # K
        iso_dict['logR'] = 0.5 * iso_dict['logL'] - 2 * (iso_dict['logTe'] - np.log10(T_sun))
        iso_dict['R_Rsun'] = 10**iso_dict['logR']
    
    # Add age in Myr
    if logage is not None:
        iso_dict['age_Myr'] = 10**(logage - 6)
    else:
        iso_dict['age_Myr'] = np.nan
    
    # Add metadata (same value for all stars in this isochrone)
    iso_dict['MH_value'] = mh if mh is not None else np.nan
    iso_dict['logAge_value'] = logage if logage is not None else np.nan
    
    return iso_dict


def save_isochrones_pickle(isochrones, output_dir='.'):
    """
    Save isochrones to pickle files.
    
    Two options:
    1. Save all isochrones in one file
    2. Save each isochrone separately
    
    Parameters:
    -----------
    isochrones : list of dict
        Isochrone data
    output_dir : str
        Directory to save pickle files
    """
    import os
    
    # Option 1: Save all isochrones in one file
    all_iso_file = os.path.join(output_dir, 'all_isochrones.pkl')
    with open(all_iso_file, 'wb') as f:
        pickle.dump(isochrones, f)
    print(f"Saved all isochrones to: {all_iso_file}")
    
    # Option 2: Save each isochrone separately
    for i, iso in enumerate(isochrones):
        mh = iso['MH_value']
        age_myr = iso['age_Myr']
        filename = f'isochrone_{i:02d}_MH{mh:+.2f}_age{age_myr:.1f}Myr.pkl'
        filepath = os.path.join(output_dir, filename)
        
        with open(filepath, 'wb') as f:
            pickle.dump(iso, f)
        print(f"Saved isochrone {i}: {filename}")
    
    # Also save a summary file
    summary = []
    for i, iso in enumerate(isochrones):
        summary.append({
            'index': i,
            'MH': iso['MH_value'],
            'logAge': iso['logAge_value'],
            'age_Myr': iso['age_Myr'],
            'n_stars': len(iso['Mass']),
            'mass_range': (np.min(iso['Mass']), np.max(iso['Mass'])),
            'filename': f'isochrone_{i:02d}_MH{iso["MH_value"]:+.2f}_age{iso["age_Myr"]:.1f}Myr.pkl'
        })
    
    summary_file = os.path.join(output_dir, 'isochrones_summary.pkl')
    with open(summary_file, 'wb') as f:
        pickle.dump(summary, f)
    print(f"\nSaved summary to: {summary_file}")
    
    return summary


def load_isochrones_pickle(filename):
    """
    Load isochrones from pickle file.
    
    Parameters:
    -----------
    filename : str
        Path to pickle file
    
    Returns:
    --------
    data : list or dict
        Loaded isochrone data
    """
    with open(filename, 'rb') as f:
        data = pickle.load(f)
    return data


def print_summary(isochrones):
    """
    Print summary of loaded isochrones.
    """
    print(f"\nLoaded {len(isochrones)} isochrones:")
    print("-" * 80)
    
    for i, iso in enumerate(isochrones):
        print(f"\nIsochrone {i}:")
        print(f"  [M/H] = {iso['MH_value']:.3f}")
        print(f"  log(Age) = {iso['logAge_value']:.3f}")
        print(f"  Age = {iso['age_Myr']:.1f} Myr")
        print(f"  Number of stars: {len(iso['Mass'])}")
        print(f"  Mass range: {np.min(iso['Mass']):.3f} - {np.max(iso['Mass']):.3f} M_sun")
        
        if 'BP_RP' in iso and not np.all(np.isnan(iso['BP_RP'])):
            bp_rp_min = np.nanmin(iso['BP_RP'])
            bp_rp_max = np.nanmax(iso['BP_RP'])
            print(f"  BP-RP range: {bp_rp_min:.3f} - {bp_rp_max:.3f}")
        
        if 'R_Rsun' in iso and not np.all(np.isnan(iso['R_Rsun'])):
            r_min = np.nanmin(iso['R_Rsun'])
            r_max = np.nanmax(iso['R_Rsun'])
            print(f"  Radius range: {r_min:.3f} - {r_max:.3f} R_sun")
        
        print(f"  Available arrays: {[k for k in iso.keys() if isinstance(iso[k], np.ndarray)]}")


# Example usage
if __name__ == "__main__":
    # Read the PARSEC file
    filename = 'output646247342069.dat.txt'  # Replace with your filename
    
    print("Reading PARSEC isochrones...")
    isochrones = read_parsec_isochrones(filename)
    
    # Print summary
    print_summary(isochrones)
    
    # Save to pickle files
    print("\n" + "="*80)
    print("Saving to pickle files...")
    summary = save_isochrones_pickle(isochrones, output_dir='.')
    
    # Example: Load back the data
    print("\n" + "="*80)
    print("Example: Loading data back...")
    loaded_isochrones = load_isochrones_pickle('all_isochrones.pkl')
    print(f"Successfully loaded {len(loaded_isochrones)} isochrones")
    
    # Example: Access specific isochrone data
    print("\n" + "="*80)
    print("Example: Accessing first isochrone data...")
    iso_0 = loaded_isochrones[0]
    print(f"First 5 masses: {iso_0['Mass'][:5]}")
    print(f"First 5 G magnitudes: {iso_0['Gmag'][:5]}")
    print(f"First 5 BP-RP colors: {iso_0['BP_RP'][:5]}")
    
    # Example: Load individual isochrone file
    if summary:
        first_file = summary[0]['filename']
        print(f"\nLoading individual file: {first_file}")
        single_iso = load_isochrones_pickle(first_file)
        print(f"Loaded isochrone with {len(single_iso['Mass'])} stars")