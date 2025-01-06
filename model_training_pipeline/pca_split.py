import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import pandas as pd
from sklearn.preprocessing import StandardScaler
from rdkit import Chem
from mordred import Calculator, descriptors
import os
# Function to calculate molecular descriptors using Mordred
def calculate_descriptors(df, smiles_column='CanonicalSMILES'):
    """
    Calculate molecular descriptors for the CanonicalSMILES column using Mordred.

    Args:
        df (pd.DataFrame): Input DataFrame containing a column of CanonicalSMILES.
        smiles_column (str): Name of the column containing SMILES strings.

    Returns:
        pd.DataFrame: DataFrame with molecular descriptors as additional columns.
    """
    # Create a Mordred calculator with all descriptors
    calc = Calculator(descriptors, ignore_3D=True)

    # Prepare to store descriptors
    descriptor_data = []

    for smiles in df[smiles_column]:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            # Calculate descriptors
            descriptor_values = calc(mol)
            descriptor_data.append(descriptor_values)
        else:
            # Handle invalid SMILES
            descriptor_data.append([None] * len(calc.descriptors))

    # Create a DataFrame for the descriptors
    descriptors_df = pd.DataFrame(descriptor_data, columns=[str(d) for d in calc.descriptors])

    # Concatenate original DataFrame with descriptor DataFrame
    result_df = pd.concat([df, descriptors_df], axis=1)

    return result_df
from sklearn.model_selection import train_test_split
def pca_split_and_save(df, output_dir='output', split_method='stratified', n_clusters=5, n_bins=10, test_size=0.2, random_state=999):
    """
    Standardize the dataset, perform PCA, split into training and test sets with similar chemical spaces,
    and save results as CSV files.

    Args:
        df (pd.DataFrame): Input DataFrame containing only features.
        output_dir (str, optional): Directory to save the output files and plot. Defaults to 'output'.
        split_method (str, optional): Method to split the dataset ('stratified' or 'grid'). Defaults to 'stratified'.
        n_clusters (int, optional): Number of clusters for stratified sampling. Defaults to 5.
        n_bins (int, optional): Number of bins for grid-based splitting. Defaults to 10.
        test_size (float, optional): Proportion of the dataset to include in the test split. Defaults to 0.2.
        random_state (int, optional): Random state for reproducibility. 

    Returns:
        tuple: Two DataFrames - train and test sets.
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Extract features
    features = df.iloc[:, :]

    # Standardize features
    scaler = StandardScaler()
    features_scaled = scaler.fit_transform(features)

    # Perform PCA
    pca = PCA()
    df_pca = pca.fit_transform(features_scaled)

    # Create DataFrame for PCA components
    df_tmp = pd.DataFrame(df_pca[:, :2], columns=['pc1', 'pc2'])
    df_tmp['index'] = df.index

    # Split the data
    print(f"Splitting the dataset using {split_method} method...")
    if split_method == 'stratified':
        # Stratified sampling based on clustering
        from sklearn.cluster import KMeans
        kmeans = KMeans(n_clusters=n_clusters, random_state=random_state)
        df_tmp['cluster'] = kmeans.fit_predict(df_tmp[['pc1', 'pc2']])
        train_indices, test_indices = train_test_split(
            df_tmp['index'],
            test_size=test_size,
            stratify=df_tmp['cluster'],
            random_state=random_state
        )
    elif split_method == 'grid':
        # Grid-based stratification
        pc1_bins = pd.cut(df_tmp['pc1'], bins=n_bins, labels=False)
        pc2_bins = pd.cut(df_tmp['pc2'], bins=n_bins, labels=False)
        df_tmp['bin'] = pc1_bins * n_bins + pc2_bins
        train_indices, test_indices = train_test_split(
            df_tmp['index'],
            test_size=test_size,
            stratify=df_tmp['bin'],
            random_state=random_state
        )
    else:
        raise ValueError("Invalid split_method. Choose 'stratified' or 'grid'.")

    # Plot PCA components
    train_data = df_tmp[df_tmp['index'].isin(train_indices)]
    test_data = df_tmp[df_tmp['index'].isin(test_indices)]

    plt.scatter(train_data['pc1'], train_data['pc2'], color='none', edgecolor='green', label='Train')
    plt.scatter(test_data['pc1'], test_data['pc2'], color='none', edgecolor='red', label='Test')
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.legend(loc="upper left")
    plt.title('PCA Train/Test Split')
    plt.savefig(f'{output_dir}/pca_split.png')
    plt.close()

    # Split the DataFrame into train and test sets
    df_train = df.loc[train_indices]
    df_test = df.loc[test_indices]
    print("saving the dataframes as csv files")
    # Save the DataFrames as CSV files
    df_train.to_csv(f'{output_dir}/train.csv', index=False)
    df_test.to_csv(f'{output_dir}/test.csv', index=False)

    # Save train and test indices
    pd.DataFrame(train_indices, columns=['Index']).to_csv(f'{output_dir}/train_indices.csv')
    pd.DataFrame(test_indices, columns=['Index']).to_csv(f'{output_dir}/test_indices.csv')

    return df_train, df_test
