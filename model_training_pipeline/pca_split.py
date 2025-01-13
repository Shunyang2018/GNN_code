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
from sklearn.cluster import KMeans
def pca_split_and_save(df, output_dir='output', n_clusters=5, test_size=0.2, random_state=999):
    """
    Standardize the dataset, perform PCA, split into training and test sets with similar chemical spaces
    using stratified sampling, and save results as standardized train/test CSV files.

    Args:
        df (pd.DataFrame): Input DataFrame containing only features.
        output_dir (str, optional): Directory to save the output files and plot. Defaults to 'output'.
        n_clusters (int, optional): Number of clusters for stratified sampling. Defaults to 5.
        test_size (float, optional): Proportion of the dataset to include in the test split. Defaults to 0.2.
        random_state (int, optional): Random state for reproducibility. 

    Returns:
        tuple: Two DataFrames - standardized train and test sets.
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

    # Stratified sampling based on clustering
    print(f"Splitting the dataset using stratified sampling with {n_clusters} clusters...")
    kmeans = KMeans(n_clusters=n_clusters, random_state=random_state)
    df_tmp['cluster'] = kmeans.fit_predict(df_tmp[['pc1', 'pc2']])

    # Perform stratified split
    train_indices, test_indices = train_test_split(
        df_tmp['index'],
        test_size=test_size,
        stratify=df_tmp['cluster'],
        random_state=random_state
    )

    # Plot PCA components
    train_data = df_tmp[df_tmp['index'].isin(train_indices)]
    test_data = df_tmp[df_tmp['index'].isin(test_indices)]

    plt.scatter(train_data['pc1'], train_data['pc2'], color='none', edgecolor='green', label='Train')
    plt.scatter(test_data['pc1'], test_data['pc2'], color='none', edgecolor='red', label='Test')
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.legend(loc="upper left")
    plt.title('PCA Train/Test Split (Stratified)')
    plt.savefig(f'{output_dir}/pca_split_stratified.png')
    plt.close()

    # Split the DataFrame into train and test sets
    df_train = df.loc[train_indices]
    df_test = df.loc[test_indices]

    # Standardize train and test datasets using the same scaler
    df_train_standardized = pd.DataFrame(
        scaler.transform(df_train), 
        columns=df.columns, 
        index=df_train.index
    )
    df_test_standardized = pd.DataFrame(
        scaler.transform(df_test), 
        columns=df.columns, 
        index=df_test.index
    )

    print("Saving the standardized dataframes as CSV files...")
    # Save the standardized DataFrames as CSV files
    df_train_standardized.to_csv(f'{output_dir}/train_standardized.csv', index=False)
    df_test_standardized.to_csv(f'{output_dir}/test_standardized.csv', index=False)

    # Save train and test indices
    pd.DataFrame(train_indices, columns=['Index']).to_csv(f'{output_dir}/train_indices.csv', index=False)
    pd.DataFrame(test_indices, columns=['Index']).to_csv(f'{output_dir}/test_indices.csv', index=False)

    return df_train_standardized, df_test_standardized
