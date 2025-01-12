## Retention Time Prediction for Small Molecules in Untargeted Metabolomics

Unidentified peaks remain a significant challenge in untargeted metabolomics using LC-MS/MS. Combining MS/MS matching with retention time significantly enhances confidence in peak annotations. Here, I demonstrate how retention times can be accurately predicted from molecular structures.
Retention times for approximately 5,000 molecules were collected using a 10-minute C18 HPLC method. The data was split into training, validation, and testing sets with a ratio of 8:1:1. Three molecular embedding approaches—fingerprint, descriptor, and graph-based methods—were evaluated. Among these, the graph-based method with an attention mechanism (AttentiveFP) demonstrated the best performance. A Bayesian optimization fine-tuning process was then employed to optimize the model parameters and achieve the best results. The final model was deployed in a Docker container, enabling both batch predictions and API-based single-molecule predictions.

### Jupyter Notebook Contents

1. **Step 1**: Training Data Input, Data Cleaning, and Quality Control  
   - Imported and cleaned raw data for training. (SMILES standardization, duplicated samples)
   - Conducted quality control to identify irreducible errors. (mean difference 0.196 min for duplicates)

2. **Step 2**: Baseline Model (Molecular Descriptors and XGBoost), PCA  
   - Built a baseline regression model using molecular descriptors, fingerprints and XGBoost.
   - Performed PCA for feature dimensionality reduction and visualization.
   - Mean Absolute Error (MAE) as loss function; robust to large deviations.

3. **Step 3**: Testing Different GNN Architectures  
   - Evaluated various GNN models including GCN, GAT, MPNN, and AttentiveFP.
   - Early stop
   - Best parameters from publications as default
   

4. **Step 4**: Fine-Tuning on Best-Performing Model  
   - Conducted grid search hyperparameter tuning on AttentiveFP, the best-performing model.
   - Achieved significant improvements in prediction accuracy.

5. **Scalability and Deployment**:  
   - [API](./GNN_api) and [docker](./GNN_docker)

### Summary
<img width="1430" alt="image" src="https://user-images.githubusercontent.com/30486093/186263575-c5e24ea6-3182-43e6-b0d0-ef3d106f2c41.png">