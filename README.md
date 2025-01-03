## ML Model for RT Prediction (2025)

I developed and trained this model from two years ago. Due to time constraints, I could not include:

1. Masking pretraining for enhanced feature learning.
2. Combining molecular descriptors with the GNN model for integrated predictions.

### Jupyter Notebook Contents

1. **Step 1**: Training Data Input, Data Cleaning, and Quality Control  
   - Imported and cleaned raw data for training.
   - Conducted quality control to identify and correct intrinsic errors.

2. **Step 2**: Baseline Model (Molecular Descriptors and XGBoost), PCA  
   - Built a baseline regression model using molecular descriptors and XGBoost.
   - Performed PCA for feature dimensionality reduction and visualization.

3. **Step 3**: Testing Different GNN Architectures  
   - Evaluated various GNN models including GCN, GAT, MPNN, and AttentiveFP.
   - Tested the integration of GNN and molecular descriptors.
    <img width="1430" alt="image" src="https://user-images.githubusercontent.com/30486093/186263575-c5e24ea6-3182-43e6-b0d0-ef3d106f2c41.png">
4. **Step 4**: Fine-Tuning on Best-Performing Model  
   - Conducted grid search hyperparameter tuning on AttentiveFP, the best-performing model.
   - Achieved significant improvements in prediction accuracy.

5. **Scalability and Deployment**:  
   - [API](./GNN_api) and [docker](./GNN_docker)