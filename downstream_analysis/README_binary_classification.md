# Binary Classification for High-Dimensional Single-Cell Data

This module provides comprehensive tools for binary classification tasks on high-dimensional single-cell RNA sequencing and spatial transcriptomics data.

## Overview

The module includes:
- **Feature selection methods** for high-dimensional data (SelectKBest, RFE, PCA)
- **Multiple classification algorithms** (Random Forest, Logistic Regression, SVM)
- **Robust evaluation** with cross-validation and multiple metrics
- **Visualization tools** for results interpretation
- **Feature importance analysis** for biological insight

## Files

- `binary_classification.py` - Main module with the `HighDimBinaryClassifier` class
- `binary_classification_demo.ipynb` - Comprehensive Jupyter notebook demonstration
- `example_usage.py` - Simple example script
- `README_binary_classification.md` - This documentation

## Quick Start

```python
from binary_classification import HighDimBinaryClassifier
import scanpy as sc

# Load your data
adata = sc.read_h5ad('your_data.h5ad')

# Initialize classifier
classifier = HighDimBinaryClassifier(random_state=42)

# Prepare data for binary classification
X, y = classifier.prepare_data(
    adata, 
    target_column='cell_type',  # Column with your labels
    positive_class='Malignant', # Positive class
    negative_class='Benign'     # Negative class
)

# Select important features
X_selected = classifier.select_features(X, y, method='selectk', n_features=500)

# Train model
classifier.train_model(X_selected, y, model_type='random_forest')

# Evaluate performance
metrics = classifier.evaluate_model(X_selected, y, cv_folds=5)

# Visualize results
classifier.plot_results(X_selected, y)
```

## Classification Scenarios

### 1. Cell Type Classification
Distinguish between different cell types (e.g., malignant vs benign cells):
```python
X, y = classifier.prepare_data(adata, 'malignant_status', 'Malignant', 'Benign')
```

### 2. Spatial Region Classification
Classify spatial regions (e.g., tumor vs normal tissue):
```python
X, y = classifier.prepare_data(adata, 'tissue_region', 'Tumor', 'Normal')
```

### 3. Treatment Response Prediction
Predict treatment response:
```python
X, y = classifier.prepare_data(adata, 'response_status', 'Responder', 'Non-responder')
```

## Feature Selection Methods

### SelectKBest
Best for: Large feature sets, interpretability
```python
X_selected = classifier.select_features(X, y, method='selectk', n_features=1000)
```

### Recursive Feature Elimination (RFE)
Best for: Model-specific feature selection
```python
X_selected = classifier.select_features(X, y, method='rfe', n_features=500)
```

### Principal Component Analysis (PCA)
Best for: Dimensionality reduction, correlated features
```python
X_selected = classifier.select_features(X, y, method='pca', n_features=100)
```

## Classification Algorithms

### Random Forest
- Good default choice
- Handles high-dimensional data well
- Provides feature importance
```python
classifier.train_model(X, y, model_type='random_forest', n_estimators=200)
```

### Logistic Regression
- Fast and interpretable
- Good for linearly separable data
- Provides probability estimates
```python
classifier.train_model(X, y, model_type='logistic', C=1.0)
```

### Support Vector Machine (SVM)
- Good for complex decision boundaries
- Effective in high dimensions
```python
classifier.train_model(X, y, model_type='svm', C=1.0, kernel='rbf')
```

## Evaluation Metrics

The module provides comprehensive evaluation:
- **Accuracy**: Overall correct predictions
- **Precision**: Positive predictive value
- **Recall**: Sensitivity/True positive rate
- **F1-score**: Harmonic mean of precision and recall
- **ROC-AUC**: Area under the receiver operating characteristic curve

## Best Practices

### Data Preparation
1. **Quality Control**: Remove low-quality cells and genes
2. **Normalization**: Use appropriate normalization (e.g., log1p transformation)
3. **Feature Selection**: Use highly variable genes or domain knowledge

### Model Training
1. **Cross-Validation**: Always use cross-validation for robust evaluation
2. **Feature Scaling**: Scale features for distance-based algorithms
3. **Class Balance**: Check for class imbalance and handle appropriately

### Interpretation
1. **Feature Importance**: Examine which genes/features are most informative
2. **Biological Validation**: Validate findings with known biology
3. **Independent Testing**: Test on independent datasets when possible

## Integration with Repository Data

This module is designed to work with the integrated scRNA-seq and MERFISH data in this repository:

```python
# Example with repository data structure
import scanpy as sc
from binary_classification import HighDimBinaryClassifier

# Load integrated data (adjust path as needed)
# adata = sc.read_h5ad('path_to_integrated_data.h5ad')

# Typical use cases for this repository:
# 1. Classify cell types between modalities
# 2. Identify spatial tumor regions
# 3. Assess integration quality
# 4. Predict biological states
```

## Troubleshooting

### Common Issues

**Error: "Insufficient data: some classes have < 10 cells"**
- Solution: Increase `min_cells` parameter or combine rare classes

**Poor performance (accuracy < 0.6)**
- Try different feature selection methods
- Increase number of selected features
- Check data quality and preprocessing
- Ensure classes are actually distinguishable

**Memory issues with large datasets**
- Reduce number of features first
- Use PCA for dimensionality reduction
- Process data in batches

### Tips for Better Performance

1. **Use highly variable genes** for feature selection
2. **Try different algorithms** - some may work better for your data
3. **Tune hyperparameters** using grid search or random search
4. **Combine multiple models** with ensemble methods
5. **Include spatial information** for spatial transcriptomics data

## Example Applications

### Cancer Biology
- Malignant vs benign cell classification
- Drug resistance prediction
- Metastatic potential assessment

### Developmental Biology
- Cell fate prediction
- Differentiation state classification
- Temporal progression analysis

### Spatial Biology
- Tissue region identification
- Cell-cell interaction prediction
- Spatial pattern recognition

## Citation

If you use this module in your research, please cite the original repository and relevant methods papers.

## License

This module is part of the HTAN Data Jamboree project and follows the same license terms.