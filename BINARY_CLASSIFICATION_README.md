# Binary Classification for High-Dimensional Single-Cell Data

This module provides a comprehensive solution for binary classification of high-dimensional single-cell genomics data. It's specifically designed to handle the challenges of single-cell RNA sequencing (scRNA-seq) and spatial transcriptomics data where the number of features (genes) often exceeds the number of samples (cells).

## Features

- **High-Dimensional Data Support**: Efficiently handles datasets with thousands of features
- **Dimensionality Reduction**: Automatic PCA for noise reduction and computational efficiency
- **Feature Selection**: Statistical feature selection to identify the most informative genes
- **Multiple Classifiers**: Support for Logistic Regression, Random Forest, and SVM
- **Cross-Validation**: Built-in cross-validation for robust performance estimation
- **AnnData Integration**: Direct support for AnnData objects (standard in single-cell analysis)
- **Comprehensive Evaluation**: Multiple metrics including ROC AUC, accuracy, and confusion matrices
- **Visualization**: Built-in plotting for performance assessment

## Quick Start

### Basic Usage

```python
from binary_classification import HighDimBinaryClassifier
from sklearn.datasets import make_classification
from sklearn.model_selection import train_test_split

# Generate high-dimensional data
X, y = make_classification(n_samples=1000, n_features=2000, n_informative=100)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)

# Initialize classifier
classifier = HighDimBinaryClassifier(
    classifier_type='random_forest',
    n_components=50,        # PCA components
    feature_selection_k=100 # Top features to select
)

# Fit and evaluate
classifier.fit(X_train, y_train)
results = classifier.evaluate(X_test, y_test)

print(f"ROC AUC: {results['roc_auc']:.3f}")
print(f"Accuracy: {results['accuracy']:.3f}")
```

### Working with AnnData Objects

```python
from binary_classification import classify_anndata
import anndata as ad

# Load your single-cell data
adata = ad.read_h5ad('your_data.h5ad')

# Classify tumor vs normal cells
classifier, results = classify_anndata(
    adata=adata,
    target_column='cell_type',
    positive_class='Tumor',
    classifier_type='random_forest'
)

print(f"Classification accuracy: {results['accuracy']:.3f}")
```

## Class Reference

### HighDimBinaryClassifier

The main classifier class for high-dimensional binary classification.

#### Parameters

- `classifier_type` (str): Type of classifier ('logistic', 'random_forest', 'svm')
- `n_components` (int, optional): Number of PCA components. Auto-determined if None
- `feature_selection_k` (int, optional): Number of top features to select
- `random_state` (int): Random state for reproducibility

#### Methods

- `fit(X, y)`: Fit the classifier to training data
- `predict(X)`: Make predictions on new data
- `predict_proba(X)`: Get prediction probabilities
- `evaluate(X_test, y_test)`: Comprehensive evaluation with multiple metrics
- `plot_performance(X_test, y_test)`: Visualize classifier performance
- `get_feature_importance()`: Get feature importance scores

## Use Cases

### 1. Cell Type Classification

Classify cells as tumor vs. normal, immune vs. non-immune, etc.

```python
# Example: Tumor vs Normal classification
classifier = HighDimBinaryClassifier(
    classifier_type='random_forest',
    n_components=100,
    feature_selection_k=200
)
```

### 2. Treatment Response Prediction

Identify cells that respond to specific treatments.

```python
# Example: Drug response classification
classifier = HighDimBinaryClassifier(
    classifier_type='logistic',
    n_components=50,
    feature_selection_k=100
)
```

### 3. Spatial Classification

Classify cells based on spatial location (e.g., tumor core vs. periphery).

```python
# Example: Spatial region classification
classifier = HighDimBinaryClassifier(
    classifier_type='svm',
    n_components=75,
    feature_selection_k=150
)
```

## Performance Considerations

### Dimensionality Reduction

For high-dimensional data (>1000 features), PCA is recommended:

- **Small datasets** (<500 cells): Use 20-50 components
- **Medium datasets** (500-2000 cells): Use 50-100 components  
- **Large datasets** (>2000 cells): Use 100-200 components

### Feature Selection

Feature selection can improve performance and interpretability:

- Start with top 10-20% of PCA components
- For gene expression: 100-500 features often work well
- Use cross-validation to optimize the number

### Classifier Choice

- **Logistic Regression**: Fast, interpretable, good baseline
- **Random Forest**: Handles non-linear relationships, provides feature importance
- **SVM**: Good for complex decision boundaries, slower on large datasets

## Examples

See `binary_classification_demo.ipynb` for comprehensive examples including:

- Synthetic high-dimensional data classification
- Simulated single-cell tumor vs. normal classification
- Comparison of different classifiers
- Feature importance analysis
- Working with AnnData objects

## Testing

Run the test suite to verify functionality:

```bash
python test_binary_classification.py
```

The tests cover:
- Basic functionality
- Different classifier types
- Feature selection
- Automatic PCA component determination
- Error handling
- AnnData integration (if available)

## Integration with Existing Workflows

This module is designed to integrate seamlessly with existing single-cell analysis pipelines:

- **Scanpy**: Works with AnnData objects from scanpy workflows
- **Seurat**: Can import data from Seurat via anndata conversion
- **Cell type annotation**: Can be used downstream of clustering/annotation
- **Spatial analysis**: Compatible with spatial transcriptomics workflows

## Best Practices

1. **Data Preprocessing**: Ensure proper normalization before classification
2. **Cross-Validation**: Always use cross-validation for performance estimation
3. **Feature Selection**: Use feature selection for interpretability
4. **Class Balance**: Consider class imbalance in your data
5. **Validation**: Test on independent datasets when possible

## Requirements

- numpy
- pandas  
- scikit-learn
- matplotlib
- seaborn
- anndata (optional, for AnnData support)

## Contributing

This module is part of the HTAN Data Jamboree project for integrating single-cell and spatial omics data. Contributions and improvements are welcome.

## License

This project is licensed under the MIT License - see the LICENSE file for details.