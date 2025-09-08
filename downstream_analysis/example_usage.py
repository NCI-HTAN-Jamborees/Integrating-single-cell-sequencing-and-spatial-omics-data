#!/usr/bin/env python3
"""
Simple example script demonstrating binary classification for high-dimensional single-cell data.

This script shows how to:
1. Create or load single-cell data
2. Perform binary classification
3. Evaluate results
4. Visualize performance

Run with: python example_usage.py
"""

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from binary_classification import HighDimBinaryClassifier, create_synthetic_classification_data


def create_example_data():
    """Create example single-cell dataset for demonstration."""
    print("Creating example single-cell dataset...")
    
    # Set random seed for reproducibility
    np.random.seed(42)
    
    # Dataset parameters
    n_cells = 1000
    n_genes = 2000
    
    # Generate synthetic expression data
    X = np.random.negative_binomial(5, 0.3, size=(n_cells, n_genes)).astype(float)
    
    # Add biological structure - two cell populations with different expression
    n_diff_genes = n_genes // 10  # 10% of genes are differentially expressed
    
    # Population A (first half) - higher expression of first gene set
    X[:n_cells//2, :n_diff_genes] *= 2.5
    
    # Population B (second half) - higher expression of second gene set
    X[n_cells//2:, n_diff_genes:2*n_diff_genes] *= 2.5
    
    # Create AnnData object
    adata = ad.AnnData(X)
    adata.var_names = [f'Gene_{i:04d}' for i in range(n_genes)]
    adata.obs_names = [f'Cell_{i:04d}' for i in range(n_cells)]
    
    # Add metadata
    adata.obs['cell_type'] = ['TypeA' if i < n_cells//2 else 'TypeB' for i in range(n_cells)]
    adata.obs['batch'] = np.random.choice(['Batch1', 'Batch2'], n_cells)
    
    # Add some spatial coordinates
    adata.obs['x'] = np.random.uniform(0, 100, n_cells)
    adata.obs['y'] = np.random.uniform(0, 100, n_cells)
    
    # Create malignant status based on cell type and some noise
    malignant_prob = np.where(adata.obs['cell_type'] == 'TypeA', 0.7, 0.3)
    adata.obs['malignant'] = np.random.binomial(1, malignant_prob, n_cells)
    adata.obs['malignant'] = adata.obs['malignant'].map({1: 'Malignant', 0: 'Benign'})
    
    print(f"Created dataset with {n_cells} cells and {n_genes} genes")
    print(f"Cell type distribution: {adata.obs['cell_type'].value_counts().to_dict()}")
    print(f"Malignant status: {adata.obs['malignant'].value_counts().to_dict()}")
    
    return adata


def run_classification_example():
    """Run a complete binary classification example."""
    print("=" * 60)
    print("Binary Classification Example for High-Dimensional Data")
    print("=" * 60)
    
    # Create example data
    adata = create_example_data()
    
    # Basic preprocessing
    print("\nPerforming basic preprocessing...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=500)
    
    # Use only highly variable genes for classification
    adata_hv = adata[:, adata.var.highly_variable].copy()
    print(f"Using {adata_hv.n_vars} highly variable genes for classification")
    
    # Example 1: Cell Type Classification
    print("\n" + "-" * 40)
    print("Example 1: Cell Type Classification")
    print("-" * 40)
    
    classifier1 = HighDimBinaryClassifier(random_state=42)
    
    # Prepare data
    X1, y1 = classifier1.prepare_data(
        adata_hv, 
        target_column='cell_type',
        positive_class='TypeA',
        negative_class='TypeB'
    )
    
    # Feature selection
    X1_selected = classifier1.select_features(X1, y1, method='selectk', n_features=100)
    
    # Train Random Forest model
    classifier1.train_model(X1_selected, y1, model_type='random_forest', n_estimators=100)
    
    # Evaluate
    metrics1 = classifier1.evaluate_model(X1_selected, y1, cv_folds=5)
    
    print(f"Best performing metric: ROC-AUC = {metrics1['roc_auc_mean']:.3f}")
    
    # Get feature importance
    top_features1 = classifier1.get_feature_importance(top_n=10)
    print(f"\nTop 5 features for cell type classification:")
    for i, row in top_features1.head(5).iterrows():
        print(f"  {row['feature']}: {row['importance']:.4f}")
    
    # Example 2: Malignant vs Benign Classification
    print("\n" + "-" * 40)
    print("Example 2: Malignant vs Benign Classification")
    print("-" * 40)
    
    classifier2 = HighDimBinaryClassifier(random_state=42)
    
    # Prepare data
    X2, y2 = classifier2.prepare_data(
        adata_hv, 
        target_column='malignant',
        positive_class='Malignant',
        negative_class='Benign'
    )
    
    # Try PCA for dimensionality reduction
    X2_selected = classifier2.select_features(X2, y2, method='pca', n_features=50)
    
    # Train Logistic Regression model
    classifier2.train_model(X2_selected, y2, model_type='logistic', C=1.0)
    
    # Evaluate
    metrics2 = classifier2.evaluate_model(X2_selected, y2, cv_folds=5)
    
    print(f"Best performing metric: ROC-AUC = {metrics2['roc_auc_mean']:.3f}")
    
    # Model Comparison
    print("\n" + "-" * 40)
    print("Model Comparison")
    print("-" * 40)
    
    models_to_test = [
        ('Random Forest', 'random_forest', {'n_estimators': 50}),
        ('Logistic Regression', 'logistic', {'C': 1.0}),
        ('SVM', 'svm', {'C': 1.0, 'kernel': 'linear'})
    ]
    
    print("Comparing models for malignant vs benign classification:")
    
    best_score = 0
    best_model = None
    
    for model_name, model_type, params in models_to_test:
        # Use a subset of data for faster comparison
        X_subset = X2_selected[:500]  # Use first 500 cells
        y_subset = y2[:500]
        
        clf = HighDimBinaryClassifier(random_state=42)
        clf.train_model(X_subset, y_subset, model_type=model_type, **params)
        metrics = clf.evaluate_model(X_subset, y_subset, cv_folds=3)
        
        roc_auc = metrics['roc_auc_mean']
        print(f"  {model_name}: ROC-AUC = {roc_auc:.3f} ± {metrics['roc_auc_std']:.3f}")
        
        if roc_auc > best_score:
            best_score = roc_auc
            best_model = model_name
    
    print(f"\nBest performing model: {best_model} (ROC-AUC = {best_score:.3f})")
    
    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print("✓ Successfully demonstrated binary classification on high-dimensional data")
    print("✓ Tested multiple feature selection methods (SelectKBest, PCA)")
    print("✓ Compared different classification algorithms")
    print("✓ Evaluated performance with cross-validation")
    print("✓ Analyzed feature importance")
    
    print(f"\nKey Results:")
    print(f"- Cell type classification accuracy: {metrics1['accuracy_mean']:.3f}")
    print(f"- Malignant classification accuracy: {metrics2['accuracy_mean']:.3f}")
    print(f"- Best overall model: {best_model}")
    
    print(f"\nNext steps:")
    print("- Apply to your actual integrated scRNA-seq/MERFISH data")
    print("- Tune hyperparameters for better performance") 
    print("- Validate results with biological knowledge")
    print("- Use spatial information for spatial transcriptomics")
    
    return classifier1, classifier2


if __name__ == "__main__":
    # Run the example
    try:
        classifier1, classifier2 = run_classification_example()
        print("\n🎉 Example completed successfully!")
        
    except Exception as e:
        print(f"\n❌ Error running example: {e}")
        print("Make sure all required packages are installed:")
        print("pip install scanpy pandas numpy scikit-learn matplotlib seaborn anndata")