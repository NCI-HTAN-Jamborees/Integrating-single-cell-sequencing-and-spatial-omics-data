"""
Binary Classification for High-Dimensional Single-Cell and Spatial Transcriptomics Data

This module provides tools for binary classification tasks on high-dimensional
single-cell RNA sequencing and spatial transcriptomics data, including:
- Cell type classification (e.g., malignant vs non-malignant)
- Spatial region classification (e.g., tumor vs normal regions)
- Integration quality assessment
- Feature selection and dimensionality reduction for high-dimensional data
"""

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from sklearn.model_selection import train_test_split, cross_val_score, StratifiedKFold
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.feature_selection import SelectKBest, f_classif, RFE
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.metrics import (
    accuracy_score, precision_score, recall_score, f1_score,
    roc_auc_score, confusion_matrix, classification_report,
    roc_curve, precision_recall_curve
)
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Tuple, Dict, List, Optional, Union
import warnings
warnings.filterwarnings('ignore')


class HighDimBinaryClassifier:
    """
    Binary classifier for high-dimensional single-cell and spatial transcriptomics data.
    
    This class provides methods for:
    - Feature selection and dimensionality reduction
    - Model training with multiple algorithms
    - Cross-validation and evaluation
    - Visualization of results
    """
    
    def __init__(self, random_state: int = 42):
        """
        Initialize the classifier.
        
        Parameters:
        -----------
        random_state : int
            Random state for reproducibility
        """
        self.random_state = random_state
        self.scaler = StandardScaler()
        self.feature_selector = None
        self.dim_reducer = None
        self.model = None
        self.label_encoder = LabelEncoder()
        self.feature_names = None
        self.selected_features = None
        
    def prepare_data(self, 
                    adata: ad.AnnData, 
                    target_column: str,
                    positive_class: str,
                    negative_class: str,
                    feature_type: str = 'X',
                    min_cells: int = 10) -> Tuple[np.ndarray, np.ndarray]:
        """
        Prepare data for binary classification.
        
        Parameters:
        -----------
        adata : AnnData
            Annotated data object containing expression data and metadata
        target_column : str
            Column name in adata.obs containing target labels
        positive_class : str
            Label for positive class
        negative_class : str
            Label for negative class
        feature_type : str
            Type of features to use ('X', 'X_pca', 'obsm_key')
        min_cells : int
            Minimum number of cells per class
            
        Returns:
        --------
        X : np.ndarray
            Feature matrix
        y : np.ndarray
            Binary target labels (0/1)
        """
        print(f"Preparing binary classification data...")
        print(f"Target column: {target_column}")
        print(f"Positive class: {positive_class}, Negative class: {negative_class}")
        
        # Filter data to include only specified classes
        mask = adata.obs[target_column].isin([positive_class, negative_class])
        adata_filtered = adata[mask].copy()
        
        # Check minimum cell counts
        class_counts = adata_filtered.obs[target_column].value_counts()
        print(f"Class distribution: {dict(class_counts)}")
        
        if any(class_counts < min_cells):
            raise ValueError(f"Insufficient data: some classes have < {min_cells} cells")
        
        # Extract features
        if feature_type == 'X':
            X = adata_filtered.X
            if hasattr(X, 'toarray'):  # Handle sparse matrices
                X = X.toarray()
            self.feature_names = adata_filtered.var_names.tolist()
        elif feature_type == 'X_pca':
            if 'X_pca' not in adata_filtered.obsm:
                sc.pp.pca(adata_filtered, n_comps=50)
            X = adata_filtered.obsm['X_pca']
            self.feature_names = [f'PC{i+1}' for i in range(X.shape[1])]
        else:
            raise ValueError(f"Unsupported feature_type: {feature_type}")
        
        # Prepare binary labels
        y_labels = adata_filtered.obs[target_column].values
        y = (y_labels == positive_class).astype(int)
        
        print(f"Data shape: {X.shape}")
        print(f"Binary label distribution: {np.bincount(y)}")
        
        return X, y
    
    def select_features(self, 
                       X: np.ndarray, 
                       y: np.ndarray, 
                       method: str = 'selectk',
                       n_features: int = 1000) -> np.ndarray:
        """
        Select most informative features for classification.
        
        Parameters:
        -----------
        X : np.ndarray
            Feature matrix
        y : np.ndarray
            Target labels
        method : str
            Feature selection method ('selectk', 'rfe', 'pca')
        n_features : int
            Number of features to select
            
        Returns:
        --------
        X_selected : np.ndarray
            Feature matrix with selected features
        """
        print(f"Selecting features using {method}...")
        
        if method == 'selectk':
            self.feature_selector = SelectKBest(score_func=f_classif, k=min(n_features, X.shape[1]))
            X_selected = self.feature_selector.fit_transform(X, y)
            
            # Get selected feature names
            if self.feature_names:
                selected_idx = self.feature_selector.get_support()
                self.selected_features = [self.feature_names[i] for i in range(len(selected_idx)) if selected_idx[i]]
                
        elif method == 'rfe':
            base_estimator = RandomForestClassifier(n_estimators=100, random_state=self.random_state)
            self.feature_selector = RFE(base_estimator, n_features_to_select=min(n_features, X.shape[1]))
            X_selected = self.feature_selector.fit_transform(X, y)
            
            # Get selected feature names
            if self.feature_names:
                selected_idx = self.feature_selector.get_support()
                self.selected_features = [self.feature_names[i] for i in range(len(selected_idx)) if selected_idx[i]]
                
        elif method == 'pca':
            self.dim_reducer = PCA(n_components=min(n_features, X.shape[1]), random_state=self.random_state)
            X_selected = self.dim_reducer.fit_transform(X)
            self.selected_features = [f'PC{i+1}' for i in range(X_selected.shape[1])]
            
        else:
            raise ValueError(f"Unsupported feature selection method: {method}")
        
        print(f"Selected {X_selected.shape[1]} features from {X.shape[1]}")
        return X_selected
    
    def train_model(self, 
                   X: np.ndarray, 
                   y: np.ndarray, 
                   model_type: str = 'random_forest',
                   scale_features: bool = True,
                   **model_params) -> None:
        """
        Train binary classification model.
        
        Parameters:
        -----------
        X : np.ndarray
            Feature matrix
        y : np.ndarray
            Target labels
        model_type : str
            Type of model ('random_forest', 'logistic', 'svm')
        scale_features : bool
            Whether to scale features
        **model_params
            Additional parameters for the model
        """
        print(f"Training {model_type} model...")
        
        # Scale features if requested
        if scale_features:
            X = self.scaler.fit_transform(X)
        
        # Initialize model
        if model_type == 'random_forest':
            default_params = {'n_estimators': 100, 'random_state': self.random_state}
            default_params.update(model_params)
            self.model = RandomForestClassifier(**default_params)
            
        elif model_type == 'logistic':
            default_params = {'random_state': self.random_state, 'max_iter': 1000}
            default_params.update(model_params)
            self.model = LogisticRegression(**default_params)
            
        elif model_type == 'svm':
            default_params = {'random_state': self.random_state, 'probability': True}
            default_params.update(model_params)
            self.model = SVC(**default_params)
            
        else:
            raise ValueError(f"Unsupported model type: {model_type}")
        
        # Train model
        self.model.fit(X, y)
        print(f"Model training completed.")
    
    def evaluate_model(self, 
                      X: np.ndarray, 
                      y: np.ndarray, 
                      cv_folds: int = 5) -> Dict[str, float]:
        """
        Evaluate model performance using cross-validation.
        
        Parameters:
        -----------
        X : np.ndarray
            Feature matrix
        y : np.ndarray
            Target labels
        cv_folds : int
            Number of cross-validation folds
            
        Returns:
        --------
        metrics : dict
            Dictionary containing evaluation metrics
        """
        print(f"Evaluating model with {cv_folds}-fold cross-validation...")
        
        # Scale features if scaler was fitted
        if hasattr(self.scaler, 'mean_'):
            X_scaled = self.scaler.transform(X)
        else:
            X_scaled = X
        
        # Cross-validation
        cv = StratifiedKFold(n_splits=cv_folds, shuffle=True, random_state=self.random_state)
        
        cv_scores = {
            'accuracy': cross_val_score(self.model, X_scaled, y, cv=cv, scoring='accuracy'),
            'precision': cross_val_score(self.model, X_scaled, y, cv=cv, scoring='precision'),
            'recall': cross_val_score(self.model, X_scaled, y, cv=cv, scoring='recall'),
            'f1': cross_val_score(self.model, X_scaled, y, cv=cv, scoring='f1'),
            'roc_auc': cross_val_score(self.model, X_scaled, y, cv=cv, scoring='roc_auc')
        }
        
        # Calculate mean and std for each metric
        metrics = {}
        for metric, scores in cv_scores.items():
            metrics[f'{metric}_mean'] = np.mean(scores)
            metrics[f'{metric}_std'] = np.std(scores)
        
        # Print results
        print("\nCross-validation results:")
        for metric in ['accuracy', 'precision', 'recall', 'f1', 'roc_auc']:
            mean_score = metrics[f'{metric}_mean']
            std_score = metrics[f'{metric}_std']
            print(f"{metric.capitalize()}: {mean_score:.3f} ± {std_score:.3f}")
        
        return metrics
    
    def get_feature_importance(self, top_n: int = 20) -> pd.DataFrame:
        """
        Get feature importance scores.
        
        Parameters:
        -----------
        top_n : int
            Number of top features to return
            
        Returns:
        --------
        importance_df : pd.DataFrame
            DataFrame with feature names and importance scores
        """
        if self.model is None:
            raise ValueError("Model must be trained first")
        
        if hasattr(self.model, 'feature_importances_'):
            importances = self.model.feature_importances_
        elif hasattr(self.model, 'coef_'):
            importances = np.abs(self.model.coef_[0])
        else:
            raise ValueError("Model does not provide feature importance")
        
        # Create DataFrame
        feature_names = self.selected_features or [f'Feature_{i}' for i in range(len(importances))]
        importance_df = pd.DataFrame({
            'feature': feature_names,
            'importance': importances
        }).sort_values('importance', ascending=False).head(top_n)
        
        return importance_df
    
    def plot_results(self, 
                    X: np.ndarray, 
                    y: np.ndarray, 
                    show_feature_importance: bool = True,
                    show_roc_curve: bool = True,
                    show_confusion_matrix: bool = True,
                    figsize: Tuple[int, int] = (15, 5)) -> None:
        """
        Plot classification results and diagnostics.
        
        Parameters:
        -----------
        X : np.ndarray
            Feature matrix
        y : np.ndarray
            True labels
        show_feature_importance : bool
            Whether to show feature importance plot
        show_roc_curve : bool
            Whether to show ROC curve
        show_confusion_matrix : bool
            Whether to show confusion matrix
        figsize : tuple
            Figure size
        """
        if self.model is None:
            raise ValueError("Model must be trained first")
        
        # Scale features if scaler was fitted
        if hasattr(self.scaler, 'mean_'):
            X_scaled = self.scaler.transform(X)
        else:
            X_scaled = X
        
        # Get predictions
        y_pred = self.model.predict(X_scaled)
        y_prob = self.model.predict_proba(X_scaled)[:, 1] if hasattr(self.model, 'predict_proba') else None
        
        # Calculate number of subplots
        n_plots = sum([show_feature_importance, show_roc_curve, show_confusion_matrix])
        
        if n_plots == 0:
            return
        
        fig, axes = plt.subplots(1, n_plots, figsize=(figsize[0], figsize[1]))
        if n_plots == 1:
            axes = [axes]
        
        plot_idx = 0
        
        # Feature importance plot
        if show_feature_importance:
            try:
                importance_df = self.get_feature_importance(top_n=15)
                ax = axes[plot_idx]
                sns.barplot(data=importance_df, x='importance', y='feature', ax=ax)
                ax.set_title('Top Feature Importances')
                ax.set_xlabel('Importance Score')
                plot_idx += 1
            except Exception as e:
                print(f"Could not plot feature importance: {e}")
                plot_idx += 1
        
        # ROC curve
        if show_roc_curve and y_prob is not None:
            ax = axes[plot_idx]
            fpr, tpr, _ = roc_curve(y, y_prob)
            auc_score = roc_auc_score(y, y_prob)
            ax.plot(fpr, tpr, label=f'ROC Curve (AUC = {auc_score:.3f})')
            ax.plot([0, 1], [0, 1], 'k--', alpha=0.5)
            ax.set_xlabel('False Positive Rate')
            ax.set_ylabel('True Positive Rate')
            ax.set_title('ROC Curve')
            ax.legend()
            plot_idx += 1
        
        # Confusion matrix
        if show_confusion_matrix:
            ax = axes[plot_idx]
            cm = confusion_matrix(y, y_pred)
            sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', ax=ax)
            ax.set_xlabel('Predicted')
            ax.set_ylabel('Actual')
            ax.set_title('Confusion Matrix')
        
        plt.tight_layout()
        plt.show()


def create_synthetic_classification_data(adata: ad.AnnData, 
                                        target_gene_sets: Dict[str, List[str]] = None) -> ad.AnnData:
    """
    Create synthetic binary classification targets based on gene expression patterns.
    
    This function is useful when ground truth labels are not available.
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data object
    target_gene_sets : dict
        Dictionary mapping class names to lists of marker genes
        
    Returns:
    --------
    adata_with_targets : AnnData
        AnnData object with synthetic classification targets
    """
    adata_syn = adata.copy()
    
    if target_gene_sets is None:
        # Create default gene sets based on high variance genes
        sc.pp.highly_variable_genes(adata_syn, n_top_genes=2000)
        highly_var_genes = adata_syn.var_names[adata_syn.var.highly_variable].tolist()
        
        # Split into two sets
        mid_point = len(highly_var_genes) // 2
        target_gene_sets = {
            'group_A': highly_var_genes[:mid_point],
            'group_B': highly_var_genes[mid_point:]
        }
    
    # Calculate mean expression for each gene set
    for group_name, genes in target_gene_sets.items():
        available_genes = [g for g in genes if g in adata_syn.var_names]
        if available_genes:
            gene_indices = [adata_syn.var_names.get_loc(g) for g in available_genes]
            if hasattr(adata_syn.X, 'toarray'):
                expression_matrix = adata_syn.X.toarray()
            else:
                expression_matrix = adata_syn.X
            mean_expression = np.mean(expression_matrix[:, gene_indices], axis=1)
            adata_syn.obs[f'{group_name}_score'] = mean_expression
    
    # Create binary classification target
    if len(target_gene_sets) == 2:
        group_names = list(target_gene_sets.keys())
        score_diff = adata_syn.obs[f'{group_names[0]}_score'] - adata_syn.obs[f'{group_names[1]}_score']
        adata_syn.obs['binary_class'] = (score_diff > np.median(score_diff)).astype(str)
        adata_syn.obs['binary_class'] = adata_syn.obs['binary_class'].map({'True': group_names[0], 'False': group_names[1]})
    
    return adata_syn


# Example usage and testing functions
def run_classification_example(adata: ad.AnnData, 
                             target_column: str = None,
                             positive_class: str = None,
                             negative_class: str = None) -> HighDimBinaryClassifier:
    """
    Run a complete binary classification example.
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data object
    target_column : str
        Column name containing target labels (if None, creates synthetic data)
    positive_class : str
        Positive class label
    negative_class : str
        Negative class label
        
    Returns:
    --------
    classifier : HighDimBinaryClassifier
        Trained classifier object
    """
    # Initialize classifier
    classifier = HighDimBinaryClassifier(random_state=42)
    
    # Create synthetic data if no target provided
    if target_column is None:
        print("No target column provided, creating synthetic classification data...")
        adata = create_synthetic_classification_data(adata)
        target_column = 'binary_class'
        unique_classes = adata.obs[target_column].unique()
        positive_class = unique_classes[0]
        negative_class = unique_classes[1]
    
    # Prepare data
    X, y = classifier.prepare_data(adata, target_column, positive_class, negative_class)
    
    # Feature selection
    X_selected = classifier.select_features(X, y, method='selectk', n_features=500)
    
    # Train model
    classifier.train_model(X_selected, y, model_type='random_forest')
    
    # Evaluate model
    metrics = classifier.evaluate_model(X_selected, y)
    
    # Plot results
    classifier.plot_results(X_selected, y)
    
    return classifier


if __name__ == "__main__":
    print("Binary Classification Module for High-Dimensional Single-Cell Data")
    print("This module provides tools for binary classification on scRNA-seq and spatial transcriptomics data.")
    print("Use run_classification_example() to test with your data.")