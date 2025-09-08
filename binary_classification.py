"""
Binary Classification Module for High-Dimensional Single-Cell Data

This module provides tools for binary classification of high-dimensional single-cell 
genomics data, including dimensionality reduction, feature selection, and cross-validation.
Designed to work with AnnData objects commonly used in single-cell analysis.
"""

import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split, cross_val_score, StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.feature_selection import SelectKBest, f_classif, mutual_info_classif
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.metrics import classification_report, confusion_matrix, roc_auc_score, roc_curve
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Union, Tuple, Optional, Dict, Any
import warnings


class HighDimBinaryClassifier:
    """
    Binary classifier optimized for high-dimensional single-cell genomics data.
    
    This class provides a complete pipeline for binary classification including:
    - Dimensionality reduction (PCA)
    - Feature selection
    - Multiple classifier options
    - Cross-validation
    - Performance evaluation
    """
    
    def __init__(self, 
                 classifier_type: str = 'logistic',
                 n_components: Optional[int] = None,
                 feature_selection_k: Optional[int] = None,
                 random_state: int = 42):
        """
        Initialize the binary classifier.
        
        Parameters:
        -----------
        classifier_type : str, default='logistic'
            Type of classifier to use. Options: 'logistic', 'random_forest', 'svm'
        n_components : int, optional
            Number of PCA components. If None, will be determined automatically
        feature_selection_k : int, optional
            Number of top features to select. If None, uses all features
        random_state : int, default=42
            Random state for reproducibility
        """
        self.classifier_type = classifier_type
        self.n_components = n_components
        self.feature_selection_k = feature_selection_k
        self.random_state = random_state
        
        # Initialize components
        self.scaler = StandardScaler()
        self.pca = None
        self.feature_selector = None
        self.classifier = None
        
        # Results storage
        self.cv_scores = None
        self.feature_importance = None
        self.is_fitted = False
        
    def _initialize_classifier(self):
        """Initialize the classifier based on the specified type."""
        if self.classifier_type == 'logistic':
            self.classifier = LogisticRegression(
                random_state=self.random_state,
                max_iter=1000,
                class_weight='balanced'
            )
        elif self.classifier_type == 'random_forest':
            self.classifier = RandomForestClassifier(
                n_estimators=100,
                random_state=self.random_state,
                class_weight='balanced'
            )
        elif self.classifier_type == 'svm':
            self.classifier = SVC(
                random_state=self.random_state,
                class_weight='balanced',
                probability=True
            )
        else:
            raise ValueError(f"Unsupported classifier type: {self.classifier_type}")
    
    def _determine_pca_components(self, X: np.ndarray, variance_threshold: float = 0.95) -> int:
        """
        Determine optimal number of PCA components based on explained variance.
        
        Parameters:
        -----------
        X : array-like
            Input features
        variance_threshold : float, default=0.95
            Cumulative variance threshold for PCA components
            
        Returns:
        --------
        int : Optimal number of components
        """
        temp_pca = PCA()
        temp_pca.fit(X)
        
        cumsum_variance = np.cumsum(temp_pca.explained_variance_ratio_)
        n_components = np.argmax(cumsum_variance >= variance_threshold) + 1
        
        # Ensure we don't exceed the maximum possible components
        max_components = min(X.shape[0] - 1, X.shape[1])
        n_components = min(n_components, max_components)
        
        return max(1, n_components)  # Ensure at least 1 component
    
    def preprocess(self, X: np.ndarray, y: np.ndarray = None) -> np.ndarray:
        """
        Preprocess the data with scaling, PCA, and feature selection.
        
        Parameters:
        -----------
        X : array-like
            Input features
        y : array-like, optional
            Target labels (required for feature selection)
            
        Returns:
        --------
        array-like : Preprocessed features
        """
        # Scale the data
        if not hasattr(self.scaler, 'scale_'):
            X_scaled = self.scaler.fit_transform(X)
        else:
            X_scaled = self.scaler.transform(X)
        
        # Apply PCA if specified or if PCA was already fitted
        if self.n_components is not None or self.pca is not None or y is not None:
            if self.pca is None:
                n_comp = self.n_components
                if n_comp is None:
                    n_comp = self._determine_pca_components(X_scaled)
                
                self.pca = PCA(n_components=n_comp, random_state=self.random_state)
                X_pca = self.pca.fit_transform(X_scaled)
            else:
                X_pca = self.pca.transform(X_scaled)
            
            X_processed = X_pca
        else:
            X_processed = X_scaled
        
        # Apply feature selection if specified and during training
        if self.feature_selection_k is not None and y is not None:
            if self.feature_selector is None:
                self.feature_selector = SelectKBest(
                    score_func=f_classif, 
                    k=min(self.feature_selection_k, X_processed.shape[1])
                )
                X_processed = self.feature_selector.fit_transform(X_processed, y)
            else:
                # During prediction, only transform if feature selector exists
                X_processed = self.feature_selector.transform(X_processed)
        elif self.feature_selection_k is not None and self.feature_selector is not None:
            # During prediction phase when feature selector was fitted during training
            X_processed = self.feature_selector.transform(X_processed)
        
        return X_processed
    
    def fit(self, X: np.ndarray, y: np.ndarray, cv_folds: int = 5) -> 'HighDimBinaryClassifier':
        """
        Fit the binary classifier.
        
        Parameters:
        -----------
        X : array-like
            Input features
        y : array-like
            Binary target labels
        cv_folds : int, default=5
            Number of cross-validation folds
            
        Returns:
        --------
        self : Returns the fitted classifier
        """
        # Validate inputs
        X = np.array(X)
        y = np.array(y)
        
        if len(np.unique(y)) != 2:
            raise ValueError("This classifier is designed for binary classification only")
        
        # Preprocess the data
        X_processed = self.preprocess(X, y)
        
        # Initialize and fit classifier
        self._initialize_classifier()
        self.classifier.fit(X_processed, y)
        
        # Perform cross-validation
        cv = StratifiedKFold(n_splits=cv_folds, shuffle=True, random_state=self.random_state)
        self.cv_scores = cross_val_score(self.classifier, X_processed, y, cv=cv, scoring='roc_auc')
        
        # Store feature importance if available
        if hasattr(self.classifier, 'feature_importances_'):
            self.feature_importance = self.classifier.feature_importances_
        elif hasattr(self.classifier, 'coef_'):
            self.feature_importance = np.abs(self.classifier.coef_[0])
        
        self.is_fitted = True
        return self
    
    def predict(self, X: np.ndarray) -> np.ndarray:
        """
        Make predictions on new data.
        
        Parameters:
        -----------
        X : array-like
            Input features
            
        Returns:
        --------
        array-like : Predicted class labels
        """
        if not self.is_fitted:
            raise ValueError("Classifier must be fitted before making predictions")
        
        X_processed = self.preprocess(X)
        return self.classifier.predict(X_processed)
    
    def predict_proba(self, X: np.ndarray) -> np.ndarray:
        """
        Predict class probabilities.
        
        Parameters:
        -----------
        X : array-like
            Input features
            
        Returns:
        --------
        array-like : Predicted probabilities for each class
        """
        if not self.is_fitted:
            raise ValueError("Classifier must be fitted before making predictions")
        
        X_processed = self.preprocess(X)
        return self.classifier.predict_proba(X_processed)
    
    def evaluate(self, X_test: np.ndarray, y_test: np.ndarray) -> Dict[str, Any]:
        """
        Evaluate the classifier performance.
        
        Parameters:
        -----------
        X_test : array-like
            Test features
        y_test : array-like
            Test labels
            
        Returns:
        --------
        dict : Dictionary containing evaluation metrics
        """
        y_pred = self.predict(X_test)
        y_proba = self.predict_proba(X_test)[:, 1]
        
        results = {
            'accuracy': np.mean(y_pred == y_test),
            'roc_auc': roc_auc_score(y_test, y_proba),
            'classification_report': classification_report(y_test, y_pred),
            'confusion_matrix': confusion_matrix(y_test, y_pred),
            'cv_mean_score': np.mean(self.cv_scores) if self.cv_scores is not None else None,
            'cv_std_score': np.std(self.cv_scores) if self.cv_scores is not None else None
        }
        
        return results
    
    def plot_performance(self, X_test: np.ndarray, y_test: np.ndarray, figsize: Tuple[int, int] = (12, 4)):
        """
        Plot performance metrics including ROC curve and confusion matrix.
        
        Parameters:
        -----------
        X_test : array-like
            Test features
        y_test : array-like
            Test labels
        figsize : tuple, default=(12, 4)
            Figure size for the plots
        """
        y_pred = self.predict(X_test)
        y_proba = self.predict_proba(X_test)[:, 1]
        
        fig, axes = plt.subplots(1, 3, figsize=figsize)
        
        # ROC Curve
        fpr, tpr, _ = roc_curve(y_test, y_proba)
        auc = roc_auc_score(y_test, y_proba)
        
        axes[0].plot(fpr, tpr, label=f'ROC Curve (AUC = {auc:.3f})')
        axes[0].plot([0, 1], [0, 1], 'k--', label='Random')
        axes[0].set_xlabel('False Positive Rate')
        axes[0].set_ylabel('True Positive Rate')
        axes[0].set_title('ROC Curve')
        axes[0].legend()
        axes[0].grid(True)
        
        # Confusion Matrix
        cm = confusion_matrix(y_test, y_pred)
        sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', ax=axes[1])
        axes[1].set_xlabel('Predicted')
        axes[1].set_ylabel('Actual')
        axes[1].set_title('Confusion Matrix')
        
        # Cross-validation scores
        if self.cv_scores is not None:
            axes[2].bar(range(len(self.cv_scores)), self.cv_scores)
            axes[2].axhline(y=np.mean(self.cv_scores), color='r', linestyle='--', 
                           label=f'Mean: {np.mean(self.cv_scores):.3f}')
            axes[2].set_xlabel('CV Fold')
            axes[2].set_ylabel('ROC AUC Score')
            axes[2].set_title('Cross-Validation Scores')
            axes[2].legend()
        
        plt.tight_layout()
        plt.show()
    
    def get_feature_importance(self, feature_names: Optional[list] = None, top_k: int = 20) -> pd.DataFrame:
        """
        Get feature importance scores.
        
        Parameters:
        -----------
        feature_names : list, optional
            Names of the features
        top_k : int, default=20
            Number of top features to return
            
        Returns:
        --------
        DataFrame : Feature importance scores
        """
        if self.feature_importance is None:
            warnings.warn("Feature importance not available for this classifier type")
            return pd.DataFrame()
        
        if feature_names is None:
            feature_names = [f'Feature_{i}' for i in range(len(self.feature_importance))]
        
        importance_df = pd.DataFrame({
            'feature': feature_names[:len(self.feature_importance)],
            'importance': self.feature_importance
        })
        
        return importance_df.nlargest(top_k, 'importance')


def classify_anndata(adata, 
                    target_column: str,
                    positive_class: str,
                    classifier_type: str = 'logistic',
                    test_size: float = 0.2,
                    n_components: Optional[int] = None,
                    feature_selection_k: Optional[int] = None,
                    random_state: int = 42) -> Tuple[HighDimBinaryClassifier, Dict[str, Any]]:
    """
    Convenience function for binary classification on AnnData objects.
    
    Parameters:
    -----------
    adata : AnnData
        Single-cell data object
    target_column : str
        Column name in adata.obs containing the target variable
    positive_class : str
        Value in target_column to be treated as positive class
    classifier_type : str, default='logistic'
        Type of classifier ('logistic', 'random_forest', 'svm')
    test_size : float, default=0.2
        Proportion of data to use for testing
    n_components : int, optional
        Number of PCA components
    feature_selection_k : int, optional
        Number of features to select
    random_state : int, default=42
        Random state for reproducibility
        
    Returns:
    --------
    tuple : (fitted_classifier, evaluation_results)
    """
    try:
        import anndata
    except ImportError:
        raise ImportError("anndata package is required for this function")
    
    # Prepare data
    X = adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X
    y = (adata.obs[target_column] == positive_class).astype(int)
    
    # Split data
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=test_size, random_state=random_state, stratify=y
    )
    
    # Initialize and fit classifier
    classifier = HighDimBinaryClassifier(
        classifier_type=classifier_type,
        n_components=n_components,
        feature_selection_k=feature_selection_k,
        random_state=random_state
    )
    
    classifier.fit(X_train, y_train)
    
    # Evaluate
    results = classifier.evaluate(X_test, y_test)
    
    # Add additional info
    results['positive_class'] = positive_class
    results['target_column'] = target_column
    results['n_samples'] = len(y)
    results['n_features'] = X.shape[1]
    results['class_distribution'] = dict(zip(*np.unique(y, return_counts=True)))
    
    return classifier, results


# Example usage and utility functions
def example_usage():
    """
    Example of how to use the HighDimBinaryClassifier.
    """
    print("Example usage of HighDimBinaryClassifier:")
    print("=========================================")
    
    # Generate synthetic high-dimensional data
    from sklearn.datasets import make_classification
    
    X, y = make_classification(
        n_samples=1000,
        n_features=2000,  # High-dimensional
        n_informative=50,
        n_redundant=100,
        n_clusters_per_class=1,
        random_state=42
    )
    
    print(f"Generated dataset: {X.shape[0]} samples, {X.shape[1]} features")
    print(f"Class distribution: {dict(zip(*np.unique(y, return_counts=True)))}")
    
    # Split data
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42, stratify=y
    )
    
    # Initialize classifier with PCA and feature selection
    classifier = HighDimBinaryClassifier(
        classifier_type='logistic',
        n_components=50,  # Reduce to 50 components
        feature_selection_k=100,  # Select top 100 features after PCA
        random_state=42
    )
    
    # Fit the classifier
    print("\nFitting classifier...")
    classifier.fit(X_train, y_train)
    
    # Evaluate performance
    print("\nEvaluating performance...")
    results = classifier.evaluate(X_test, y_test)
    
    print(f"Accuracy: {results['accuracy']:.3f}")
    print(f"ROC AUC: {results['roc_auc']:.3f}")
    print(f"CV Score (mean ± std): {results['cv_mean_score']:.3f} ± {results['cv_std_score']:.3f}")
    
    return classifier, results


if __name__ == "__main__":
    example_usage()