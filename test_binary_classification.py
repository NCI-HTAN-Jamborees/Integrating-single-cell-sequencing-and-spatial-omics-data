"""
Test suite for the binary classification module.

This module contains basic tests to ensure the binary classification functionality
works correctly with different types of input data and configurations.
"""

import numpy as np
import pandas as pd
from sklearn.datasets import make_classification
from sklearn.model_selection import train_test_split
import warnings

# Import the module to test
from binary_classification import HighDimBinaryClassifier, classify_anndata


def test_basic_functionality():
    """Test basic classifier functionality."""
    # Generate test data
    X, y = make_classification(
        n_samples=200,
        n_features=100,
        n_informative=20,
        n_redundant=10,
        random_state=42
    )
    
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.3, random_state=42
    )
    
    # Test logistic regression
    classifier = HighDimBinaryClassifier(
        classifier_type='logistic',
        n_components=10,
        random_state=42
    )
    
    # Fit and predict
    classifier.fit(X_train, y_train)
    predictions = classifier.predict(X_test)
    probabilities = classifier.predict_proba(X_test)
    
    # Basic assertions
    assert len(predictions) == len(y_test)
    assert probabilities.shape == (len(y_test), 2)
    assert classifier.is_fitted
    assert classifier.cv_scores is not None
    
    # Evaluate
    results = classifier.evaluate(X_test, y_test)
    assert 'accuracy' in results
    assert 'roc_auc' in results
    assert 0 <= results['accuracy'] <= 1
    assert 0 <= results['roc_auc'] <= 1
    
    print("✓ Basic functionality test passed")


def test_different_classifiers():
    """Test different classifier types."""
    X, y = make_classification(
        n_samples=150,
        n_features=50,
        n_informative=15,
        random_state=42
    )
    
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.3, random_state=42
    )
    
    classifiers = ['logistic', 'random_forest', 'svm']
    
    for clf_type in classifiers:
        classifier = HighDimBinaryClassifier(
            classifier_type=clf_type,
            n_components=10,
            random_state=42
        )
        
        classifier.fit(X_train, y_train)
        results = classifier.evaluate(X_test, y_test)
        
        assert results['accuracy'] > 0.5  # Should be better than random
        print(f"✓ {clf_type} classifier test passed")


def test_feature_selection():
    """Test feature selection functionality."""
    X, y = make_classification(
        n_samples=200,
        n_features=100,
        n_informative=20,
        random_state=42
    )
    
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.3, random_state=42
    )
    
    # With feature selection
    classifier = HighDimBinaryClassifier(
        classifier_type='random_forest',
        n_components=20,
        feature_selection_k=10,
        random_state=42
    )
    
    classifier.fit(X_train, y_train)
    results = classifier.evaluate(X_test, y_test)
    
    assert classifier.feature_selector is not None
    assert results['accuracy'] > 0.5
    
    # Test feature importance
    feature_importance = classifier.get_feature_importance(top_k=5)
    assert len(feature_importance) <= 5
    assert 'feature' in feature_importance.columns
    assert 'importance' in feature_importance.columns
    
    print("✓ Feature selection test passed")


def test_pca_auto_components():
    """Test automatic PCA component determination."""
    X, y = make_classification(
        n_samples=100,
        n_features=50,
        n_informative=10,
        random_state=42
    )
    
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.3, random_state=42
    )
    
    # Let PCA components be determined automatically
    classifier = HighDimBinaryClassifier(
        classifier_type='logistic',
        n_components=None,  # Auto-determine
        random_state=42
    )
    
    classifier.fit(X_train, y_train)
    
    assert classifier.pca is not None
    assert classifier.pca.n_components_ > 0
    assert classifier.pca.n_components_ <= min(X_train.shape[0] - 1, X_train.shape[1])
    
    print("✓ Auto PCA components test passed")


def test_edge_cases():
    """Test edge cases and error handling."""
    X, y = make_classification(
        n_samples=100,
        n_features=20,
        n_informative=5,
        random_state=42
    )
    
    # Test with invalid classifier type
    try:
        classifier = HighDimBinaryClassifier(classifier_type='invalid')
        classifier.fit(X, y)
        assert False, "Should have raised ValueError"
    except ValueError:
        print("✓ Invalid classifier type error handling passed")
    
    # Test with non-binary labels
    y_multi = np.array([0, 1, 2] * (len(y) // 3) + [0] * (len(y) % 3))
    try:
        classifier = HighDimBinaryClassifier()
        classifier.fit(X, y_multi)
        assert False, "Should have raised ValueError for non-binary labels"
    except ValueError:
        print("✓ Non-binary labels error handling passed")
    
    # Test prediction before fitting
    try:
        classifier = HighDimBinaryClassifier()
        classifier.predict(X)
        assert False, "Should have raised ValueError"
    except ValueError:
        print("✓ Prediction before fitting error handling passed")


def test_anndata_integration():
    """Test AnnData integration if available."""
    try:
        import anndata as ad
        
        # Create test data
        X, y = make_classification(
            n_samples=100,
            n_features=50,
            n_informative=10,
            random_state=42
        )
        
        # Create AnnData object
        adata = ad.AnnData(X=X)
        adata.obs['cell_type'] = ['Class_A' if label == 0 else 'Class_B' for label in y]
        
        # Test classify_anndata function
        classifier, results = classify_anndata(
            adata=adata,
            target_column='cell_type',
            positive_class='Class_B',
            classifier_type='logistic',
            test_size=0.3,
            random_state=42
        )
        
        assert classifier.is_fitted
        assert 'accuracy' in results
        assert 'positive_class' in results
        assert results['positive_class'] == 'Class_B'
        
        print("✓ AnnData integration test passed")
        
    except ImportError:
        print("✓ AnnData not available, skipping integration test")


def run_all_tests():
    """Run all tests."""
    print("Running binary classification tests...")
    print("=" * 50)
    
    # Suppress warnings for cleaner output
    warnings.filterwarnings('ignore')
    
    test_basic_functionality()
    test_different_classifiers()
    test_feature_selection()
    test_pca_auto_components()
    test_edge_cases()
    test_anndata_integration()
    
    print("=" * 50)
    print("All tests completed successfully! ✓")


if __name__ == "__main__":
    run_all_tests()