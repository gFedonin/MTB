import numpy as np
from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.utils import check_X_y, check_array
from sklearn.utils.multiclass import unique_labels
from sklearn.utils.validation import check_is_fitted


class VocabularyClf(BaseEstimator, ClassifierMixin):

    def fit(self, X, y):
        # Check that X and y have correct shape
        X, y = check_X_y(X, y, accept_sparse=True)
        # Store the classes seen during fit
        self.classes_ = unique_labels(y)
        self.X_ = X
        self.y_ = y
        # Return the classifier
        return self

    def decision_function(self, X):
        res = np.zeros((X.shape[0], 2))
        nnz = X.getnnz(axis=1)
        for i in range(X.shape[0]):
            if nnz[i] > 0:
                res[i, 1] = 1
            else:
                res[i, 0] = 1
        return res

    def predict(self, X):
        # Check is fit had been called
        check_is_fitted(self, ['X_', 'y_'])

        # Input validation
        X = check_array(X, accept_sparse=True)

        D = self.decision_function(X)
        return self.classes_[np.argmax(D, axis=1)]

    def predict_proba(self, X):
        return self.decision_function(X)