import numpy as np
from sklearn.feature_selection import SelectFromModel, SelectKBest, mutual_info_classif
from sklearn.linear_model import LogisticRegression

path_to_selected_features = './res/MIFSLR_features.txt'


class MIFSLR(LogisticRegression):

    def __init__(self, penalty='l2', dual=False, tol=1e-4, C=1.0,
                 fit_intercept=True, intercept_scaling=1, class_weight=None,
                 random_state=None, solver='liblinear', max_iter=100,
                 multi_class='ovr', verbose=0, warm_start=False, n_jobs=1, k=1000):

        super(MIFSLR, self).__init__(penalty, dual, tol, C, fit_intercept, intercept_scaling, class_weight, random_state,
                                     solver, max_iter, multi_class, verbose, warm_start, n_jobs)
        self.mifs = SelectKBest(mutual_info_classif, k=k)
        # self.dummy = DummyClassifier()

    def fit(self, X, y, sample_weight=None):
        self.mifs.fit(X, y)
        # with open(path_to_selected_features, 'a') as f:
        #     f_list = []
        #     for i in self.model.get_support(True):
        #         f_list.append(str(i))
        #     f.write(';'.join(f_list) + '\n')
        X_new = self.mifs.transform(X)
        # if X_new.shape[1] == 0:
        #     X_new = np.zeros(shape=(X_new.shape[0], 1))
        super(MIFSLR, self).fit(X_new, y, sample_weight)
        return self

    def predict_proba(self, X):
        X_new = self.mifs.transform(X)
        # if X_new.shape[1] == 0:
        #     X_new = np.zeros(shape=(X_new.shape[0], 1))
        return super(MIFSLR, self).predict_proba(X_new)

    def decision_function(self, X):
        X_new = self.mifs.transform(X)
        # if X_new.shape[1] == 0:
        #     X_new = np.zeros(shape=(X_new.shape[0], 1))
        return super(MIFSLR, self).decision_function(X_new)
