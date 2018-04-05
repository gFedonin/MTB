from SFS import SequentialFeatureSelector as SFS
from sklearn.feature_selection import SelectKBest, mutual_info_classif
from sklearn.linear_model import LogisticRegression

# path_to_selected_features = './res/MIGFFSLR_features.txt'
# print_features = False


class MIGFFSLR(LogisticRegression):

    def __init__(self, penalty='l2', dual=False, tol=1e-4, C=1.0,
                 fit_intercept=True, intercept_scaling=1, class_weight=None,
                 random_state=None, solver='liblinear', max_iter=100,
                 multi_class='ovr', verbose=0, warm_start=False, n_jobs=1, mi_k=10000, fs_k=1000, cv=3,
                 forward=True, floating=True, scoring='f1'):

        super(MIGFFSLR, self).__init__(penalty, dual, tol, C, fit_intercept, intercept_scaling, class_weight,
                                     random_state, solver, max_iter, multi_class, verbose, warm_start, n_jobs)
        self.mifs = SelectKBest(mutual_info_classif, k=mi_k)
        self.sfs = SFS(LogisticRegression(penalty, dual, tol, C, fit_intercept, intercept_scaling, class_weight,
                                          random_state, solver, max_iter, multi_class, verbose, warm_start, n_jobs=1),
               k_features=(1, fs_k),
               forward=forward,
               floating=floating,
               # verbose=2,
               scoring=scoring,
               cv=cv, n_jobs=n_jobs)

    def fit(self, X, y, sample_weight=None):
        self.mifs.fit(X, y)
        X_mi = self.mifs.transform(X)
        self.sfs.fit(X_mi, y)
        X_sfs = self.sfs.transform(X_mi)
        super(MIGFFSLR, self).fit(X_sfs, y, sample_weight)
        return self

    def predict_proba(self, X):
        X_mi = self.mifs.transform(X)
        X_sfs = self.sfs.transform(X_mi)
        return super(MIGFFSLR, self).predict_proba(X_sfs)

    def decision_function(self, X):
        X_mi = self.mifs.transform(X)
        X_sfs = self.sfs.transform(X_mi)
        return super(MIGFFSLR, self).decision_function(X_sfs)
