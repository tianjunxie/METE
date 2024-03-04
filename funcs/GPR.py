### GPR class
import os
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C
# from sklearn.metrics import r2_score
# from sklearn.model_selection import train_test_split
import pickle

class GPRModel:
    def __init__(self, autosave, filename):
        self.kernel = C(1.0, (1e-4, 1e4)) * RBF(0.001, (1e-5, 1e5))
        self.alpha = 0.01
        self.random_state = 0
        self.filename = filename
        self.autosave = autosave
        self.model = None
        if os.path.isfile(self.filename):
            self.load(self.filename)

    ### train if given input
    def fit(self, X, y):
        self.model = GaussianProcessRegressor(kernel=self.kernel,\
                                              alpha=self.alpha,\
                                              n_restarts_optimizer=5,\
                                              random_state=self.random_state)
        r2 = self.model.fit(X, y).score(X, y)
        return(r2)
    
    ### make predictions
    def predict(self, X):
        return self.model.predict(X, return_std=True)
    
    ### save the GPR model
    def save(self, filename):
        with open(filename, 'wb') as f:
            pickle.dump(self.model, f)
    
    ### load from cpickle file
    def load(self, filename):
        with open(filename, 'rb') as f:
            self.model = pickle.load(f)


# print(y_pred, std)


