from sklearn.datasets import load_iris
from sklearn.datasets import load_wine

iris = load_iris().data

for i in range(0, iris.shape[0]):
    for j in range(0, iris.shape[1]):
        iris[i][j] *= 10
        
wine = load_wine().data

for i in range(0, wine.shape[0]):
    for j in range(0, wine.shape[1]):
        wine[i][j] *= 1000

ip = 'localhost';
port = 7766;