import numpy as np 
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LogisticRegression
from matplotlib import pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D

def generate_sphere_datapoints(r=1, N=20):
    theta_vec = np.linspace(0, np.pi, N)
    phi_vec = np.linspace(0, 2*np.pi, 2*N)
    [th, phi] = np.meshgrid(theta_vec, phi_vec)
    R = r*np.ones(th.shape)
    
    x = R*np.sin(th)*np.cos(phi)
    y = R*np.sin(th)*np.sin(phi)
    z = R*np.cos(th)

    return x, y, z

if __name__=="__main__":
    r = 1;
    N = 20;
    x, y, z = generate_sphere_datapoints(r, N)
    X = np.vstack([x.flatten(), y.flatten(), z.flatten()])

    xp = x - z;
    yp = y;
    zp = -z;
    Xp = np.vstack([xp.flatten(), yp.flatten(), zp.flatten()])


    normXp = np.linalg.norm(Xp, axis=0)
    notX = X[:, normXp >= r**2]
    X = X[:, normXp < r**2]
    Xp = Xp[:, normXp < r**2]
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.scatter(X[0,:], X[1,:], X[2,:], c='r', marker='o')
    ax.scatter(Xp[0,:], Xp[1,:], Xp[2,:], c='b', marker='o')

    ax.set_xlabel('x_cm')
    ax.set_ylabel('v_cm')
    ax.set_zlabel('x_leg')

    plt.show()

    t = np.ones((X.shape[1], 1))
    X = np.vstack([X.T, notX.T])
    t = np.vstack([t, np.zeros((notX.shape[1], 1))])

    poly = PolynomialFeatures(4)
    X = poly.fit_transform(X) 

    classifier = LogisticRegression().fit(X, t)
    coef = classifier.coef_[0]
    coef = np.atleast_2d(coef)
    coef = coef
    print(coef.shape)
    print(coef@X.T > 0)
    print(sum((classifier.predict(X) - t.flatten())**2))
