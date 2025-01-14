from sklearn.gaussian_process import GaussianProcessClassifier, GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, DotProduct, WhiteKernel, Product, Sum, ConstantKernel, RBF, Matern, ExpSineSquared
from sklearn.metrics.pairwise import polynomial_kernel, rbf_kernel
import numpy as np

def gpc_poly_prod(x,y,degree):
    kernel1 = (DotProduct(sigma_0=0.1) ** degree)
    kernel2 = (DotProduct(sigma_0=0.1) ** degree)
    kernel = Product(kernel1,kernel2)
    gpc = GaussianProcessClassifier(kernel=kernel).fit(x,y)
    return gpc


def gpc_poly(x,y,degree):
    kernel = (DotProduct() ** degree)
    gpc = GaussianProcessClassifier(kernel=kernel).fit(x,y)
    return gpc


def gp_test(x,y,degree1):
    k = (DotProduct() ** degree1)
    gpr =  GaussianProcessRegressor(kernel=k).fit(x,y)
    GramM = polynomial_kernel(x,degree=degree1,coef0=(gpr.kernel_.theta))
    GramM2 = gpr.kernel_(x)

    return {'gpr': gpr, 'GramMatrix': GramM, 'GramMatrix2': GramM2}



def gp(x,y,degree1,degree2,k1,k2,operator):
  if operator == "prod":
    #print("PROD")
    if k1 == "poly" and k2 == "poly":
      k     =  ConstantKernel()*( (DotProduct() ** degree1)*(DotProduct() ** degree2) )  + WhiteKernel()
      #k     =  ( (DotProduct() ** degree1)*(DotProduct() ** degree2) )  + WhiteKernel()
      #k     =  ( (DotProduct(sigma_0=0.1) ** degree1)*(DotProduct(sigma_0=0.1) ** degree2) ) + WhiteKernel(noise_level=1.0)
      gpr   = GaussianProcessRegressor(kernel=k).fit(x,y)
      #GramM = np.exp((gpr.kernel_.theta[1]))*(polynomial_kernel(x,degree=degree1,coef0=np.exp(gpr.kernel_.theta[2]))*polynomial_kernel(x,degree=degree2,coef0=np.exp(gpr.kernel_.theta[3])))

    if k1 == "rbf" and k2 == "rbf":
      k     = ConstantKernel()*( (RBF())*(RBF()) ) + WhiteKernel()
      k     = ( (RBF())*(RBF()) ) + WhiteKernel()
      gpr   = GaussianProcessRegressor(kernel=k).fit(x,y)
      #GramM = rbf_kernel(x,gamma=gpr.kernel_.theta[2])*rbf_kernel(x,gamma=gpr.kernel_.theta[3])


    if k1 == "poly" and k2 == "rbf":
      k     = ConstantKernel()*( (DotProduct() ** degree1)*(RBF()) ) + WhiteKernel()
      #k     = ( (DotProduct() ** degree1)*(RBF()) ) + WhiteKernel()
      gpr   = GaussianProcessRegressor(kernel=k).fit(x,y)
      #GramM = polynomial_kernel(x,degree=degree1,coef0=gpr.kernel_.theta[2])*rbf_kernel(x,gamma=gpr.kernel_.theta[3])


    if k1 == "rbf" and k2 == "poly":
      k     = ConstantKernel()*( ( RBF() )*( DotProduct() ** degree2  ) ) + WhiteKernel()
      gpr   = GaussianProcessRegressor(kernel=k).fit(x,y)
      #GramM = rbf_kernel(x,gamma=gpr.kernel_.theta[2])*polynomial_kernel(x,degree=degree1,coef0=gpr.kernel_.theta[3])
      
    if k1 == "matern1.5" and k2 == "matern1.5":
      k     = ConstantKernel()*( ( Matern(nu=1.5) )*( Matern(nu=1.5)  ) ) + WhiteKernel()
      k     = ( ( Matern(nu=1.5) )*( Matern(nu=1.5)  ) ) + WhiteKernel()
      gpr   = GaussianProcessRegressor(kernel=k).fit(x,y)
      
    if k1 == "matern2.5" and k2 == "matern2.5":
      k     = ConstantKernel()*( ( Matern(nu=2.5) )*( Matern(nu=2.5)  ) ) + WhiteKernel()
      k     = ( ( Matern(nu=2.5) )*( Matern(nu=2.5)  ) ) + WhiteKernel()
      gpr   = GaussianProcessRegressor(kernel=k).fit(x,y)
      
    if k1 == "expsine2" and k2 == "expsine2":
      k     = ConstantKernel()*( ( ExpSineSquared() )*( ExpSineSquared()  ) ) + WhiteKernel()
      gpr   = GaussianProcessRegressor(kernel=k).fit(x,y)  
      
      
    if k1 == "poly" and k2 == "matern1.5":
      k     = ConstantKernel()*( ( DotProduct() ** degree1 )*( Matern(nu=1.5)  ) ) + WhiteKernel()
      gpr   = GaussianProcessRegressor(kernel=k).fit(x,y)  
      
    if k1 == "poly" and k2 == "matern2.5":
      k     = ConstantKernel()*( ( DotProduct() ** degree1 )*( Matern(nu=2.5)  ) ) + WhiteKernel()
      gpr   = GaussianProcessRegressor(kernel=k).fit(x,y) 
      
    if k1 == "matern1.5" and k2 == "poly":
      k     = ConstantKernel()*( ( Matern(nu=1.5) )*( DotProduct() ** degree2  ) ) + WhiteKernel()
      gpr   = GaussianProcessRegressor(kernel=k).fit(x,y)  
      
    if k1 == "matern2.5" and k2 == "poly":
      k     = ConstantKernel()*( ( Matern(nu=2.5) )*( DotProduct() ** degree2  ) ) + WhiteKernel()
      gpr   = GaussianProcessRegressor(kernel=k).fit(x,y)  
      
    
    
    if k1 == "poly" and k2 == "rbf":
      k     = ConstantKernel()*( ( DotProduct() ** degree1 )*( RBF()  ) ) + WhiteKernel()
      gpr   = GaussianProcessRegressor(kernel=k).fit(x,y)  
      
    if k1 == "rbf" and k2 == "poly":
      k     = ConstantKernel()*( ( RBF() )*( DotProduct() ** degree2  ) ) + WhiteKernel()
      gpr   = GaussianProcessRegressor(kernel=k).fit(x,y)   
      
      
      
  if operator == "sum":
      #print("SUM")
      if k1 == "poly" and k2 == "poly":
        k     =  ConstantKernel()*( (DotProduct() ** degree1) + (DotProduct() ** degree2) ) + WhiteKernel()
        gpr   = GaussianProcessRegressor(kernel=k).fit(x,y)
        #GramM = np.exp((gpr.kernel_.theta[1]))*(polynomial_kernel(x,degree=degree1,coef0=np.exp(gpr.kernel_.theta[2]))*polynomial_kernel(x,degree=degree2,coef0=np.exp(gpr.kernel_.theta[3])))
  
      if k1 == "rbf" and k2 == "rbf":
        k     = ConstantKernel()*( (RBF()) + (RBF()) ) + WhiteKernel()
        gpr   = GaussianProcessRegressor(kernel=k).fit(x,y)
        #GramM = rbf_kernel(x,gamma=gpr.kernel_.theta[2])*rbf_kernel(x,gamma=gpr.kernel_.theta[3])
  
  
      if k1 == "poly" and k2 == "rbf":
        k     = ConstantKernel()*( (DotProduct() ** degree1) + (RBF()) ) + WhiteKernel()
        gpr   = GaussianProcessRegressor(kernel=k).fit(x,y)
        #GramM = polynomial_kernel(x,degree=degree1,coef0=gpr.kernel_.theta[2])*rbf_kernel(x,gamma=gpr.kernel_.theta[3])
  
  
      if k1 == "rbf" and k2 == "poly":
        k     = ConstantKernel()*( ( RBF() ) + ( DotProduct() ** degree2  ) ) + WhiteKernel()
        gpr   = GaussianProcessRegressor(kernel=k).fit(x,y)
        #GramM = rbf_kernel(x,gamma=gpr.kernel_.theta[2])*polynomial_kernel(x,degree=degree1,coef0=gpr.kernel_.theta[3])  


  if operator == "single":
      #print(" -- SINGLE -- ")
      if k1 == "poly":
        k     = ConstantKernel()*( ( DotProduct()**degree1  ) ) + WhiteKernel()
        #k     = ( ( DotProduct(sigma_0=1)**degree1  ) ) + WhiteKernel()
        #k     = ConstantKernel()*( ( DotProduct(sigma_0=1)**degree1  ) ) + WhiteKernel()
        gpr   = GaussianProcessRegressor(kernel=k).fit(x,y)
      if k1 == "rbf":
        k     = ConstantKernel()*( ( RBF()  ) ) + WhiteKernel()
        gpr   = GaussianProcessRegressor(kernel=k).fit(x,y)
      if k1 == "matern1.5":
        k     = ConstantKernel()*( ( Matern(nu=1.5)  ) ) + WhiteKernel()
        gpr   = GaussianProcessRegressor(kernel=k).fit(x,y)   
      if k1 == "matern2.5":
        k     = ConstantKernel()*( ( Matern(nu=2.5)  ) ) + WhiteKernel()
        gpr   = GaussianProcessRegressor(kernel=k).fit(x,y) 

  GramM = gpr.kernel_(x)



  return {'gpr': gpr, 'GramMatrix': GramM}
  
  
  
  
  
  
  
  
  
