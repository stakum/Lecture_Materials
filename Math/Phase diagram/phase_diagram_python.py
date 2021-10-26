#!/usr/bin/env python
# coding: utf-8

# #  相図のプロット
# Original code made by Prof. Kinefuchi
# $$
# \begin{cases}
#     \dot{x}_{1}=x_1-2x_2+1 \\
#     \dot{x}_{2}=x_1-x_2+3
# \end{cases}
# $$

# In[1]:


import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


# ## 微分方程式の定義

# In[2]:


# Ridht-hand side of the differential equation
def rhs(t,x):
    dx1 = x[0] - 2.0 * x[1] + 1.0
    dx2 = x[0] - x[1] + 3.0
    dxdt = [dx1,dx2]
    return dxdt


# ## 初期値の設定

# In[3]:


# initial condition
x_ini = [[-5,-7],[-5,-6],[-5,-5],[-5,-4],[-5,-3],[-5,-2],[-
5,-1],[-5,0]]


# In[4]:


# Interval of integration
t_span=[0,100]


# ## 微分方程式解法

# In[5]:


# Integration for each initial condition
for xi in x_ini:
    sol = solve_ivp(rhs,t_span,xi,rtol=1e-10,atol=1e-10)
    x1 = sol.y[0,:]
    x2 = sol.y[1,:]
    plt.plot(x1,x2,'b-')

    # Figure
plt.xlabel('x1')
plt.ylabel('x2')
plt.axis('equal')
plt.show()


# In[ ]:




