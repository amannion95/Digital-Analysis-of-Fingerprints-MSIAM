import numpy as np
import matplotlib.pyplot as plt
import os.path
from mpl_toolkits.mplot3d import Axes3D

import sys


x=[]
y=[]

x_opti=[]
y_opti=[]

x_opti_cov_like=[]
y_opti_cov_like=[]


x2=[]
y2=[]
z2=[]

x2_opti=[]
y2_opti=[]
z2_opti=[]

x2_opti_cov_like=[]
y2_opti_cov_like=[]
z2_opti_cov_like=[]


if(os.path.isfile('topython_1D.txt') ):
	i=0
	with open('topython_1D.txt', 'r') as my_file:
		for line in my_file:
		
			a=line
			if i%2==0:
				x.append(float(a))
			else:
				y.append(float(a))
			i+=1
	plt.plot(x,y)
	plt.title('1D classic method')
	print("for 1D classic method :\n","px = ",x[y.index(min(y))])
plt.show()

if(os.path.isfile('topython_1D_opti.txt') ):
	

	i=0
	with open('topython_1D_opti.txt', 'r') as my_file:
		for line in my_file:
		
			a=line
			if i%2==0:
				x_opti.append(float(a))
			else:
				y_opti.append(float(a))
			i+=1
	plt.plot(x_opti,y_opti)
	plt.title('1D opti method ')
	print("for 1D opti method :\n","px = ",x_opti[y_opti.index(min(y_opti))])


if(os.path.isfile('topython_2D.txt') ):
	i=0
	with open('topython_2D.txt', 'r') as my_file:
		for line in my_file:
			
			a=line
			if i%3==0:
				x2.append(float(a))
			elif i%3==1:
				y2.append(float(a))
			else:
				z2.append(float(a))
			i+=1
	
	fig = plt.figure()
	ax3 = plt.axes(projection='3d')
	ax3.plot(x2,y2,z2)	
	plt.title('2D classic method')
	print("for 2D classic method :\n","px =",x2[z2.index(min(z2))],"  py = ",y2[z2.index(min(z2))])
	
if(os.path.isfile('topython_2D_opti.txt') ):
	i=0
	with open('topython_2D_opti.txt', 'r') as my_file:
		for line in my_file:
			
			a=line
			if i%3==0:
				x2_opti.append(float(a))
			elif i%3==1:
				y2_opti.append(float(a))
			else:
				z2_opti.append(float(a))
			i+=1

	fig = plt.figure()
	ax4 = plt.axes(projection='3d')
	ax4.plot(x2_opti,y2_opti,z2_opti)
	
	plt.title('2D opti method')
	print("for 2D opti method :\n","px =",x2_opti[z2_opti.index(min(z2_opti))],"  py = ",y2_opti[z2_opti.index(min(z2_opti))])


if(os.path.isfile('topython_2D_opti_covariance.txt') ):
	i=0
	with open('topython_2D_opti_covariance.txt', 'r') as my_file:
		for line in my_file:
			
			a=line
			if i%3==0:
				x2_opti_cov_like.append(float(a))
			elif i%3==1:
				y2_opti_cov_like.append(float(a))
			else:
				z2_opti_cov_like.append(float(a))
			i+=1
	
	fig = plt.figure()
	ax5 = plt.axes(projection='3d')
	ax5.plot(x2_opti_cov_like,y2_opti_cov_like,z2_opti_cov_like)
	plt.title('2D opti method with cov like error')
	print("for 2D opti method with cov like error :\n","px =",x2_opti_cov_like[z2_opti_cov_like.index(max(z2_opti_cov_like))],"  py = ",y2_opti_cov_like[z2_opti_cov_like.index(max(z2_opti_cov_like))])

plt.show()
