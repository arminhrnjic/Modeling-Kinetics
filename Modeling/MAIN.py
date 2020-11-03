import numpy as np
import numanal as na
import matplotlib.pyplot as plt

def func_example1(x): 
    #Example 6e on the page 75
    return np.exp(x)-3*x*x

def func_example2(x):
    return np.cos(x)

def func_example3(x):
    #Example 5b on the page 177
    return np.exp(x/3)+x*x
"""
print('First task: ')
na.solve_newton(func_example1,1,0.001,1000)
na.solve_newton(func_example1,3,0.001,1000)
na.solve_bisection(func_example1,0,1,0.00001,100) # Example 6e on the page 75 range [0,1]
na.solve_bisection(func_example1,3,5,0.00001,100) # Example 6e on the page 75 range [3,5]
"""


#Exam 5e data
def exam_5e_data():

    exam5e_data=[[0,1,2,3,4,5,6,7,8,9,10],
                [2.0000,1.9926,1.9853,1.9780,1.9707,1.9634,1.9562,1.9490,1.9418,1.9346,1.9275],
                [2.0000,1.9463,1.8941,1.8433,1.7938,1.7457,1.6988,1.6533,1.6089,1.5657,1.5237],
                [2.0000,1.7904,1.6028,1.4349,1.2845,1.1499,1.0294,0.9215,0.8250,0.7385,0.6611]]

    ex5edata=np.array(exam5e_data)

    plt.xlabel('Time (min)')
    plt.ylabel('C (mol/L)')
    plt.plot(ex5edata[0],ex5edata[1], label='200°C')

    #plt.figure(2)
    plt.plot(ex5edata[0],ex5edata[2],label='300°C')

    #plt.figure(3)
    plt.plot(ex5edata[0],ex5edata[3],label='400°C')
    plt.legend()
    plt.show()

def exam_3a_data():
    e3a_data=[  [0,     1,      2,      3,      4,      5,      6,      7,      8,      9,      10],
                [1.0000,0.9732,	0.9471,	0.9216,	0.8969,	0.8728,	0.8494,	0.8266,	0.8044,	0.7829,	0.7619],
                [2.0000,1.9463,	1.8941,	1.8433,	1.7938,	1.7457,	1.6988,	1.6533,	1.6089,	1.5657,	1.5237],
                [4.0000,3.8927,	3.7882,	3.6866,	3.5876,	3.4914,	3.3977,	3.3065,	3.2178,	3.1314,	3.0474]]

    e3a=np.array(e3a_data)
    
    plt.figure(1)
    plt.xlabel('Time (min)')
    plt.ylabel('C (mol/L)')
    plt.plot(e3a[0],e3a[1], label='1.00 mol/L')
       
    plt.plot(e3a[0],e3a[2],label='2.00 mol/L')
    plt.plot(e3a[0],e3a[3],label='3.00 mol/L')
    plt.legend()

    plt.figure(2)
    plt.xlabel('Time (min)')
    plt.ylabel('-dC/dt (mol/Ls)')
    time1=na.diff_desc_fd(e3a[0],e3a[1])[0]
    time2=na.diff_desc_fd(e3a[0],e3a[2])[0]
    time3=na.diff_desc_fd(e3a[0],e3a[1])[0]
    rate1=[-i for i in na.diff_desc_fd(e3a[0],e3a[1])[1]]
    rate2=[-i for i in na.diff_desc_fd(e3a[0],e3a[2])[1]]
    rate3=[-i for i in na.diff_desc_fd(e3a[0],e3a[3])[1]]
    plt.plot( na.diff_desc_fd(e3a[0],e3a[1])[0], rate1 ,label='1.00 mol/L') 
    plt.plot( na.diff_desc_fd(e3a[0],e3a[2])[0], rate2 ,label='2.00 mol/L') 
    plt.plot( na.diff_desc_fd(e3a[0],e3a[3])[0], rate3 ,label='3.00 mol/L') 
    print(rate1[0],rate2[0],rate3[0])

    plt.figure(3)
    plt.xlabel('t (min)')
    plt.ylabel('ln(rate)')
    ln_of_rate1=[np.log(i) for i in rate1]
    ln_of_rate2=[np.log(i) for i in rate2]
    ln_of_rate3=[np.log(i) for i in rate3]
    plt.plot(time1, ln_of_rate1, label='1.00 mol/L')
    plt.plot(time2, ln_of_rate2, label='2.00 mol/L')
    plt.plot(time3, ln_of_rate3, label='3.00 mol/L')

    suma= na.diff_desc_fd(time1,ln_of_rate1)[1] + na.diff_desc_fd(time2,ln_of_rate2)[1] +na.diff_desc_fd(time3,ln_of_rate3)[1]
    rate_constant=sum(suma)/len(suma)

    time_analytical=[i for i in range(181)]
    concentration_reactant_analytical1=[1*np.exp(-0.0272*i) for i in range(181)]
    concentration_reactant_analytical2=[1*np.exp(-0.0272*i) for i in range(181)]
    concentration_reactant_analytical3=[1*np.exp(-0.0272*i) for i in range(181)]

    concentration_product_analytical1=[1-i for i in concentration_reactant_analytical1]

    plt.figure(4)
    plt.xlabel('Time (min)')
    plt.ylabel('C (mol/L)')
    plt.plot(time_analytical,concentration_reactant_analytical1,label='Reactant 1.00 mol/L')
    plt.plot(time_analytical,concentration_product_analytical1,label='Prodcuts 1.00 mol/L')





    plt.legend()
    plt.show()

exam_3a_data()
