"""
Armin Hrnjić, Kemijski inštitut armin.hrnjic@ki.si armin.hrnjic@outlook.com

Python script for a course 'Multi-scale Materials Modelling and Engineering' at UNG 

To do:

- Handle errors
- Make a GUI 

"""


import numpy as np
 
def func_example1(x):
    return np.exp(x)-3*x*x

def func_example2(x):
    return np.cos(x)

def func_example3(x):
    return np.exp(x/3)+x*x

def func_example4(x):
    return x
    
def diff_fd(func, x, dx):
    """
    Forward difference method for differentiation. 
    
    args:
    func - Function for the differentiation
    x - Value of the fanction to estimate differential 
    dx - Value of the forward difference 

    """
    return (func(x+dx)-func(x))/dx

def diff_bd(func, x, dx):
    """
    Backward difference method for differentiation. 
    
    args:
    func - Function for the differentiation
    x - Value of the fanction to estimate differential 
    dx - Value of the backward difference 

    """
    return (func(x-dx)-func(x))/dx

def diff_tpf(func,x,dx):
    """
    Not really sure about Three point formula differentiation method

    args:
    func - Function for the differentiation
    x - Value of the fanction to estimate differential 
    dx - Value of the difference 
    """
    return (0.5/dx)*(func(x+dx)-func(x-dx)) #-((dx*dx/6)*(func(x)*func(x)*func(x)))

def diff_fpf(func,x,dx):
    """
    Differentiation method using a five point function.

    args:
    func - Function for the differentiation
    x - Value of the fanction to estimate differential 
    dx - Value of the forward difference 
    """
    return (1/(12*dx))*(func(x-2*dx)-8*func(x-dx)+8*func(x+dx)-func(x+2*dx))

def diff_desc_fd(x,y):
    """
    Calculates a slope between each points of a descreate function. (Front difference)

    x - time (type: list)
    y - concentration (type: list)

    returns two lists bufx and buffer -- bufx - time without the last point -- buffer - slope between each point.
    """

    buffer=[(y[i+1]-y[i])/(x[i+1]-x[i]) for i in range(len(x)-1)]
    
    
    
    bufx=x[0:-1]
    
    return bufx,buffer


def diff_dess_fpf(x,y):

    dx=x[1]-x[0]
    diff=[(1/(12*dx))*(y[i]-8*y[i+1]+8*y[i+2]-y[i+3]) for i in range(len(x)-3)]

    

    return diff
        
        

def solve_newton(func,initial,toc,iter_limit):
    """
    Solves a function from func(). 
    The method used is Newtons method based.

    args:
    func - takes any function as an argument (This is added so the function can be modulated for future work!)
    initial - initial aproximation of the result
    toc - relative error (tolerance)
    iter_limit - limit of number of iterations
    
    To do:

    - HANDLE ERRORS IDIOT!!

    
    """
    
    
    rnd=-np.log10(toc)# Finds number of digits to round a result

    for i in range(iter_limit): # Loops for the iter_limit number of iterations
        p=initial-(func(initial)/diff_fpf(func,initial,toc)) #Newton method formula for a solution aproximation
        if np.absolute(p-initial)<toc: # Checks if the result is in the range of relative error
            print('The result is '+ str(round(p,int(rnd-1)))) 
            break

        initial=p
    print('Newton method ended')



def solve_bisection(func,lower_limit,upper_limit,toc,iter_limit):
    """
    Solves a function from func(). 
    The method used is bisection method based on the Intermediate Value Theorem.

    args:
    func - takes any function as an argument (This is added so the function can be modulated for future work!)
    lower_limit - lower estimate of the solution range
    higher_limit - higher estimate of the solution range
    toc - relative error (tolerance)
    iter_limit - limit of number of iterations
    
    To do:

    - Find a way to find multiple zeros for any order of polinomial
    - Write an algorithm that calculates iter_limit
    """
    
    rnd=-np.log10(toc) # Finds number of digits to round a result
    
    FA=func(lower_limit) # Finds a value of function for a lower limit as parameter
    for i in range(iter_limit):  # Loops for the iter_limit number of iterations
        TOC=((upper_limit-lower_limit)/2) 
        p=lower_limit+TOC # Calculates middle value of the range
        
        if func(p)==0 or TOC<toc: # Checks the result and prints it to the expected number of digits
            print('The result is '+ str(round(p,int(rnd-1))))
            break
        
        if FA*func(p)>0: # Sets a middle value as new lower limit or higer limit (f(lower_limit)*f(higher_limit) needs to be negative if the zero is between them)
            lower_limit=p
        else:
            upper_limit=p

    print('Bisectional method ended')       


def solve_secant(func, initial0, initial1, toc, iter_limit):
    """
    Function finds a solution to the func() using Newton method with secant algorithm. 
    args:
    func - takes any function as an argument (This is added so the function can be modulated for future work!)
    initial0 - first initial guess of the result
    initial1 - second initial guess of the result (These can be viewed as a range)
    toc - relative error (tolerance)
    iter_limit - limit of number of iterations

    """
    rnd=-np.log10(toc) # Finds number of digits to round a result
    q0=func(initial0) # Finds a value of function for a first initial guess as parameter
    q1=func(initial1) # Finds a value of function for a second initial guess as parameter
    

    for i in range(2,iter_limit): # Loops for the iter_limit number of iterations
        p=initial1-(q1*(initial1-initial0)/(q1-q0))  # Calculates using secant method formula 
        if np.absolute(p-initial1)<toc: # Checks the result and prints it to the expected number of digits
            print('The result is '+ str(round(p,int(rnd-1))))
            break
        initial0=initial1 # Sets new values for the next itteration
        q0=q1
        initial1=p
        q1=func(p)
    
    print('Secant method ended')


def integrate_trapezoid(func, a,b, n):
    """
    Function integration with trapezoid method 
    args:
    func - function for the integration
    a - lower limit for definitive integral
    b - upper limit for definitive integral
    n - step size (dx)

    returns value of the integral (type: float)
    """

    f_range=b-a
    step=f_range/n
    b1=a+step
    a1=a
    s=0
    for i in range(n):
        s=s+(step/2)*(func(a1)+func(b1))
        a1=b1
        b1=b1+step
        
    
    return s

def integrate_simpson(func, a,b, n):
    """
    Function integration with simpson method 
    args:
    func - function for the integration
    a - lower limit for definitive integral
    b - upper limit for definitive integral
    n - step size (dx)

    returns value of the integral (type: float)
    """

    f_range=b-a
    step=f_range/n
    b1=a+step
    a1=a
    mid=(a1+b1)/2
    s=0
    for i in range(n):
        s=s+((step/6)*(func(a1)+4*func(mid)+func(b1)))
        a1=b1
        b1=b1+step
        mid=(a1+b1)/2
        
    
    return s



if __name__ == "__main__":
    
    solve_bisection(func_example1,0,1,0.00001,100) # Example 6e on the page 75 range [0,1]
    solve_bisection(func_example1,3,5,0.00001,100) # Example 6e on the page 75 range [3,5]
    solve_secant(func_example1,0.5,0.75,0.001,100) # Example 6e on the page 75 range [0,1]
    solve_secant(func_example1,3,5,0.00001,100) # Example 6e on the page 75 range [3,5]
    solve_newton(func_example1,1,0.00001,100) # Example 6e on the page 75 with initial guess of 1

    print('Function 3 is e^x/3+x*x')
    print('Forward difference method result for function 3:')
    print(diff_fd(func_example3,0.5,0.01))
    print('Triple point formula differentiation result for function 3:')
    print(diff_tpf(func_example3,0.5,0.01))
    print('Five point formula differentiation result for function 3:')
    print( diff_fpf(func_example3,0.5,0.01))
    print('Analytical solution for differentiation of function 3:')
    print((1/3)*(6*0.5+ np.exp(0.5/3)))
    print('Integral value of function 3 in the range of 0 to 15 is:')
    print(integrate_trapezoid(func_example3,0,1,1))
    print('Integral value of function 3 in the range of 0 to 15 is:')
    print(integrate_simpson(func_example3,0,1,1))