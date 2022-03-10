'''
Description: Forward Euler Approximation code for PC5210 HW2
'''

import numpy as np
from matplotlib import pyplot as plt

# Global parameters
omega = 1

# ==============================
#  Forward euler approximation
# ==============================
def forward_euler_approx(q_initial, q_final, p_initial, delta):
    # Initalize parameters
    q0 = q_initial
    qf = q_final
    p0 = p_initial
    n = int((qf-q0)//delta) # Number of points, floor div and set to int
    
    # Initialize arrays
    q_array = np.linspace(q0, qf, n)
    p_array = np.zeros(n)

    # Forward Euler Approx
    p_array[0] = p0
    for i in range(1,n):
        q_array[i] = delta*p_array[i-1] + q_array[i-1]
        p_array[i] = -delta*omega**2*q_array[i-1] + p_array[i-1]

    return(q_array, p_array)

# Data for plotting
deltas = [0.2, 0.1, 0.05, 0.01]
q_start = 5
q_end = 20
q1, p1 = forward_euler_approx(5, 20, 0, deltas[0])
q2, p2 = forward_euler_approx(5, 20, 0, deltas[1])
q3, p3 = forward_euler_approx(5, 20, 0, deltas[2])
q4, p4 = forward_euler_approx(5, 20, 0, deltas[3])

# Multiplot q(t)
def multiplot_qt():
    fig, axs = plt.subplots(2, 2)
    fig.suptitle(r'q(t) for various $\delta$')
    axs[0, 0].plot(np.linspace(0, int(q_end-q_start), int((q_end-q_start)//deltas[0])), q1)
    axs[0, 0].set_title(r'$\delta$ = ' + str(deltas[0]))
    axs[0, 1].plot(np.linspace(0, int(q_end-q_start), int((q_end-q_start)//deltas[1])), q2)
    axs[0, 1].set_title(r'$\delta$ = ' + str(deltas[1]))
    axs[1, 0].plot(np.linspace(0, int(q_end-q_start), int((q_end-q_start)//deltas[2])), q3)
    axs[1, 0].set_title(r'$\delta$ = ' + str(deltas[2]))
    axs[1, 1].plot(np.linspace(0, int(q_end-q_start), int((q_end-q_start)//deltas[3])), q4)
    axs[1, 1].set_title(r'$\delta$ = ' + str(deltas[3]))
    for ax in axs.flat:
        ax.set(xlabel='t', ylabel='q(t)')
    # for ax in axs.flat:
    #     ax.label_outer()
    plt.tight_layout()
    plt.show()

# Multiplot p(q)
def multiplot_pq():
    fig, axs = plt.subplots(2, 2)
    fig.suptitle(r'p(q) for various $\delta$')
    axs[0, 0].plot(q1, p1)
    axs[0, 0].set_title(r'$\delta$ = ' + str(deltas[0]))
    axs[0, 1].plot(q2, p2)
    axs[0, 1].set_title(r'$\delta$ = ' + str(deltas[1]))
    axs[1, 0].plot(q3, p3)
    axs[1, 0].set_title(r'$\delta$ = ' + str(deltas[2]))
    axs[1, 1].plot(q4, p4)
    axs[1, 1].set_title(r'$\delta$ = ' + str(deltas[3]))
    for ax in axs.flat:
        ax.set(xlabel='q(t)', ylabel='p(t)')
    # for ax in axs.flat:
    #     ax.label_outer()
    plt.tight_layout()
    plt.show()

# Phase portrait p(q)
def phase_portrait(delta):
    q1, p1 = forward_euler_approx(4.9, 20, 0, delta)
    q2, p2 = forward_euler_approx(5.0, 20, 0, delta)
    q3, p3 = forward_euler_approx(5.1, 20, 0, delta)
    q4, p4 = forward_euler_approx(5.2, 20, 0, delta)
    plt.plot(q1, p1, color='r', label = '4.9')
    plt.plot(q2, p2, color='g', label = '5.0')
    plt.plot(q3, p3, color='b', label = '5.1')
    plt.plot(q4, p4, color='k', label = '5.2')
    plt.xlabel('q(t)')
    plt.ylabel('p(t)')
    plt.title('Forward Euler Approximation - Phase portrait')
    plt.legend()
    plt.show()

# Total error calculation
def total_error(step_sizes: list):
    actual = -5*np.sin(np.pi/2)
    print(f'actual: {actual}')
    errors = []
    for delta in step_sizes:
        q, p = forward_euler_approx(5, 20, 0, delta)
        q_final = q[-1]
        print(f'delta: {delta}, q_final: {q_final}')
        total_error = q_final-actual
        errors.append(total_error)
    return step_sizes, errors

#Total error plot
def total_error_plot(step_sizes: list):
    x, y = total_error(step_sizes)
    plt.plot(x, y, 'b', label = 'Total Error')
    x_quad = np.linspace(0,0.2,100)
    y_quad = x_quad**2
    plt.plot(x_quad, y_quad, 'r', linestyle='dashed', label = r'$+x^2$')
    plt.plot(x_quad, -y_quad, 'r', linestyle='dotted', label = r'$-x^2$')
    plt.xlabel(r'$\delta$')
    plt.ylabel('Total Error')
    plt.legend()
    plt.title(r'Graph of total error vs step size $\delta$')
    plt.show()

# Stanalone function for comparing numerical approx at step-size
# 0.1 with the analytical solution. Mainly used to make sure
# I have the correct start and end points.
def comparison_plot():
    q, p = forward_euler_approx(5, 20, 0, deltas[2])
    t = np.linspace(0, int(q_end-q_start), int((q_end-q_start)//deltas[2]))
    print(f't: {t}')
    plt.plot(t, q)
    x = np.arange(0,15,deltas[2])   # start,stop,step
    y = np.cos(x)
    plt.plot(x,y,'r')
    plt.show()

# Comment these in or out depending on which plot you want to see
# multiplot_qt()
# multiplot_pq()
# phase_portrait(deltas[2])
total_error_plot(np.arange(0.01, 0.2, 0.002))
# comparison_plot()
