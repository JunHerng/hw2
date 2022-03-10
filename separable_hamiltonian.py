'''
Description: Forward Euler Approximation code for PC5210 HW2
with seperable Hamiltonian.
'''

import numpy as np
from matplotlib import pyplot as plt

# Global parameters
omega = 1

# ==============================
#  Forward euler approximation
# ==============================
def separable_forward_euler_approx(q_initial, q_final, p_initial, delta):
    # Initalize parameters
    q0 = q_initial
    qf = q_final
    p0 = p_initial
    n = int((qf-q0)//delta*2) # Number of points, floor div and set to int
    
    # Initialize arrays
    q_array = np.linspace(q0, qf, n)
    p_array = np.zeros(n)
    p_array[0] = p0

    # Forward Euler Approx
    # First Half
    for i in range(1,n):
        if i%2 == 1: # Odd
            q_array[i] = delta/2*p_array[i-1] + q_array[i-1]
            p_array[i] = p_array[i-1]
    # Second Half
        elif i%2 == 0: # Even
            q_array[i] = q_array[i-1]
            p_array[i] = p_array[i-1] - delta/2*omega**2*q_array[i-1]
        else: # Debug
            print('Math error???')

    return(q_array, p_array)

# Data for plotting
deltas = [0.2, 0.1, 0.05, 0.01]
q_start = 5
q_end = 35
q1, p1 = separable_forward_euler_approx(q_start, q_end, 0, deltas[0])
q2, p2 = separable_forward_euler_approx(q_start, q_end, 0, deltas[1])
q3, p3 = separable_forward_euler_approx(q_start, q_end, 0, deltas[2])
q4, p4 = separable_forward_euler_approx(q_start, q_end, 0, deltas[3])

# Multiplot q(t)
def multiplot_qt():
    fig, axs = plt.subplots(2, 2)
    fig.suptitle(r'q(t) for various $\delta$')
    axs[0, 0].plot(np.linspace(0, int(q_end-q_start), int((q_end-q_start)//deltas[0])*2), q1)
    axs[0, 0].set_title(r'$\delta$ = ' + str(deltas[0]))
    axs[0, 1].plot(np.linspace(0, int(q_end-q_start), int((q_end-q_start)//deltas[1])*2), q2)
    axs[0, 1].set_title(r'$\delta$ = ' + str(deltas[1]))
    axs[1, 0].plot(np.linspace(0, int(q_end-q_start), int((q_end-q_start)//deltas[2])*2), q3)
    axs[1, 0].set_title(r'$\delta$ = ' + str(deltas[2]))
    axs[1, 1].plot(np.linspace(0, int(q_end-q_start), int((q_end-q_start)//deltas[3])*2), q4)
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

# Plotting p^2+q^2 as (scaled) energy
def plot_energy():
    # Plot 1
    q, p = separable_forward_euler_approx(5,35,0,0.01)
    q_plot = q[::2]
    p_plot = p[::2]
    t = np.linspace(0,15,np.size(q_plot))
    plt.plot(t, q_plot**2+p_plot**2, 'r', label='0.01')
    # Plot 2
    q1, p1 = separable_forward_euler_approx(5,35,0,0.02)
    q_plot1 = q1[::2]
    p_plot1 = p1[::2]
    t1 = np.linspace(0,15,np.size(q_plot1))
    plt.plot(t1, q_plot1**2+p_plot1**2, 'g', label='0.02')
    # Plot 3
    q2, p2 = separable_forward_euler_approx(5,35,0,0.03)
    q_plot2 = q2[::2]
    p_plot2 = p2[::2]
    t2 = np.linspace(0,15,np.size(q_plot2))
    plt.plot(t2, q_plot2**2+p_plot2**2, 'b', label='0.03')
    # Labels and Aesthetics
    plt.title(r'Plots of $(p^2+q^2)$ for various $\delta$')
    plt.legend()
    plt.show()

# Phase portrait p(q)
def phase_portrait(delta):
    q1, p1 = separable_forward_euler_approx(4.9, 20, 0, delta)
    q2, p2 = separable_forward_euler_approx(5.0, 20, 0, delta)
    q3, p3 = separable_forward_euler_approx(5.1, 20, 0, delta)
    q4, p4 = separable_forward_euler_approx(5.2, 20, 0, delta)
    plt.plot(q1, p1, color='r', label = '4.9')
    plt.plot(q2, p2, color='g', label = '5.0')
    plt.plot(q3, p3, color='b', label = '5.1')
    plt.plot(q4, p4, color='k', label = '5.2')
    plt.xlabel('q(t)')
    plt.ylabel('p(t)')
    plt.title('Forward Euler Approximation - Phase portrait')
    plt.legend()
    plt.show()

# multiplot_qt()
# multiplot_pq()
# plot_energy()
# phase_portrait(deltas[1])
# phase_portrait(deltas[2])
phase_portrait(deltas[3])