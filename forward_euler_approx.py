'''
Description: Forward Euler Approximation code for PC5210 HW2
'''

import numpy as np
from matplotlib import pyplot as plt

# Global parameters
omega = 1

# Forward euler approximation function
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
fig, axs = plt.subplots(2, 2)
fig.suptitle(r'q(t) for various $\delta$')
axs[0, 0].plot(np.linspace(1, 100, int((q_end-q_start)//deltas[0])), q1)
axs[0, 0].set_title(r'$\delta$ = ' + str(deltas[0]))
axs[0, 1].plot(np.linspace(1, 100, int((q_end-q_start)//deltas[1])), q2)
axs[0, 1].set_title(r'$\delta$ = ' + str(deltas[1]))
axs[1, 0].plot(np.linspace(1, 100, int((q_end-q_start)//deltas[2])), q3)
axs[1, 0].set_title(r'$\delta$ = ' + str(deltas[2]))
axs[1, 1].plot(np.linspace(1, 100, int((q_end-q_start)//deltas[3])), q4)
axs[1, 1].set_title(r'$\delta$ = ' + str(deltas[3]))
for ax in axs.flat:
    ax.set(xlabel='t', ylabel='q(t)')
# for ax in axs.flat:
#     ax.label_outer()
plt.tight_layout()
plt.show()

# Multiplot p(q)
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