import numpy as np
import scipy as sp
import scipy.integrate
import matplotlib.pyplot as plt

N = 100  # Number of spatial points
width_fuel = 50.0 # Width of the fuel region in cm
width_reflector = 20.0 # Width of the reflector region in cm
total_width = width_fuel + 2 * width_reflector # The absolute total width of your simulated reactor universe.
dx = total_width / N  # Spatial step size
# define geometry and grid parameters
d_fuel = 0.35 # Diffusion coefficient in the fuel region (cm)
sigma_a_fuel = 0.065 # Absorption cross-section in the fuel region (cm^-1)
sigma_f_fuel = 0.052 # Fission cross-section in the fuel region (cm^-1)
# reflector properties
d_reflector = 0.16 # Diffusion coefficient in the reflector region (cm)
sigma_a_reflector = 0.012 # Absorption cross-section in the reflector region (cm^-1)
sigma_f_reflector = 0.0 # Fission cross-section in the reflector region (cm^-1)

nu = 2.43 # Average number of neutrons produced per fission
w_per_fission = 3.2e-11 # Energy released per fission in Joules
# Create a "Blueprint" array that represents what material is present at each point across the reactor. 
x = np.linspace(0, total_width, N) # Create a spatial grid from 0 to total_width with N points.
D = np.zeros(N)
sigma_a = np.zeros(N)
sigma_f = np.zeros(N)
fuel_start = width_reflector # Starting point of the fuel region
fuel_end = width_reflector + width_fuel # Ending point of the fuel region

for i in range(N):
    if fuel_start <= x[i] <= fuel_end:
        # If the point is within the fuel region, assign fuel properties
        D[i] = d_fuel
        sigma_a[i] = sigma_a_fuel
        sigma_f[i] = sigma_f_fuel
    else:
        # If the point is within the reflector region, assign reflector properties
        D[i] = d_reflector
        sigma_a[i] = sigma_a_reflector
        sigma_f[i] = sigma_f_reflector

# size of internal system 
M = N - 2
A = np.zeros((M, M))
B = np.zeros((M, M))
for i in range(M):
    # map the internal points to the global grid 
    #Matrix A represents the diffusion and absorption temrs.
    g = i + 1
    A[i, i] = (2 *  D[g] / dx**2) + sigma_a[g] # Diagonal term includes diffusion and absorption
    if i > 0:
        A[i, i - 1] = -D[g] / dx**2
    if i < M - 1:
        A[i, i + 1] = -D[g] / dx**2

    # Matrix B represents the fission term.
    for i in range(M):
        g = i + 1
        B[i, i] = nu * sigma_f[g] # Fission term only contributes on the diagonal

# Solve the eigenvalue problem A * phi = (1/k) * B * phi
import scipy.linalg
eigenvalues, eigenvectors = sp.linalg.eig(A, B)
real_eigenvalues = np.real(eigenvalues) # Extract the real part of the eigenvalue
physical_mask = real_eigenvalues > 1e-10 # Filter out non-physical eigenvalues
filtered_lambdas = real_eigenvalues[physical_mask] # Keep only the physical eigenvalues
filtered_vectors = np.real(eigenvectors[:, physical_mask])
k_modes = 1.0 / filtered_lambdas
fundamental_index = np.argmax(k_modes)
k_eff = k_modes[fundamental_index]
phi_internal = filtered_vectors[:, fundamental_index]

# Enforce a physical rule: if the math flipped the graph upside down, flip it back!
if np.max(phi_internal) < 0:
    phi_internal = -phi_internal

print(f"Effective Multiplication Factor (k_eff): {k_eff:.4f}")

import matplotlib.pyplot as plt
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
ax1.plot(x, sigma_a, color = 'red', linewidth = 2, label='Absorption Cross-Section')
ax1.plot(x, sigma_f, color = 'blue', linewidth = 2, label='Fission Cross-Section')
ax1.set_xlabel('Position (cm)')
ax1.set_ylabel('Cross-Section (cm^-1)')
ax1.set_title('Cross-Sections Across the Reactor')
ax1.legend()
ax1.grid()

# Plot the neutron flux distribution
phi_full = np.zeros(N)
phi_full[1:-1] = phi_internal
ax2.plot(x, phi_full, color = 'green', linewidth = 2, label ='Neutron Flux')
ax2.set_xlabel('Position (cm)')
ax2.set_ylabel('Neutron Flux (arbitrary units)')
ax2.set_title('Neutron Flux Distribution')
ax2.legend()
ax2.grid()
plt.tight_layout()
plt.show()
