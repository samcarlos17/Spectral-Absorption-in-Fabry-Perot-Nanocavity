import numpy as np
import matplotlib.pyplot as plt
import meep as meep

# Physical constants
c = 3e8  # Speed of light in m/s

# --- Drude parameters in Meep units (1/μm) ---
eps_inf = 12
f_p = 0.1167   # Meep units: 1/μm
gamma = 0.0011  # Meep units: 1/μm

# Convert to angular frequency (omega = 2πf)
omega_p = 2 * np.pi * f_p
gamma_w = 2 * np.pi * gamma

# Frequency range (Meep units, 1/μm)
wavelengths=np.linspace(0.45,0.65,50)
freqs = 1/wavelengths  # avoid 0 to prevent div-by-zero
omegas = 2 * np.pi * freqs

# Compute complex permittivity ε(ω)
eps_real = []
eps_imag = []

for omega in omegas:
    eps = eps_inf - (omega_p**2) / (omega**2 + 1j * gamma_w * omega)
    eps_real.append(np.real(eps))
    eps_imag.append(np.imag(eps))

eps_real = np.array(eps_real)
eps_imag = np.array(eps_imag)

# --- Plotting ---
plt.figure(figsize=(10, 5))

plt.plot(wavelengths, eps_real, label='Re(ε)', color='blue')
plt.plot(wavelengths, eps_imag, label='Im(ε)', color='red')
plt.xlim(0.45, 0.65)  # Adjust as needed for your spectrum
plt.xlabel('Wavelength (μm)')
plt.ylabel('Permittivity ε(ω)')
plt.title('Drude Model ε(ω) for InAs Material')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.gca().invert_xaxis()  # So longer wavelengths (low freq) are on the left
plt.show()
