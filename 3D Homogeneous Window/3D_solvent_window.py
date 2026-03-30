import numpy as np
import plotly.graph_objects as go

# =============================================================================
# 1. PHYSICAL & CHEMICAL PARAMETERS (Hansen Solubility Parameters)
# =============================================================================
# R: Universal gas constant [J/(mol·K)], T: Temperature [K]
R = 8.314
T = 298.15
RT = R * T

# Hansen Solubility Parameters (HSP) for components: 
# V: Molar Volume, d: Dispersion, p: Polar, h: Hydrogen bonding
# Component X (e.g., Water), Y (e.g., Ethanol), Z (e.g., Ethyl Acetate)
V_x, d_x, p_x, h_x = 18.1, 15.5, 16.0, 42.3
V_y, d_y, p_y, h_y = 58.5, 15.5, 8.8, 19.4
V_z, d_z, p_z, h_z = 98.5, 15.8, 5.3, 7.2

# Range for Acetalization Degree (AD) - Vertical axis of the prism
AD_min = 20
AD_max = 100

# =============================================================================
# 2. GRID GENERATION & COORDINATE TRANSFORMATION
# =============================================================================
grid_res = 120
# Define Cartesian grid for the 3D space
X_cart = np.linspace(0, 1, grid_res)
Y_cart = np.linspace(0, np.sqrt(3)/2, int(grid_res * np.sqrt(3)/2))
Z_AD = np.linspace(AD_min, AD_max, grid_res)

X, Y, Z = np.meshgrid(X_cart, Y_cart, Z_AD, indexing='ij')

# Transformation: Cartesian (X, Y) to Barycentric/Ternary coordinates (x_frac, y_frac, z_frac)
# This maps the rectangular grid onto a triangular base
z_frac = Y / (np.sqrt(3) / 2)
y_frac = X - (np.sqrt(3) / 3) * Y
x_frac = 1.0 - y_frac - z_frac

# =============================================================================
# 3. THERMODYNAMIC CALCULATIONS (Flory-Huggins Chi Parameter)
# =============================================================================
# Signed Distance Function (SDF) to identify points outside the ternary triangle
sd = np.maximum(-x_frac, np.maximum(-y_frac, -z_frac))

# Clipping and normalizing fractions to ensure valid mixture properties
x_c = np.clip(x_frac, 0, 1)
y_c = np.clip(y_frac, 0, 1)
z_c = np.clip(z_frac, 0, 1)
sum_c = x_c + y_c + z_c
x_c /= sum_c
y_c /= sum_c
z_c /= sum_c

# Calculating mixture properties based on volume fractions (phi)
V_mix = V_x*x_c + V_y*y_c + V_z*z_c
phi_x = (V_x * x_c) / V_mix
phi_y = (V_y * y_c) / V_mix
phi_z = (V_z * z_c) / V_mix

# Linear mixing rule for HSP of the solvent blend
d_mix = d_x*phi_x + d_y*phi_y + d_z*phi_z
p_mix = p_x*phi_x + p_y*phi_y + p_z*phi_z
h_mix = h_x*phi_x + h_y*phi_y + h_z*phi_z

# Polymer HSP as a function of Acetalization Degree (AD) - Empirical regression
poly_d = -0.02112 * Z + 18.82287
poly_p = -0.09456 * Z + 15.21316
poly_h = -0.11071 * Z + 20.24478

# Calculate Hansen distance squared and the interaction parameter (Chi)
# Distance formula: Ra^2 = 4(dd1-dd2)^2 + (dp1-dp2)^2 + (dh1-dh2)^2
distance_sq = 4*(d_mix - poly_d)**2 + (p_mix - poly_p)**2 + (h_mix - poly_h)**2
chi_base = 0.6 * (V_mix / RT) * distance_sq

# =============================================================================
# 4. BOUNDARY CLIPPING (SDF Barrier)
# =============================================================================
# Apply a steep penalty (barrier) to points outside the ternary boundaries 
# to ensure the isosurface is clipped cleanly within the prism.
barrier = 0.5 + 500.0 * sd
Chi_grid = np.maximum(chi_base, barrier)

# =============================================================================
# 5. 3D VISUALIZATION & RENDERING (Publication Quality)
# =============================================================================
fig = go.Figure()

# Isosurface rendering for the solubility volume
fig.add_trace(go.Isosurface(
    x=X.flatten(),
    y=Y.flatten(),
    z=Z.flatten(),
    value=Chi_grid.flatten(),
    isomin=0.0,
    isomax=0.5,           # Threshold for solubility (usually Chi < 0.5)
    surface_count=2,
    colorscale='Viridis',
    reversescale=False, 
    caps=dict(x_show=False, y_show=False, z_show=False),
    opacity=0.9,
  
    # Advanced lighting for enhanced 3D depth perception (important for SCI figures)
    lighting=dict(
        ambient=0.3,      # Shadows in non-illuminated areas
        diffuse=0.8,      # Smooth gradients on curved surfaces
        specular=0.5,     # Glossy highlights
        roughness=0.4,    # Matte-finish texture
        fresnel=0.2       # Rim lighting for better edge definition
    ),
    name='Solubility Volume',
  
    # Colorbar configuration
    colorbar=dict(
        title=dict(
            text="<b>Interaction Parameter (χ)</b>", 
            font=dict(family="Arial", size=22, color="black"), 
            side="top"
        ),
        x=0.75, y=0.4, len=0.75, thickness=20,  
        tickfont=dict(family="Arial", size=16, color="black"), 
        outlinewidth=1, outlinecolor="black"
    )
))

# --- Prism Wireframe Reconstruction ---
pt_water = [0, 0]
pt_ethanol = [1, 0]
pt_ester = [0.5, np.sqrt(3)/2]

# Vertical and base boundary lines
lines_x = [pt_water[0], pt_ethanol[0], pt_ester[0], pt_water[0], pt_water[0], pt_ethanol[0], pt_ethanol[0], pt_ester[0], pt_ester[0]]
lines_y = [pt_water[1], pt_ethanol[1], pt_ester[1], pt_water[1], pt_water[1], pt_ethanol[1], pt_ethanol[1], pt_ester[1], pt_ester[1]]
lines_z = [AD_min, AD_min, AD_min, AD_min, AD_max, AD_max, AD_min, AD_min, AD_max]

fig.add_trace(go.Scatter3d(
    x=lines_x, y=lines_y, z=lines_z,
    mode='lines',
    line=dict(color='black', width=4),
    name='Phase Boundaries',
    hoverinfo='skip'
))

# Top triangle boundary
fig.add_trace(go.Scatter3d(
    x=[pt_water[0], pt_ethanol[0], pt_ester[0], pt_water[0]],
    y=[pt_water[1], pt_ethanol[1], pt_ester[1], pt_water[1]],
    z=[AD_max, AD_max, AD_max, AD_max],
    mode='lines',
    line=dict(color='black', width=4),
    showlegend=False,
    hoverinfo='skip'
))

# --- Vertex Labels ---
for z_pos in [AD_min, AD_max]:
    fig.add_trace(go.Scatter3d(
        x=[pt_water[0]-0.06, pt_ethanol[0]+0.06, pt_ester[0]],
        y=[pt_water[1]-0.06, pt_ethanol[1]-0.06, pt_ester[1]+0.06],
        z=[z_pos, z_pos, z_pos],
        mode='text',
        text=['Water', 'Ethanol', 'Ethyl Acetate'],
        textfont=dict(family="Arial", size=18, color='black'),
        showlegend=False,
        hoverinfo='skip'
    ))

# --- Global Layout & Aesthetic Refinement ---
fig.update_layout(
    font=dict(family="Arial", color="black"),
    paper_bgcolor='white',
    plot_bgcolor='white',
    margin=dict(l=0, r=0, b=0, t=20),
  
    scene=dict(
        # Hide 2D grid axes for a cleaner ternary prism look
        xaxis=dict(visible=False, showbackground=False),
        yaxis=dict(visible=False, showbackground=False),
      
        # Customize Z-axis (Acetalization Degree)
        zaxis=dict(
            visible=True,
            showbackground=False, 
            showgrid=False,        
            zeroline=False,        
            showline=True,         
            linecolor='black',     
            linewidth=2,           
            ticks='outside',       
            tickcolor='black',     
            tickwidth=2,
            ticklen=6,
            title=dict(
                text='<b>Acetalization Degree (AD)</b>', 
                font=dict(family="Arial", size=20, color='black')
            ),
            tickfont=dict(family="Arial", size=16, color='black')
        ),
      
        aspectmode='manual',
        aspectratio=dict(x=1, y=np.sqrt(3)/2, z=1),
        camera=dict(eye=dict(x=1.5, y=1.5, z=0.5))
    )
)

fig.show()

# Export high-resolution image for publication (Requires 'kaleido' library)
# fig.write_image("Phase_Diagram_3D.png", width=1200, height=1000, scale=3)
