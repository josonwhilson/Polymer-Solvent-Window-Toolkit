import numpy as np
import pandas as pd
import plotly.express as px

# --- 1. Define basic parameters ---
R = 8.314
T = 298.15
RT = R * T

# Solvent parameters: x=Water, y=Ethanol, z=Ethyl Acetate
V_x, d_x, p_x, h_x = 18.1, 15.5, 16.0, 42.3
V_y, d_y, p_y, h_y = 58.5, 15.8, 8.8, 19.4
V_z, d_z, p_z, h_z = 98.5, 15.8,  5.3,  7.2

# --- 2. New: Function to calculate PVB solubility parameters based on Degree of Acetalization (AD) ---
def get_pvb_hsp(AD):
    """Calculate d, p, h values of PVB based on degree of acetalization"""
    poly_d = -0.02112 * AD + 18.82287
    poly_p = -0.09456 * AD + 15.21316
    poly_h = -0.11071 * AD + 20.24478
    return poly_d, poly_p, poly_h

# --- 3. Encapsulation: Core function to find matching solvents ---
def find_solvents_for_pvb(AD, step=0.05, chi_threshold=0.5):
    """
    Given degree of acetalization AD, find all Water/Ethanol/Ethyl Acetate ratios satisfying Chi <= threshold
    """
    # 3.1 Get polymer parameters for current AD
    poly_d, poly_p, poly_h = get_pvb_hsp(AD)
    
    # 3.2 Generate grid (vectorized calculation)
    x_vals = np.arange(0, 1 + step, step)
    y_vals = np.arange(0, 1 + step, step)
    X_grid, Y_grid = np.meshgrid(x_vals, y_vals)
    
    X_flat = X_grid.flatten()
    Y_flat = Y_grid.flatten()
    
    valid_mask = (X_flat + Y_flat) <= (1.0 + 1e-6)
    x = X_flat[valid_mask]
    y = Y_flat[valid_mask]
    z = 1.0 - x - y
    z[z < 0] = 0  # Correct extremely small negative numbers
    
    # 3.3 Substitute into formula for calculation
    V_mix = V_x*x + V_y*y + V_z*z
    
    phi_x = (V_x * x) / V_mix
    phi_y = (V_y * y) / V_mix
    phi_z = (V_z * z) / V_mix
    
    d_mix = d_x*phi_x + d_y*phi_y + d_z*phi_z
    p_mix = p_x*phi_x + p_y*phi_y + p_z*phi_z
    h_mix = h_x*phi_x + h_y*phi_y + h_z*phi_z
    
    # Calculate Hansen distance squared and Chi parameter
    distance_sq = 4*(d_mix - poly_d)**2 + (p_mix - poly_p)**2 + (h_mix - poly_h)**2
    chi = 0.6 * (V_mix / RT) * distance_sq
    
    # 3.4 Filter results
    solution_mask = chi <= chi_threshold
    
    results = pd.DataFrame({
        'x_Water': x[solution_mask].round(3),
        'y_Ethanol': y[solution_mask].round(3),
        'z_EthylAcetate': z[solution_mask].round(3),
        'Chi': chi[solution_mask].round(4)
    })
    
    return results, poly_d, poly_p, poly_h

# ==========================================
# --- 4. Actual usage and testing (🌟 SCI-level plot beautification 🌟) ---
# ==========================================

target_AD = 100
 # <-- You can modify this value at any time
# Note: 0.0051 generates a large number of points; overlapping scatter points can create an area-like effect
step_size = 0.0021 # Extract step size as a variable for easier area calculation

# Call function
df_solutions, current_d, current_p, current_h = find_solvents_for_pvb(target_AD, step=step_size)

print(f"[Current Setting] PVB Degree of Acetalization AD = {target_AD}")
print(f"[Polymer HSP] d={current_d:.2f}, p={current_p:.2f}, h={current_h:.2f}")
print(f"[Calculation Result] Found {len(df_solutions)} ratio combinations meeting the condition (Chi<=0.5).\n")
# ---------------------------------------------------------
# 📐 New: Solubility window area calculation module
# ---------------------------------------------------------
if len(df_solutions) > 0:
    # 1. Calculate the number of soluble points
    count_soluble = len(df_solutions)
    
    # 2. Calculate the representative area of a single grid point (dx * dy)
    area_per_point = step_size * step_size
    
    # 3. Calculate the absolute area of the solubility region on the (x, y) projection plane
    area_soluble = count_soluble * area_per_point
    
    # 4. Calculate the total phase diagram area (x + y <= 1 is a right triangle with base 1 and height 1, area is 0.5)
    area_total = 0.5
    
    # 5. Calculate percentage
    area_ratio = (area_soluble / area_total) * 100
    
    print(f"[Area Calculation] Step size: {step_size}")
    print(f"[Area Calculation] Soluble points: {count_soluble}")
    print(f"[Area Calculation] Solubility window ratio: {area_ratio:.2f}% (Relative to the entire phase diagram)")
    
    # ---------------------------------------------------------
    # 🎨 Plotting section (Add area to the title)
    # ---------------------------------------------------------
    print("Top 5 feasible ratios (x: Water, y: Ethanol, z: Ethyl Acetate):")
    print(df_solutions.head(5))
    
    fig = px.scatter_ternary(
        df_solutions, 
        a="x_Water", b="y_Ethanol", c="z_EthylAcetate",
        color="Chi", 
        color_continuous_scale="Viridis",
        labels={
            "x_Water": "<b>Water</b>",
            "y_Ethanol": "<b>Ethanol</b>",
            "z_EthylAcetate": "<b>Ethyl Acetate</b>",
            "Chi": "<b>Interaction Parameter (<i>&chi;</i>)</b>"
        }
    )
    
    fig.update_traces(marker=dict(size=4, opacity=1.0, line=dict(width=0)))

    
    # 🌟 3. Colorbar fine customization 🌟
    fig.update_coloraxes(
        colorbar=dict(
            title=dict(
                text="<b>Chi Parameter (χ)</b>", 
                font=dict(family="Arial", size=18, color="black"),
                side="top"
            ),
            thickness=40,    # Thickness parameter
            len=0.75,        # Length parameter
            x=0.80,          # X position slightly to the right to avoid blocking the right corner
            y=0.5,           # Center symmetrically
            tickfont=dict(family="Arial", size=14, color="black"), 
            outlinewidth=1,  # Add a black outline to the colorbar
            outlinecolor="black"
        )
    )

    # 🌟 4. Global font and ternary axis refinement settings 🌟
    # Extract common axis configuration (black outline, black ticks, light internal auxiliary grid)
    axis_config = dict(
        title=dict(font=dict(family="Arial", size=20, color="black")),
        tickfont=dict(family="Arial", size=24, color="black"),
        showgrid=True, gridcolor="rgba(211,211,211,0.5)", gridwidth=1, # Light gray auxiliary grid
        linecolor="black", linewidth=2, showline=True,  # Bold black outline
        ticks="outside", ticklen=6, tickwidth=2, tickcolor="black" # Professional outer ticks
    )

    fig.update_layout(
        font=dict(family="Arial", color="black"), # Global font settings
        paper_bgcolor='white', # Outer pure white background
        plot_bgcolor='white',  # Inner pure white background
        
        # SCI figure title setting (optional, can also add Figure caption directly in LaTeX/Word)
        title=dict(
             
            font=dict(family="Arial", size=22, color="black"),
            x=0.5, xanchor='center', y=0.95
        ),
        
        # Deeply rewrite ternary coordinate system styles
        ternary=dict(
            bgcolor="rgba(0,0,0,0)", # Force the central area to be fully transparent
            aaxis=axis_config,
            baxis=axis_config,
            caxis=axis_config
        ),
        margin=dict(l=40, r=40, b=80, t=80) # Optimize margins
    )
    
    fig.show()
    # Recommended saving method:
    # fig.write_image("SCI_Ternary_Phase_Diagram.png", scale=3)
else:
    print("No mixed solvent ratio meeting the condition was found at this degree of acetalization.")
