# Polymer Solubility Toolkit: Solvent Window Mapping & Visualization

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Plotly](https://img.shields.io/badge/Visualization-Plotly-orange.svg)](https://plotly.com/)
[![Pandas](https://img.shields.io/badge/Data-Pandas-green.svg)](https://pandas.pydata.org/)

This repository provides a comprehensive computational suite for identifying and visualizing the **Solvent Window** in polymer systems, specifically optimized for the homogeneous synthesis of **Polyvinyl Butyral (PVB)**. 

By integrating **Hansen Solubility Parameters (HSP)** and **Flory-Huggins Theory**, this toolkit allows researchers to map the solubility boundaries across various chemical states (Acetalization Degree) and solvent compositions.

---

## 🌟 Key Components

The toolkit consists of two primary modules designed for different analytical scales:

### 1. 3D Global Evolution Mapper (`3D_solvent_window.py`)
*   **Purpose**: Visualizes the "Solubility Volume" across the entire range of Acetalization Degree (AD 20–100).
*   **Feature**: Uses **Signed Distance Functions (SDF)** to render a smooth, continuous prism volume.
*   **Benefit**: Ideal for observing how the solvent window shifts and shrinks as the polymer structure evolves.

### 2. 2D Quantitative Traversal Tool (`2D_ternary_traversal.py`)
*   **Purpose**: Performs an **exhaustive grid search** at a user-defined AD value.
*   **Feature**: Calculates the **Solubility Window Area Ratio (%)** relative to the entire ternary phase space.
*   **Aesthetics**: Generates SCI-level ternary plots with high-density scatter mapping to simulate area-based regions.
*   **Benefit**: Provides precise solvent ratios (Water/Ethanol/Ethyl Acetate) and quantitative metrics for process optimization.

---

## 🧪 Scientific Methodology

The core engine calculates the interaction parameter ($\chi$) between the polymer and the solvent mixture:

1.  **Mixing Rules**: Calculates the effective HSP ($d_D, d_P, d_H$) and molar volume ($V_{mix}$) of the ternary mixture.
2.  **Polymer State Modeling**: Dynamically updates the polymer's HSP based on its **Acetalization Degree (AD)** using empirically derived linear regression models.
3.  **Boundary Criterion**: Defines the homogeneous window at $\chi \leq 0.5$.
4.  **Area Quantization**: For 2D slices, the solubility window area is calculated as:
    $$\text{Area Ratio} = \frac{\sum \text{Soluble Grid Points} \times \text{Step}^2}{0.5 (\text{Total Phase Space Area})} \times 100\%$$

---

## 🚀 Getting Started

### Prerequisites
```bash
pip install numpy pandas plotly kaleido
```

### Running the Tools
**To generate the 3D Prism Volume:**
```bash
python 3D_solvent_window.py
```

**To perform 2D quantitative analysis at a specific AD:**
```bash
# Open the script and set your target_AD (e.g., 100)
python 2D_ternary_traversal.py
```

---

## 📊 Visual Gallery

| 3D Global Window | 2D Quantitative Slice |
| :---: | :---: |
| ![3D Preview](3D%20Homogeneous%20Window/img_example/Solubility%20Window%2020-100.png) | ![2D Preview](https://via.placeholder.com/400x300.png?text=2D+Ternary+Plot) |
| *Evolution of solubility from AD 20 to 100* | *Solvent window area and $\chi$ distribution* |

---

## 🛠 Transferability
While pre-configured for PVB (Water/Ethanol/Ethyl Acetate), this method is **highly transferable**. To adapt it to other systems:
1.  Update the **Solvent HSP constants** in the "Basic Parameters" section.
2.  Modify the **`get_pvb_hsp(AD)`** function to reflect the solubility parameters of your specific polymer.

## 📂 Project Structure
- `3D_solvent_window.py`: Continuous 3D volume rendering using Plotly Isosurface.
- `2D_ternary_traversal.py`: Exhaustive grid search, area calculation, and ternary scatter plotting.
- `requirements.txt`: Environment dependencies.

## 📝 Citation
If this methodology or toolkit assists in your research, please consider citing it as:
> *Methodology for Identifying and Visualizing Solvent Windows in Polymer Synthesis: A Generalizable Approach.*

## 📄 License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---
