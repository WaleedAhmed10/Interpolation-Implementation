# Numerical Interpolation

This suite comprises three interactive Python graphical user interface (GUI) applications built using **Tkinter**. Each lab demonstrates a specific set of **Numerical Analysis and Interpolation methods** applied to equally spaced data points, featuring difference table generation, mathematical verification, multi-formula evaluation, and Matplotlib-based geometric visualization.

---

## Suite Overview

| Script Name | Main Focus | Implemented Methods | Features |
| :--- | :--- | :--- | :--- |
| **`Lab #5.py`** | Polynomial Interpolation | Newton Forward, Newton Backward | Table formatting, closed-form polynomial expansion ($P(x)$), responsive UI |
| **`Lab #6.py`** | Central Differences | Stirling, Bessel, Everett, Gauss Forward, Gauss Backward | Central difference grid layout, multi-formula comparative table, error calculation |
| **`Lab #7.py`** | Applied Geometry / Visualization | Hybrid Engine (Newton Forward + Stirling + Newton Backward) | Face outline reconstruction, coordinate mirroring ($x \to -x$), Matplotlib plot export |

---

## Detailed Lab Breakdown

### 1. Lab #5: Newton Forward & Backward Interpolation (`Lab #5.py`)

#### Overview
Lab #5 provides an interactive calculator and difference table engine for **Newton's Forward and Backward Difference Interpolation**. It validates input spacing, builds higher-order forward difference matrices, generates the exact polynomial representation of degree $n-1$, and handles target value queries.

#### Key Functions & Logic
* `difference_table(y)`: Computes the full forward difference table up to the maximum order available for the dataset.
* `interpolate(x, table, query)`: Automatically determines whether to apply Newton Forward (lower half of range) or Newton Backward (upper half of range) based on the target $x$-location.
* `generate_polynomial(x, table)`: Converts tabular difference coefficients into explicit standard polynomial coefficients $a_n x^n + \dots + a_1 x + a_0$.
* `ScrollableFrame`: Custom Tkinter canvas-backed frame that handles dynamic resizing and scroll wheel behavior across operating systems.

---

### 2. Lab #6: Central Difference Interpolation (`Lab #6.py`)

#### Overview
Lab #6 focuses on **Central Difference Interpolation** techniques. It constructs a central difference arrangement and presents side-by-side evaluation of all major central interpolation formulas for any valid query point.

#### Implemented Formulae
1. **Stirling's Formula:** Optimal when $|p| \le 0.25$.
2. **Bessel's Formula:** Optimal when $p \approx 0.5$.
3. **Everett's Formula:** Expressed using $p$ and $q = 1 - p$.
4. **Gauss Forward Formula:** Uses forward central steps.
5. **Gauss Backward Formula:** Uses backward central steps.

#### Key Features
* **Central Difference Table Visualizer:** Renders a clean visual central grid ($\delta^k f$) using staggered line spacing.
* **Comparative Multi-Formula Output:** Computes and displays $x_0$, $p$, $f(x)$, absolute error (when ground truth is provided), and accuracy alerts for all five methods simultaneously.
* **Statistical Summary:** Evaluates the mean and spread ($\max - \min$) across all five methods to measure model stability.

---

### 3. Lab #7: Symmetrical Face Outline Reconstruction (`Lab #7.py`)

#### Overview
Lab #7 applies numerical interpolation to a practical 2D geometric reconstruction problem: generating a smooth, symmetrical human face contour from sparse half-outline data points $(y, x)$.

#### Architecture & Features
* **Hybrid `Interpolator` Class:**
  * **Newton Forward:** Used for the bottom region of the outline ($y < \text{low threshold}$).
  * **Stirling Central:** Used for the middle section of the face ($y$ between low and high thresholds).
  * **Newton Backward:** Used for the upper forehead region ($y > \text{high threshold}$).
* **3-Step Reconstruction Process:**
  1. **Reconstruct Left Half:** Computes dense interpolated coordinates along the input curve.
  2. **Mirror Right Half:** Applies spatial reflection ($x \to -x$) across the vertical midline.
  3. **Plot Full Outline:** Fills and renders the complete face contour using Matplotlib.
* **Plot Export:** Built-in dialog to export generated figures as high-resolution PNG images.

---

## Installation & Setup

### Prerequisites
* **Python Version:** Python 3.8 or higher.
* Standard Python modules used: `tkinter`, `math`, `re`.

### External Dependencies
`Lab #7.py` requires external libraries for numerical array manipulations and plot rendering. Install them via `pip`:

```bash
pip install numpy matplotlib
```

---

## Usage

To launch any module, execute the desired Python file from your terminal:

```bash
# Run Lab 5 (Newton Forward/Backward & Polynomial Expansion)
python "Lab #5.py"

# Run Lab 6 (Central Difference Comparative Analysis)
python "Lab #6.py"

# Run Lab 7 (Symmetrical Face Outline Reconstruction)
python "Lab #7.py"
```