# Lab 8 - Vaccine Cold-Chain Temperature Excursion

Numerical Methods, Open Ended Lab. One MATLAB script (`lab8.m`) that solves every task in the lab,
prints a full report in the Command Window and opens a tabbed window with tables, results and plots.

## What you need

- `lab8.m` (the only file)
- **MATLAB desktop app** (R2016b or newer) if you have it installed. No toolboxes are required.
- **MATLAB Online** (https://matlab.mathworks.com) if you do not have the desktop app. Sign in with a free MathWorks account.

Use the desktop app when it is available. Otherwise use MATLAB Online. The script runs the same in both.

## How to run

### Option A: MATLAB desktop app
1. Put `lab8.m` in a folder, for example `C:\NumLab`.
2. Open MATLAB and make that folder the current folder (Current Folder panel, or `cd C:\NumLab`).
3. Type `lab8` in the Command Window and press Enter. (Or open the file and click **Run**.)

### Option B: MATLAB Online
1. Go to https://matlab.mathworks.com and sign in.
2. Upload `lab8.m`: click **Upload** on the Home tab, or drag the file into the Current Folder panel.
3. Check that the Current Folder is the one holding `lab8.m` (usually MATLAB Drive).
4. Type `lab8` in the Command Window and press Enter. (Or open the file and click **Run**.)
5. If the "Cold chain lab" window does not appear, allow pop-ups for the MATLAB Online site in your browser.

Alternative for either option: **New Script**, paste the code, save it as exactly `lab8.m`, then press **Run**.

## What you get

The report prints in the Command Window and the same text is in the window tabs:

| Tab | Content |
|---|---|
| Difference table | Forward differences up to the 6th, spacing check (h = 2 h), degree supported |
| Estimates | T at 5.0, 6.4 and 11.0 h with formula, origin, p, degree, every term, second-method check, accuracy comment. Quick estimator box at the bottom for any t |
| Maximum | Interpolated maximum, degree used, no-larger-value check, spline cross-check |
| Crossings | Times T crosses 8.0 deg C, total time above the limit, inverse-interpolation check, release decision |
| Validation | T(6.4 h) by Stirling, Bessel, backward and natural cubic spline, spread, where they agree and disagree |
| Plots | Readings with all four curves, the 8.0 deg C limit, crossings, maximum, and the spread of the curves |

Maximize the window if the quick estimator at the bottom of the Estimates tab is cut off.

## Methods used

| Task | Method | Why |
|---|---|---|
| T(5.0 h) | Bessel, origin 4 h, p = 0.5, degree 5 | Middle of an interval: odd terms vanish |
| T(6.4 h) | Stirling, origin 6 h, p = 0.2, degree 6 | Near the centre node, symmetric nodes |
| T(11.0 h) | Newton backward, origin 12 h, p = -0.5, degree 6 | Near the end of the table |
| Maximum and crossings | Degree-6 polynomial through all 7 readings | Unique polynomial through the data |
| Second methods | Lagrange on the same nodes, natural cubic spline | Independent checks |

## Results (4 decimal places)

| Item | Value |
|---|---|
| T(5.0 h) | 8.7539 deg C |
| T(6.4 h) | 9.7283 deg C |
| T(11.0 h) | 8.5193 deg C |
| Maximum | 10.1052 deg C at t = 7.8329 h |
| Crossings of 8.0 deg C | 4.1150 h and 11.5798 h |
| Time above 8.0 deg C | 7.4648 h |
| Validation at 6.4 h | Stirling 9.7283, Bessel 9.7241, Backward 9.7258, Spline 9.7275 (spread 0.0041) |
| Decision | Do not release |

## Settings you can change (top of `lab8.m`, section 0)

- `MAX_ALLOWED_H = 0` : allowed time above 8.0 deg C. The lab gives no number, so 0 h is assumed
  (any excursion means no automatic release). Set it to your course or protocol value if different.
- `DEG_STIR`, `DEG_BESS`, `DEG_BACK` : degrees used in the validation at 6.4 h (defaults 6, 5, 4).
  The backward formula at degree 6 is identical to Stirling, which is why degree 4 is used.

## Notes

- Interpolated values are estimates, not measurements. The logger resolves 0.1 deg C, so the
  4th decimal is shown only because the lab asks for it.
- Values near t = 0 and t = 12 are the least reliable. The backward and Bessel curves extrapolate
  on 0-4 h, which is why the curves part company there.
- The crossing times are roots of the direct interpolant. The reversed-variable (inverse) result is
  only a check and differs by 0.0492 h and 0.2249 h.
- If the tabbed window cannot be built, the script shows a warning and draws a plain figure.
  The Command Window report still prints in full.

## Troubleshooting

| Problem | Fix |
|---|---|
| `Undefined function or variable 'lab8'` | The current folder is not the one holding `lab8.m`. Change folder or upload the file there |
| Error about functions in a script | Your MATLAB is older than R2016b. Use MATLAB Online |
| Window does not show in MATLAB Online | Allow pop-ups for the site, or look for the figure window behind the browser |
| Text looks cut off | Maximize the figure window |
