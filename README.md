# INVDFE - Finite Difference Inversion of the Eikonal Equation

using the 2D eikonal solver of Podvin & Lecomte (1991) and a modified version of the simplex downhill methof of the Numerical Recipies.

## Overview
INVDFE is a program for **finite difference inversion** of the Eikonal equation. The program has undergone several updates and improvements over time, as outlined below.
But its kernel is still 25 years old, still written in Fortran77 with all its limitations. 

## Version History

### Version 2.5.2 (20/07/09)
- The maximum dimensions are defined through `statement parameter` in `maxval.com`.
- Maximum dimensions for **stations and shots** increased to **1000**.

### Version 2.5.1 (15/07/04)
- Fixed an issue in `getprec` to allow the user to omit the seed when using `mode=start` and `mode=joint`.
- Now supports both `mode=start` and `mode=start,##` syntaxes.
- **(29/10/2006)** Fixed a **bug in dimension settings** (750 → 1600).

### Version 2.5.0 (19/01/04)
- Added **resolution mode**, which calculates the matrix of partial derivatives around a given model.
- Works like `mode=start`, with the reference model provided via the `apriori` command in the `precision` file.
- The **resolution calculation** introduced in version **2.2.0** is now **permanently removed**.

### Version 2.4.1 (03/12/03)
- Fixed a **bug in back ray tracing**.

### Version 2.4.0 (24/02/03)
- Added a **water option** to introduce a mask below the topography with constant velocity and variable thickness.

### Version 2.3.1 (23/11/00)
- Fixed various bugs.
- The `misfit` function is now set to **infinity** if:
  - The interface is **not fully enclosed** between the topography and the model bottom.
  - The **reflection impact point** is on one of the model’s borders.
- Fixed a **bug in `comp_matvit2`** affecting non-constant velocity joint inversion.

### Version 2.3.0 (01/06/00)
- Added **interface migration** within a known medium (**special case of `mode=joint`**).
- Activated when `nxt + nzt = 0`.
- Reference medium provided via `apriori` in `precision` file.
- Modified: `getdata`, `compmisfit`, `comp_matvit2`.

### Version 2.2.1 (20/03/00)
- Changed output format of `hwd.dat` file.
- Now generates **one file per shot** (`hwd??.dat`).
- Internal format changed from `offset, time` → `offset, relative depth to topography, time`.

### Version 2.2.0 (21/10/99)
- Introduced **resolution mode**, activated by a switch in the `precision` file.

### Version 2.1.1 (22/09/99)
- Bug fixes.

### Version 2.1.0 (20/09/99)
- **Decoupled control point discretization** between velocity media and interfaces.
- `nxt` and `nzt` are now **two-element arrays**.
- Number of control points for the interface specified by `ni` in `geometry.dat`.
- Introduced variable `nm = mode + 1` to simplify **testing between smooth and interface modes**.

### Version 2.0.1 (11/08/99)
- Optimized `comp_misfit` for **6× faster performance**.

### Version 2.0.0 (10/07/99)
- Introduced **joint inversion** with **reflected waves**.
- The medium is now discretized as **two velocity layers** (above/below interface) + the interface itself.
- New `mode` parameter options:
  - `mode=start(0)`: Inversion with a **smooth medium**.
  - `mode=joint(1)`: Inversion with an **interface**.
  - `mode=continue(2)`: Continue an inversion, **regardless of type**.
  - `mode=direct(3)`: Forward modeling in a **smooth medium**.
  - `mode=dirref(4)`: Forward modeling in a **medium with an interface**.
- **Major code overhaul!**

## Contributors
- André Herrero: Code development
- Luigi Improta, Fabio Villani: Testers

## Contact
For any issues or contributions, please submit a pull request or open an issue in this repository.
