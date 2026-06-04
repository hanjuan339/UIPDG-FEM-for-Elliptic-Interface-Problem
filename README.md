# An Unfitted Interface Penalty DG–FE Method for Elliptic Interface Problems

This repository contains the MATLAB implementation to reproduce the numerical results presented in the manuscript:
> **An Unfitted Interface Penalty DG–FE Method for Elliptic Interface Problems**

## Authors
- **Authors**: Juan Han, Haijun Wu, and Yuanming Xiao

*If you find this code useful for your research, please cite our manuscript or this repository as a reference.*

---

## 1. Prerequisites & Environment

To run the code and reproduce the results, your system needs to meet the following requirements:

- **MATLAB**: R2018a or later is recommended.
- **Penalty Parameter**: The default penalty parameter is set to **100** across the main numerical experiments to ensure stability.
- **C Compiler**: A working C compiler (e.g., Xcode/Clang, GCC, or MinGW) is required *only if* you need to recompile the MEX components.

---

## 2. Compilation & Troubleshooting (MEX Files)

This repository already includes pre-compiled binaries for the C-based MEX acceleration function (`aij`). However, due to differences in operating systems, hardware architectures, or MATLAB versions, these pre-compiled files might not run out of the box on your machine.

**If you encounter an error stating that the MEX file cannot be executed, please recompile the source file yourself in the MATLAB Command Window:**

```matlab
mex aij.c
```

## 3. Code Structure & Description

Here is the structure of the repository and the description of the core files:

```text
├── ex1a.m                  # Main execution script (corresponds to Fig. 11-14)
├── ex1e.m                  # Main execution script (corresponds to Fig. 15)
├── exifun.m                # Input data configuration
├── uipdgfem.m              # Core solver implementation
└── aij.c                   # C source code for fast assembly
