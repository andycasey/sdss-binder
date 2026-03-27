# CubeGen Notebook — Data Cube Reconstruction Pipeline

This notebook provides an example of how to use 3DCubeGen within a notebook environment using BinderHub and provides a complete workflow to:

- Generate IFU data cubes  
- Produce emission-line maps  
- Create RGB composite images  
- Apply 2D deconvolution  
- Generate final science-quality products  

This workflow is designed to work with **LVM DRP outputs** and the **CubeGen / 3DCubeGen** package.

---

# 📚 Documentation

For full documentation of **3DCubeGen**, please visit:

https://hjibarram.github.io/3DCubeGen/

GitHub Repository:

https://github.com/hjibarram/3DCubeGen

---

# Requirements

```bash
git clone https://github.com/hjibarram/3DCubeGen.git
cd 3DCubeGen
git switch dev
pip install .
```

You also need:

- LVM DRP outputs

---

# Imports

```python
from CubeGen.tools.cubegen import map_ifu
from CubeGen.tools.mapgen import gen_map
from CubeGen.tools.images import get_jpg
from CubeGen.tools.deconvolve import deconvolve_2dfile
import CubeGen.tools.tools as tools
```

---

# Step 1 — Generate Data Cube

The `map_ifu` function is used to reconstruct the IFU data cube.

## Main Parameters

| Parameter | Description |
|-----------|-------------|
| redux_dir | Directory containing reduced data |
| redux_ver | Reduction version |
| out_path | Output directory |
| sigm_s | Spatial smoothing kernel |
| pix_s | Pixel size |
| alph_s | Kernel shape parameter |
| flu16 | Output flux normalization |
| use_slitmap | Use slitmap astrometry |
| pbars | Show progress bars |
| spec_range | Spectral range (optional) |
| nameF | Output cube name |

Example:

```python
map_ifu(...)
```

Output:

```
lvmCube-*.fits
```

---

# Step 2 — Generate Emission Line Maps

he `gen_map` function generates emission-line maps from the reconstructed cube.

## Available Emission Lines

| Option | Emission Line |
|--------|---------------|
| 1 | [OIII] |
| 2 | Hα |
| 3 | [SII] |
| 4 | Continuum |
| 5 | Custom spectral range |

Example:


```python
gen_map(...)
```

Creates:

- OIII maps
- Hα maps
- SII maps

Output:

```
lvmMap-*.fits
```


---

# Step 3 — Create RGB Composite Image

Create color composite images:

```python
get_jpg(...)
```

Output:

```
lvmMap-*_OHS.jpeg
```

---

# Step 4 — Apply Deconvolution

Apply PSF deconvolution (optional):

```python
deconvolve_2dfile(...)
```

---

# Step 5 — Final RGB

```python
get_jpg(...)
```

---

# Workflow

Raw LVM DRP  
↓  
map_ifu  
↓  
gen_map  
↓  
get_jpg  
↓  
deconvolve_2dfile  
↓  
get_jpg  

---

# Output Files

- lvmCube-*.fits  
- lvmMap-*_OIII.fits  
- lvmMap-*_HI.fits  
- lvmMap-*_SII.fits  
- lvmMap-*_OHS.jpeg  
- lvmMap-*_decv.fits  
- lvmMap-*_OHS_decv.jpeg  

---

# Documentation

https://hjibarram.github.io/3DCubeGen/

---

# GitHub

https://github.com/hjibarram/3DCubeGen
