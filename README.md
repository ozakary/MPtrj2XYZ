# MPtrj2XYZ

---

📄 Author: **Ouail Zakary**  
- 📧 Email: [Ouail.Zakary@oulu.fi](mailto:Ouail.Zakary@oulu.fi)  
- 🔗 ORCID: [0000-0002-7793-3306](https://orcid.org/0000-0002-7793-3306)  
- 🌐 Website: [Personal Webpage](https://cc.oulu.fi/~nmrwww/members/Ouail_Zakary.html)  
- 📁 Portfolio: [GitHub Portfolio](https://ozakary.github.io/)

---

Fast and memory-efficient converter for the Materials Project Trajectory (MPtrj) dataset to extended XYZ format with element-based filtering.

## Overview

MPtrj2XYZ extracts structures from the large MPtrj dataset (12+ GB JSON file containing 1.5M+ structures with their corresponding DFT data) and converts them to extended XYZ format with forces, energies, and stress tensors. The tool supports element-based filtering with three strategies and uses streaming JSON parsing for minimal memory footprint.

## Features

- **Memory Efficient**: Processes 12+ GB JSON files using streaming parser (constant ~500 MB memory usage)
- **Fast**: Processes 4,000-6,000 structures per second
- **Element Filtering**: Three filtering strategies (any, only, all) or convert entire dataset
- **Complete Data Preservation**: Extracts lattice, positions, forces, energies, and stress tensors
- **Extended XYZ Format**: Output compatible with ASE and standard molecular dynamics tools
- **Progress Tracking**: Real-time progress bar with statistics

## Installation

```bash
pip install ase ijson tqdm numpy
```

## Quick Start

```bash
# Extract structures containing only H, C, O, S, Xe
python extract_mptrj_structures.py MPtrj_2022.9_full.json -e H C O S Xe -s only

# Convert entire dataset to XYZ (no filtering)
python extract_mptrj_structures.py MPtrj_2022.9_full.json -e all

# Extract with multiple strategies
python extract_mptrj_structures.py MPtrj_2022.9_full.json -e H C N O -s any only all
```

## Usage

```
python extract_mptrj_structures.py [JSON_FILE] [OPTIONS]

Required:
  JSON_FILE              Path to MPtrj JSON file

Options:
  -e, --elements         Target elements (default: H C N O F S Cl Xe)
                        Use "all" to convert entire dataset
  -s, --strategy         Filter strategies: any, only, all (default: all three)
  -o, --output           Output directory (default: mptrj_extracted)
  -c, --chunk-size       Progress update frequency (default: 1000)
```

## Filtering Strategies

### `any` - At least one target element
Structures containing at least one of the specified elements (can have others).

**Example**: Target = {H, C, O}
- ✅ H₂O (contains H and O)
- ✅ CH₄ (contains C and H)  
- ✅ NaCl·H₂O (contains H and O, plus Na and Cl)
- ❌ NaCl (no target elements)

### `only` - Only target elements
Structures containing only the specified elements (no other elements allowed).

**Example**: Target = {H, C, O}
- ✅ H₂O (only H and O)
- ✅ CH₄ (only C and H)
- ✅ C₆H₁₂O₆ (only C, H, and O)
- ❌ NaCl·H₂O (contains Na and Cl)

### `all` - All target elements present
Structures must contain all specified elements (can have others).

**Example**: Target = {H, C, O}
- ✅ C₆H₁₂O₆ (has all three: C, H, O)
- ✅ NaHCO₃ (has all three, plus Na)
- ❌ H₂O (missing C)
- ❌ CH₄ (missing O)

## Output Format

Extended XYZ files with the following structure:

```
[number_of_atoms]
Lattice="..." Properties=species:S:1:pos:R:3:forces:R:3 stress="..." free_energy=... energy=... pbc="T T T"
Element1  x1  y1  z1  fx1  fy1  fz1
Element2  x2  y2  z2  fx2  fy2  fz2
...
```

**Included data:**
- Lattice vectors (3×3 matrix)
- Atomic species and Cartesian positions
- Atomic forces (eV/Å)
- Stress tensor (kBar, 9 components as flattened 3×3 matrix)
- Energies: `free_energy` (corrected), `energy` (uncorrected)
- Periodic boundary conditions

## Examples

### Organic molecules
```bash
python extract_mptrj_structures.py MPtrj_2022.9_full.json -e H C N O -s only
```

### Water-containing systems
```bash
python extract_mptrj_structures.py MPtrj_2022.9_full.json -e H O -s any
```

### Noble gas compounds
```bash
python extract_mptrj_structures.py MPtrj_2022.9_full.json -e Xe F O -s all
```

### Complete dataset conversion
```bash
python extract_mptrj_structures.py MPtrj_2022.9_full.json -e all -o full_dataset
```

## Output Statistics

The tool provides detailed extraction statistics:

```
======================================================================
EXTRACTION STATISTICS
======================================================================
Total structures processed: 1,580,395
Total time: 372.4 seconds (6.2 minutes)
Processing rate: 4244.7 structures/second

Extracted structures by strategy:
----------------------------------------------------------------------
  only  :     15,234 structures ( 0.96%) |    45.23 MB
======================================================================
```

## Performance

- **Processing speed**: 4,000-6,000 structures/second (CPU-dependent)
- **Memory usage**: ~500 MB (independent of dataset size)
- **Full dataset conversion**: ~30-60 minutes for 1.5M structures

## Dataset Information

The MPtrj dataset contains DFT relaxation trajectories from the Materials Project:
- **Structures**: 1,580,395 trajectory frames
- **Materials**: 154,719 unique materials
- **Source**: [MPtrj dataset](https://doi.org/10.6084/m9.figshare.23713842)

## Technical Details

### Streaming Architecture
Uses `ijson` to parse JSON incrementally:
- Processes one structure at a time
- Never loads entire 12 GB file into memory
- Constant memory footprint regardless of file size

### Stress Tensor Handling
- Input: 3×3 symmetric matrix from DFT (may have small numerical asymmetries)
- Output: Flattened 9-component vector preserving original values
- Format: [σ_xx, σ_xy, σ_xz, σ_yx, σ_yy, σ_yz, σ_zx, σ_zy, σ_zz]
- Units: kBar (matches VASP OUTCAR stress units)

## Loading Converted Data

```python
from ase.io import read

# Read all structures
structures = read('mptrj_extracted/structures_only.xyz', index=':')

# Read single structure
structure = read('mptrj_extracted/structures_only.xyz', index=0)

# Access data
atoms = structures[0]
forces = atoms.get_forces()
stress = atoms.info['stress']  # 9-component vector as string
energy = atoms.info['energy']
```

## Requirements

- Python ≥ 3.7
- ase ≥ 3.22.0
- ijson ≥ 3.2.0
- tqdm ≥ 4.65.0
- numpy ≥ 1.24.0

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.

## Acknowledgments

- Materials Project for the MPtrj dataset
- ASE (Atomic Simulation Environment) developers
