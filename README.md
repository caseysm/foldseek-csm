# Foldseek Align

A package for running Foldseek analysis and creating distance matrices.

## Installation

### Step 1: Install Foldseek

Foldseek is a required dependency and is best installed using conda:

```bash
conda install -c bioconda foldseek
```

### Step 2: Install foldseek_align

After installing Foldseek, you can install foldseek_align using pip:

```bash
pip install foldseek_align
```

Alternatively, if you want pip to be aware of the Foldseek requirement (though it won't install Foldseek):

```bash
pip install foldseek_align[foldseek]
```

## Usage

```python
from foldseek_align import foldseek_main

foldseek_main(input_dir="path/to/input", output_dir="path/to/output")
```

Or from the command line:

```bash
foldseek_align path/to/input path/to/output
```

## Development

For development, clone the repository and install in editable mode:

```bash
git clone https://github.com/yourusername/foldseek_align.git
cd foldseek_align
conda install -c bioconda foldseek
pip install -e .[foldseek]
```

This will install the package in editable mode with the Foldseek requirement specified.