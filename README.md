# Sirius MS2 vs MS1+MS2 Annotations
[Sirius](https://github.com/sirius-ms/sirius) uses MS1 isotope patterns to improve molecular formula and structure annotations. But how much do MS1 data actually help? This repository contains code and data to compare Sirius annotations with and without MS1 information.

## Download the data
WARNING: Each file is about 30 GB in size, so make sure you have enough disk space before downloading.

```bash
mkdir archive
cd archive
wget https://zenodo.org/records/19681322/files/emi-positive-archive.tar.zst
wget https://zenodo.org/records/19692818/files/emi-negative-archive.tar.zst 
cd ..
```

## Requirements
### Sirius CLI Tool
You will need to have Sirius installed and available on your PATH. You can download Sirius from [here](https://github.com/sirius-ms/sirius/releases) and follow the installation instructions for your platform.


### Metabolite Annotator
```bash
uv tool install git+https://github.com/earth-metabolome-initiative/metabolite-annotator
```


## Run
```bash
uv run main.py
````

Or if you want to build the summary of the results:
```bash
uv run summary.py
```

## Results
See here: [`summary.md`](./summary.md)