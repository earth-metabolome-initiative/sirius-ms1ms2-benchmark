# Sirius MS2 vs MS1+MS2 Annotations
[Sirius](https://github.com/sirius-ms/sirius) uses MS1 isotope patterns to improve molecular formula and structure annotations. But how much do MS1 data actually help? This repository contains code and data to compare Sirius annotations with and without MS1 information.

## Download the data
```bash
mkdir archive
cd archive
wget TODO
wget TODO
cd ..
```

## Requirements
### Sirius CLI Tool
You will need to have Sirius installed and available on your PATH. You can download Sirius from [here](https://github.com/sirius-ms/sirius/releases) and follow the installation instructions for your platform.


### Metabolite Annotator
```bash
uv tool install git+https://github.com/earth-metabolome-initiative/metabolite-annotator
```

## Results
See here: [`summary.md`](./summary.md)