# Mineral Interchange Project
This GitHub repository is associated with a research paper titled "Trusted Data Sharing for Mineral Exploration and Mining Tenements" which is currently under consideration for publication in the journal Computers & Geosciences.


## Description
The project introduces an approach to accelerate the exploration-discovery-mining process by facilitating secure sharing of mineral data while maintaining confidentiality. To elaborate, we explore the integration of data confidentiality/privacy risk assessment, data privatization techniques, and obfuscation technologies. This integration aims to enable the sharing of sensitive data that would otherwise remain private. We've established a set of metrics to quantify the trade-off between data confidentiality loss and utility gain when sharing data. Additionally, we've devised various approaches to obfuscate data, including value removal, binning, and sampling, before sharing. We've examined and empirically validated the impact of different obfuscation methods on both confidentiality and utility using an actual mineral dataset provided by an Australian minerals company.

**It's important to note that the genesis of this project lies in the notebook labeled [paper_notebook_code/Input_Generator.ipynb](file://paper_notebook_code/Input_Generator.ipynb).**

**In light of confidentiality considerations, we regretfully cannot disclose two files:** 
* the file that compiles tenements along with their respective sizes and coordinates, 
* and the file that compiles Surface Auger Sampling.
## Installation

Please make sure you have **Python 3.10** installed.

When setting up your project or running experiments, consider creating a **conda environment** to isolate your project's dependencies. This helps prevent conflicts between packages and ensures reproducibility. 

```
conda create --name my_project_env python=3.10
```
The necessary dependencies and packages are listed in the "requirements.txt" file. 
```commandline
pip install -r requirements.txt

```

This repository is organized into two main directories:

* **utils**: This directory contains various utility functions that are useful for the experiments conducted in the paper.

* **paper_notebook_code**: In this directory, you will find all the code used for running the experiments mentioned in the paper. These notebooks and scripts are related to the experiments discussed in the paper and provide insights and results based on the conducted research.


## License

This project is licensed under the CSIRO BSD / MIT licence  - see the [LICENSE](CSIRO BSD MIT Licence v2.0-5.txt) file for details.

## Contact
If you have any questions or suggestions, feel free to contact me at sirine.mrabet@data61.csiro.au






_Last updated 17/08/2023 by Sirine M'rabet_
