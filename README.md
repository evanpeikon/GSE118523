# Analysis from CNS cell-type specific gene profiling of aging P301S tau transgenic mice
1. Download dockerfile (Dockerfile) and python analysis script (analysis_script.py)
2. Build docker image: ```docker build -t gene-analysis .```
3. Run the analysis: ```docker run --rm -v $(pwd):/app gene-analysis```

# Results
The code above analysis will produce the following images:

1. Differential expression analysis (young WT vs. yound transgenic mice)
![young_mice_volcano_plot](https://github.com/user-attachments/assets/89569c0f-4dac-4ec0-b205-b905d832e1e6)

2. Top 10 enriched GO terms for young mice
![young_mice_go_terms](https://github.com/user-attachments/assets/e1ee4c45-cfe0-4295-b201-264959354329)

3. Differential expression analysis (old WT vs. old transgenic mice)
![old_mice_volcano_plot](https://github.com/user-attachments/assets/86ad8950-b227-477b-bddd-afea40d92a3e)

4. 2. Top 10 enriched GO terms for old mice
![old_mice_go_terms](https://github.com/user-attachments/assets/adc1bb89-c529-4d06-97cc-8a9eaed07c74)
