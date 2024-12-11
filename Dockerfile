# Dockerfile
FROM python:3.9-slim

# Set up working directory
WORKDIR /app

# Install system dependencies (wget or curl)
RUN apt-get update && apt-get install -y wget

# Install Python dependencies
RUN pip install --no-cache-dir gprofiler-official pandas matplotlib

# Copy analysis script
COPY analysis_script.py .

# Optional: If you want to add the CSV files directly into the container
# COPY GSE118523_20161109_old_wt_tg.csv /app/
# COPY GSE118523_20161109_young_wt_tg.csv /app/

# Set the entry point
ENTRYPOINT ["python", "analysis_script.py"]

