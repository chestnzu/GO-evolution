FROM python:3.9-slim

WORKDIR /app
COPY ../example_data ../basic* ../covi* ../requirements.txt ../output /app/

RUN apt-get update && apt-get install -y python3-pip && \
pip install --upgrade pip && \
pip3 install -r requirements.txt && mkdir -p /app/output && \
    mkdir -p /app/example_data && \
    mv /app/covidGroups /app/example_data/covidGroups && \
    mv /app/data /app/output/data && \
    mv /app/figure /app/output/figure

# Volumes are only available for a container, not while building an image.
# Therefore when mounting, use /app/output, instead of /output
VOLUME ['/output']   


ENTRYPOINT ["python3","covid_enrichment.py"]
CMD ['2017-12-15','2021-12-15','./example_data/covidGroups']
