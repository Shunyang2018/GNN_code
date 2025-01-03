FROM conda/miniconda3

RUN conda install -c conda-forge rdkit==2018.09.3


RUN conda update -n base -c defaults conda

RUN pip install torch==1.11.0\
    && pip install tqdm

RUN pip install dgl

RUN pip install dgllife
RUN pip install Flask
RUN pip install gunicorn

SHELL ["/bin/bash", "-c"]

COPY api.py .
COPY utils.py .
COPY c18/ ./C18/
COPY pfp/ ./PFP/

COPY run_server.sh ./
RUN chmod +x run_server.sh

EXPOSE 9002

ENTRYPOINT ["./run_server.sh"]
