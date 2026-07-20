FROM mambaorg/micromamba:2.0.5

WORKDIR /app
COPY --chown=$MAMBA_USER:$MAMBA_USER . /app
RUN micromamba install -y -n base -f environment.yml && micromamba clean --all --yes

EXPOSE 8501
CMD ["streamlit", "run", "app/streamlit_app.py", "--server.address=0.0.0.0"]

