FROM cannin/jupyter-r:enrichmentmap

COPY * /workspace/
WORKDIR /workspace

# https://github.com/ipython/ipython/issues/7062
CMD sh -c "ipython notebook --ip=*"
