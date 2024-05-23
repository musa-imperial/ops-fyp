sudo docker run --gpus all --rm -it -v `pwd`:`pwd` -w `pwd` -u $(id -u):$(id -g) devitocodes/devito:nvidia-nvc-latest python 2d_diffusion.py
