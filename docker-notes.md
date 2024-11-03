## Using `Oligopool Calculator` with `docker`

On non-`Linux` and `Linux` systems both, you can use `oligopool` via [`docker`](https://docs.docker.com/desktop/).
```bash
$ git clone https://github.com/ayaanhossain/oligopool.git # Or, manually download repository
$ cd oligopool # Navigate into the downloaded repository
$ docker build -t oligopool-docker . # Build docker image and save it as oligopool-docker
```
> **Note** that the name of the image (`oligopool-docker`) is flexible, and can be anything we want.

You can see a list of all `docker` images built on your system.
```bash
$ docker images
REPOSITORY         TAG       IMAGE ID       CREATED              SIZE
oligopool-docker   latest    414f88c5e29c   About a minute ago   1.76GB
```

You can now mount your project directory and start using `Oligopool Calculator`.
```bash
$ cd path/to/your/project # Navigate to your project directory
$ # For MacOS and Linux
$ docker run -it -v $(pwd):/op-workspace oligopool-docker # Loads your project directory
$ # For Windows
$ docker run -it -v ${PWD}:/op-workspace oligopool-docker # Loads your project directory
$ # After Loading
$ ll # Shows you your files! You are in a Linux environment now.
```
> **Note** `op-workspace` is just a name for your project directory inside docker, could be anything.