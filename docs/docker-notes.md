# Using Oligopool Calculator with Docker

On non-Linux and Linux systems both, you can use `oligopool` via [Docker](https://docs.docker.com/desktop/).

**Related documentation**:
- [User Guide](docs.md) - Tutorials, examples, and workflows
- [API Reference](api.md) - Complete parameter documentation for all modules
- [AI Agent Guide](agent-skills.md) - Decision trees, best practices, and gotchas

To follow this mini guide you will need to have [`git`](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) and `docker` installed.


### Building `Oligopool Calculator` `docker` image.

The first step in the process is to build an `oligopool` `docker` image.

1. Either manually or via terminal download or clone this repository.
```bash
$ git clone https://github.com/ayaanhossain/oligopool.git # Or, manually download repository
```

2. Navigate into the repo.
```bash
$ cd oligopool
```

3. And then, build the `docker` image from the repo.
```bash
$ docker build -t oligopool-docker . # Note the dot at the end
```
> **Note** that the name of the image (`oligopool-docker`) is flexible, and can be anything we want.

4. You can see a list of all `docker` images built on your system, if you need to recall for later.
```bash
$ docker images
REPOSITORY         TAG       IMAGE ID       CREATED              SIZE
oligopool-docker   latest    414f88c5e29c   About a minute ago   1.76GB
```

### Attaching your project directory to a `docker` container.

With the `docker` image created, we can now start a container within which `oligopool` is available. You can then mount your project directory and start using `Oligopool Calculator`.

1. First, navigate to your project directory.
```bash
$ cd path/to/your/project
```

2. If you are on `MacOS` or `Linux`, you can spin a container from `oligopool-docker` image and load your project directory (`pwd` path) as `op-workspace` inside the container.
```bash
$ docker run -it -v $(pwd):/op-workspace --name op-container oligopool-docker # Loads your project directory
```
> **Note** `op-workspace` is just a name for your project directory inside the container we spun up; it can be anything. Similarly `op-container` is the name of the container we spun up from the `oligopool-docker` image.

3. If you are on `Windows`, the syntax for this is slightly different.
```powershell
$ docker run -it -v ${PWD}:/op-workspace --name op-container oligopool-docker # Loads your project directory
```

4. Once the container is up and running, you can access your directory content.
```bash
$ ls -la # for example, will show you your project directory content
```
The terminal within the `docker` container will now be a `Linux` (Ubuntu) `bash` prompt. You can now execute any `oligopool` function or script from the terminal.
```bash
$ python
Python 3.12.3 (main, Sep 11 2024, 14:17:37) [GCC 13.2.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> import oligopool as op
>>> op.__version__
```
You can also use the CLI entry points:
```bash
$ op
$ op manual
$ oligopool
```

5. When you are done with the container work, simply exit from the terminal.
```bash
$ exit
```

6. You can see a list of all running `docker` containers (not images) using
```bash
$ docker ps -a
```
which will even show you stopped containers you spinned earlier. If you want to remove stopped `docker` containers, check the [documentation](https://docs.docker.com/reference/cli/docker/container/rm/).


### Using `Oligopool Calculator` within `docker` via `jupyter`.

`Oligopool Calculator` for `Design Mode` is best used interactively via `jupyter`. `Analysis Mode` functions could need to run for longer so perhaps better to be used via scripts.

Containers from the `oligopool-docker` image we just created can connect to `jupyter`. We will need to first map the port `8080` within the container to a port on our host machine.

1. From your host machine terminal (this is outside the `container` prompt, so exit if you are within the container), map the port `8080` within the container to port `8888` on your host and start an interactive `bash` session from your project directory.
```bash
$ docker run -p 8888:8080 -it -v $(pwd):/op-workspace --name op-container oligopool-docker
```
> **Note** The port `8080` within the `docker` image is exposed by design. If you want to expose different ports, modify the `EXPOSE` line in `Dockerfile` in the repo and rebuild the image.
> Similarly, feel free to map to any port other than `8888` on your local host. You can even map `8080` within container to `8080` on your host.

If you are on `Windows`, use:
```powershell
$ docker run -p 8888:8080 -it -v ${PWD}:/op-workspace --name op-container oligopool-docker
```

2. From the `bash` terminal within the container simply use the following
```bash
$ jupyter notebook --ip=0.0.0.0 --port=8080 --no-browser --allow-root
```
which will start a `jupyter` notebook within the container, but now accessible at
```
http://localhost:8888
```
on your local browser. This is very similar to how we would start `jupyter` notebooks from the terminal on our host.


> **Note** If you want more packages within this docker image, you can modify the `pip install` lines in our `Dockerfile`.

---

## Next Steps

Once your Docker container is running, refer to:
- [User Guide](docs.md) - Learn the design and analysis workflows
- [API Reference](api.md) - Look up specific function parameters
- [AI Agent Guide](agent-skills.md) - Best practices and troubleshooting tips
