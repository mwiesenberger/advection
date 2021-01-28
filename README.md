# Advection project

Testing and exploration ground for advection and regularization schemes.
Generate the site https://mwiesenberger.github.io/advection

## Install
In order to generate the static website we use
[jupyter-book](https://jupyterbook.org).

```bash
# install jupyter-book
pip install -U jupyter-book
# clone this repository
git clone https://github.com/mwiesenberger/advection
# build the book
jupyter-book build path/to/advection
# we use ghp-import to publish changes on the github pages
pip install ghp-import
```
In order to locally generate the simulation data you will need the
[Feltor](https://github.com/feltor-dev/feltor) code repository.  Follow the
quick-start guide to install.  It is recommended to keep Feltor and this
repository next to each other.  If you prefer not to, you need to set the
`FELTOR_PATH` environment variable in order for the `execute.sh` script to
compile and execute the `path/to/feltor/src/lamb_dipole/shu_b` code.

## Usage
Build the book with
```bash
jupyter-book build path/to/advection
```
To publish changes after the book was built
```bash
cd path/to/advection
ghp-import -n -f -p _build/html
```

## How to reproduce the results
The data that the simulation generates is too large to host alongside the
notebooks so unfortunately we will not be able to run an interactive notebook
online. You will have to download and run this repository and its dependencies
described above locally on your machine in order to reproduce the data and
interact with the notebooks.

## Author
Matthias Wiesenberger
