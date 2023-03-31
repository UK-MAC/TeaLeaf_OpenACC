# TeaLeaf OpenACC

This is an OpenACC port of Tealeaf reference

## Compling

In most cases one needs to load the appropriate modules (including those for the accelerator in question) and then:

```
make clean; make COMPILER=CRAY
```

or,

```
make clean; make COMPILER=PGI
```
