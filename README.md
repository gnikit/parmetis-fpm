# ParMETIS for fpm

This provides a Fortran API and a source repackaging of the
[ParMETIS](https://github.com/KarypisLab/ParMETIS) library originally developed by Karypis Lab

## TODO

- [ ] Run tests in CI
- [ ] Make CI open a PR for changes if the repository is dirty
- [ ] Add Fortran API procedural and OOP interfaces

## Usage

```sh
fpm build --flag "-DIDXTYPEWIDTH=64 -DREALTYPEWIDTH=64" --compiler mpif90 --c-flag "-DIDXTYPEWIDTH=64 -DREALTYPEWIDTH=64" --c-compiler mpicc
```

To use `parmetis` as a dependency in your `fpm` project, add the following to your `fpm.toml` file:

```toml
[dependencies]
parmetis = { git = "https://github.com/gnikit/parmetis-fpm.git" }
```

## License

> MIT License

[metis-fpm](https://github.com/gnikit/metis-fpm) is distributed under the MIT License. See [LICENSE](LICENSE) for more information.
[METIS](https://github.com/KarypisLab/METIS) is distributed under the terms of the Apache License, Version 2.0
