# OpenBLAS Linear Algebra

A collection of commonly used functions for linear algebra for Fortran using [OpenBLAS](https://github.com/xianyi/OpenBLAS).

*Suggestions and code contributions are welcome.*

```fortran
use open_linalg_m, only: det, inv, matmul, operator(.i.), operator(.x.), solve
```

## Build with [fortran-lang/fpm](https://github.com/fortran-lang/fpm)

```sh
fpm run --example --list
```

```toml
[dependencies]
open_linalg = { git = "https://github.com/zoziha/open_linalg.git" }
```

## Reference

- [fortran-lang/stdlib](https://github.com/fortran-lang/stdlib)
- [fortran-fans/forlab](https://github.com/fortran-fans/forlab)
- [QcmPlab/scifortran](https://github.com/QcmPlab/SciFortran)
