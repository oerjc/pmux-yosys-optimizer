# Third-party license index

Original copyright and license notices continue to apply. This index locates
notices; it does not replace or relicense any third-party component. Paths below
are relative to this `licenses` directory unless specified otherwise.

| Component | Notice location and source |
|---|---|
| Yosys v0.69 | `COPYING`, from upstream commit `9f75ca1f9834a39a863915b5dae0c7b1e33533bc` |
| ABC | `abc/copyright.txt` and notices under `abc/src`; source commit `a51bf4a7bb38ce06ea3a544b1463b40832213917` |
| Boost.Regex | `libs/boost_regex/LICENSE_1_0.txt`; Boost Software License 1.0 |
| json11 | `libs/json11/LICENSE.txt`; Dropbox MIT notice extracted from upstream `libs/json11/json11.hpp` |
| ezSAT | `libs/ezsat/LICENSE.txt`; ISC notice extracted from upstream `libs/ezsat/ezsat.h` |
| cxxopts | `libs/cxxopts/LICENSE` |
| fmt | `libs/fmt/LICENSE`, plus `libs/fmt/doc/LICENSE-exception` |
| MiniSat | `libs/minisat/LICENSE` |
| symfpu | `libs/symfpu/LICENSE`, `LICENSE-BSD`, and `LICENSE-GPL3` |
| toml++ | `libs/tomlplusplus/LICENSE`; its additional image notice is retained separately |
| GCC runtime libraries | `gcc-libs/README`, `COPYING3`, `COPYING.RUNTIME`, and `COPYING.LIB`; includes the GCC Runtime Library Exception |
| winpthreads | `libwinpthread/COPYING` |
| Big Integer Library | Public-domain statement in vendored Yosys `libs/bigint/README`, section “Legal” |
| SHA-1 | Public-domain statement in vendored Yosys `libs/sha1/sha1.h` |
| Portable CMake | The vendored CMake wheel retains `cmake/data/doc/cmake/LICENSE.rst`, dependency notices, and `cmake-4.4.4.dist-info/licenses/` |
| Portable Ninja | The vendored Ninja wheel retains `ninja-1.13.2.dist-info/licenses/LICENSE_Apache_20` |

The source bundle's `vendor-sha256.json` identifies the exact archives containing
upstream source and embedded notices. Windows system libraries are supplied by
Windows and are not redistributed in this package. Official competition examples
retain their authors' rights and are not covered by the project's ISC grant.

## Provenance of the three added license texts

Boost.Regex's vendored `include/boost/regex/config.hpp` explicitly identifies the
Boost Software License, Version 1.0, and points to
<https://www.boost.org/LICENSE_1_0.txt>. The archive did not contain that separate
file. The complete standard text supplied here was transcribed from the installed
Microsoft OneDrive 26.183.0920.0005 `ThirdPartyNotices.txt`, section “Boost”, on
2026-10-19; only its four-space indentation was removed. The local source file's
SHA-256 was
`89ae55bf169ed1b3d610576f6ebb243cd8048ffa274bc6414fdab8b00367386c`.
The copied license text's SHA-256 is
`c9bff75738922193e67fa726fa225535870d2aa1059f91452c411736284ad566`.
The upstream URL is a verification reference; this task did not retrieve it.

The json11 notice was extracted from the vendored Yosys source header
`libs/json11/json11.hpp`, SHA-256
`d653ab8c6e132d97d41fe95a75f9a6f089af456631e8ad541d3adfc00238b8d2`.
The standalone notice's SHA-256 is
`950e288a2548cb33f6a90d23958806b0803f510ab17c600d661619a73ae3a6c8`.

The ezSAT notice was extracted from the vendored Yosys source header
`libs/ezsat/ezsat.h`, SHA-256
`eb839cf23aad3754ec9e6bb31e550148aa8aea33f9bd9eda33a8a9d441708387`.
The standalone notice's SHA-256 is
`1e36dceee80781e78c7e52760cd66d45114303eb3d20935cd8a1c2cec6ed51b0`.
For both header extracts, comment delimiters were removed and the original
copyright, grant, condition, and disclaimer text were retained.
