# Security Policy

## Supported versions

Metabolopan is distributed as prebuilt binaries and source on the
[Releases](https://github.com/SCArcherKuo/metabolopan/releases) page. Security
fixes are made against the latest release; please reproduce any issue on the
most recent version before reporting.

## Reporting a vulnerability

Please report security issues **privately** — do not open a public issue for a
vulnerability.

- Preferred: use GitHub's private vulnerability reporting on this repository
  (the **Security** tab -> **Report a vulnerability**).
- Alternatively, email **archerkuo9006@gmail.com** with a description, steps to
  reproduce, and the affected version.

You can expect an initial acknowledgement within about a week. Once a fix is
available it will ship in a new release, and the report will be credited unless
you prefer to remain anonymous.

## Scope

Metabolopan runs entirely on the user's machine and reaches only the public
KEGG and PubChem REST APIs over HTTPS to fetch reference data; it does not run a
server or transmit your input data anywhere. Reports about the desktop
application, its release artifacts, and its handling of fetched/cached data are
all in scope.
