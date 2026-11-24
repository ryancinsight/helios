# ADR 0022: Align the published Helios physics provider boundary

- Status: Accepted
- Date: 2026-11-20
- Class: `[major]`

## Context

The 0.1.0 release exposes Aequitas 0.4 quantities and Leto 0.51 array and
geometry types. The registry gives the corresponding 0.5 and 0.52 types
distinct Rust identities. Applications must be able to exchange these types
with RITK and the other Atlas providers under one dependency graph.

## Decision

Move the 11 affected crates to 0.2.0 with Aequitas 0.5, Leto 0.52, Coeus
0.14, Gaia Mesh 0.8, RITK DICOM 0.6, and RITK Registration 0.57. Keep the
allocator at 0.1.0 and Rust 1.95. The 0.2.0 cohort uses Rust 1.97, the
compiler required by its Gaia dependency. Select only the 11 matching
package versions for automatic publication.

Pinning separate old and new provider graphs would compile some isolated
paths while rejecting nominal types at API boundaries. A single published
source allows Cargo to resolve compatible feature combinations and exposes
missing or stale archives as explicit resolution failures.

## Verification

Run the workspace checks and verify the 11 normalized archives. A fresh
registry consumer must resolve one Aequitas, Eunomia, and Leto provider set.
The native CT round trip must preserve the grid pose and calibrated values.
