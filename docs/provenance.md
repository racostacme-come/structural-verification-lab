# Portable snapshot provenance

The original repositories and their pinned commits are unchanged. The initial import used Git archives on Windows. MeshBeacon's `text=auto` attributes caused eleven archived text files to have CRLF bytes, although their original Git blobs have LF bytes. The import manifest hashed those archive bytes. The integration commit stored LF blobs, so Linux verification failed while a Windows checkout passed.

The 21 October 2026 repair replaces those eleven manifest hashes with hashes of the original pinned Git blobs and normalizes the local copies to those exact bytes. `provenance-repair.json` retains both the previous archive hash and the canonical source hash for every affected path. There are no numerical source changes.

One intentional imported-file patch changes MeshBeacon's first `.gitattributes` rule from `* text=auto` to `* text=auto eol=lf`. This prevents the same platform-dependent checkout behavior. `suite.json` records the original and patched hashes under `integration_patches`; its ordinary file manifest verifies the patched bytes. All other imported files must match their original pinned Git blobs exactly.

The integrity gate remains byte-exact. It does not silently normalize data, ignore a failed hash, or rewrite original history. The CI workflow checks it on Windows and Linux. Root tests also demonstrate that changed numerical source and missing files fail the gate.
