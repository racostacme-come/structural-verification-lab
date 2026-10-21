# Contributing

Keep model changes separate from presentation changes. Use a focused feature
branch and a descriptive commit; run the README's formatting, Python, native,
and packaging checks before merging. Rebuild the extension after C++ changes.
Run the benchmark and inspect its JSON and PNG before updating tracked results.

Tests should use independently derived expected values, meaningful invariants,
or failure cases. Do not replace failed numerical checks with looser thresholds
without explaining and justifying the change. Preserve failure history and state
the compiler/platform on which a numerical claim was measured.

Do not add course handouts, exam material, exercise solutions, private data,
credentials, binary build products, or virtual environments. The project has no
private data dependency. Use MIT-compatible original contributions.
