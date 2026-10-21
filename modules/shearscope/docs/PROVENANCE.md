# Scope and provenance

Project selected from the local AFEM topic directories for Timoshenko beam
elements and transverse shear locking, and the Structural Dynamics subject.
Only folder and file names were used for topic discovery. No course document,
exam, exercise implementation, or private data was opened or incorporated.

All numerical source, tests, prose and example data were created for this project.
Standard beam theory informs the equations. The numerical implementation was
written independently; no public tutorial code was copied. The public references
in the README explain the underlying mechanics and provide further reading.

Development was assisted by OpenAI Codex. Commits retain the user's configured Git
identity and actual timestamps. No external review or experimental validation is
claimed. The sample CSV and PNG files are simulations, not measured data.

Initial workflow: scaffold on main; feature/beam-kernel for element and solvers;
feature/validation-portfolio for campaign, CLI, docs and CI. Feature merges follow
passing local checks. GitHub CI is a separate check after publication.
