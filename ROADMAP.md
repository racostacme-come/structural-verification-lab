# Seven-day major-project cycle

Thermoelastic transfer reached day-7 acceptance on 2 November 2026, after the
cycle begun 27 October. All numerical methods, independent studies and the
original four-page synthesis are complete. The release candidate includes a
[nine-criterion acceptance matrix](validation/release-acceptance.md), measured
[JSON evidence](validation/release-acceptance.json) and [changelog](CHANGELOG.md).

The final candidate must pass Linux, Windows and manuscript CI before merge;
the merge tree must match the tested tree before v0.2.0 publication. The v0.1.0
beam study and imported source provenance remain preserved. No contribution
is backdated to 28 October.

After release, rotate to Transport Verification Lab and record its next
original seven-day research question in the series portfolio. Retain this
Structural cycle if any release gate fails; 3-4 November are available only
for justified repair, with an explicit blocker/recovery plan on day 9.

Future substantive extensions use a seven-day target, with a permitted 5â€“9-day window:
1. State research question, acceptance criteria and shared-model assumptions.
2. Implement the first numerical improvement and independent oracle.
3. Integrate the second component and test interface failure cases.
4. Complete the end-to-end experiment and parameter sensitivity.
5. Verify convergence, regression tests and cross-platform builds.
6. Draft and compile the <=5-page paper, reproduce figures and inspect pages.
7. Release after evidence review; document limitations and next questions.

Days 8â€“9 are reserved for justified fixes. A missed acceptance criterion blocks a completed-release claim. Each day must record a meaningful contribution and push the corresponding branch. The series portfolio tracks the active cycle.
