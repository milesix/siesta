List of software and documentation tasks required to release Siesta (version):

1. [ ] Ensure correct version and date in the top entry of `ReleaseNotes.md`.

2. [ ] Ensure copyright year in `Docs/Contributors.txt` is up to date.

3. [ ] Bump `version.info` to release value (if it is still in use).

----

4. [ ] Run `release.sh` script to create manuals, tarball for distribution, checksums and signatures.

5. [ ] Check that everything is OK with those assets:
    - [ ] Code compiles.
    - [ ] Code runs correctly.
    - [ ] Manuals look OK.
    - [ ] Checksums are correct.
    - [ ] Signatures are valid.

----

6. [ ] Tag the release commit (with leading `v` for versions `4.*` and lower, without `v` if not).

7. [ ] Run `release.py` script to create GitLab release.

----

8. [ ] Bump `version.info` to post-release value (if it is still in use).

9. [ ] Bump version in `ReleaseNotes.md`.
