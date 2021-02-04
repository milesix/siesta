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

7. [ ] Create a page on the GitLab wiki with the Release notes for this version..

8. [ ] Update the [Guide to Siesta versions](https://gitlab.com/siesta-project/siesta/-/wikis/Guide-to-Siesta-versions) on GitLab.
       Check that the new Release Notes wiki page is linked from it.

9. [ ] Run `release.py` script to create GitLab release.
       GitLab will send a notification of the release to users that track the repository.

10. [ ] Update the [Guide to Siesta versions](https://gitlab.com/siesta-project/siesta/-/wikis/Guide-to-Siesta-versions) on GitLab
        by adding a link to the GitLab release that was just created.

----

11. [ ] Bump `version.info` to post-release value (if it is still in use).

12. [ ] Bump version in `ReleaseNotes.md`.
