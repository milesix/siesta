List of software and documentation tasks required to release Siesta (version):

1. [ ] Ensure correct version and date in the top entry of `ReleaseNotes.md`.

2. [ ] Ensure copyright year in `Docs/Contributors.txt` is up to date.

3. [ ] Check test outputs and update them if needed.

4. [ ] Bump `version.info` to release value (if it is still in use).

----

5. [ ] Run `release.sh` script to create manuals, tarball for distribution, checksums and signatures.

6. [ ] Check that everything is OK with those assets:
    - [ ] Code compiles.
    - [ ] Code runs correctly.
    - [ ] Manuals look OK.
    - [ ] Checksums are correct.
    - [ ] Signatures are valid.

----

7. [ ] Tag the release commit (with leading `v` for versions `4.*` and lower, without `v` if not).

8. [ ] Create a page on the GitLab wiki with the Release notes for this version.

9. [ ] Update the [Guide to Siesta versions](https://gitlab.com/siesta-project/siesta/-/wikis/Guide-to-Siesta-versions) on GitLab.
       Check that the new Release Notes wiki page is linked from it.

10. [ ] Run `release.py` script to create GitLab release.
       GitLab will automatically send an email notification of the release to users that track the repository.

11. [ ] Update the [Guide to Siesta versions](https://gitlab.com/siesta-project/siesta/-/wikis/Guide-to-Siesta-versions) on GitLab
        by adding a link to the GitLab release that was just created.

----

12. [ ] Bump `version.info` to post-release value (if it is still in use).

13. [ ] Bump version in `ReleaseNotes.md`.

14. [ ] Create new milestone for next release to enable tracking issues and merge requests.
