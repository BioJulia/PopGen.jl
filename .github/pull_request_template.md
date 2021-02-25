# A clear and descriptive title (No issue numbers please)

> _This template is extensive, so fill out all that you can. If are a new contributor or unsure about any section, leave it empty and a reviewer will help you_ :smile:. If anything in this template isn't relevant, just delete or ignore it._

## Types of changes

This PR implements the following changes:
_(tick any that apply)_

* [ ] :sparkles: New feature (A non-breaking change which adds functionality).
* [ ] :bug: Bug fix (A non-breaking change, which fixes an issue).
* [ ] :boom: Breaking change (fix or feature that would cause existing functionality to change).

## :clipboard: Additional details
>_replace this block of text with your information_
- If you have implemented new features or behaviour
  - **Provide a description of the addition** in as many details as possible.

  - **Provide justification of the addition**.

  - **Provide a runnable example of use of your addition**. This lets reviewers
    and others try out the feature before it is merged or makes it's way to release.

- If you have changed current behaviour...
  - **Describe the behaviour prior to you changes**

  - **Describe the behaviour after your changes** and justify why you have made the changes,
    Please describe any breakages you anticipate as a result of these changes.

  - **Does your change alter APIs or existing exposed methods/types?**
    If so, this may cause dependency issues and breakages, so the maintainer
    will need to consider this when versioning the next release.

  - If you are implementing changes that are intended to increase performance, you
    should provide the results of a simple performance benchmark exercise
    demonstrating the improvement. Especially if the changes make code less legible.

## :ballot_box_with_check: Checklist
>_it's ok if not all the boxes are checked :smile:_
- [ ] :art: The changes implemented is consistent with the [julia style guide](https://docs.julialang.org/en/v1/manual/style-guide/).
- [ ] :blue_book: I have updated and added relevant docstrings, in a manner consistent with the [documentation styleguide](https://docs.julialang.org/en/v1/manual/documentation/).
- [ ] :blue_book: I have added or updated relevant user and developer manuals/documentation in `docs/src/`.
- [ ] :ok: There are unit tests that cover the code changes I have made.
- [ ] :ok: The unit tests cover my code changes AND they pass.
- [ ] :ok: All changes should be compatible with the latest stable version of Julia.
- [ ] :thought_balloon: I have commented liberally for any complex pieces of internal code.
