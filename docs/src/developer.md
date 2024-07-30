# [Developer documentation](@id dev_docs)

!!! note "Contributing guidelines"
    If you haven't, please read the
    [Contributing guidelines](contributing.md) first.

If you want to make contributions to this package that involves code,
then this guide is for you.

## First time clone

!!! tip "If you have writing rights"

    If you have writing rights, you don't have to fork. Instead,
    simply clone and skip ahead. Whenever **upstream** is mentioned,
    use **origin** instead.

If this is the first time you work with this repository, follow the
instructions below to clone the repository.

1. Fork this repo
2. Clone your repo (this will create a `git remote` called `origin`)
3. Add this repo as a remote:

   ```bash
   git remote add upstream https://github.com/moble/PostNewtonian.jl
   ```

This will ensure that you have two remotes in your git: `origin` and
`upstream`.  You will create branches and push to `origin`, and you
will fetch and update your local `main` branch from `upstream`.

## Linting and formatting

Install a plugin on your editor to use
[EditorConfig](https://editorconfig.org).  This will ensure that your
editor is configured with important formatting settings.

We use [https://pre-commit.com](https://pre-commit.com) to run the
linters and formatters.  In particular, the Julia code is formatted
using
[JuliaFormatter.jl](https://github.com/domluna/JuliaFormatter.jl), so
please install it globally first:

```julia-repl
julia> # Press ]
pkg> activate
pkg> add JuliaFormatter
```

To install `pre-commit`, use your standard python method — for
example:

```bash
python -m pip install pre-commit
```

With `pre-commit` installed, activate it as a pre-commit hook by
running this command in the repository folder:

```bash
pre-commit install
```

To run the linting and formatting manually, enter the command below:

```bash
pre-commit run -a
```

**Now, you can only commit if all the pre-commit tests pass**.

## Testing

As with most Julia packages, you can just open Julia in the repository
folder, activate the environment, and run `test`:

```julia-repl
julia> # press ]
pkg> activate .
pkg> test
```

## Working on a new issue

We try to keep a linear history in this repo, so it is important to
keep your branches up-to-date.

1. Fetch from the remote and fast-forward your local main

   ```bash
   git fetch upstream
   git switch main
   git merge --ff-only upstream/main
   ```

2. Branch from `main` to address the issue (see below for naming)

   ```bash
   git switch -c 42-add-answer-universe
   ```

3. Push the new local branch to your personal remote repository

   ```bash
   git push -u origin 42-add-answer-universe
   ```

4. Create a pull request to merge your remote branch into the org main.

### Branch naming

- If there is an associated issue, add the issue number.
- If there is no associated issue, **and the changes are small**, add
  a prefix such as "typo", "hotfix", "small-refactor", according to
  the type of update.
- If the changes are not small and there is no associated issue, then
  create the issue first, so we can properly discuss the changes.
- Use dash separated imperative wording related to the issue (e.g.,
  `14-add-tests`, `15-fix-model`, `16-remove-obsolete-files`).

### Commit message

- Use imperative or present tense, for instance: *Add feature* or *Fix
  bug*.
- Have informative titles.
- When necessary, add a body with details.
- If there are breaking changes, add the information to the commit
  message.

### Before creating a pull request

!!! tip "Atomic git commits"
    Try to create "atomic git commits" (recommended reading: [The Utopic
    Git
    History](https://blog.esciencecenter.nl/the-utopic-git-history-d44b81c09593)).

- Make sure the tests pass.
- Make sure the pre-commit tests pass.
- Fetch any `main` updates from upstream and rebase your branch, if
  necessary:

  ```bash
  git fetch upstream
  git rebase upstream/main BRANCH_NAME
  ```

- Then you can open a pull request and work with the reviewer to
  address any issues.

## Building and viewing the documentation locally

Following the latest suggestions, we recommend using `LiveServer` to
build the documentation.  Here is how you do it:

1. Run `julia --project=docs` to open Julia in the environment of the
   docs.
1. If this is the first time building the docs
   1. Press `]` to enter `pkg` mode
   1. Run `pkg> dev .` to use the development version of your package
   1. Press backspace to leave `pkg` mode
1. Run `julia> using LiveServer`
1. Run `julia> servedocs()`

## Making a new release

To create a new release, once the code is merged into `main`, just go
to [the `register`
action](https://github.com/moble/PostNewtonian.jl/actions/workflows/register.yml),
click the gray "Run workflow" button, enter the version component you
want to bump, and click the green "Run workflow" button.

After that, you only need to wait and verify:

- Wait for the bot to comment (should take < 1m) with a link to a RP
  to the registry
- Follow the link and wait for a comment on the auto-merge
- The comment should said all is well and auto-merge should occur
  shortly
- After the merge happens, TagBot will trigger and create a new GitHub
  tag. Check on <https://github.com/moble/PostNewtonian.jl/releases>
- After the release is create, a "docs" GitHub action will start for
  the tag.
- After it passes, a deploy action will run.
- After that runs, the [stable
  docs](https://moble.github.io/PostNewtonian.jl/stable) should be
  updated. Check them and look for the version number.
