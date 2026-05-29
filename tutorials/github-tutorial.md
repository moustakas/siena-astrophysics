# Git and GitHub: A Hands-On Tutorial

A ~1-hour introduction to version control for undergraduate science
students.  No prior experience assumed.

## Contents

- [Before You Begin](#before-you-begin)
- [1. Why Version Control?](#1-why-version-control)
- [2. Configuration](#2-configuration)
- [3. Core Concepts](#3-core-concepts)
- [4. Your First Repository](#4-your-first-repository)
- [5. History and Differences](#5-history-and-differences)
- [6. Pushing to GitHub](#6-pushing-to-github)
- [7. Cloning and Forking](#7-cloning-and-forking)
- [8. Branches](#8-branches)
- [9. Pull Requests and Merging](#9-pull-requests-and-merging)
- [10. Daily Workflow](#10-daily-workflow)
- [Going Further](#going-further)
- [Additional Resources](#additional-resources)

---

## Before You Begin

Please complete the following **before** the tutorial
session.

### 1. Create a GitHub account

Go to [github.com](https://github.com) and sign up for a free
account. If you use your `.edu` email address you can create free
private repositories via the [GitHub
Education](https://education.github.com) program.

### 2. Install Git

**macOS** — Open Terminal (Applications → Utilities → Terminal) and run:
```bash
xcode-select --install
```
This installs Apple's command-line developer tools, which include Git. Click
*Install* when the dialog appears.

**Windows** — Download and install **Git for Windows** from
[gitforwindows.org](https://gitforwindows.org). Accept all defaults during
installation. This installs Git, **Git Bash** (the terminal you will use
throughout this tutorial), and Git Credential Manager (which handles
authentication automatically).

> **Windows users:** use **Git Bash** for all terminal commands in this
> tutorial — not PowerShell and not Command Prompt. After installation, find
> it by searching "Git Bash" in the Start menu.

**Linux (Ubuntu/Debian):**
```bash
sudo apt install git
```

Verify the installation:
```bash
git --version
```

### 3. Install VS Code

Download and install Visual Studio Code from
[code.visualstudio.com](https://code.visualstudio.com). It runs on macOS,
Windows, and Linux.

VS Code is a free text editor with an integrated terminal. Any text editor
works — `nano` is a simpler terminal-based alternative — but this tutorial
assumes VS Code. To open a terminal inside VS Code, use **View → Terminal**
(or `` Ctrl+` ``).

### 4. Authenticate with GitHub

The **GitHub CLI** (`gh`) handles authentication so that `git push` works
without passwords or tokens.

**macOS** — First install Homebrew from [brew.sh](https://brew.sh) if you
don't already have it, then:
```bash
brew install gh
```

**Windows** — Download the installer from [cli.github.com](https://cli.github.com).

**Linux** — See installation instructions at [cli.github.com](https://cli.github.com).

Once installed, run:
```bash
gh auth login
```
Choose *GitHub.com* → *HTTPS* → *Login with a web browser* and follow the prompts.

> **Windows users:** Git Credential Manager (bundled with Git for Windows) may
> authenticate you automatically the first time you push — a browser window will
> open asking you to sign in to GitHub. If that happens, you don't need `gh`.

---

## 1. Why Version Control?

You've probably ended up with something like this:

```
spectrum_fit.py
spectrum_fit_v2.py
spectrum_fit_final.py
spectrum_fit_final_FIXED.py
spectrum_fit_final_FIXED_v2.py
```

This is version control by hand — and it fails. You can't see *what* changed
between versions, *why* it changed, or recover the exact state from three weeks
ago.

**Git** tracks every change you make to a set of files. It stores a complete
history, lets you label and revisit any past state, and lets multiple people
work on the same code without overwriting each other.

**GitHub** hosts Git repositories online, making it easy to back up your work,
share it, and collaborate.

---

## 2. Configuration

Open a terminal. On macOS/Linux: Terminal or the VS Code integrated terminal.
On Windows: **Git Bash**.

Tell Git who you are:
```bash
git config --global user.name "Your Name"
git config --global user.email "you@example.com"
```
Use the same email address as your GitHub account.

Set VS Code as your default editor (used for commit messages and merge conflicts):
```bash
git config --global core.editor "code --wait"
```

If you prefer `nano` instead:
```bash
git config --global core.editor "nano"
```

Set the default branch name to `main` to match GitHub:
```bash
git config --global init.defaultBranch main
```

Verify your settings:
```bash
git config --list
```

---

## 3. Core Concepts

Three terms before anything makes sense:

**Repository (repo):** A directory that Git is tracking. Every file change
inside it can be recorded.

**Commit:** A snapshot of the repository at a point in time — a permanent save
point with a short message explaining what changed and why.

**Staging area:** A holding zone where you specify *which* changes to include
in the next commit. This lets you split unrelated edits into separate commits.

The flow is always:

```
make changes  →  git add (stage)  →  git commit (snapshot)
```

```
Working directory         Staging area           Repository
─────────────────         ────────────           ──────────
  edit files         →    git add <file>    →    git commit
                          (mark for commit)       (saved to history)
```

---

## 4. Your First Repository

Create a new directory and initialize a Git repository inside it:

```bash
mkdir astro-notes
cd astro-notes
git init
```

Git creates a hidden `.git/` folder inside `astro-notes/`. This is where all
history is stored — don't modify it.

### Create a file

Open VS Code in the current directory:
```bash
code .
```

Create a new file called `notes.txt`. Add a line or two — anything works, for
example:

```
Wien's displacement law: lambda_max * T = 2.898e6 nm K
Stefan-Boltzmann law: L = 4 * pi * R^2 * sigma * T^4
```

Save the file.

### Stage and commit

Check what Git sees:
```bash
git status
```

`notes.txt` appears as an **untracked file** — Git knows it exists but is not
yet recording changes to it.

Stage the file:
```bash
git add notes.txt
```

Run `git status` again. The file is now listed under *"Changes to be
committed"* — it is in the staging area.

Make the first commit:
```bash
git commit -m "Add notes file"
```

The `-m` flag lets you write the message inline. Good commit messages are short
and describe *what changed and why*, not just what the code does.

---

## 5. History and Differences

### Viewing history

```bash
git log
```

You'll see the commit hash (a long hex string that uniquely identifies the
commit), your name, the date, and the message. For a compact view:

```bash
git log --oneline
```

### Making a second commit

Edit `notes.txt` — add or change a line. Save the file.

```bash
git status
```

The file is now *modified*. Before staging it, see exactly what changed:

```bash
git diff
```

Lines starting with `+` are additions; lines starting with `-` are deletions.

Stage and commit:
```bash
git add notes.txt
git commit -m "Add Planck function note"
```

Run `git log --oneline` again — you now have two commits in the history.

---

## 6. Pushing to GitHub

Your repository exists only on your local machine. Pushing it to GitHub backs
it up and makes it shareable.

### Create a repository on GitHub

1. Go to [github.com](https://github.com) and sign in.
2. Click **+** (top right) → **New repository**.
3. Name it `astro-notes`.
4. Leave all options at their defaults — **do not** add a README, `.gitignore`,
   or license. If you add any of these, the initial histories will conflict.
5. Click **Create repository**.

GitHub shows a setup page with your repository URL:
```
https://github.com/YOUR_USERNAME/astro-notes.git
```

### Connect your local repo and push

```bash
git remote add origin https://github.com/YOUR_USERNAME/astro-notes.git
git branch -M main
git push -u origin main
```

- `git remote add origin ...` — registers the GitHub URL under the name
  `origin`. This name is the conventional label for the primary remote.
- `git branch -M main` — renames the local branch to `main`.
- `git push -u origin main` — uploads your commits. The `-u` flag sets the
  default, so future pushes only require `git push`.

Reload your GitHub page — your file and both commits are there.

---

## 7. Cloning and Forking

### Cloning

**Cloning** downloads a complete copy of a repository from GitHub, including
its full history.

```bash
cd ~
git clone https://github.com/moustakas/github-tutorial.git
cd github-tutorial
ls
git log --oneline
```

Run the script:
```bash
python blackbody.py                       # the Sun (5778 K, default)
python blackbody.py --temperature 3000    # cool M-dwarf star
python blackbody.py --temperature 30000   # hot O-type star
python blackbody.py --help
```

### Forking

You cannot push changes back to `moustakas/github-tutorial` — you don't have
write access to it. To contribute, you first **fork** the repository: this
creates your own copy on GitHub that you fully control.

1. Go to [github.com/moustakas/github-tutorial](https://github.com/moustakas/github-tutorial).
2. Click **Fork** (top right) → **Create fork**.
3. Now clone *your fork* (replace `YOUR_USERNAME`):

```bash
cd ~
git clone https://github.com/YOUR_USERNAME/github-tutorial.git
cd github-tutorial
```

Verify that the remote points to your fork:
```bash
git remote -v
```

---

## 8. Branches

A **branch** is an independent line of development. The default branch is
`main`. When you want to add a feature or fix something, you create a new
branch, do the work there, and merge it back into `main` when it's ready.
This keeps `main` stable and lets multiple people work in parallel without
interfering with each other.

### Create a branch

```bash
git checkout -b add-savefig
```

`-b` creates the branch and switches to it. Confirm which branch you're on:
```bash
git branch
```

The asterisk marks the active branch.

### Make a change

Open `blackbody.py` in VS Code:
```bash
code blackbody.py
```

Add a `--savefig` argument so the plot can be saved to a file instead of
displayed on screen.

In the `argparse` block, after the `--wavemax` argument, add:
```python
parser.add_argument('--savefig', type=str, default=None,
                    help='Save the figure to this filename (e.g., spectrum.png)')
```

Replace the `plt.show()` line at the end of `main()` with:
```python
if args.savefig:
    plt.savefig(args.savefig)
    print(f'Figure saved to {args.savefig}')
else:
    plt.show()
```

Test it:
```bash
python blackbody.py --savefig sun.png
ls *.png
```

### Commit and push the branch

```bash
git add blackbody.py
git commit -m "Add --savefig argument to save figure to file"
git push origin add-savefig
```

---

## 9. Pull Requests and Merging

Your change lives on your fork, on the `add-savefig` branch. A **pull request
(PR)** is how you formally propose that those changes be merged into another
repository — in this case, back into `moustakas/github-tutorial`.

### Open a pull request

1. Go to your fork on GitHub: `github.com/YOUR_USERNAME/github-tutorial`.
2. GitHub displays a banner: *"add-savefig had recent pushes — Compare & pull
   request."* Click it.
3. Confirm the **base repository** is `moustakas/github-tutorial`, base branch
   `main`; the **head repository** is your fork, branch `add-savefig`.
4. Write a short title (e.g., *"Add --savefig argument"*) and a sentence or two
   describing what the PR does and why.
5. Click **Create pull request**.

### Review and merge

The repository owner can comment on specific lines of your diff, request
changes, or approve the PR. Once approved, clicking **Merge pull request**
incorporates your changes into `main`.

### Sync after a merge

After the merge, bring your local `main` up to date:
```bash
git checkout main
git pull origin main
```

---

## 10. Daily Workflow

The commands you will use in nearly every session:

```bash
git status                    # what has changed?
git diff                      # show exact line-by-line changes
git add <file>                # stage a file
git commit -m "message"       # snapshot staged changes
git push                      # upload commits to GitHub
git pull                      # download commits from GitHub
git log --oneline             # compact history
git checkout -b <branchname>  # create and switch to a new branch
```

A typical session looks like:

```bash
git pull                          # start by syncing with GitHub
# ... edit files ...
git status
git add <changed files>
git commit -m "Describe the change"
git push
```

---

## Going Further

Topics not covered today, in rough order of usefulness:

- **`.gitignore`** — tell Git to ignore generated files (e.g., `*.png`,
  `__pycache__/`)
- **`git stash`** — temporarily set aside uncommitted changes
- **Merge conflicts** — what happens when two people edit the same line, and
  how to resolve it
- **`git rebase`** — an alternative to merging that produces a cleaner linear
  history
- **SSH authentication** — a faster, token-free alternative to HTTPS
- **`git blame`** — see who last changed each line of a file, and when
- **GitHub Actions** — automated workflows triggered by commits (running tests,
  building documentation)
- **Protected branches and code review workflows** — how larger teams use PRs
  to enforce quality

---

## Additional Resources

- [git - the simple guide](https://rogerdudler.github.io/git-guide/) — a concise, no-frills one-page reference
- [Version Control with Git](https://swcarpentry.github.io/git-novice/) — a thorough lesson from Software Carpentry
- [Pro Git](https://git-scm.com/book/en/v2) — the complete reference book, freely available online
- [GitHub Docs](https://docs.github.com) — official documentation for all GitHub features
