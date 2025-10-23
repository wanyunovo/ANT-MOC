# Contributing to ANT-MOC

Thanks for taking the time to contribute!

## Indentation

- Use UNIX line endings when committing any change.

- Do not introduce TABs into source files.

- Each line of the code has an indentation of 2 spaces. If you are familiar with Vim, you can type the
following command to enable the local settings for Vim

```bash
$ echo "set secure exrc" >> $HOME/.vimrc
```

## Git convention

Each commit must have a readable message formatted according to
[Commit Message Guidlines for AngularJS](https://github.com/angular/angular/blob/master/CONTRIBUTING.md#commit).

For example, a multiline commit message seems like

```
feat(geometry): add hexagonal lattices to CSG

Now the code can read and construct hexagonal lattices
```

The commonly used message types are

- **build**: Changes that affect the build system or external dependencies (example scopes: cmake, hdf5, hip)
- **ci**: Changes to our CI configuration files and scripts (example scopes: gitlab, codeclimate)
- **docs**: Documentation only changes
- **feat**: A new feature
- **fix**: A bug fix
- **perf**: A code change that improves performance
- **refactor**: A code change that neither fixes a bug nor adds a feature
- **style**: Changes that do not affect the meaning of the code (white-space, formatting, missing semi-colons, etc)
- **test**: Adding missing tests or correcting existing tests

## Multiline messages

Basically, there are two ways to commit multiline messages

- Use bash here-doc;
- Use external tools or plugins integrated with your IDEs.

For the latter, we recommend [Commitizen](https://www.npmjs.com/package/commitizen). You can access it from either npm or VSCode extensions.

For here-doc, the options `-F` with `-` will make `git commit` read the message from the standard input.

```bash
$ git commit -F- << END
> feat(geometry): add hexagonal lattices to CSG
>
> Now the code can read and structured hexagonal lattices
> END
```

This command will generate a here-doc to the standard input and then feed the input to `commit`.
For more information about here-docs for bash, see [Here Documents](http://www.gnu.org/software/bash/manual/bash.html#Here-Documents).

