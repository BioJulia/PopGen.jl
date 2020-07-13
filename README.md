![error_cactus](/static/img/logo_banner.png)
![build_status](https://img.shields.io/github/workflow/status/pdimens/PopGen.jl/documentation?label=Build%20Status&logo=GitHub&style=for-the-badge)
# Website
The documentation website is built using [Docusaurus 2](https://v2.docusaurus.io/), a modern static website generator.

There are two easy ways to edit this documentation:

## 1. The Easier Way - Submit Pull Requests
We have GitHub Actions configured to auto-deploy the documentation site when changes are pushed to the `documentation` branch. That means as long as we edit the source content, building and deployment will be handled automatically!

## 2. The More Technical Way - Clone locally
1. Clone the respository onto your system with `git clone https://github.com/pdimens/PopGen.jl.git` and switch to the `documentation` branch.
2. Install the correct NodeJS modules by navigating to the repository folder and using the command `yarn install`, which will parse the `package.json` file and install all the necessary NodeJS modules into the directory.
    1. You will need NodeJS and Yarn installed on your system for this to work.
3. Use `yarn start` to start a local live-reloading development server in a browser window. Most changes are reflected live without having to restart the server.
4. Make your changes.
5. Since GitHub Actions takes care of building and deploying, you just need to submit your changes as a Pull Request and things magically work. **this is why we generally recommend "The Easier Way".**
